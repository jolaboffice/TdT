#!/usr/bin/env python3
"""
GUI Point Picker for SEM Image Stitching

Displays consecutive image pairs side by side.
Click the same physical point in each image to define overlap.
Results are saved to config.json for use with stitching.py.

Usage:
    python pick_points.py /path/to/tiles/
    python pick_points.py /path/to/tiles/ --ribbon 70

Controls:
    Left click       : Select overlap point (A first, then B)
    Right click drag  : Pan image
    Scroll wheel      : Zoom in/out at cursor position
    Enter             : Confirm pair and move to next
    Left/Right arrow  : Navigate between pairs
    'g'               : Go to specific pair number
    'z'               : Undo last point
    'r'               : Reset zoom to full view
    'c'               : Reset contrast/brightness
    'q'               : Quit and save progress

Dependencies:
    pip install numpy matplotlib tifffile
"""

import numpy as np
import json
import argparse
from pathlib import Path
import tifffile
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Slider


class PanZoomHandler:
    """Handles right-click pan and scroll zoom for a matplotlib axes."""

    def __init__(self, ax, zoom_factor=1.3):
        self.ax = ax
        self.zoom_factor = zoom_factor
        self._dragging = False
        self._drag_start = None  # (x_data, y_data) at drag start
        self._xlim_start = None
        self._ylim_start = None

    def on_press(self, event):
        if event.inaxes != self.ax:
            return False
        if event.button == 3:  # right click
            self._dragging = True
            self._drag_start = (event.x, event.y)  # pixel coords
            self._xlim_start = self.ax.get_xlim()
            self._ylim_start = self.ax.get_ylim()
            return True
        return False

    def on_release(self, event):
        if event.button == 3:
            self._dragging = False
            return True
        return False

    def on_motion(self, event):
        if not self._dragging or event.x is None or event.y is None:
            return False

        # Compute pixel delta
        dx_px = event.x - self._drag_start[0]
        dy_px = event.y - self._drag_start[1]

        # Convert pixel delta to data coordinates
        xlim = self._xlim_start
        ylim = self._ylim_start
        ax_bbox = self.ax.get_window_extent()
        ax_w = ax_bbox.width
        ax_h = ax_bbox.height

        dx_data = -(dx_px / ax_w) * (xlim[1] - xlim[0])
        dy_data = (dy_px / ax_h) * (ylim[1] - ylim[0])  # y is inverted for images

        self.ax.set_xlim(xlim[0] + dx_data, xlim[1] + dx_data)
        self.ax.set_ylim(ylim[0] + dy_data, ylim[1] + dy_data)
        return True

    def on_scroll(self, event):
        if event.inaxes != self.ax:
            return False
        if event.xdata is None or event.ydata is None:
            return False

        if event.button == "up":
            scale = 1.0 / self.zoom_factor
        elif event.button == "down":
            scale = self.zoom_factor
        else:
            return False

        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        # Zoom centered on cursor position
        cx, cy = event.xdata, event.ydata
        new_w = (xlim[1] - xlim[0]) * scale
        new_h = (ylim[1] - ylim[0]) * scale

        # Fraction of cursor within current view
        fx = (cx - xlim[0]) / (xlim[1] - xlim[0])
        fy = (cy - ylim[0]) / (ylim[1] - ylim[0])

        self.ax.set_xlim(cx - fx * new_w, cx + (1 - fx) * new_w)
        self.ax.set_ylim(cy - fy * new_h, cy + (1 - fy) * new_h)
        return True


class PointPicker:
    def __init__(self, files, ribbon_height=70, output="config.json",
                 existing_pairs=None):
        self.files = files
        self.names = [f.name for f in files]
        self.ribbon_height = ribbon_height
        self.output = output

        self.total_pairs = len(files) - 1
        self.current_pair = 0
        self.clicks = []
        self.markers = []

        # Index-based pairs list: None = not yet defined
        self.pairs = [None] * self.total_pairs
        self._load_existing_pairs(existing_pairs)

        # Find first undefined pair to start at
        start_idx = self._first_undefined_pair()

        # Contrast/brightness state
        self.contrast_a = 1.0
        self.brightness_a = 0.0
        self.contrast_b = 1.0
        self.brightness_b = 0.0
        self.img_display_a = None
        self.img_display_b = None
        self.raw_img_a = None
        self.raw_img_b = None

        self.fig, (self.ax_a, self.ax_b) = plt.subplots(1, 2, figsize=(16, 8))
        self.fig.subplots_adjust(bottom=0.15)
        self.fig.canvas.manager.set_window_title("SEM Stitching Point Picker")

        # Sliders for contrast/brightness
        # Image A sliders (left side)
        ax_contrast_a = plt.axes([0.08, 0.06, 0.35, 0.025])
        ax_brightness_a = plt.axes([0.08, 0.025, 0.35, 0.025])
        # Image B sliders (right side)
        ax_contrast_b = plt.axes([0.57, 0.06, 0.35, 0.025])
        ax_brightness_b = plt.axes([0.57, 0.025, 0.35, 0.025])

        self.slider_contrast_a = Slider(ax_contrast_a, 'A Contrast', 0.1, 5.0,
                                         valinit=1.0, valstep=0.1)
        self.slider_brightness_a = Slider(ax_brightness_a, 'A Bright', -1.0, 1.0,
                                           valinit=0.0, valstep=0.05)
        self.slider_contrast_b = Slider(ax_contrast_b, 'B Contrast', 0.1, 5.0,
                                         valinit=1.0, valstep=0.1)
        self.slider_brightness_b = Slider(ax_brightness_b, 'B Bright', -1.0, 1.0,
                                           valinit=0.0, valstep=0.05)

        self.slider_contrast_a.on_changed(self._on_contrast_a_change)
        self.slider_brightness_a.on_changed(self._on_brightness_a_change)
        self.slider_contrast_b.on_changed(self._on_contrast_b_change)
        self.slider_brightness_b.on_changed(self._on_brightness_b_change)

        # Pan/zoom handlers for each axes
        self.pz_a = PanZoomHandler(self.ax_a)
        self.pz_b = PanZoomHandler(self.ax_b)

        self.fig.canvas.mpl_connect("button_press_event", self._on_press)
        self.fig.canvas.mpl_connect("button_release_event", self._on_release)
        self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.fig.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        self._load_pair(start_idx)
        plt.show()

    def _load_existing_pairs(self, existing_pairs):
        """Map loaded pairs from config into self.pairs by matching image names."""
        if not existing_pairs:
            return
        # Build lookup: (image_a, image_b) -> pair dict
        for pair in existing_pairs:
            img_a = pair.get("image_a")
            img_b = pair.get("image_b")
            for i in range(self.total_pairs):
                if self.names[i] == img_a and self.names[i + 1] == img_b:
                    self.pairs[i] = pair
                    break
        loaded = sum(1 for p in self.pairs if p is not None)
        print(f"  Loaded {loaded}/{self.total_pairs} existing pairs from config.")

    def _first_undefined_pair(self):
        """Return index of first pair without coordinates, or 0."""
        for i, p in enumerate(self.pairs):
            if p is None:
                return i
        return 0

    # --- Contrast/Brightness ---

    def _apply_contrast_brightness(self, img, contrast, brightness):
        """Apply contrast and brightness to image, return 0-1 normalized."""
        img_norm = (img - img.min()) / (img.max() - img.min() + 1e-10)
        adjusted = (img_norm - 0.5) * contrast + 0.5 + brightness
        return np.clip(adjusted, 0, 1)

    def _on_contrast_a_change(self, val):
        self.contrast_a = val
        if self.raw_img_a is not None and self.img_display_a is not None:
            adjusted = self._apply_contrast_brightness(
                self.raw_img_a, self.contrast_a, self.brightness_a)
            self.img_display_a.set_data(adjusted)
            self.fig.canvas.draw_idle()

    def _on_brightness_a_change(self, val):
        self.brightness_a = val
        if self.raw_img_a is not None and self.img_display_a is not None:
            adjusted = self._apply_contrast_brightness(
                self.raw_img_a, self.contrast_a, self.brightness_a)
            self.img_display_a.set_data(adjusted)
            self.fig.canvas.draw_idle()

    def _on_contrast_b_change(self, val):
        self.contrast_b = val
        if self.raw_img_b is not None and self.img_display_b is not None:
            adjusted = self._apply_contrast_brightness(
                self.raw_img_b, self.contrast_b, self.brightness_b)
            self.img_display_b.set_data(adjusted)
            self.fig.canvas.draw_idle()

    def _on_brightness_b_change(self, val):
        self.brightness_b = val
        if self.raw_img_b is not None and self.img_display_b is not None:
            adjusted = self._apply_contrast_brightness(
                self.raw_img_b, self.contrast_b, self.brightness_b)
            self.img_display_b.set_data(adjusted)
            self.fig.canvas.draw_idle()

    # --- Event routing ---

    def _on_press(self, event):
        # Right click -> pan
        if event.button == 3:
            self.pz_a.on_press(event)
            self.pz_b.on_press(event)
            return

        # Left click -> point selection
        if event.button == 1:
            self._pick_point(event)

    def _on_release(self, event):
        if event.button == 3:
            redraw = self.pz_a.on_release(event) or self.pz_b.on_release(event)
            if redraw:
                self.fig.canvas.draw_idle()

    def _on_motion(self, event):
        redraw = self.pz_a.on_motion(event) or self.pz_b.on_motion(event)
        if redraw:
            self.fig.canvas.draw_idle()

    def _on_scroll(self, event):
        redraw = self.pz_a.on_scroll(event) or self.pz_b.on_scroll(event)
        if redraw:
            self.fig.canvas.draw_idle()

    # --- Pair loading ---

    def _load_pair(self, idx):
        if idx < 0 or idx >= self.total_pairs:
            return

        self.current_pair = idx
        self.clicks = []
        self._clear_markers()

        self.ax_a.clear()
        self.ax_b.clear()

        img_a = tifffile.imread(str(self.files[idx]))
        img_b = tifffile.imread(str(self.files[idx + 1]))

        # Convert to grayscale if RGB
        if img_a.ndim == 3:
            img_a = np.mean(img_a[:, :, :3], axis=2)
        if img_b.ndim == 3:
            img_b = np.mean(img_b[:, :, :3], axis=2)

        # Store raw images
        self.raw_img_a = img_a.astype(np.float64)
        self.raw_img_b = img_b.astype(np.float64)

        # Reset sliders to default
        self.contrast_a = 1.0
        self.brightness_a = 0.0
        self.contrast_b = 1.0
        self.brightness_b = 0.0
        self.slider_contrast_a.set_val(1.0)
        self.slider_brightness_a.set_val(0.0)
        self.slider_contrast_b.set_val(1.0)
        self.slider_brightness_b.set_val(0.0)

        # Display with contrast/brightness applied
        adjusted_a = self._apply_contrast_brightness(
            self.raw_img_a, self.contrast_a, self.brightness_a)
        adjusted_b = self._apply_contrast_brightness(
            self.raw_img_b, self.contrast_b, self.brightness_b)

        self.img_display_a = self.ax_a.imshow(adjusted_a, cmap="gray", vmin=0, vmax=1)
        self.img_display_b = self.ax_b.imshow(adjusted_b, cmap="gray", vmin=0, vmax=1)

        # Ribbon boundary
        h = img_a.shape[0]
        ribbon_y = h - self.ribbon_height
        self.ax_a.axhline(ribbon_y, color="red", linewidth=0.5,
                          linestyle="--", alpha=0.5)
        self.ax_b.axhline(ribbon_y, color="red", linewidth=0.5,
                          linestyle="--", alpha=0.5)

        # Show existing points if this pair was already defined
        saved = self.pairs[idx]
        if saved is not None:
            pa = saved["point_a"]
            pb = saved["point_b"]
            self.clicks = [tuple(pa), tuple(pb)]
            self._add_marker(self.ax_a, pa[0], pa[1], "yellow",
                             f"A: ({pa[0]}, {pa[1]})")
            self._add_marker(self.ax_b, pb[0], pb[1], "orange",
                             f"B: ({pb[0]}, {pb[1]})")
            self.ax_a.set_title(
                f"Image A: {self.names[idx]}\n"
                f"Saved A: ({pa[0]}, {pa[1]})  [click to re-pick]",
                fontsize=9, color="orange")
            self.ax_b.set_title(
                f"Image B: {self.names[idx+1]}\n"
                f"Saved B: ({pb[0]}, {pb[1]})  [click to re-pick]",
                fontsize=9, color="orange")
        else:
            self.ax_a.set_title(
                f"Image A: {self.names[idx]}\n(click overlap point)",
                fontsize=9)
            self.ax_b.set_title(
                f"Image B: {self.names[idx+1]}\n(click same point)",
                fontsize=9)

        defined = sum(1 for p in self.pairs if p is not None)
        status = "SAVED" if saved is not None else "empty"
        self.fig.suptitle(
            f"Pair {idx+1}/{self.total_pairs} [{status}]  |  "
            f"Done: {defined}/{self.total_pairs}\n"
            f"Enter=confirm  "
            f"\u2190\u2192=navigate  g=goto  z=undo  r=reset  q=quit+save",
            fontsize=10)

        self.fig.canvas.draw()

    # --- Point selection ---

    def _pick_point(self, event):
        if event.inaxes not in (self.ax_a, self.ax_b):
            return
        if event.xdata is None or event.ydata is None:
            return

        x, y = round(event.xdata), round(event.ydata)

        # If editing a saved pair, clicking A resets both points
        if event.inaxes == self.ax_a:
            if len(self.clicks) >= 1:
                # Reset: remove all existing markers and start fresh
                self._clear_markers()
                self.clicks = []
            self.clicks.append((x, y))
            self._add_marker(self.ax_a, x, y, "lime",
                             f"A: ({x}, {y})")
            self.ax_a.set_title(
                f"Image A: {self.names[self.current_pair]}\n"
                f"Point A: ({x}, {y})", fontsize=9, color="green")
            self.ax_b.set_title(
                f"Image B: {self.names[self.current_pair+1]}\n"
                f"(click same point)", fontsize=9)
            self.fig.canvas.draw_idle()

        elif event.inaxes == self.ax_b and len(self.clicks) == 1:
            self.clicks.append((x, y))
            self._add_marker(self.ax_b, x, y, "cyan",
                             f"B: ({x}, {y})")
            self.ax_b.set_title(
                f"Image B: {self.names[self.current_pair+1]}\n"
                f"Point B: ({x}, {y}) -- press Enter",
                fontsize=9, color="green")
            self.fig.canvas.draw_idle()

    def _add_marker(self, ax, x, y, color, label):
        # Scale marker size to current zoom level
        xlim = ax.get_xlim()
        view_w = abs(xlim[1] - xlim[0])
        radius = max(3, view_w * 0.008)

        circle = Circle((x, y), radius, fill=False,
                         edgecolor=color, linewidth=2)
        ax.add_patch(circle)
        txt = ax.text(x + radius * 1.5, y, label, color=color,
                      fontsize=8, fontweight="bold",
                      bbox=dict(boxstyle="round,pad=0.2",
                                facecolor="black", alpha=0.7))
        # Crosshair lines for precise visibility
        line_len = radius * 3
        lh = ax.plot([x - line_len, x + line_len], [y, y],
                     color=color, linewidth=0.8, alpha=0.6)[0]
        lv = ax.plot([x, x], [y - line_len, y + line_len],
                     color=color, linewidth=0.8, alpha=0.6)[0]
        self.markers.extend([circle, txt, lh, lv])

    # --- Key events ---

    def _on_key(self, event):
        if event.key == "enter":
            if len(self.clicks) == 2:
                self.pairs[self.current_pair] = {
                    "image_a": self.names[self.current_pair],
                    "image_b": self.names[self.current_pair + 1],
                    "point_a": list(self.clicks[0]),
                    "point_b": list(self.clicks[1])
                }
                print(f"  Pair {self.current_pair+1}: "
                      f"A({self.clicks[0][0]},{self.clicks[0][1]}) "
                      f"B({self.clicks[1][0]},{self.clicks[1][1]})")
                # Move to next undefined pair, or next pair
                nxt = self._first_undefined_pair()
                if self.pairs[nxt] is not None:
                    # All defined; move to next sequentially
                    nxt = self.current_pair + 1
                if nxt >= self.total_pairs:
                    self._save_and_exit()
                else:
                    self._load_pair(nxt)

        elif event.key == "left":
            self._load_pair(self.current_pair - 1)

        elif event.key == "right":
            self._load_pair(self.current_pair + 1)

        elif event.key == "g":
            # Jump to specific pair number (terminal input)
            try:
                n = int(input(f"  Go to pair (1-{self.total_pairs}): "))
                if 1 <= n <= self.total_pairs:
                    self._load_pair(n - 1)
                else:
                    print(f"  Invalid: must be 1-{self.total_pairs}")
            except (ValueError, EOFError):
                print("  Cancelled.")

        elif event.key == "z":
            if self.clicks:
                self.clicks.pop()
                # Remove last 4 artists (circle, text, hline, vline)
                for _ in range(min(4, len(self.markers))):
                    m = self.markers.pop()
                    m.remove()
                self.fig.canvas.draw_idle()

        elif event.key == "r":
            self.ax_a.autoscale()
            self.ax_b.autoscale()
            self.fig.canvas.draw_idle()

        elif event.key == "q":
            self._save_and_exit()

        elif event.key == "c":
            # Reset contrast/brightness to default
            self.slider_contrast_a.set_val(1.0)
            self.slider_brightness_a.set_val(0.0)
            self.slider_contrast_b.set_val(1.0)
            self.slider_brightness_b.set_val(0.0)

    # --- Cleanup ---

    def _clear_markers(self):
        for m in self.markers:
            try:
                m.remove()
            except Exception:
                pass
        self.markers = []

    def _save_and_exit(self):
        defined = [p for p in self.pairs if p is not None]
        if not defined:
            print("No pairs selected. Nothing to save.")
            plt.close(self.fig)
            return

        config = {
            "images": self.names,
            "pairs": defined,
            "parameters": {
                "search_range_px": 8,
                "ribbon_height_px": self.ribbon_height
            },
            "output": "stitched.tif"
        }

        with open(self.output, "w") as f:
            json.dump(config, f, indent=4)
        print(f"\nConfig saved: {self.output}")
        print(f"  {len(defined)}/{self.total_pairs} pairs defined")
        print("\nRun stitching:")
        print(f"  python stitching.py {self.output}")
        plt.close(self.fig)


def main():
    parser = argparse.ArgumentParser(
        description="GUI point picker for SEM stitching")
    parser.add_argument("tiles_dir", help="Directory containing .tif tiles")
    parser.add_argument("--ribbon", type=int, default=70,
                        help="Ribbon height in pixels (default: 70)")
    args = parser.parse_args()

    tiles_dir = Path(args.tiles_dir).resolve()
    files = sorted(tiles_dir.glob("*.tif"))
    if len(files) < 2:
        print(f"Need at least 2 .tif files in {tiles_dir}, found {len(files)}")
        return

    # Save config inside the tiles directory, named config_<foldername>.json
    folder_name = tiles_dir.name
    output_path = str(tiles_dir / f"config_{folder_name}.json")

    # Load existing config if present
    existing_pairs = None
    if Path(output_path).exists():
        try:
            with open(output_path) as f:
                existing = json.load(f)
            existing_pairs = existing.get("pairs", [])
            print(f"Found existing config: {output_path}")
        except (json.JSONDecodeError, KeyError):
            print("Warning: Could not parse existing config, starting fresh.")

    print(f"Found {len(files)} tiles. {len(files)-1} pairs to define.")
    print(f"Config will be saved to: {output_path}")
    print("Controls:")
    print("  Left click       : Select point (A first, then B)")
    print("  Right click drag : Pan image")
    print("  Scroll wheel     : Zoom in/out at cursor")
    print("  Enter            : Confirm pair, next")
    print("  Left/Right arrow : Navigate between pairs")
    print("  g                : Go to specific pair number")
    print("  z                : Undo last point")
    print("  r                : Reset zoom")
    print("  c                : Reset contrast/brightness")
    print("  q                : Quit and save\n")

    PointPicker(files, ribbon_height=args.ribbon, output=output_path,
                existing_pairs=existing_pairs)


if __name__ == "__main__":
    main()
