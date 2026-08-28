#!/usr/bin/env python3
"""
SEM Image Stitching Tool

Stitches multiple TIF images using gradient-based normalized cross-correlation.
Designed for negative-stained SEM images with beam damage in overlap regions.

Features:
  - Gradient-based NCC (robust to beam damage brightness changes)
  - Manual overlap point specification per image pair
  - Fine alignment: +/-6px translation, +/-0.5 degree rotation
  - Subpixel refinement (parabolic fitting)
  - Bottom ribbon exclusion from alignment computation
  - Lossless pixel placement (earlier-scanned image prioritized in overlap)

Usage:
    python stitching.py --create-config           # Generate template config.json
    python stitching.py config.json               # Run stitching

Dependencies:
    pip install numpy scipy tifffile
"""

import numpy as np
import json
import sys
import argparse
from pathlib import Path
from collections import deque
from scipy import ndimage  # for Sobel gradient
import tifffile


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def create_template_config(output_path):
    """Generate a template config.json for the user to fill in."""
    config = {
        "images": [
            "image_001.tif",
            "image_002.tif",
            "image_003.tif"
        ],
        "pairs": [
            {
                "image_a": "image_001.tif",
                "image_b": "image_002.tif",
                "point_a": [1000, 500],
                "point_b": [50, 500],
                "_comment": "point_a/b are [x, y] coordinates of the same physical point in each image"
            },
            {
                "image_a": "image_001.tif",
                "image_b": "image_003.tif",
                "point_a": [500, 900],
                "point_b": [500, 50]
            }
        ],
        "parameters": {
            "search_range_px": 8,
            "ribbon_height_px": 70
        },
        "output": "stitched.tif"
    }
    with open(output_path, "w") as f:
        json.dump(config, f, indent=4)
    print(f"Template config saved to: {output_path}")
    print("Edit the file with your image paths and overlap points, then run:")
    print(f"  python stitching.py {output_path}")


def load_config(config_path):
    """Load and validate configuration."""
    with open(config_path) as f:
        config = json.load(f)

    required_keys = ["images", "pairs"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Config missing required key: '{key}'")

    if not config["images"]:
        raise ValueError("No images specified in config")
    if not config["pairs"]:
        raise ValueError("No pairs specified in config")

    for i, pair in enumerate(config["pairs"]):
        for field in ("image_a", "image_b", "point_a", "point_b"):
            if field not in pair:
                raise ValueError(f"Pair {i} missing field: '{field}'")
        if pair["image_a"] not in config["images"]:
            raise ValueError(f"Pair {i}: image_a '{pair['image_a']}' not in images list")
        if pair["image_b"] not in config["images"]:
            raise ValueError(f"Pair {i}: image_b '{pair['image_b']}' not in images list")

    return config


# ---------------------------------------------------------------------------
# Image Processing
# ---------------------------------------------------------------------------

def compute_gradient_magnitude(image):
    """Compute Sobel gradient magnitude. Robust to brightness changes from beam damage."""
    img = image.astype(np.float64)
    gx = ndimage.sobel(img, axis=1)
    gy = ndimage.sobel(img, axis=0)
    return np.sqrt(gx ** 2 + gy ** 2)


def compute_ncc(template, patch):
    """Normalized cross-correlation between two identically-sized arrays."""
    if template.shape != patch.shape:
        return -1.0
    t = template.ravel().astype(np.float64)
    p = patch.ravel().astype(np.float64)
    t_mean = t.mean()
    p_mean = p.mean()
    tc = t - t_mean
    pc = p - p_mean
    t_std = np.sqrt(np.dot(tc, tc))
    p_std = np.sqrt(np.dot(pc, pc))
    if t_std < 1e-10 or p_std < 1e-10:
        return 0.0
    return np.dot(tc, pc) / (t_std * p_std)


def subpixel_peak_1d(v0, v1, v2):
    """Parabolic subpixel refinement from three consecutive NCC values.
    Returns fractional offset from the centre value (clamped to -0.5 .. +0.5)."""
    denom = 2.0 * (2.0 * v1 - v0 - v2)
    if abs(denom) < 1e-12:
        return 0.0
    return float(np.clip((v0 - v2) / denom, -0.5, 0.5))


# ---------------------------------------------------------------------------
# Pair-wise Alignment
# ---------------------------------------------------------------------------

def find_optimal_offset(img_a, img_b, approx_offset,
                        search_range=8, ribbon_height=75):
    """
    Find optimal alignment of img_b relative to img_a using gradient-based NCC.
    Translation-only search within ±search_range pixels.

    Parameters
    ----------
    img_a, img_b : 2-D arrays (grayscale)
    approx_offset : (ox, oy) approximate position of B's origin in A's frame
    search_range : half-width of translation search in pixels
    ribbon_height: height (px) of the bottom info ribbon to exclude

    Returns
    -------
    (offset_x, offset_y, ncc_score)
    """
    ox_approx, oy_approx = approx_offset
    h_a, w_a = img_a.shape[:2]
    h_b, w_b = img_b.shape[:2]

    eh_a = h_a - ribbon_height
    eh_b = h_b - ribbon_height

    iox = int(round(ox_approx))
    ioy = int(round(oy_approx))

    # Overlap region in A's pixel coords (excluding ribbon)
    oa_x0 = max(0, iox)
    oa_x1 = min(w_a, iox + w_b)
    oa_y0 = max(0, ioy)
    oa_y1 = min(eh_a, ioy + eh_b)

    overlap_w = oa_x1 - oa_x0
    overlap_h = oa_y1 - oa_y0
    margin = search_range

    if overlap_w < margin * 4 or overlap_h < margin * 4:
        print(f"  WARNING: Overlap region too small ({overlap_w}x{overlap_h}px). "
              "Using approximate offset.")
        return (ox_approx, oy_approx, 0.0)

    # Template from A — inset by margin so B has room to shift
    tmpl_x0 = oa_x0 + margin
    tmpl_x1 = oa_x1 - margin
    tmpl_y0 = oa_y0 + margin
    tmpl_y1 = oa_y1 - margin

    if tmpl_x1 <= tmpl_x0 or tmpl_y1 <= tmpl_y0:
        print("  WARNING: Template region empty after margin. Using approximate offset.")
        return (ox_approx, oy_approx, 0.0)

    template = img_a[tmpl_y0:tmpl_y1, tmpl_x0:tmpl_x1]
    th, tw = template.shape

    # Search region from B (corresponding area + margins for shifting)
    sb_x0 = max(0, tmpl_x0 - iox - margin)
    sb_x1 = min(w_b, tmpl_x1 - iox + margin)
    sb_y0 = max(0, tmpl_y0 - ioy - margin)
    sb_y1 = min(eh_b, tmpl_y1 - ioy + margin)

    search_region = img_b[sb_y0:sb_y1, sb_x0:sb_x1]
    sh, sw = search_region.shape

    if sh < th or sw < tw:
        print(f"  WARNING: Search region ({sw}x{sh}) < template ({tw}x{th}). "
              "Using approximate offset.")
        return (ox_approx, oy_approx, 0.0)

    # Gradient magnitudes
    template_grad = compute_gradient_magnitude(template)
    search_grad = compute_gradient_magnitude(search_region)

    # Expected template position inside search region (at zero correction)
    exp_dx = (tmpl_x0 - iox) - sb_x0
    exp_dy = (tmpl_y0 - ioy) - sb_y0

    best_ncc = -2.0
    best_dx = best_dy = 0
    ncc_map = {}

    n_trans = (2 * search_range + 1) ** 2
    print(f"  Searching {n_trans} translations "
          f"(template {tw}x{th}, search {sw}x{sh}) ...")

    for dy in range(-search_range, search_range + 1):
        for dx in range(-search_range, search_range + 1):
            py = int(round(exp_dy)) + dy
            px = int(round(exp_dx)) + dx
            if py < 0 or px < 0 or py + th > sh or px + tw > sw:
                continue
            patch = search_grad[py:py + th, px:px + tw]
            val = compute_ncc(template_grad, patch)
            ncc_map[(dx, dy)] = val
            if val > best_ncc:
                best_ncc = val
                best_dx, best_dy = dx, dy

    # --- Subpixel refinement ---
    key_c = (best_dx, best_dy)
    sub_dx = 0.0
    kx_m = (best_dx - 1, best_dy)
    kx_p = (best_dx + 1, best_dy)
    if kx_m in ncc_map and kx_p in ncc_map:
        sub_dx = subpixel_peak_1d(ncc_map[kx_m], ncc_map[key_c], ncc_map[kx_p])

    sub_dy = 0.0
    ky_m = (best_dx, best_dy - 1)
    ky_p = (best_dx, best_dy + 1)
    if ky_m in ncc_map and ky_p in ncc_map:
        sub_dy = subpixel_peak_1d(ncc_map[ky_m], ncc_map[key_c], ncc_map[ky_p])

    final_ox = ox_approx + best_dx + sub_dx
    final_oy = oy_approx + best_dy + sub_dy

    print(f"  Result: offset=({final_ox:.2f}, {final_oy:.2f}), NCC={best_ncc:.4f}")
    print(f"  Correction from approx: dx={best_dx + sub_dx:+.2f}, "
          f"dy={best_dy + sub_dy:+.2f}")

    return (final_ox, final_oy, best_ncc)


# ---------------------------------------------------------------------------
# Global Position Graph
# ---------------------------------------------------------------------------

def compute_global_positions(image_names, optimized_offsets):
    """
    BFS from the first image to assign a global (x, y, angle) to every image.

    Parameters
    ----------
    image_names : list of str
    optimized_offsets : dict  {(name_a, name_b): (ox, oy, angle, ncc), ...}

    Returns
    -------
    dict  {image_name: (global_x, global_y)}
    """
    adj = {name: [] for name in image_names}
    for (na, nb), (ox, oy, _ncc) in optimized_offsets.items():
        adj[na].append((nb, ox, oy))
        adj[nb].append((na, -ox, -oy))

    positions = {}
    root = image_names[0]
    positions[root] = (0.0, 0.0)

    queue = deque([root])
    visited = {root}

    while queue:
        cur = queue.popleft()
        cx, cy = positions[cur]
        for neighbor, ox, oy in adj[cur]:
            if neighbor not in visited:
                positions[neighbor] = (cx + ox, cy + oy)
                visited.add(neighbor)
                queue.append(neighbor)

    for name in image_names:
        if name not in positions:
            print(f"  WARNING: '{name}' is not connected to any pair. Placing at origin.")
            positions[name] = (0.0, 0.0)

    return positions


# ---------------------------------------------------------------------------
# Canvas Assembly
# ---------------------------------------------------------------------------

def _gray_mean(img):
    """Mean intensity of an image region (handles RGB and grayscale)."""
    if img.ndim == 3:
        return np.mean(img[:, :, :3].astype(np.float64))
    return np.mean(img.astype(np.float64))


def create_stitched_image(images_dict, positions, image_order, ribbon_height=75):
    """
    Assemble the final stitched image.

    - Images are placed at integer pixel positions without interpolation.
    - In overlap regions, the brighter image's pixels are used
      (beam damage darkens the previously-scanned area, so the brighter
      version is the undamaged one).
    - Supports both 2-D grayscale and 3-D RGB arrays.
    - Original dtype is preserved. Empty areas are 0 (black).
    """
    min_x = min_y = float("inf")
    max_x = max_y = float("-inf")

    for name in image_order:
        gx, gy = positions[name]
        h, w = images_dict[name].shape[:2]
        ix, iy = int(round(gx)), int(round(gy))
        min_x, min_y = min(min_x, ix), min(min_y, iy)
        max_x, max_y = max(max_x, ix + w), max(max_y, iy + h)

    canvas_w = max_x - min_x
    canvas_h = max_y - min_y
    print(f"  Canvas size: {canvas_w} x {canvas_h}")

    first = images_dict[image_order[0]]
    dtype = first.dtype
    if first.ndim == 3:
        canvas = np.zeros((canvas_h, canvas_w, first.shape[2]), dtype=dtype)
    else:
        canvas = np.zeros((canvas_h, canvas_w), dtype=dtype)
    # Track which image occupies each pixel (by index, -1 = empty)
    owner = np.full((canvas_h, canvas_w), -1, dtype=np.int32)

    # Pre-compute canvas placement info for each image
    placements = {}
    for i, name in enumerate(image_order):
        gx, gy = positions[name]
        img = images_dict[name]
        h, w = img.shape[:2]
        ix = int(round(gx)) - min_x
        iy = int(round(gy)) - min_y

        src_y0 = max(0, -iy)
        src_y1 = min(h, canvas_h - iy)
        src_x0 = max(0, -ix)
        src_x1 = min(w, canvas_w - ix)
        dst_y0 = iy + src_y0
        dst_y1 = iy + src_y1
        dst_x0 = ix + src_x0
        dst_x1 = ix + src_x1
        placements[i] = (name, img, src_y0, src_y1, src_x0, src_x1,
                         dst_y0, dst_y1, dst_x0, dst_x1)

    # First pass: place all images, keeping track of owners
    for i, name in enumerate(image_order):
        name, img, sy0, sy1, sx0, sx1, dy0, dy1, dx0, dx1 = placements[i]
        empty = owner[dy0:dy1, dx0:dx1] == -1

        if img.ndim == 3:
            canvas[dy0:dy1, dx0:dx1][empty] = img[sy0:sy1, sx0:sx1][empty]
        else:
            canvas[dy0:dy1, dx0:dx1][empty] = img[sy0:sy1, sx0:sx1][empty]
        owner[dy0:dy1, dx0:dx1][empty] = i

    # Second pass: resolve overlaps by brightness
    # For each pair of overlapping images, compare mean intensity in the
    # overlap region (excluding ribbon) and let the brighter one win.
    for i in range(len(image_order)):
        for j in range(i + 1, len(image_order)):
            ni, imgi, syi0, syi1, sxi0, sxi1, dyi0, dyi1, dxi0, dxi1 = placements[i]
            nj, imgj, syj0, syj1, sxj0, sxj1, dyj0, dyj1, dxj0, dxj1 = placements[j]

            # Overlap rectangle on canvas
            ov_y0 = max(dyi0, dyj0)
            ov_y1 = min(dyi1, dyj1)
            ov_x0 = max(dxi0, dxj0)
            ov_x1 = min(dxi1, dxj1)

            if ov_y0 >= ov_y1 or ov_x0 >= ov_x1:
                continue  # no overlap

            # Corresponding source regions (exclude ribbon from comparison)
            hi = images_dict[ni].shape[0]
            hj = images_dict[nj].shape[0]

            # Overlap in image i's source coords
            si_y0 = syi0 + (ov_y0 - dyi0)
            si_y1 = syi0 + (ov_y1 - dyi0)
            si_x0 = sxi0 + (ov_x0 - dxi0)
            si_x1 = sxi0 + (ov_x1 - dxi0)

            # Overlap in image j's source coords
            sj_y0 = syj0 + (ov_y0 - dyj0)
            sj_y1 = syj0 + (ov_y1 - dyj0)
            sj_x0 = sxj0 + (ov_x0 - dxj0)
            sj_x1 = sxj0 + (ov_x1 - dxj0)

            # Exclude ribbon rows from brightness comparison
            ribbon_limit_i = hi - ribbon_height
            ribbon_limit_j = hj - ribbon_height
            comp_si_y1 = min(si_y1, ribbon_limit_i)
            comp_sj_y1 = min(sj_y1, ribbon_limit_j)

            if comp_si_y1 <= si_y0 or comp_sj_y1 <= sj_y0:
                continue  # overlap is entirely in ribbon area

            patch_i = imgi[si_y0:comp_si_y1, si_x0:si_x1]
            patch_j = imgj[sj_y0:comp_sj_y1, sj_x0:sj_x1]

            mean_i = _gray_mean(patch_i)
            mean_j = _gray_mean(patch_j)

            # Brighter = undamaged = should be on top
            if mean_i >= mean_j:
                winner = i
            else:
                winner = j

            wname, wimg, wsy0, _, wsx0, _, wdy0, _, wdx0, _ = placements[winner]

            # Overwrite overlap region with winner's pixels
            ov_src_y0 = wsy0 + (ov_y0 - wdy0)
            ov_src_y1 = wsy0 + (ov_y1 - wdy0)
            ov_src_x0 = wsx0 + (ov_x0 - wdx0)
            ov_src_x1 = wsx0 + (ov_x1 - wdx0)

            canvas[ov_y0:ov_y1, ov_x0:ov_x1] = \
                wimg[ov_src_y0:ov_src_y1, ov_src_x0:ov_src_x1]
            owner[ov_y0:ov_y1, ov_x0:ov_x1] = winner

            print(f"  Overlap {image_order[i]} vs {image_order[j]}: "
                  f"mean {mean_i:.1f} vs {mean_j:.1f} -> {image_order[winner]} wins")

    for i, name in enumerate(image_order):
        _, _, _, _, _, _, dy0, dy1, dx0, dx1 = placements[i]
        ix = dx0
        iy = dy0
        print(f"  Placed {name} at ({ix}, {iy})")

    return canvas


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="SEM Image Stitching Tool — gradient-based NCC alignment")
    parser.add_argument("config", nargs="?",
                        help="Path to config JSON file")
    parser.add_argument("--create-config", action="store_true",
                        help="Create a template config.json and exit")
    args = parser.parse_args()

    if args.create_config:
        create_template_config("config.json")
        return

    if not args.config:
        parser.print_help()
        sys.exit(1)

    # ---- Load config ----
    config_path = Path(args.config)
    config = load_config(config_path)
    config_dir = config_path.parent

    params = config.get("parameters", {})
    search_range = params.get("search_range_px", 8)
    ribbon_height = params.get("ribbon_height_px", 70)

    # ---- Load images ----
    print("Loading images...")
    images = {}
    for img_name in config["images"]:
        img_path = config_dir / img_name
        if not img_path.exists():
            raise FileNotFoundError(f"Image not found: {img_path}")
        img = tifffile.imread(str(img_path))
        # Convert to 2-D grayscale if needed
        if img.ndim == 3:
            img = np.mean(img[:, :, :3], axis=2).astype(img.dtype)
        images[img_name] = img
        print(f"  {img_name}: {img.shape[1]}x{img.shape[0]}, dtype={img.dtype}")

    # ---- Pair-wise alignment ----
    print("\n--- Pair-wise alignment ---")
    optimized = {}
    for pair in config["pairs"]:
        name_a = pair["image_a"]
        name_b = pair["image_b"]
        pt_a = np.array(pair["point_a"], dtype=np.float64)  # [x, y]
        pt_b = np.array(pair["point_b"], dtype=np.float64)  # [x, y]

        approx_offset = (pt_a[0] - pt_b[0], pt_a[1] - pt_b[1])

        print(f"\n[{name_a}] <-> [{name_b}]")
        print(f"  Approximate offset (x, y): ({approx_offset[0]:.0f}, {approx_offset[1]:.0f})")

        result = find_optimal_offset(
            images[name_a], images[name_b],
            approx_offset, search_range, ribbon_height)

        optimized[(name_a, name_b)] = result

    # ---- Global positions ----
    print("\n--- Computing global positions ---")
    positions = compute_global_positions(config["images"], optimized)
    for name in config["images"]:
        gx, gy = positions[name]
        print(f"  {name}: x={gx:.2f}, y={gy:.2f}")

    # ---- Stitch ----
    print("\n--- Assembling stitched image ---")
    result = create_stitched_image(
        images, positions, config["images"], ribbon_height)

    # ---- Save ----
    folder_name = config_dir.name
    output_path = config_dir / f"{folder_name}.tif"
    tifffile.imwrite(str(output_path), result)
    print(f"\nSaved: {output_path}")
    print(f"  Size: {result.shape[1]}x{result.shape[0]}, dtype={result.dtype}")
    print("Done.")


if __name__ == "__main__":
    main()
