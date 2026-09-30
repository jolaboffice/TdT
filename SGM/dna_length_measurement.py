"""SEM DNA path tracing and length measurement.

Traces the DNA backbone along a least-cost path (within a corridor) between the
user-clicked points, measures per-segment length in nm/bp (applying the DNA stretch),
and saves .csv / .maps / _points.json. DNALengthMeasurer=headless, InteractiveDNAMeasurer=matplotlib GUI.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons, Slider
from scipy import ndimage
from skimage import filters, io
from skimage.graph import route_through_array
import json
import os
from pathlib import Path


class DNALengthMeasurer:
    """Base class for DNA length measurement"""

    def __init__(self, image_path=None, image_array=None):
        if image_path is not None:
            self.image_path = image_path
            img = io.imread(image_path)
            if len(img.shape) == 3:
                self.original_image = np.mean(img, axis=2).astype(np.float64)
            else:
                self.original_image = img.astype(np.float64)
            print(f"Image loaded: {image_path}")
            print(f"Image size: {self.original_image.shape}")
        elif image_array is not None:
            self.image_path = None
            self.original_image = image_array.astype(np.float64)
        else:
            raise ValueError("provide either image_path or image_array")

        self.image = self.original_image.copy()
        self.mask = None
        self.points = []
        self.path_coords = []
        self.segment_lengths = []
        self.length_nm = 0
        self.length_bp = 0
        self.pixel_to_nm = 1.0
        self.nm_per_bp = 0.338  # B-form DNA: 0.338 nm/bp
        self.dna_stretch = 1.0  # DNA stretch default 100%
        self.corridor_width = 70  # corridor width (pixels)

    def set_scale(self, pixel_to_nm):
        """Set the pixel-to-nm scale"""
        self.pixel_to_nm = pixel_to_nm
        print(f"Scale set: 1 pixel = {pixel_to_nm} nm")

    def set_dna_stretch(self, stretch_percentage):
        """
        Set the DNA stretch

        Parameters:
        -----------
        stretch_percentage : float
            DNA stretch ratio (0-100, e.g. 100 = 100% stretched)
        """
        self.dna_stretch = stretch_percentage / 100.0
        effective_nm_per_bp = self.nm_per_bp * self.dna_stretch
        print(f"DNA stretch set: {stretch_percentage}%")
        print(f"Effective length: {effective_nm_per_bp:.4f} nm/bp")

    def set_corridor_width(self, pixels):
        """
        Set the path-search corridor width

        Parameters:
        -----------
        pixels : int
            max distance from the straight line between two points (pixels)
        """
        self.corridor_width = pixels
        print(f"Corridor width set: {pixels} pixels")

    def nm_to_bp(self, length_nm):
        """Convert nm to bp"""
        effective_nm_per_bp = self.nm_per_bp * self.dna_stretch
        return length_nm / effective_nm_per_bp

    def create_extreme_intensity_mask(self, white_threshold=0.95, black_threshold=0.05):
        """Build the extreme-intensity region mask"""
        normalized = (self.original_image - self.original_image.min()) / \
                     (self.original_image.max() - self.original_image.min())

        white_regions = normalized > white_threshold
        black_regions = normalized < black_threshold
        extreme_regions = white_regions | black_regions

        self.mask = ~extreme_regions

        white_count = np.sum(white_regions)
        black_count = np.sum(black_regions)
        total_pixels = extreme_regions.size
        valid_percentage = (np.sum(self.mask) / total_pixels) * 100

        print("\nExtreme-intensity region detection:")
        print(f"  - white region: {white_count:,} pixels ({white_count/total_pixels*100:.2f}%)")
        print(f"  - black region: {black_count:,} pixels ({black_count/total_pixels*100:.2f}%)")
        print(f"  - valid region: {np.sum(self.mask):,} pixels ({valid_percentage:.2f}%)")

        return self.mask

    def preprocess_image(self, sigma=1.5):
        """Image preprocessing"""
        self.image = filters.gaussian(self.original_image, sigma=sigma, preserve_range=True)
        return self.image

    def create_cost_map(self, dark_bias=2.0, gradient_weight=0.5, mask_penalty=1000.0):
        """Build the cost map"""
        normalized = (self.image - self.image.min()) / (self.image.max() - self.image.min())
        intensity_cost = normalized ** dark_bias

        gradient_y = ndimage.sobel(self.image, axis=0)
        gradient_x = ndimage.sobel(self.image, axis=1)
        gradient = np.sqrt(gradient_x**2 + gradient_y**2)
        gradient_norm = (gradient - gradient.min()) / (gradient.max() - gradient.min() + 1e-10)

        cost_map = intensity_cost * (1 - gradient_weight) + gradient_norm * gradient_weight

        if self.mask is not None:
            cost_map[~self.mask] = mask_penalty

        cost_map = cost_map + 1e-6
        self.cost_map = cost_map
        return cost_map

    def find_path_between_points(self, start_point, end_point):
        """Find the minimum-cost path between two points (with corridor constraint)"""
        cost_map_local = self.cost_map.copy()

        # corridor mask: only within corridor_width of the line between the two points
        h, w = cost_map_local.shape
        y1, x1 = start_point
        y2, x2 = end_point

        # all pixel coordinates
        yy, xx = np.ogrid[:h, :w]

        # direction/length of the line through the two points
        dx = x2 - x1
        dy = y2 - y1
        line_length = np.sqrt(dx**2 + dy**2)

        if line_length > 0:
            # point -> line perpendicular distance
            distance = np.abs((yy - y1) * dx - (xx - x1) * dy) / line_length

            # segment projection parameter t (bounds the region beyond the endpoints)
            t = ((yy - y1) * dy + (xx - x1) * dx) / (line_length ** 2)

            # slightly extend the segment range (margin around start/end)
            margin = self.corridor_width / line_length
            outside_segment = (t < -margin) | (t > 1 + margin)

            # high penalty outside the corridor
            outside_corridor = (distance > self.corridor_width) | outside_segment
            cost_map_local[outside_corridor] = 1000.0

        indices, weight = route_through_array(
            cost_map_local, start_point, end_point, fully_connected=True
        )
        return indices

    def calculate_segment_length(self, path_coords):
        """Compute the length of a path segment"""
        length_pixels = 0
        for i in range(len(path_coords) - 1):
            dy = path_coords[i+1][0] - path_coords[i][0]
            dx = path_coords[i+1][1] - path_coords[i][1]
            length_pixels += np.sqrt(dx**2 + dy**2)

        length_nm = length_pixels * self.pixel_to_nm
        length_bp = self.nm_to_bp(length_nm)
        return length_pixels, length_nm, length_bp

    def trace_dna_path(self, user_points, point_types=None, point_segments=None,
                      dark_bias=2.0, gradient_weight=0.5, sigma=1.5,
                      white_threshold=0.95, black_threshold=0.05, auto_mask=True):
        """Trace the DNA path along the user-specified points

        Parameters:
        -----------
        user_points : list of (x, y)
            all points (measure + guide)
        point_types : list of str, optional
            type of each point ('measure' or 'guide').
            if None, all are treated as 'measure'.
        point_segments : list of (int or None), optional
            segment number each point belongs to (1-indexed).
            None for measure, segment number for guide.
            if None, use the default automatic assignment.
        """
        print(f"\n{'='*70}")
        print("DNA path tracing started")
        print(f"{'='*70}")

        if point_types is None:
            point_types = ['measure'] * len(user_points)

        # extract measure points
        measure_pts = [(pt, i) for i, (pt, t) in enumerate(zip(user_points, point_types)) if t == 'measure']
        num_segs = max(1, len(measure_pts) - 1)
        seg_guides = [[] for _ in range(num_segs)]

        # assign guide points to segments
        # user-specified mode only if point_segments exists and any guide has a segment set
        has_assigned_segments = (point_segments is not None and
                                 any(seg is not None for seg, ptype in zip(point_segments, point_types) if ptype == 'guide'))
        if has_assigned_segments:
            # use user-specified segments
            for i, (pt, ptype, seg) in enumerate(zip(user_points, point_types, point_segments)):
                if ptype == 'guide' and seg is not None:
                    seg_idx = seg - 1  # 1-indexed -> 0-indexed
                    if 0 <= seg_idx < num_segs:
                        seg_guides[seg_idx].append((pt, i))  # (coord, input order)
        else:
            # default automatic assignment (nearest segment)
            guide_pts = [pt for pt, t in zip(user_points, point_types) if t == 'guide']
            for gx, gy in guide_pts:
                best_seg = 0
                best_dist = float('inf')
                for s in range(num_segs):
                    mx1, my1 = measure_pts[s][0]
                    mx2, my2 = measure_pts[s + 1][0]
                    dx, dy = mx2 - mx1, my2 - my1
                    seg_len_sq = dx**2 + dy**2
                    if seg_len_sq > 0:
                        t = ((gx - mx1) * dx + (gy - my1) * dy) / seg_len_sq
                        t_c = max(0, min(1, t))
                        cx, cy = mx1 + t_c * dx, my1 + t_c * dy
                    else:
                        cx, cy = mx1, my1
                    dist = np.sqrt((gx - cx)**2 + (gy - cy)**2)
                    if dist < best_dist:
                        best_dist = dist
                        best_seg = s
                seg_guides[best_seg].append(((gx, gy), 0))  # 0 instead of input order

        # sort each segment's guides by input order (when using point_segments)
        # or by t value (automatic assignment)
        if has_assigned_segments:
            for s in range(num_segs):
                seg_guides[s].sort(key=lambda x: x[1])  # sort by input order
        else:
            # compute t then sort
            for s in range(num_segs):
                mx1, my1 = measure_pts[s][0]
                mx2, my2 = measure_pts[s + 1][0]
                dx, dy = mx2 - mx1, my2 - my1
                seg_len_sq = dx**2 + dy**2
                new_list = []
                for (gpt, _) in seg_guides[s]:
                    gx, gy = gpt
                    t_val = ((gx - mx1) * dx + (gy - my1) * dy) / seg_len_sq if seg_len_sq > 0 else 0
                    new_list.append((gpt, t_val))
                new_list.sort(key=lambda x: x[1])
                seg_guides[s] = new_list

        # build the sorted point list
        ordered_points = []
        ordered_types = []
        for s in range(num_segs):
            ordered_points.append(measure_pts[s][0])
            ordered_types.append('measure')
            for (gpt, _) in seg_guides[s]:
                ordered_points.append(gpt)
                ordered_types.append('guide')
        # add the last measure point
        ordered_points.append(measure_pts[-1][0])
        ordered_types.append('measure')

        user_points = ordered_points
        point_types = ordered_types

        self.points = user_points
        self.point_types = point_types

        measure_count = sum(1 for t in point_types if t == 'measure')
        guide_count = sum(1 for t in point_types if t == 'guide')
        print(f"Points: {measure_count} measure, {guide_count} guide")
        print(f"DNA stretch: {self.dna_stretch*100:.1f}%")

        # preprocessing
        self.preprocess_image(sigma=sigma)

        # extreme-intensity region mask
        if auto_mask:
            self.create_extreme_intensity_mask(
                white_threshold=white_threshold,
                black_threshold=black_threshold
            )

        # build the cost map
        self.create_cost_map(dark_bias=dark_bias, gradient_weight=gradient_weight)

        # trace the path between adjacent points (through all points)
        print("\nSearching path...")
        full_path = []
        sub_paths = []  # path for each adjacent pair

        for i in range(len(user_points) - 1):
            start = (user_points[i][1], user_points[i][0])
            end = (user_points[i+1][1], user_points[i+1][0])

            sub_path = self.find_path_between_points(start, end)
            sub_paths.append(sub_path)

            if i == 0:
                full_path.extend(sub_path)
            else:
                full_path.extend(sub_path[1:])

        self.path_coords = full_path

        # extract measure-point indices
        measure_indices = [i for i, t in enumerate(point_types) if t == 'measure']

        # display path: only the first-measure to last-measure span
        first_m = measure_indices[0]
        last_m = measure_indices[-1]
        display_path = []
        for j in range(first_m, last_m):
            if j == first_m:
                display_path.extend(sub_paths[j])
            else:
                display_path.extend(sub_paths[j][1:])
        self.display_path = display_path

        # compute segment lengths between measure points (summing guide sub-paths)
        self.segment_lengths = []
        seg_id = 0
        for m in range(len(measure_indices) - 1):
            seg_id += 1
            start_idx = measure_indices[m]
            end_idx = measure_indices[m + 1]

            # sum the sub_paths between start_idx and end_idx
            merged_path = []
            for j in range(start_idx, end_idx):
                if j == start_idx:
                    merged_path.extend(sub_paths[j])
                else:
                    merged_path.extend(sub_paths[j][1:])

            seg_pixels, seg_nm, seg_bp = self.calculate_segment_length(merged_path)
            self.segment_lengths.append({
                'segment_id': seg_id,
                'start_point': user_points[start_idx],
                'end_point': user_points[end_idx],
                'length_pixels': seg_pixels,
                'length_nm': seg_nm,
                'length_bp': seg_bp
            })

            guide_via = end_idx - start_idx - 1
            via_str = f" (via {guide_via} guide)" if guide_via > 0 else ""
            print(f"  Segment {seg_id}: {seg_pixels:.1f} px = {seg_nm:.2f} nm = {seg_bp:.0f} bp{via_str}")

        # compute total length
        total_pixels, total_nm, total_bp = self.calculate_segment_length(full_path)
        self.length_nm = total_nm
        self.length_bp = total_bp

        print(f"\n{'='*70}")
        print("Tracing complete!")
        print(f"{'='*70}")
        print(f"Total length: {total_nm:.2f} nm = {total_bp:.0f} bp")
        print(f"Segment count: {len(self.segment_lengths)}")
        print(f"{'='*70}\n")

        return full_path, total_nm, total_bp

    def visualize_result(self, save_path=None, figsize=(24, 20), show_segment_labels=True):
        """Visualize the result (vertical layout)"""
        fig, axes = plt.subplots(3, 1, figsize=figsize)

        # original image
        axes[0].imshow(self.original_image, cmap='gray')
        axes[0].set_title('Original Image', fontsize=12)
        axes[0].axis('off')

        # Cost map
        im1 = axes[1].imshow(self.cost_map, cmap='hot', vmin=0, vmax=10)
        axes[1].set_title('Cost Map', fontsize=12)
        axes[1].axis('off')
        plt.colorbar(im1, ax=axes[1], fraction=0.046, orientation='horizontal', pad=0.02)

        # result image
        axes[2].imshow(self.original_image, cmap='gray')

        # draw the path (only the first-measure to last-measure span)
        draw_path = getattr(self, 'display_path', self.path_coords)
        if len(draw_path) > 0:
            path_array = np.array(draw_path)
            axes[2].plot(path_array[:, 1], path_array[:, 0], 'r-',
                        linewidth=2, label='Traced Path', alpha=0.8)

        # show measure points only
        if len(self.points) > 0:
            point_types = getattr(self, 'point_types', ['measure'] * len(self.points))
            m_count = 0
            for i, (x, y) in enumerate(self.points):
                ptype = point_types[i] if i < len(point_types) else 'measure'
                if ptype != 'measure':
                    continue
                m_count += 1
                axes[2].plot(x, y, 'go', markersize=10,
                            markeredgecolor='yellow', markeredgewidth=2)
                axes[2].text(x, y, str(m_count), color='yellow', fontsize=10,
                           ha='center', va='center', weight='bold')

        axes[2].set_title(f'Result: {self.length_nm:.2f} nm ({self.length_bp:.0f} bp)',
                         fontsize=12)
        axes[2].legend(loc='upper right')
        axes[2].axis('off')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Result image saved: {save_path}")

        plt.show()

        self.print_segment_summary()

        return fig

    def print_segment_summary(self):
        """Print the per-segment length summary"""
        if len(self.segment_lengths) == 0:
            return

        print("\n" + "="*85)
        print("Per-segment length summary")
        print("="*85)
        print(f"{'Seg':<5} {'From':<10} {'To':<10} {'Length (nm)':<15} {'Length (bp)':<15} {'Length (px)':<12}")
        print("-" * 85)

        for seg in self.segment_lengths:
            seg_id = seg['segment_id']
            from_pt = f"Point {seg_id}"
            to_pt = f"Point {seg_id + 1}"

            print(f"{seg_id:<5} {from_pt:<10} {to_pt:<10} {seg['length_nm']:<15.2f} "
                  f"{seg['length_bp']:<15.0f} {seg['length_pixels']:<12.2f}")

        print("-" * 85)
        print(f"{'TOTAL':<5} {'':<10} {'':<10} {self.length_nm:<15.2f} "
              f"{self.length_bp:<15.0f} {'':<12}")
        print("="*85 + "\n")

    def save_segment_data(self, csv_path=None):
        """Save segment data to a CSV file"""
        import csv

        # build CSV filename (saved in the image folder)
        if csv_path is None:
            if self.image_path:
                img_dir = Path(self.image_path).parent
                base_name = Path(self.image_path).stem
                csv_path = str(img_dir / f"{base_name}.csv")
            else:
                csv_path = "dna_segment_lengths.csv"

        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)

            # header
            writer.writerow(['# DNA Length Measurement Results'])
            writer.writerow(['# Image:', self.image_path if self.image_path else 'N/A'])
            writer.writerow(['# Scale:', f'{self.pixel_to_nm} nm/pixel'])
            writer.writerow(['# DNA Stretch:', f'{self.dna_stretch*100:.1f}%'])
            writer.writerow(['# Effective:', f'{self.nm_per_bp * self.dna_stretch:.4f} nm/bp'])
            writer.writerow([])

            # data
            writer.writerow(['Segment_ID', 'From_Point', 'To_Point', 
                           'Length_nm', 'Length_bp', 'Length_pixels'])

            for seg in self.segment_lengths:
                writer.writerow([
                    seg['segment_id'],
                    f"Point {seg['segment_id']}",
                    f"Point {seg['segment_id'] + 1}",
                    f"{seg['length_nm']:.2f}",
                    f"{seg['length_bp']:.0f}",
                    f"{seg['length_pixels']:.2f}"
                ])

            # total
            total_px = sum(s['length_pixels'] for s in self.segment_lengths)
            writer.writerow(['TOTAL', '', '', 
                           f"{self.length_nm:.2f}", 
                           f"{self.length_bp:.0f}", 
                           f"{total_px:.2f}"])

        print(f"Segment data saved: {csv_path}")
        return csv_path

    def save_maps_data(self, maps_path=None):
        """Save to a .maps file (bp, rounded integers)"""
        if maps_path is None:
            if self.image_path:
                img_dir = Path(self.image_path).parent
                base_name = Path(self.image_path).stem
                maps_path = str(img_dir / f"{base_name}.maps")
            else:
                maps_path = "dna_length.maps"

        # filename (without extension)
        file_name = Path(self.image_path).stem if self.image_path else "unknown"

        # each fragment length (bp, rounded)
        fragment_lengths = [round(seg['length_bp']) for seg in self.segment_lengths]

        # total length (sum of fragments to avoid rounding error)
        total_bp = sum(fragment_lengths)

        # fragment count
        num_fragments = len(self.segment_lengths)

        # save tab-separated
        with open(maps_path, 'w', encoding='utf-8') as f:
            row = [file_name, str(total_bp), str(num_fragments)]
            row.extend([str(length) for length in fragment_lengths])
            f.write('\t'.join(row) + '\n')

        print(f"Maps file saved: {maps_path}")
        return maps_path

    def save_image_with_points(self, save_path=None, save_adjusted=True):
        """Save at original quality with only the points and path overlaid on the original image

        Parameters:
        -----------
        save_path : str, optional
            save path; if None, auto-generated in the image folder
        save_adjusted : bool
            if True, also save contrast/brightness-adjusted versions (with and without trace)
        """
        from PIL import Image, ImageDraw, ImageFont

        # determine save path (in the image's folder)
        if save_path is None:
            if self.image_path:
                img_dir = Path(self.image_path).parent
                base_name = Path(self.image_path).stem
                save_path = str(img_dir / f"{base_name}_result.tif")
            else:
                save_path = "dna_result.tif"

        save_dir = Path(save_path).parent
        base_name = Path(save_path).stem.replace('_result', '')

        # normalize the original image to uint8
        img_normalized = (self.original_image - self.original_image.min()) / \
                         (self.original_image.max() - self.original_image.min())
        img_uint8 = (img_normalized * 255).astype(np.uint8)

        # convert grayscale to RGB
        img_rgb = Image.fromarray(img_uint8).convert('RGB')
        draw = ImageDraw.Draw(img_rgb)

        # marker size (proportional to image size)
        img_min_dim = min(self.original_image.shape)
        marker_radius = max(8, img_min_dim // 150)
        line_width = max(2, img_min_dim // 500)

        # draw the path (red, only the first-measure to last-measure span)
        draw_path = getattr(self, 'display_path', self.path_coords)
        if len(draw_path) > 1:
            path_points = [(int(coord[1]), int(coord[0])) for coord in draw_path]
            draw.line(path_points, fill='red', width=line_width)

        # draw measure points only
        point_types = getattr(self, 'point_types', ['measure'] * len(self.points))
        m_count = 0
        font_size = max(12, marker_radius * 2)
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", font_size)
        except:
            try:
                font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
            except:
                font = ImageFont.load_default()

        for i, (x, y) in enumerate(self.points):
            ptype = point_types[i] if i < len(point_types) else 'measure'
            if ptype != 'measure':
                continue
            m_count += 1

            bbox = [x - marker_radius, y - marker_radius,
                    x + marker_radius, y + marker_radius]
            draw.ellipse(bbox, fill='green', outline='yellow', width=max(2, line_width))

            label = str(m_count)
            text_bbox = draw.textbbox((x, y), label, font=font)
            text_w = text_bbox[2] - text_bbox[0]
            text_h = text_bbox[3] - text_bbox[1]
            draw.text((x - text_w//2, y - text_h//2), label, fill='yellow', font=font)

        img_rgb.save(save_path)
        print(f"Result image saved (original quality): {save_path}")

        # save contrast/brightness-adjusted versions
        if save_adjusted and hasattr(self, 'contrast') and hasattr(self, 'brightness'):
            # create the adjusted image
            adjusted = self.apply_contrast_brightness()
            adjusted_uint8 = (adjusted * 255).astype(np.uint8)

            # 1. adjusted image only (no trace)
            adjusted_only_path = str(save_dir / f"{base_name}_adjusted.tif")
            Image.fromarray(adjusted_uint8).save(adjusted_only_path)
            print(f"Adjusted image saved (no trace): {adjusted_only_path}")

            # 2. adjusted image + trace
            adjusted_rgb = Image.fromarray(adjusted_uint8).convert('RGB')
            draw_adj = ImageDraw.Draw(adjusted_rgb)

            # draw the path
            if len(draw_path) > 1:
                draw_adj.line(path_points, fill='red', width=line_width)

            # draw the points
            m_count = 0
            for i, (x, y) in enumerate(self.points):
                ptype = point_types[i] if i < len(point_types) else 'measure'
                if ptype != 'measure':
                    continue
                m_count += 1
                bbox = [x - marker_radius, y - marker_radius,
                        x + marker_radius, y + marker_radius]
                draw_adj.ellipse(bbox, fill='green', outline='yellow', width=max(2, line_width))
                label = str(m_count)
                text_bbox = draw_adj.textbbox((x, y), label, font=font)
                text_w = text_bbox[2] - text_bbox[0]
                text_h = text_bbox[3] - text_bbox[1]
                draw_adj.text((x - text_w//2, y - text_h//2), label, fill='yellow', font=font)

            adjusted_traced_path = str(save_dir / f"{base_name}_adjusted_traced.tif")
            adjusted_rgb.save(adjusted_traced_path)
            print(f"Adjusted image saved (with trace): {adjusted_traced_path}")

        return save_path

    def save_points(self, json_path=None):
        """Save points"""
        if json_path is None:
            if self.image_path:
                img_dir = Path(self.image_path).parent
                base_name = Path(self.image_path).stem
                json_path = str(img_dir / f"{base_name}_points.json")
            else:
                json_path = "dna_points.json"

        # save the original click order (clicked_*) instead of the reordered self.*
        points_to_save = getattr(self, 'clicked_points', self.points)
        types_to_save = getattr(self, 'clicked_point_types', getattr(self, 'point_types', ['measure'] * len(points_to_save)))
        segments_to_save = getattr(self, 'clicked_point_segments', getattr(self, 'point_segments', None))

        data = {
            'image_path': self.image_path,
            'image_shape': list(self.original_image.shape),
            'points': points_to_save,
            'point_types': types_to_save,
            'point_segments': segments_to_save,
            'pixel_to_nm': self.pixel_to_nm,
            'dna_stretch': self.dna_stretch * 100,  # saved as a percentage
            'contrast': getattr(self, 'contrast', 1.0),
            'brightness': getattr(self, 'brightness', 0.0)
        }

        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)

        print(f"Points saved: {json_path}")
        return json_path

    def load_points(self, json_path):
        """Load points"""
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # check image size
        if list(self.original_image.shape) != data['image_shape']:
            print("Warning: image size differs!")
            print(f"  current: {self.original_image.shape}")
            print(f"  saved info: {data['image_shape']}")
            response = input("Continue? (y/n): ")
            if response.lower() != 'y':
                return None

        self.points = [tuple(p) for p in data['points']]
        self.point_types = data.get('point_types', ['measure'] * len(self.points))

        # validate point_segments length (recover corrupted files)
        loaded_segments = data.get('point_segments')
        if loaded_segments is None or len(loaded_segments) != len(self.points):
            if loaded_segments is not None:
                print(f"  Warning: point_segments length mismatch ({len(loaded_segments)} vs {len(self.points)}), reset")
            self.point_segments = [None] * len(self.points)
        else:
            self.point_segments = loaded_segments
        self.pixel_to_nm = data.get('pixel_to_nm', 1.0)
        self.dna_stretch = data.get('dna_stretch', 100.0) / 100.0

        # load contrast/brightness (for InteractiveDNAMeasurer)
        loaded_contrast = data.get('contrast', 1.0)
        loaded_brightness = data.get('brightness', 0.0)
        if hasattr(self, 'contrast'):
            self.contrast = loaded_contrast
        if hasattr(self, 'brightness'):
            self.brightness = loaded_brightness

        measure_count = sum(1 for t in self.point_types if t == 'measure')
        guide_count = sum(1 for t in self.point_types if t == 'guide')
        print(f"Points loaded: {json_path}")
        print(f"  - Measure: {measure_count}, Guide: {guide_count}")
        print(f"  - scale: {self.pixel_to_nm} nm/pixel")
        print(f"  - DNA stretch: {self.dna_stretch*100:.1f}%")
        print(f"  - contrast: {loaded_contrast:.1f}, brightness: {loaded_brightness:.2f}")

        return self.points


# ============================================================================
# Enhanced Interactive GUI
# ============================================================================

class InteractiveDNAMeasurer(DNALengthMeasurer):
    """Interactive GUI adding zoom/pan and point save/load"""

    def __init__(self, image_path):
        super().__init__(image_path)
        self.fig = None
        self.ax = None
        self.img_display = None  # image display object
        self.clicked_points = []
        self.clicked_point_types = []  # 'measure' or 'guide'
        self.clicked_point_segments = []  # segment number for guides (None for measure)
        self.point_plots = []
        self.text_plots = []
        self.mode = 'navigate'
        self.home_view = None
        self.selected_point_idx = None  # index of the point selected in Move mode
        self.insert_after_idx = None  # reference-point index in Insert After mode
        self.contrast = 1.0  # contrast (0.1 ~ 3.0)
        self.brightness = 0.0  # brightness (-1.0 ~ 1.0)
        self.pending_guide_point = None  # guide coordinate awaiting segment selection
        self.segment_buttons = []  # segment-selection buttons
        self.segment_button_axes = []  # button axes

    def apply_contrast_brightness(self):
        """Return the image with contrast and brightness applied"""
        # normalize to 0~1
        img_norm = (self.original_image - self.original_image.min()) / \
                   (self.original_image.max() - self.original_image.min() + 1e-10)
        # apply contrast (around the median)
        adjusted = (img_norm - 0.5) * self.contrast + 0.5
        # apply brightness
        adjusted = adjusted + self.brightness
        # clipping
        adjusted = np.clip(adjusted, 0, 1)
        return adjusted

    def on_contrast_change(self, val):
        """Contrast slider callback"""
        self.contrast = val
        adjusted = self.apply_contrast_brightness()
        self.img_display.set_data(adjusted)
        self.fig.canvas.draw_idle()

    def on_brightness_change(self, val):
        """Brightness slider callback"""
        self.brightness = val
        adjusted = self.apply_contrast_brightness()
        self.img_display.set_data(adjusted)
        self.fig.canvas.draw_idle()

    def _get_point_label(self, idx):
        """Build the display label for a point index (per-type ordinal)"""
        ptype = self.clicked_point_types[idx]
        if ptype == 'measure':
            m_count = sum(1 for t in self.clicked_point_types[:idx+1] if t == 'measure')
            return str(m_count)
        else:
            g_count = sum(1 for t in self.clicked_point_types[:idx+1] if t == 'guide')
            return str(g_count)

    def _add_point_plot(self, x, y, idx):
        """Add a point to the plot"""
        ptype = self.clicked_point_types[idx]
        label = self._get_point_label(idx)

        if ptype == 'measure':
            color, edge_color = 'go', 'yellow'
            text_color = 'yellow'
        else:
            color, edge_color = 'bo', 'cyan'
            text_color = 'cyan'

        point_plot, = self.ax.plot(x, y, color, markersize=10,
                    markeredgecolor=edge_color, markeredgewidth=2)
        text_plot = self.ax.text(x, y, label,
                    color=text_color, fontsize=10,
                    ha='center', va='center', weight='bold')

        self.point_plots.append(point_plot)
        self.text_plots.append(text_plot)

    def onclick(self, event):
        """Mouse click event"""
        if event.inaxes != self.ax:
            return

        if self.mode == 'add' and event.button == 1:
            x, y = int(event.xdata), int(event.ydata)
            self.clicked_points.append((x, y))
            self.clicked_point_types.append('measure')
            self.clicked_point_segments.append(None)

            idx = len(self.clicked_points) - 1
            self._add_point_plot(x, y, idx)
            self.fig.canvas.draw_idle()

            label = self._get_point_label(idx)
            print(f"Add Measure {label}: ({x}, {y})")

        elif self.mode == 'add_guide' and event.button == 1:
            x, y = int(event.xdata), int(event.ydata)

            # check measure count
            measure_count = sum(1 for t in self.clicked_point_types if t == 'measure')
            if measure_count < 2:
                print("Warning: at least 2 measure points are required!")
                return

            num_segs = measure_count - 1
            if num_segs == 1:
                # if only one segment, add directly
                self._add_guide_to_segment(x, y, 1)
            else:
                # if multiple segments, choose via buttons
                self.pending_guide_point = (x, y)
                self._show_segment_buttons(num_segs)
                print(f"Select a segment (1-{num_segs})")

        elif self.mode == 'delete' and event.button == 1:
            if len(self.clicked_points) == 0:
                return

            x, y = event.xdata, event.ydata
            distances = [np.sqrt((px - x)**2 + (py - y)**2)
                        for px, py in self.clicked_points]
            nearest_pt_idx = int(np.argmin(distances))

            # if a measure, compute its index before deletion
            removed_type = self.clicked_point_types[nearest_pt_idx]
            removed_seg = self.clicked_point_segments[nearest_pt_idx]

            if removed_type == 'measure':
                # measure index of the measure to delete (0-indexed)
                m_idx = sum(1 for t in self.clicked_point_types[:nearest_pt_idx] if t == 'measure')

                # adjust guide segment numbers
                for i in range(len(self.clicked_point_segments)):
                    if self.clicked_point_types[i] == 'guide':
                        seg = self.clicked_point_segments[i]
                        if seg is not None:
                            # segments involving m_idx: m_idx (before), m_idx+1 (after)
                            # 1-indexed: segments m_idx and m_idx+1
                            if seg == m_idx + 1:
                                # merge the deleted measure's following segment into the preceding one
                                self.clicked_point_segments[i] = max(1, m_idx)
                            elif seg > m_idx + 1:
                                # later segments shift down by 1
                                self.clicked_point_segments[i] = seg - 1

            # delete the point
            removed_point = self.clicked_points.pop(nearest_pt_idx)
            self.clicked_point_types.pop(nearest_pt_idx)
            self.clicked_point_segments.pop(nearest_pt_idx)
            self.point_plots[nearest_pt_idx].remove()
            self.text_plots[nearest_pt_idx].remove()
            self.point_plots.pop(nearest_pt_idx)
            self.text_plots.pop(nearest_pt_idx)

            for i in range(len(self.text_plots)):
                self.text_plots[i].set_text(self._get_point_label(i))
            self.fig.canvas.draw_idle()
            seg_info = f" (segment {removed_seg})" if removed_seg else ""
            print(f"{removed_type} point deleted{seg_info}: {removed_point}")

        elif self.mode == 'move' and event.button == 1:
            if len(self.clicked_points) == 0:
                return

            x, y = event.xdata, event.ydata

            if self.selected_point_idx is None:
                # first click: select the point to move
                distances = [np.sqrt((px - x)**2 + (py - y)**2)
                            for px, py in self.clicked_points]
                self.selected_point_idx = np.argmin(distances)

                # highlight the selected point
                self.point_plots[self.selected_point_idx].set_markeredgecolor('red')
                self.point_plots[self.selected_point_idx].set_markersize(14)
                self.fig.canvas.draw_idle()

                label = self._get_point_label(self.selected_point_idx)
                print(f"{label} selected - click the new position")
            else:
                # second click: move the selected point to the new position
                new_x, new_y = int(x), int(y)
                old_point = self.clicked_points[self.selected_point_idx]
                ptype = self.clicked_point_types[self.selected_point_idx]
                self.clicked_points[self.selected_point_idx] = (new_x, new_y)

                # update the plot
                edge_color = 'yellow' if ptype == 'measure' else 'cyan'
                self.point_plots[self.selected_point_idx].set_data([new_x], [new_y])
                self.point_plots[self.selected_point_idx].set_markeredgecolor(edge_color)
                self.point_plots[self.selected_point_idx].set_markersize(10)
                self.text_plots[self.selected_point_idx].set_position((new_x, new_y))

                self.fig.canvas.draw_idle()
                label = self._get_point_label(self.selected_point_idx)
                print(f"{label} moved: {old_point} -> ({new_x}, {new_y})")

                self.selected_point_idx = None

        elif self.mode == 'insert_after' and event.button == 1:
            x, y = event.xdata, event.ydata

            if self.insert_after_idx is None:
                # step 1: select a reference point (click an existing measure point)
                nearest_idx = self._find_nearest_measure_point(x, y)
                if nearest_idx is not None:
                    self.insert_after_idx = nearest_idx
                    # visual marker (green outline)
                    self.point_plots[nearest_idx].set_markeredgecolor('lime')
                    self.point_plots[nearest_idx].set_markeredgewidth(3)
                    self.fig.canvas.draw_idle()
                    label = self._get_point_label(nearest_idx)
                    print(f"Reference point selected: {label} - a new point will be inserted at the next click")
                else:
                    print("No measure point nearby")
            else:
                # step 2: insert the new point
                insert_idx = self.insert_after_idx + 1
                self._insert_measure_point(int(x), int(y), insert_idx)
                # clear the reference-point marker
                self._reset_insert_reference()

    def on_scroll(self, event):
        """Mouse-wheel zoom"""
        if event.inaxes != self.ax or self.mode != 'navigate':
            return

        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        xdata = event.xdata
        ydata = event.ydata

        if event.button == 'up':
            scale_factor = 0.8
        elif event.button == 'down':
            scale_factor = 1.25
        else:
            return

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

        relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])

        self.ax.set_xlim([xdata - new_width * (1 - relx), xdata + new_width * relx])
        self.ax.set_ylim([ydata - new_height * (1 - rely), ydata + new_height * rely])

        self.fig.canvas.draw_idle()

    def on_pan_press(self, event):
        if event.inaxes != self.ax or self.mode != 'navigate':
            return
        if event.button == 3:
            self.pan_start = (event.xdata, event.ydata)

    def on_pan_motion(self, event):
        if not hasattr(self, 'pan_start') or self.pan_start is None:
            return
        if event.inaxes != self.ax or self.mode != 'navigate':
            return

        dx = event.xdata - self.pan_start[0]
        dy = event.ydata - self.pan_start[1]

        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()

        self.ax.set_xlim([cur_xlim[0] - dx, cur_xlim[1] - dx])
        self.ax.set_ylim([cur_ylim[0] - dy, cur_ylim[1] - dy])

        self.fig.canvas.draw_idle()

    def on_pan_release(self, event):
        self.pan_start = None

    def _add_guide_to_segment(self, x, y, seg):
        """Add a guide point to the given segment"""
        self.clicked_points.append((x, y))
        self.clicked_point_types.append('guide')
        self.clicked_point_segments.append(seg)

        idx = len(self.clicked_points) - 1
        self._add_point_plot(x, y, idx)
        self.fig.canvas.draw_idle()

        label = self._get_point_label(idx)
        print(f"Add Guide {label} (segment {seg}): ({x}, {y})")
        self.pending_guide_point = None

    def _show_segment_buttons(self, num_segs):
        """Show segment-selection buttons (horizontal, above the slider and below the image)"""
        self._hide_segment_buttons()

        # button area: above the slider (y=0.04), below the image (y=0.14)
        btn_height = 0.015
        btn_y = 0.09  # between slider and image

        # compute button width (total width / number of segments)
        total_width = 0.70  # available width
        gap = 0.003
        btn_width = min(0.0125, (total_width - gap * (num_segs - 1)) / num_segs)

        # center the buttons
        total_btn_width = btn_width * num_segs + gap * (num_segs - 1)
        start_x = 0.05 + (0.75 - total_btn_width) / 2  # center of the image area

        for i in range(num_segs):
            seg_num = i + 1
            x = start_x + i * (btn_width + gap)

            ax_btn = plt.axes([x, btn_y, btn_width, btn_height])
            btn = Button(ax_btn, str(seg_num), color='lightblue', hovercolor='deepskyblue')

            # create the callback explicitly
            def make_callback(s):
                def callback(event):
                    self._on_segment_button_click(s)
                return callback
            btn.on_clicked(make_callback(seg_num))

            self.segment_button_axes.append(ax_btn)
            self.segment_buttons.append(btn)

        self.fig.canvas.draw_idle()

    def _hide_segment_buttons(self):
        """Hide the segment-selection buttons"""
        for ax_btn in self.segment_button_axes:
            ax_btn.remove()
        self.segment_button_axes = []
        self.segment_buttons = []

    def _on_segment_button_click(self, seg_num):
        """Segment-button click callback"""
        print(f"  segment {seg_num} button clicked")
        if self.pending_guide_point is not None:
            x, y = self.pending_guide_point
            self._add_guide_to_segment(x, y, seg_num)
            self._hide_segment_buttons()
            self.fig.canvas.draw_idle()
        else:
            print("  Warning: pending_guide_point is None")

    def _find_nearest_measure_point(self, x, y, threshold=20):
        """Return the index of the nearest measure point"""
        min_dist = float('inf')
        nearest_idx = None
        for i, (px, py) in enumerate(self.clicked_points):
            if self.clicked_point_types[i] == 'measure':
                dist = np.sqrt((x - px)**2 + (y - py)**2)
                if dist < min_dist and dist < threshold:
                    min_dist = dist
                    nearest_idx = i
        return nearest_idx

    def _insert_measure_point(self, x, y, insert_idx):
        """Insert a measure point at the given position"""
        # compute measure index before insertion (0-indexed)
        m_idx = sum(1 for t in self.clicked_point_types[:insert_idx] if t == 'measure')

        # adjust guide segment numbers: segment m_idx (1-indexed) is split
        for i in range(len(self.clicked_point_segments)):
            if self.clicked_point_types[i] == 'guide':
                seg = self.clicked_point_segments[i]
                if seg is not None:
                    if seg > m_idx:
                        # segments after the insertion -> +1
                        self.clicked_point_segments[i] = seg + 1
                    elif seg == m_idx and i >= insert_idx:
                        # guides in the same segment but after the insertion point -> new segment
                        self.clicked_point_segments[i] = seg + 1

        self.clicked_points.insert(insert_idx, (x, y))
        self.clicked_point_types.insert(insert_idx, 'measure')
        self.clicked_point_segments.insert(insert_idx, None)

        # redraw all points (refresh labels)
        self.redraw_points()
        print(f"New Measure point inserted: position M{self._get_measure_index(insert_idx) + 1}")

    def _get_measure_index(self, idx):
        """Compute the measure-point index from the overall index"""
        count = 0
        for i in range(idx):
            if self.clicked_point_types[i] == 'measure':
                count += 1
        return count

    def _reset_insert_reference(self):
        """Clear the reference-point selection in Insert After mode"""
        if self.insert_after_idx is not None:
            idx = self.insert_after_idx
            if idx < len(self.point_plots):
                self.point_plots[idx].set_markeredgecolor('yellow')
                self.point_plots[idx].set_markeredgewidth(2)
            self.insert_after_idx = None
            self.fig.canvas.draw_idle()

    def on_mode_change(self, label):
        mode_map = {
            'Navigate': 'navigate',
            'Add Measure': 'add',
            'Insert After': 'insert_after',
            'Add Guide': 'add_guide',
            'Delete': 'delete',
            'Move Point': 'move'
        }
        self.mode = mode_map[label]

        # hide segment buttons on mode change
        self._hide_segment_buttons()
        self.pending_guide_point = None  # reset on mode change

        # reset the selected point on mode change
        if self.selected_point_idx is not None:
            ptype = self.clicked_point_types[self.selected_point_idx]
            edge_color = 'yellow' if ptype == 'measure' else 'cyan'
            self.point_plots[self.selected_point_idx].set_markeredgecolor(edge_color)
            self.point_plots[self.selected_point_idx].set_markersize(10)
            self.selected_point_idx = None

        # reset the Insert After reference point
        if self.mode != 'insert_after':
            self._reset_insert_reference()

        titles = {
            'navigate': ('Navigation Mode: Wheel=Zoom, Right-drag=Pan', 'blue'),
            'add': ('Add Measure Point: Left-click to add (green)', 'green'),
            'insert_after': ('Insert After: Click measure point, then click to insert', 'lime'),
            'add_guide': ('Add Guide Point: Left-click to add (blue)', 'deepskyblue'),
            'delete': ('Delete Mode: Left-click nearest point', 'red'),
            'move': ('Move Point Mode: Click point, then click new position', 'orange'),
        }
        title, color = titles[self.mode]
        self.ax.set_title(title, fontsize=14, color=color)

        self.fig.canvas.draw_idle()

    def redraw_points(self):
        """Redraw the points"""
        # remove existing plots
        for plot in self.point_plots:
            plot.remove()
        for text in self.text_plots:
            text.remove()
        self.point_plots = []
        self.text_plots = []

        # redraw
        for i, (x, y) in enumerate(self.clicked_points):
            self._add_point_plot(x, y, i)

        self.fig.canvas.draw_idle()

    def preprocess_image(self, sigma=1.5):
        """Image preprocessing - apply contrast/brightness"""
        # use the contrast/brightness-adjusted image
        adjusted = self.apply_contrast_brightness()
        # convert to the 0~255 range
        adjusted = adjusted * 255.0
        self.image = filters.gaussian(adjusted, sigma=sigma, preserve_range=True)
        print(f"Applied contrast/brightness: contrast={self.contrast:.1f}, brightness={self.brightness:.2f}")
        return self.image

    def on_trace(self, event):
        """Trace button"""
        measure_count = sum(1 for t in self.clicked_point_types if t == 'measure')
        if measure_count < 2:
            print("\nWarning: please add at least 2 measure points!")
            return
        if len(self.clicked_points) < 2:
            print("\nWarning: please add at least 2 points!")
            return

        path, length_nm, length_bp = self.trace_dna_path(
            self.clicked_points,
            point_types=self.clicked_point_types,
            point_segments=self.clicked_point_segments,
            dark_bias=2.5,
            gradient_weight=0.3,
            sigma=1.5,
            white_threshold=0.95,
            black_threshold=0.05,
            auto_mask=True
        )

        # show the result on screen (no save)
        self.visualize_result(save_path=None, show_segment_labels=True)
        # save the point/path image at original quality (in the image folder)
        # save_adjusted=True also saves the adjusted image
        self.save_image_with_points(save_adjusted=True)
        self.save_segment_data()
        self.save_maps_data()
        self.save_points()

    def on_clear(self, event):
        """Clear button"""
        self.clicked_points = []
        self.clicked_point_types = []
        self.clicked_point_segments = []
        for plot in self.point_plots:
            plot.remove()
        for text in self.text_plots:
            text.remove()
        self.point_plots = []
        self.text_plots = []
        self.selected_point_idx = None
        self.insert_after_idx = None

        self.fig.canvas.draw_idle()
        print("\nAll points have been reset.")

    def on_save_points(self, event):
        """Save Points button"""
        if len(self.clicked_points) < 2:
            print("\nWarning: at least 2 points are required!")
            return

        self.points = self.clicked_points
        self.point_types = self.clicked_point_types
        self.point_segments = self.clicked_point_segments
        self.save_points()

    def on_load_points(self, event):
        """Load Points button"""
        img_dir = Path(self.image_path).parent
        base_name = Path(self.image_path).stem
        json_path = str(img_dir / f"{base_name}_points.json")

        if not os.path.exists(json_path):
            print(f"\nSaved points file not found: {json_path}")
            return

        points = self.load_points(json_path)
        if points:
            self.clicked_points = points
            self.clicked_point_types = self.point_types
            self.clicked_point_segments = self.point_segments
            self.redraw_points()

            # update slider values and refresh the image
            self.slider_contrast.set_val(self.contrast)
            self.slider_brightness.set_val(self.brightness)
            adjusted = self.apply_contrast_brightness()
            self.img_display.set_data(adjusted)
            self.fig.canvas.draw_idle()

            print(f"\n{len(points)} points loaded")

    def on_reset_view(self, event):
        """Reset View button"""
        if self.home_view is not None:
            self.ax.set_xlim(self.home_view[0])
            self.ax.set_ylim(self.home_view[1])
            self.fig.canvas.draw_idle()

    def run_interactive(self):
        """Run interactive mode"""
        self.fig = plt.figure(figsize=(14, 10))

        ax_left = 0.05
        ax_bottom = 0.14
        ax_width = 0.75
        ax_height = 0.82
        self.ax = plt.axes([ax_left, ax_bottom, ax_width, ax_height])
        # display the normalized image
        img_norm = (self.original_image - self.original_image.min()) / \
                   (self.original_image.max() - self.original_image.min() + 1e-10)
        self.img_display = self.ax.imshow(img_norm, cmap='gray', vmin=0, vmax=1)

        # expand ylim to fill the axes area
        fig_w, fig_h = self.fig.get_size_inches()
        axes_w_inch = fig_w * ax_width
        axes_h_inch = fig_h * ax_height
        axes_aspect = axes_h_inch / axes_w_inch

        img_h, img_w = self.original_image.shape
        total_y_range = img_w * axes_aspect
        pad = (total_y_range - img_h) / 2 * 0.8  # reduce padding by 20%
        self.ax.set_ylim(img_h + pad, -pad)

        self.ax.set_title('Navigation Mode: Wheel=Zoom, Right-drag=Pan',
                         fontsize=14, color='blue')
        self.ax.axis('off')

        self.home_view = (self.ax.get_xlim(), self.ax.get_ylim())

        # contrast/brightness sliders
        slider_contrast_ax = plt.axes([0.15, 0.04, 0.25, 0.025])
        slider_brightness_ax = plt.axes([0.55, 0.04, 0.25, 0.025])

        self.slider_contrast = Slider(slider_contrast_ax, 'Contrast', 0.1, 5.0,
                                       valinit=1.0, valstep=0.1)
        self.slider_brightness = Slider(slider_brightness_ax, 'Brightness', -1.0, 1.0,
                                          valinit=0.0, valstep=0.05)

        self.slider_contrast.on_changed(self.on_contrast_change)
        self.slider_brightness.on_changed(self.on_brightness_change)

        # radio buttons
        radio_ax = plt.axes([0.82, 0.55, 0.15, 0.28])
        radio_ax.set_title('Mode', fontsize=12, weight='bold')
        self.radio = RadioButtons(radio_ax,
                                   ('Navigate', 'Add Measure', 'Insert After',
                                    'Add Guide', 'Delete', 'Move Point'),
                                   active=0)
        self.radio.on_clicked(self.on_mode_change)

        # buttons
        btn_trace_ax = plt.axes([0.82, 0.48, 0.15, 0.05])
        btn_save_ax = plt.axes([0.82, 0.41, 0.15, 0.05])
        btn_load_ax = plt.axes([0.82, 0.34, 0.15, 0.05])
        btn_clear_ax = plt.axes([0.82, 0.27, 0.15, 0.05])
        btn_reset_ax = plt.axes([0.82, 0.20, 0.15, 0.05])

        btn_trace = Button(btn_trace_ax, 'Trace DNA', color='lightgreen')
        btn_save = Button(btn_save_ax, 'Save Points', color='lightyellow')
        btn_load = Button(btn_load_ax, 'Load Points', color='lightcyan')
        btn_clear = Button(btn_clear_ax, 'Clear All', color='lightcoral')
        btn_reset = Button(btn_reset_ax, 'Reset View', color='lightblue')

        btn_trace.on_clicked(self.on_trace)
        btn_save.on_clicked(self.on_save_points)
        btn_load.on_clicked(self.on_load_points)
        btn_clear.on_clicked(self.on_clear)
        btn_reset.on_clicked(self.on_reset_view)

        # guide text
        guide_text = (
            "=== User guide ===\n\n"
            "Measure (green): measurement reference point\n"
            "Guide (blue): waypoint\n"
            "Delete: delete point\n"
            "Move: move point"
        )

        guide_ax = plt.axes([0.82, 0.02, 0.15, 0.16])
        guide_ax.axis('off')
        guide_ax.text(0, 1, guide_text, fontsize=7.5,
                     verticalalignment='top', family='monospace')

        # connect events
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.fig.canvas.mpl_connect('button_press_event', self.on_pan_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_pan_motion)
        self.fig.canvas.mpl_connect('button_release_event', self.on_pan_release)

        # draw preloaded points if any
        if len(self.clicked_points) > 0:
            self.redraw_points()
            # update slider values
            self.slider_contrast.set_val(self.contrast)
            self.slider_brightness.set_val(self.brightness)
            adjusted = self.apply_contrast_brightness()
            self.img_display.set_data(adjusted)

        plt.show()


if __name__ == "__main__":
    print("\nThis is a library file.")
    print("Run run_interactive.py!")