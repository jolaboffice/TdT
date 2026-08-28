#!/usr/bin/env python3
"""
hDOM_tif2html.py
Creates optimized HTML reports for individual hDOM analysis results.
Generates individual HTML reports, summary index page, chromosome coverage,
and webserver symlink.
"""

import os
import sys
import argparse
import io
import base64
import re
import glob
import pandas as pd
import numpy as np
from PIL import Image, ImageEnhance
from DOM_lib import tif2RGimage

# Constants
BPP = 200.0
BRIGHTNESS_ENHANCE = 1.5
REF_BRIGHTNESS_ENHANCE = 1.2

# Poly-N overlay configuration
POLYN_MIN_BP = 100  # minimum N-run length to mark
POLYN_OVERLAY_COLOR = (0, 0, 255)  # solid blue (R, G, B)

def image_to_base64(image):
    """Convert PIL Image to base64 string."""
    buf = io.BytesIO()
    image.save(buf, format='PNG')
    return base64.b64encode(buf.getvalue()).decode('utf-8')

def process_query_image(tif_path):
    """Load, enhance, and pad the query TIF image."""
    if not os.path.exists(tif_path):
        print(f"Warning: .tif file not found: {tif_path}")
        return None, None
    try:
        img, _, _ = tif2RGimage(tif_path)
        img = ImageEnhance.Brightness(img).enhance(BRIGHTNESS_ENHANCE)

        # Add 10px black padding
        w, h = img.size
        padded_w = w + 20
        new_img = Image.new('RGB', (padded_w, h), (0, 0, 0))
        new_img.paste(img, (10, 0))

        return new_img, image_to_base64(new_img)
    except Exception as e:
        print(f"Error processing query image: {e}")
        return None, None

def _detect_ref_bpp(ref_map_folder):
    """Detect bpp from ref_map folder name, e.g. 'ref_map_<N>bpp' → N.
    Falls back to module-level BPP (200) if no bpp in folder name."""
    m = re.search(r'(\d+)bpp', os.path.basename(ref_map_folder or ''))
    if m:
        return float(m.group(1))
    return BPP


# Cache: {(abs_ref_map_folder, ref_bpp, min_n_bp): {chr_name: [(start_px, end_px), ...]}}
_n_regions_cache = {}


def _find_n_regions_for_ref(ref_map_folder, ref_bpp, min_n_bp=POLYN_MIN_BP):
    """Scan *.fasta files in ref_map_folder and return poly-N pixel ranges
    (in ref-image pixel units) keyed by chr name (FASTA basename without ext)."""
    n_regions = {}
    if not ref_map_folder or not os.path.isdir(ref_map_folder):
        return n_regions
    fasta_files = glob.glob(os.path.join(ref_map_folder, "*.fasta"))
    if not fasta_files:
        return n_regions
    pat = re.compile(r'N{' + str(int(min_n_bp)) + r',}')
    for fasta_file in fasta_files:
        chr_name = os.path.splitext(os.path.basename(fasta_file))[0]
        seq_parts = []
        with open(fasta_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                seq_parts.append(line.strip().upper())
        seq = ''.join(seq_parts)
        regions = [(m.start() / ref_bpp, m.end() / ref_bpp) for m in pat.finditer(seq)]
        if regions:
            n_regions[chr_name] = regions
    return n_regions


def _get_n_regions_cached(ref_map_folder, ref_bpp, min_n_bp=POLYN_MIN_BP):
    if not ref_map_folder:
        return {}
    key = (os.path.abspath(ref_map_folder), float(ref_bpp), int(min_n_bp))
    if key not in _n_regions_cache:
        _n_regions_cache[key] = _find_n_regions_for_ref(ref_map_folder, ref_bpp, min_n_bp)
    return _n_regions_cache[key]


def _apply_polyn_overlay(arr, raw_r, raw_g, ref_name, ref_map_folder, ref_bpp,
                         crop_start_ref_px, crop_end_ref_px, direction,
                         final_width):
    """Paint blue on pixels that are (a) inside a poly-N ref column AND
    (b) carry no original R/G signal — i.e. the BPP bin is purely N.

    raw_r, raw_g: per-pixel R, G channels of `arr` BEFORE channel adjustment.
    Coordinates: crop_start/end_ref_px are in ref-image pixel units (pre-crop, pre-flip).
    `arr` is already cropped, flipped (if direction==-1), and resized to final_width.
    """
    n_regions = _get_n_regions_cached(ref_map_folder, ref_bpp)
    ref_regions = n_regions.get(ref_name)
    if not ref_regions:
        return
    crop_w_ref = max(crop_end_ref_px - crop_start_ref_px, 1)
    sx = final_width / crop_w_ref
    region_cols = np.zeros(final_width, dtype=bool)
    for n_start_ref, n_end_ref in ref_regions:
        rs = max(n_start_ref, crop_start_ref_px) - crop_start_ref_px
        re_ = min(n_end_ref, crop_end_ref_px) - crop_start_ref_px
        if re_ <= rs:
            continue
        if direction == -1:
            rs, re_ = crop_w_ref - re_, crop_w_ref - rs
        x0 = int(round(rs * sx))
        x1 = int(round(re_ * sx))
        x0 = max(0, x0); x1 = min(final_width, x1)
        if x1 > x0:
            region_cols[x0:x1] = True
    if not region_cols.any():
        return
    # "No signal" = pixel sits at the per-crop baseline of both channels.
    # Sim refs use a non-zero floor (e.g. uint8 5), so == 0 misses everything.
    r_base = int(raw_r.min())
    g_base = int(raw_g.min())
    pure_n = (raw_r == r_base) & (raw_g == g_base) & region_cols[None, :]
    arr[pure_n] = POLYN_OVERLAY_COLOR


def process_ref_match(match_info, ref_map_folder, query_final_width):
    """Crop, enhance, and flip reference matching region based on direction."""
    ref_name = match_info['ref_name']
    aligned_bp = match_info['aligned_bp']
    scale = match_info['scale']
    direction = match_info.get('dir', 1)

    ref_path = os.path.join(ref_map_folder, f"{ref_name}.tif")
    if not os.path.exists(ref_path):
        return None

    ref_bpp = _detect_ref_bpp(ref_map_folder)

    try:
        ref_img, ref_w, ref_h = tif2RGimage(ref_path)
        shift_px = int(aligned_bp / ref_bpp)

        # Calculate region width
        # query_w is the query's display width in px at query bpp (BPP module constant).
        # The ref region in ref pixels must cover the same bp span scaled by `scale`,
        # converted from query bpp to ref bpp.
        query_w = max(query_final_width - 20, 1)
        ref_region_w = int(query_w * scale * BPP / ref_bpp)

        # Crop with 10px buffer
        start_x = max(0, shift_px - 10)
        end_x = min(ref_w, shift_px + ref_region_w + 10)

        if start_x >= end_x:
            return None

        crop = ref_img.crop((start_x, 0, end_x, ref_h))

        # Handle reverse direction
        if direction == -1:
            crop = crop.transpose(Image.Transpose.FLIP_LEFT_RIGHT)

        crop = crop.resize((query_final_width, ref_h))

        # Channel adjustments
        arr = np.array(crop)

        # Snapshot raw R/G *before* adjustment, for poly-N "pure N pixel" mask
        raw_r = arr[:, :, 0].copy()
        raw_g = arr[:, :, 1].copy()

        # Normalize Red
        red = arr[:, :, 0].astype(np.float32)
        r_min, r_max = np.min(red), np.max(red)
        if r_max > r_min:
            red = ((red - r_min) / (r_max - r_min)) * 255.0
        arr[:, :, 0] = np.clip(red, 0, 255).astype(np.uint8)

        # Increase Green
        arr[:, :, 1] = np.clip((arr[:, :, 1].astype(np.float32) * 2.0), 0, 255).astype(np.uint8)

        # Poly-N blue overlay: only pixels with no original R/G signal,
        # within fasta-derived N-region columns
        _apply_polyn_overlay(
            arr, raw_r, raw_g, ref_name, ref_map_folder, ref_bpp,
            crop_start_ref_px=start_x, crop_end_ref_px=end_x,
            direction=direction, final_width=query_final_width,
        )

        enhanced = Image.fromarray(arr)
        enhanced = ImageEnhance.Brightness(enhanced).enhance(REF_BRIGHTNESS_ENHANCE)

        return image_to_base64(enhanced)
    except Exception as e:
        print(f"Error processing ref {ref_name}: {e}")
        return None

def get_table_html(df, folder_name="", file_name="", ref_map_path=None):
    """Format and generate HTML table from matching results with hyperlinks on rank."""
    if df is None or df.empty:
        return "<p>No matching results available.</p>"

    df_fmt = df.copy()
    cols = df_fmt.columns

    # Formatting
    if 'aligned_bp' in cols:
        df_fmt['aligned_bp'] = df_fmt['aligned_bp'].apply(lambda x: f"{int(x):,}" if pd.notna(x) else "")
    if 'shift_bp' in cols:
        df_fmt['shift_bp'] = df_fmt['shift_bp'].apply(lambda x: f"{int(x):,}" if pd.notna(x) else "")

    for c in ['scale', 'cc_r', 'cc_g']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(lambda x: f"{x:.2f}" if pd.notna(x) and isinstance(x, (int, float)) else x)

    for c in ['cc_rg', 'cc_rg2']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(lambda x: f"{x:.1e}" if pd.notna(x) and isinstance(x, (int, float)) else x)

    # Merge log(P) (PCC score) and misc (maligner rescaled m_score) into a single
    # "score" column: PCC-method rows show log(P), maligner rows show the m-score.
    if 'method' in cols and ('log(P)' in cols or 'misc' in cols):
        methods = df['method'].tolist()
        score_vals = []
        for i in range(len(df_fmt)):
            m = str(methods[i]).lower() if i < len(methods) else ''
            if 'mal' in m:
                v = df.iloc[i].get('misc')
                score_vals.append(f"{v:.2f}" if pd.notna(v) and pd.api.types.is_number(v) else "")
            else:
                v = df.iloc[i].get('log(P)')
                score_vals.append(f"{v:.2f}" if pd.notna(v) and pd.api.types.is_number(v) else "")
        df_fmt = df_fmt.drop(columns=[c for c in ('log(P)', 'misc') if c in df_fmt.columns])
        df_fmt['score'] = score_vals

    # Hyperlink for rank. 'sel' stays the table row position (1,2,...) so view.php
    # selects the right row, but the visible label uses the per-method rank:
    # pcc-rank{n} / mal-rank{n} / KU-rank{n} (logical -> KU).
    _method_prefix = {'pcc': 'pcc', 'mal': 'mal', 'logical': 'KU'}
    if 'rank' in cols:
        methods = df['method'].tolist() if 'method' in df.columns else [''] * len(df_fmt)
        ranks = df['rank'].tolist() if 'rank' in df.columns else list(range(1, len(df_fmt) + 1))
        new_ranks = []
        for i in range(len(df_fmt)):
            sel_idx = i + 1
            method = str(methods[i]).lower() if i < len(methods) else ''
            rank_raw = str(ranks[i])
            if 'poly N' in rank_raw:
                label = rank_raw
            else:
                try:
                    rk = int(ranks[i])
                except (TypeError, ValueError):
                    rk = sel_idx
                label = f"{_method_prefix.get(method, method)}-rank{rk}"
            ref_param = f'&ref_map={ref_map_path}' if ref_map_path else ''
            new_ranks.append(f'<a href="../view.php?folder={folder_name}&file={file_name}&sel={sel_idx}{ref_param}">{label}</a>')
        df_fmt['rank'] = new_ranks

    # Display-only method name: 'logical' -> 'KU' (pcc/mal unchanged)
    if 'method' in cols:
        df_fmt['method'] = df_fmt['method'].apply(
            lambda x: 'KU' if str(x).lower() == 'logical' else x)

    if 'dir' in cols:
        df_fmt['dir'] = df_fmt['dir'].apply(lambda x: f"{int(x)}" if pd.notna(x) else "")

    # Display-only header relabel: 'ref_name' -> "chromosome".
    df_fmt = df_fmt.rename(columns={'ref_name': 'chromosome'})

    html = df_fmt.to_html(index=False, classes='matching-table', border=0, escape=False)

    # Add conditional row classes (PCC/MAL coloring and Rank highlighting)
    if 'method' in cols:
        parts = html.split('<tbody>')
        if len(parts) == 2:
            header, body = parts[0] + '<tbody>', parts[1]
            rows = re.findall(r'<tr>.*?</tr>', body, re.DOTALL)
            new_rows = []
            for i, row in enumerate(rows):
                rank_raw = str(df.iloc[i].get('rank', ''))
                m_val = str(df.iloc[i]['method']).lower() if pd.notna(df.iloc[i]['method']) else ''

                cls = []
                if 'poly N' in rank_raw: cls.append('method-polyn')
                elif 'pcc' in m_val: cls.append('method-pcc')
                elif 'mal' in m_val: cls.append('method-mal')

                try:
                    # Highlight Top 2 (Rank 1 and 2 in CSV)
                    if int(rank_raw) in (1, 2):
                        cls.append('rank-top')
                except (TypeError, ValueError):
                    pass

                if cls:
                    row = row.replace('<tr>', f'<tr class="{" ".join(cls)}">', 1)
                new_rows.append(row)

            # Reconstruct with closing tags
            footer_match = re.search(r'(</tbody>.*?</table>)', body, re.DOTALL)
            footer = footer_match.group(1) if footer_match else '</tbody></table>'
            html = header + "".join(new_rows) + footer

    return html

def get_html_style():
    """Return the CSS styles for the report."""
    return """
    body { font-family: 'Segoe UI', sans-serif; line-height: 1.6; color: #333; max-width: 1200px; margin: 0 auto; padding: 40px 20px; background-color: #f9f9f9; }
    .container { background-color: #fff; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #2980b9; margin-top: 30px; }
    .info-box { background-color: #f8f9fa; border-left: 5px solid #3498db; padding: 15px; margin-bottom: 20px; white-space: pre-wrap; font-family: monospace; font-size: 0.95em; overflow-x: auto; }
    .image-container { margin: 30px 0; text-align: center; background-color: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #ddd; }
    .comparison-box { margin: 30px 0; padding: 20px; background-color: #000; border-radius: 5px; }
    .comparison-images { display: flex; flex-direction: column; gap: 0; }
    .comparison-item { display: flex; align-items: center; width: 100%; }
    .image-label { min-width: 60px; padding: 0 15px; font-weight: bold; color: #fff; text-align: right; }
    .comparison-item .image-container { margin: 0; padding: 0; border: none; flex: 1; background-color: #000; }
    .comparison-item img { display: block; margin: 0 auto; max-width: 100%; height: auto; }
    .rank-separator { height: 1px; background-color: #444; margin: 4px 0; }
    .table-container { overflow-x: auto; margin-top: 15px; }
    .matching-table { width: 100%; border-collapse: collapse; font-size: 0.9em; }
    .matching-table th { background-color: #3498db; color: #fff; text-align: right; padding: 12px; position: sticky; top: 0; }
    .matching-table th:first-child { text-align: left; }
    .matching-table td { padding: 10px 12px; border-bottom: 1px solid #eee; text-align: right; }
    .matching-table td:first-child { text-align: left; }
    .matching-table tr:hover { background-color: #f1f7ff; }
    .method-pcc { background-color: #e8f4f8; }
    .method-pcc:nth-child(even) { background-color: #d4e9f0; }
    .method-mal { background-color: #fff4e6; }
    .method-mal:nth-child(even) { background-color: #ffe8cc; }
    .method-polyn { background-color: #f0e0f0; }
    .method-polyn:nth-child(even) { background-color: #e0d0e0; }
    .rank-top { font-weight: bold; }
    .summary-table { width: 100%; border-collapse: collapse; margin-top: 20px; }
    .summary-table th { background-color: #2c3e50; color: #fff; padding: 12px; text-align: left; cursor: pointer; user-select: none; }
    .summary-table td { padding: 10px 12px; border-bottom: 1px solid #ddd; }
    .summary-table tr:hover { background-color: #f5f5f5; }
    .text-right { text-align: right !important; }
    .text-center { text-align: center !important; }
    footer { margin-top: 50px; font-size: 0.8em; color: #888; text-align: center; }
    """

def find_ref_map_folder(data_dir):
    """Find ref_map folder: prefer bpp-specific folder if BPP != 200."""
    parent_dir = os.path.dirname(data_dir) or os.getcwd()
    # If BPP != 200, check for bpp-specific ref_map first
    if BPP != 200.0:
        bpp_name = f"ref_map_{int(BPP)}bpp"
        for candidate in [
            os.path.join(data_dir, bpp_name),
            os.path.join(parent_dir, bpp_name),
        ]:
            if os.path.isdir(candidate):
                return candidate
    for candidate in [
        os.path.join(data_dir, "ref_map"),
        os.path.join(parent_dir, "ref_map"),
        "ref_map",
    ]:
        if os.path.isdir(candidate):
            return candidate
    return "ref_map"

def read_ref_map_from_info(tif_path):
    """Read ref_map path from _info.txt file."""
    info_path = os.path.splitext(tif_path)[0] + "_info.txt"
    if not os.path.exists(info_path):
        return None
    try:
        with open(info_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('ref_map:'):
                    path = line.split(':', 1)[1].strip()
                    if os.path.isdir(path):
                        return path
    except Exception:
        pass
    return None

def tif2html(tif_path):
    """Generates an HTML report and returns summary information."""
    ext = os.path.splitext(tif_path)[1].lower()
    if ext not in (".tif", ".tiff"):
        print(f"Error: Input must be a .tif or .tiff file. Got: {tif_path}")
        return None

    base_path = os.path.splitext(tif_path)[0]
    info_path = f"{base_path}_info.txt"
    csv_path = f"{base_path}_matching_results.csv"
    html_path = f"{base_path}.html"

    # 1. Load Info
    info_content = f"Information file not found: {info_path}"
    if os.path.exists(info_path):
        with open(info_path, 'r', encoding='utf-8') as f:
            info_content = f.read()

    # 2. Process Images
    query_img, query_b64 = process_query_image(tif_path)
    if not query_img:
        return None

    query_final_w = query_img.size[0]
    molecule_length_bp = int((query_final_w - 20) * BPP)

    # 3. Read Results & Generate Comparisons
    df = pd.read_csv(csv_path) if os.path.exists(csv_path) else None
    comp_html = ""
    summary_info = {
        'tif_name': os.path.basename(tif_path),
        'html_link': os.path.basename(html_path),
        'ref_map': 'N/A',
        'aligned_bp': 'N/A',
        'length_bp': molecule_length_bp,
        'direction': '+',
        'log_p': 'N/A'
    }

    if df is not None and not df.empty:
        top_hit = df.iloc[0]
        summary_info['ref_map'] = top_hit.get('ref_name', 'N/A')
        summary_info['log_p'] = top_hit.get('log(P)', 'N/A')
        aligned = top_hit.get('aligned_bp')
        if pd.isna(aligned):
            aligned = top_hit.get('shift_bp', 0)
        summary_info['aligned_bp'] = int(aligned)

        # Handle direction in summary
        direction = top_hit.get('dir', 1)
        summary_info['direction'] = '+' if direction == 1 else '-'

        # Find ref_map folder (prefer path from _info.txt)
        tif_dir = os.path.dirname(os.path.abspath(tif_path))
        ref_map_folder = read_ref_map_from_info(tif_path) or find_ref_map_folder(tif_dir)

        # Show top-2 hits from each of the three CSV rank blocks:
        # ranks 1-2, 11-12, and 21-22 (blocks of 10 rows each).
        display_ranks = [1, 2, 11, 12, 21, 22]
        separator_after = {2, 12}

        ref_htmls = ""
        for r in display_ranks:
            idx = r - 1
            if idx >= len(df):
                continue
            row = df.iloc[idx]
            ref_name = row.get('ref_name')
            if pd.isna(ref_name) or not ref_name:
                continue
            aligned_bp = row.get('aligned_bp')
            if pd.isna(aligned_bp):
                aligned_bp = row.get('shift_bp', 0)
            scale = row.get('scale', 1.0)
            if pd.isna(scale):
                scale = 1.0
            dir_val = row.get('dir', 1)
            ref_b64 = process_ref_match(
                {'ref_name': ref_name, 'aligned_bp': aligned_bp, 'scale': scale, 'dir': dir_val},
                ref_map_folder, query_final_w)
            if ref_b64:
                ref_htmls += f'<div class="comparison-item"><div class="image-label">rank{r}</div>' \
                             f'<div class="image-container"><img src="data:image/png;base64,{ref_b64}"></div></div>'
            if r in separator_after:
                ref_htmls += '<div class="rank-separator"></div>'

        comp_html = f'<div class="comparison-box"><div class="comparison-images">' \
                    f'<div class="comparison-item"><div class="image-label">data</div>' \
                    f'<div class="image-container"><img src="data:image/png;base64,{query_b64}"></div></div>' \
                    f'{ref_htmls}</div></div>'
    else:
        comp_html = f'<div class="image-container"><img src="data:image/png;base64,{query_b64}"></div>'

    # 4. Generate Final HTML
    tif_abs = os.path.abspath(tif_path)
    folder_name = os.path.basename(os.path.dirname(tif_abs))
    file_name = os.path.basename(tif_abs)
    table_html = get_table_html(df, folder_name=folder_name, file_name=file_name, ref_map_path=ref_map_folder)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>hDOM Analysis - {os.path.basename(tif_path)}</title>
    <style>{get_html_style()}</style>
</head>
<body>
    <div class="container">
        <h1>hDOM Analysis Report</h1>
        <p><strong>File:</strong> {os.path.basename(tif_path)} ({molecule_length_bp / 1000:.1f} kbp)</p>
        {comp_html}
        <h2>Matching Results</h2>
        <div class="table-container">{table_html}</div>
        <h2>Information</h2>
        <div class="info-box">{info_content}</div>
    </div>
    <footer>Generated by hDOM_tif2html.py</footer>
</body>
</html>"""

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created: {html_path}")
    return summary_info

# --- Summary & Index ---

def generate_index_html(summaries, output_path):
    """Creates a summary index HTML page with sortable columns."""
    rows_html = ""
    for idx, s in enumerate(summaries, start=1):
        aligned_bp_fmt = f"{s['aligned_bp']:,}" if isinstance(s['aligned_bp'], (int, float)) else s['aligned_bp']
        length_bp_fmt = f"{s['length_bp']:,}" if isinstance(s['length_bp'], (int, float)) else s['length_bp']
        log_p_fmt = f"{s['log_p']:.2f}" if isinstance(s['log_p'], (int, float)) else s['log_p']
        rows_html += f"""
        <tr>
            <td class="row-num text-right">{idx}</td>
            <td data-value="{s['tif_name'].lower()}"><a href="{s['html_link']}">{s['tif_name']}</a></td>
            <td data-value="{s['ref_map'].lower()}">{s['ref_map']}</td>
            <td class="text-right" data-value="{s['aligned_bp']}">{aligned_bp_fmt}</td>
            <td class="text-right" data-value="{s['length_bp']}">{length_bp_fmt}</td>
            <td class="text-right" data-value="{s['log_p']}">{log_p_fmt}</td>
            <td class="text-center" data-value="{s['direction']}">{s['direction']}</td>
        </tr>"""

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>hDOM Analysis Summary Index</title>
    <style>{get_html_style()}</style>
    <script>
    function sortTable(n) {{
        var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
        table = document.getElementById("summaryTable");
        switching = true;
        dir = "asc";
        while (switching) {{
            switching = false;
            rows = table.rows;
            for (i = 1; i < (rows.length - 1); i++) {{
                shouldSwitch = false;
                x = rows[i].getElementsByTagName("TD")[n];
                y = rows[i + 1].getElementsByTagName("TD")[n];

                var xVal = x.getAttribute("data-value");
                var yVal = y.getAttribute("data-value");

                // Multi-level sort for Aligned Position (column 3):
                // Primary: Reference Map (column 2), Secondary: Aligned Position (column 3)
                if (n === 3) {{
                    var xRef = rows[i].getElementsByTagName("TD")[2].getAttribute("data-value").toLowerCase();
                    var yRef = rows[i + 1].getElementsByTagName("TD")[2].getAttribute("data-value").toLowerCase();
                    if (xRef !== yRef) {{
                        if (dir == "asc") {{
                            if (xRef > yRef) {{ shouldSwitch = true; break; }}
                        }} else {{
                            if (xRef < yRef) {{ shouldSwitch = true; break; }}
                        }}
                        continue;
                    }}
                }}

                // Try numeric comparison if both are numbers
                var xNum = parseFloat(xVal);
                var yNum = parseFloat(yVal);

                if (!isNaN(xNum) && !isNaN(yNum)) {{
                    if (dir == "asc") {{
                        if (xNum > yNum) {{ shouldSwitch = true; break; }}
                    }} else if (dir == "desc") {{
                        if (xNum < yNum) {{ shouldSwitch = true; break; }}
                    }}
                }} else {{
                    if (dir == "asc") {{
                        if (xVal.toLowerCase() > yVal.toLowerCase()) {{ shouldSwitch = true; break; }}
                    }} else if (dir == "desc") {{
                        if (xVal.toLowerCase() < yVal.toLowerCase()) {{ shouldSwitch = true; break; }}
                    }}
                }}
            }}
            if (shouldSwitch) {{
                rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                switching = true;
                switchcount ++;
            }} else {{
                if (switchcount == 0 && dir == "asc") {{
                    dir = "desc";
                    switching = true;
                }}
            }}
        }}
        // Re-number the # column
        rows = table.rows;
        for (i = 1; i < rows.length; i++) {{
            rows[i].getElementsByTagName("TD")[0].innerHTML = i;
        }}
    }}
    </script>
</head>
<body>
    <div class="container">
        <h1>hDOM Analysis Summary Index</h1>
        <p><a href="chr_cover.html">Chromosome Coverage</a></p>
        <table class="summary-table" id="summaryTable">
            <thead>
                <tr>
                    <th class="text-right">#</th>
                    <th onclick="sortTable(1)">TIF File</th>
                    <th onclick="sortTable(2)">Reference Map</th>
                    <th class="text-right" onclick="sortTable(3)">Aligned Position (bp)</th>
                    <th class="text-right" onclick="sortTable(4)">Molecule Length (bp)</th>
                    <th class="text-right" onclick="sortTable(5)">log(P)</th>
                    <th class="text-center" onclick="sortTable(6)">Dir</th>
                </tr>
            </thead>
            <tbody>
                {rows_html}
            </tbody>
        </table>
    </div>
    <footer>Generated by hDOM_tif2html.py</footer>
</body>
</html>"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created Summary Index: {output_path}")


def create_webserver_symlink(source_dir):
    """Create a symbolic link in the webserver directory (~/html/hDOM) and link ref_map folder."""
    try:
        source_dir = os.path.realpath(source_dir)
        webserver_base = os.path.realpath(os.path.expanduser("~/html/hDOM"))
        if not os.path.exists(webserver_base):
            os.makedirs(webserver_base, exist_ok=True)

        dir_name = os.path.basename(source_dir)
        link_path = os.path.join(webserver_base, dir_name)

        # Skip if source_dir is already under the webserver directory
        if source_dir.startswith(webserver_base):
            return link_path

        # 1. Link data folder to webserver path
        if os.path.islink(link_path) or os.path.exists(link_path):
            if os.path.islink(link_path):
                os.unlink(link_path)
            else:
                import shutil
                if os.path.isdir(link_path):
                    shutil.rmtree(link_path)
                else:
                    os.remove(link_path)

        os.symlink(source_dir, link_path)
        print(f"Created symbolic link: {link_path} -> {source_dir}")

        # 2. Link ref_map folder into the data folder
        parent_dir = os.path.dirname(source_dir)
        ref_map_src = os.path.join(parent_dir, "ref_map")
        if os.path.isdir(ref_map_src):
            ref_map_dest = os.path.join(source_dir, "ref_map")

            # Remove if exists and is not pointing to the correct place
            if os.path.exists(ref_map_dest) or os.path.islink(ref_map_dest):
                if os.path.islink(ref_map_dest):
                    if os.path.abspath(os.readlink(ref_map_dest)) != ref_map_src:
                        os.unlink(ref_map_dest)
                    else:
                        return link_path
                else:
                    print(f"Warning: {ref_map_dest} exists and is not a symbolic link. Skipping.")
                    return link_path

            if not os.path.exists(ref_map_dest):
                os.symlink(ref_map_src, ref_map_dest)
                print(f"Created internal link: {ref_map_dest} -> {ref_map_src}")

        return link_path
    except Exception as e:
        print(f"Error creating symbolic link: {e}")
        return None


def generate_reports(tif_files):
    """Generate individual HTML reports and summary index.html.
    Returns the index.html path, or None if no reports were generated.
    """
    summaries = []
    for tif_file in tif_files:
        summary = tif2html(tif_file)
        if summary:
            summaries.append(summary)

    if not summaries:
        print("No reports generated.")
        return None

    first_file_dir = os.path.dirname(os.path.abspath(tif_files[0]))

    # Generate chromosome coverage
    import hDOM_chr_cover
    ref_map_folder = read_ref_map_from_info(tif_files[0]) or find_ref_map_folder(first_file_dir)
    chr_cover_path = os.path.join(first_file_dir, "chr_cover.html")
    hDOM_chr_cover.generate_chr_cover(first_file_dir, ref_map_folder, chr_cover_path)

    index_path = os.path.join(first_file_dir, "index.html")
    generate_index_html(summaries, index_path)
    create_webserver_symlink(first_file_dir)
    return index_path


def main():
    parser = argparse.ArgumentParser(description="Convert hDOM .tif results to HTML reports.")
    parser.add_argument("tif_files", nargs='+', help="One or more .tif files or a directory containing .tif files")
    args = parser.parse_args()

    all_tif_files = []
    for path in args.tif_files:
        if os.path.isdir(path):
            files = [os.path.join(path, f) for f in os.listdir(path) if f.lower().endswith(('.tif', '.tiff'))]
            all_tif_files.extend(sorted(files))
        else:
            all_tif_files.append(path)

    if not all_tif_files:
        print("No .tif files found to process.")
        sys.exit(1)

    index_path = generate_reports(all_tif_files)
    if index_path:
        print(f"\nHTML reports generated successfully.")
        print(f"Index page: {index_path}")


if __name__ == "__main__":
    main()
