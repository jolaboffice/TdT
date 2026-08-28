#!/usr/bin/env python3
"""
hDOM_chr_cover.py
Generate chromosome coverage visualization (chr_cover.html).
Shows all chromosomes as horizontal bars with matched molecule coverage highlighted.
"""

import os
import re
import pandas as pd
from DOM_lib import tif2RGimage

BPP = 200.0


def get_chromosome_lengths(ref_map_folder):
    """Get chromosome lengths (in bp) from reference TIF file widths."""
    chr_lengths = {}
    if not os.path.isdir(ref_map_folder):
        return chr_lengths
    for f in os.listdir(ref_map_folder):
        if f.lower().endswith(('.tif', '.tiff')):
            ref_name = os.path.splitext(f)[0]
            ref_path = os.path.join(ref_map_folder, f)
            try:
                _, w, _ = tif2RGimage(ref_path)
                chr_lengths[ref_name] = int(w * BPP)
            except Exception:
                pass
    return chr_lengths


def chromosome_sort_key(name):
    """Sort chromosomes naturally: chr01_p, chr01_m, chr02_p, ..., chrX_m, chrY_p."""
    match = re.match(r'chr(\d+|X|Y)_(p|m)', name, re.IGNORECASE)
    if match:
        num_str, suffix = match.groups()
        if num_str.upper() == 'X':
            num = 23
        elif num_str.upper() == 'Y':
            num = 24
        else:
            num = int(num_str)
        suffix_order = 0 if suffix == 'p' else 1
        return (num, suffix_order)
    return (999, 0)


def collect_coverage(data_folder):
    """Collect coverage data from all matching results CSVs (top hit only)."""
    coverage = []
    for f in sorted(os.listdir(data_folder)):
        if not f.endswith('_matching_results.csv'):
            continue
        csv_path = os.path.join(data_folder, f)
        tif_name = f.replace('_matching_results.csv', '.tif')
        tif_path = os.path.join(data_folder, tif_name)
        try:
            df = pd.read_csv(csv_path)
            if df.empty:
                continue
            top = df.iloc[0]
            ref_name = top.get('ref_name')
            if pd.isna(ref_name):
                continue
            aligned_bp = top.get('aligned_bp')
            if pd.isna(aligned_bp):
                aligned_bp = top.get('shift_bp', 0)
            aligned_bp = float(aligned_bp)
            scale = top.get('scale', 1.0)
            if pd.isna(scale):
                scale = 1.0

            # Get molecule length from tif
            if os.path.exists(tif_path):
                _, w, _ = tif2RGimage(tif_path)
                mol_length_bp = w * BPP * scale
            else:
                continue

            log_p = top.get('log(P)', 0)
            if pd.isna(log_p):
                log_p = 0
            log_p = float(log_p)

            if mol_length_bp > 0:
                coverage.append({
                    'ref_name': ref_name,
                    'start_bp': int(aligned_bp),
                    'end_bp': int(aligned_bp + mol_length_bp),
                    'tif_name': tif_name,
                    'html_link': tif_name.replace('.tif', '.html'),
                    'log_p': log_p,
                })
        except Exception:
            pass
    return coverage


def generate_chr_cover(data_folder, ref_map_folder, output_path):
    """Generate chromosome coverage HTML."""
    chr_lengths = get_chromosome_lengths(ref_map_folder)
    if not chr_lengths:
        print("Warning: No chromosome data found for coverage visualization.")
        return None

    coverage = collect_coverage(data_folder)
    sorted_chrs = sorted(chr_lengths.keys(), key=chromosome_sort_key)
    max_length = max(chr_lengths.values())

    # Group coverage by chromosome
    chr_coverage = {}
    for c in coverage:
        ref = c['ref_name']
        if ref not in chr_coverage:
            chr_coverage[ref] = []
        chr_coverage[ref].append(c)

    total_molecules = len(coverage)
    covered_chrs = len(chr_coverage)
    total_data_bp = sum(c['end_bp'] - c['start_bp'] for c in coverage)
    genome_bp = sum(chr_lengths.values())

    # Build chromosome bars
    bars_html = ""
    for chr_name in sorted_chrs:
        chr_len = chr_lengths[chr_name]
        bar_width_pct = (chr_len / max_length) * 100

        segments_html = ""
        count = 0
        if chr_name in chr_coverage:
            count = len(chr_coverage[chr_name])
            for cov in chr_coverage[chr_name]:
                start_pct = (cov['start_bp'] / chr_len) * 100
                width_pct = ((cov['end_bp'] - cov['start_bp']) / chr_len) * 100
                if start_pct < 0:
                    start_pct = 0
                if start_pct + width_pct > 100:
                    width_pct = 100 - start_pct
                log_p = cov.get('log_p', 0)
                if log_p > -4:
                    color_cls = "cov-red"
                elif log_p > -10:
                    color_cls = "cov-green"
                else:
                    color_cls = "cov-black"
                tooltip = f"{cov['tif_name']}: {cov['start_bp']:,}-{cov['end_bp']:,} bp, log(P)={log_p:.1f}"
                link = cov['html_link']
                segments_html += (
                    f'<a href="{link}" class="cov-segment {color_cls}" '
                    f'style="left:{start_pct:.2f}%;width:{width_pct:.2f}%" '
                    f'title="{tooltip}"></a>\n'
                )

        chr_len_mbp = chr_len / 1_000_000

        bars_html += f"""
        <div class="chr-row">
            <div class="chr-label">{chr_name}</div>
            <div class="chr-count">{count if count > 0 else ''}</div>
            <div class="chr-bar-wrapper">
                <div class="chr-bar" style="width:{bar_width_pct:.1f}%">
                    {segments_html}
                </div>
            </div>
            <div class="chr-length">{chr_len_mbp:.1f} Mbp</div>
        </div>"""

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>Chromosome Coverage</title>
    <style>
        body {{ font-family: 'Segoe UI', sans-serif; line-height: 1.6; color: #333; margin: 0 auto; padding: 40px 20px; background-color: #f9f9f9; }}
        .container {{ background-color: #fff; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        .stats {{ margin: 15px 0 25px 0; color: #555; }}
        .chr-row {{ display: flex; align-items: center; margin: 3px 0; }}
        .chr-label {{ width: 80px; min-width: 80px; font-size: 0.85em; font-weight: bold; text-align: right; padding-right: 5px; white-space: nowrap; }}
        .chr-count {{ width: 30px; min-width: 30px; font-size: 0.8em; color: #666; text-align: right; padding-right: 8px; }}
        .chr-bar-wrapper {{ flex: 1; }}
        .chr-bar {{ height: 16px; background-color: #e0e0e0; border-radius: 3px; position: relative; overflow: hidden; border: 1px solid #ccc; }}
        .cov-segment {{ position: absolute; top: 0; height: 100%; display: block; }}
        .cov-red {{ background-color: #e74c3c; }}
        .cov-green {{ background-color: #27ae60; }}
        .cov-black {{ background-color: #000000; }}
        .cov-segment:hover {{ opacity: 0.5; z-index: 1; }}
        .chr-length {{ font-size: 0.75em; color: #888; padding-left: 8px; white-space: nowrap; width: 70px; }}
        .back-link {{ margin-bottom: 15px; }}
        .back-link a {{ color: #3498db; text-decoration: none; }}
        .back-link a:hover {{ text-decoration: underline; }}
        footer {{ margin-top: 50px; font-size: 0.8em; color: #888; text-align: center; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Chromosome Coverage</h1>
        <div class="back-link"><a href="index.html">Back to Summary Index</a></div>
        <div class="stats">{total_molecules} molecules mapped to {covered_chrs} chromosomes, {total_data_bp / 1e9:.2f} Gb / {genome_bp / 1e9:.2f} Gb ({total_data_bp / genome_bp:.2f}x)</div>
        {bars_html}
    </div>
    <footer>Generated by hDOM_chr_cover.py</footer>
</body>
</html>"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created: {output_path}")
    return output_path
