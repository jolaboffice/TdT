#!/usr/bin/env python3
"""
hDOM_evaluate.py
Generates an evaluation HTML report based on hDOM analysis results.
Reuses helpers from hDOM_tif2html.py.

Differences from hDOM_tif2html.py:
- Shows rank 1, 2 and rank 11, 12 (with thin dark-gray separator between them)
- Provides a user-selectable rank slot at the bottom (starts black; user
  enters a rank number to display that molecule)

Output: {base}_eval.html next to each input .tif, and index_eval.html
"""

import os
import sys
import argparse
import json
import pandas as pd

# Ensure this script's directory is on sys.path so sibling modules
# (hDOM_tif2html, DOM_lib, hDOM_chr_cover) import regardless of CWD.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import hDOM_tif2html as t2h
from hDOM_tif2html import (
    BPP,
    process_query_image,
    process_ref_match,
    get_table_html,
    get_html_style,
    find_ref_map_folder,
)
from DOM_lib import tif2RGimage
from hDOM_chr_cover import get_chromosome_lengths, chromosome_sort_key

# Ranks (1-based, by row index in the CSV/table) to show before the selectable slot
FIXED_RANKS = [1, 2, 11, 12]
# Where to insert separators (after these positions within FIXED_RANKS)
SEPARATOR_AFTER = {2}  # after rank 2; another separator goes before the selectable slot


def _extra_style():
    """Additional CSS on top of get_html_style()."""
    return """
    .rank-separator { height: 1px; background-color: #444; margin: 4px 0; }
    .selector-box { margin: 10px 0 0 0; padding: 10px; background-color: #222;
                    border-radius: 4px; color: #ddd; font-family: 'Segoe UI', sans-serif; }
    .selector-box label { margin-right: 8px; font-weight: bold; }
    .selector-box input { width: 60px; padding: 4px 6px; font-size: 0.95em;
                          border: 1px solid #555; background-color: #111; color: #eee;
                          border-radius: 3px; }
    .selector-box .hint { margin-left: 10px; color: #888; font-size: 0.85em; }
    .placeholder-black { width: 100%; height: 30px; background-color: #000; }
    .eval-toolbar { margin: 10px 0; padding: 10px; background-color: #222;
                    border-radius: 4px; color: #ddd; font-family: 'Segoe UI', sans-serif;
                    display: flex; align-items: center; gap: 14px; flex-wrap: wrap; }
    .eval-toolbar .group-label { font-weight: bold; color: #ccc; }
    .class-btn { width: 32px; height: 32px; padding: 0; border: 1px solid #555;
                 background-color: #111; color: #eee; border-radius: 3px; cursor: pointer;
                 font-weight: bold; font-size: 0.95em; display: inline-flex;
                 align-items: center; justify-content: center; }
    .class-btn:hover { border-color: #888; }
    .class-btn.active[data-value="S"],
    .class-btn.active[data-value="A"],
    .class-btn.active[data-value="B"],
    .class-btn.active[data-value="M"],
    .class-btn.active[data-value="C"] { background-color: #00b050; border-color: #00b050; color: #fff; }
    .class-btn.active[data-value="D"] { background-color: #ffc000; border-color: #ffc000; color: #fff; }
    .class-btn.active[data-value="F"] { background-color: #e74c3c; border-color: #e74c3c; color: #fff; }
    .star-btn { width: 32px; height: 32px; padding: 0; border: 1px solid #555;
                background-color: #111; color: #eee; border-radius: 3px; cursor: pointer;
                font-size: 1.1em; display: inline-flex; align-items: center;
                justify-content: center; }
    .star-btn.active { background-color: #ffd400; border-color: #ffd400; color: #fff; }
    .chosen-info { color: #ff0000; font-weight: bold; margin-left: 8px; }
    table.matching-table tbody tr { cursor: pointer; }
    table.matching-table tbody tr.chosen-row { background-color: #a20000 !important; }
    table.matching-table tbody tr.chosen-row td { color: #fff; font-weight: bold; }
    """


def _row_to_match_info(row):
    """Convert a DataFrame row to the dict expected by process_ref_match."""
    aligned_bp = row.get('aligned_bp')
    if pd.isna(aligned_bp):
        aligned_bp = row.get('shift_bp', 0)
    scale = row.get('scale', 1.0)
    if pd.isna(scale):
        scale = 1.0
    return {
        'ref_name': row.get('ref_name'),
        'aligned_bp': aligned_bp,
        'scale': scale,
        'dir': row.get('dir', 1),
    }


def _render_rank_item(label, b64):
    """Render a comparison-item row with a given label and base64 image."""
    if b64 is None:
        # Fallback to black placeholder if image cannot be produced
        return (f'<div class="comparison-item"><div class="image-label">{label}</div>'
                f'<div class="image-container"><div class="placeholder-black"></div></div></div>')
    return (f'<div class="comparison-item"><div class="image-label">{label}</div>'
            f'<div class="image-container">'
            f'<img src="data:image/png;base64,{b64}"></div></div>')


# Rank labels render as "{prefix}-rank{n}" (pcc-rank1, mal-rank1, KU-rank1).
METHOD_PREFIX = {'pcc': 'pcc', 'mal': 'mal', 'logical': 'KU'}


def _slot_label(df, r):
    """Label for overall rank r: '{pcc|mal|KU}-rank{per-method-rank}'.

    Uses the df's 'method'/'rank' columns when available; otherwise falls back
    to the positional convention (1-10 = pcc, 11+ = mal)."""
    try:
        row = df.iloc[r - 1]
        method = str(row.get('method', '')).lower()
        rk = int(row.get('rank'))
        if method and not pd.isna(row.get('rank')):
            return f"{METHOD_PREFIX.get(method, method)}-rank{rk}"
    except Exception:
        pass
    return ('pcc-rank%d' % r) if r <= 10 else ('mal-rank%d' % (r - 10))


def tif2eval_html(tif_path, ref_map_override=None):
    """Generate an _eval.html report for a single tif. Returns summary dict or None."""
    ext = os.path.splitext(tif_path)[1].lower()
    if ext not in ('.tif', '.tiff'):
        print(f"Error: Input must be a .tif/.tiff file. Got: {tif_path}")
        return None

    base_path = os.path.splitext(tif_path)[0]
    info_path = f"{base_path}_info.txt"
    csv_path = f"{base_path}_matching_results.csv"
    html_path = f"{base_path}_eval.html"

    # Info
    info_content = f"Information file not found: {info_path}"
    if os.path.exists(info_path):
        with open(info_path, 'r', encoding='utf-8') as f:
            info_content = f.read()

    # Query image
    query_img, query_b64 = process_query_image(tif_path)
    if not query_img:
        return None
    query_final_w = query_img.size[0]
    molecule_length_bp = int((query_final_w - 20) * BPP)

    # Results
    df = pd.read_csv(csv_path) if os.path.exists(csv_path) else None

    summary_info = {
        'tif_name': os.path.basename(tif_path),
        'html_link': os.path.basename(html_path),
        'ref_map': 'N/A',
        'aligned_bp': 'N/A',
        'length_bp': molecule_length_bp,
        'direction': '+',
        'log_p': 'N/A',
    }

    tif_dir = os.path.dirname(os.path.abspath(tif_path))
    ref_map_folder = ref_map_override if ref_map_override else find_ref_map_folder(tif_dir)

    # Raw tif width (matches hDOM_chr_cover.py: coverage-length = raw_w * BPP * scale)
    try:
        _, raw_w, _ = tif2RGimage(tif_path)
    except Exception:
        raw_w = max(query_final_w - 20, 1)

    # Pre-render all available ref images so the selectable slot has data
    all_ref_b64 = {}  # rank (1-based row idx) -> base64 string
    ranks_data = []   # per-rank info (index 0 = rank 1) for index page overrides
    if df is not None and not df.empty:
        # summary comes from top hit
        top_hit = df.iloc[0]
        summary_info['ref_map'] = top_hit.get('ref_name', 'N/A')
        summary_info['log_p'] = top_hit.get('log(P)', 'N/A')
        aligned = top_hit.get('aligned_bp')
        if pd.isna(aligned):
            aligned = top_hit.get('shift_bp', 0)
        summary_info['aligned_bp'] = int(aligned)
        summary_info['direction'] = '+' if top_hit.get('dir', 1) == 1 else '-'

        for i in range(len(df)):
            row = df.iloc[i]
            ref_name = row.get('ref_name')
            # Collect rank data for every row (used by index override logic)
            r_aligned = row.get('aligned_bp')
            if pd.isna(r_aligned):
                r_aligned = row.get('shift_bp', 0)
            r_logp = row.get('log(P)')
            r_scale = row.get('scale', 1.0)
            if pd.isna(r_scale):
                r_scale = 1.0
            r_scale = float(r_scale)
            # unified score: pcc -> log(P); mal rows have no log(P) -> use m_score (misc)
            r_misc = row.get('misc')
            r_method = row.get('method')
            if not pd.isna(r_logp):
                r_score = float(r_logp)
            elif r_method == 'mal' and not pd.isna(r_misc):
                r_score = float(r_misc)
            else:
                r_score = None
            ranks_data.append({
                'ref_name': ('' if pd.isna(ref_name) else str(ref_name)) or 'N/A',
                'aligned_bp': int(r_aligned) if not pd.isna(r_aligned) else 0,
                'log_p': float(r_logp) if not pd.isna(r_logp) else None,
                'dir': '+' if row.get('dir', 1) == 1 else '-',
                'scale': r_scale,
                'cov_length_bp': int(raw_w * BPP * r_scale),
                'score': round(r_score, 2) if r_score is not None else None,
                'label': _slot_label(df, i + 1),
            })
            if pd.isna(ref_name) or not ref_name:
                continue
            b64 = process_ref_match(_row_to_match_info(row), ref_map_folder, query_final_w)
            if b64:
                all_ref_b64[i + 1] = b64  # 1-based rank = row index + 1
    summary_info['ranks'] = ranks_data
    summary_info['ref_map_folder'] = ref_map_folder

    # Build comparison block
    # data row
    comp_items = [_render_rank_item('data', query_b64)]

    for idx, r in enumerate(FIXED_RANKS):
        b64 = all_ref_b64.get(r)
        comp_items.append(_render_rank_item(_slot_label(df, r), b64))
        if r in SEPARATOR_AFTER:
            comp_items.append('<div class="rank-separator"></div>')

    # Second separator before the selectable slot
    comp_items.append('<div class="rank-separator"></div>')

    # Selectable slot: starts black, swapped via JS
    comp_items.append(
        '<div class="comparison-item" id="sel-item">'
        '<div class="image-label"><span id="sel-label">rank ?</span></div>'
        '<div class="image-container">'
        '<img id="sel-img" style="display:none;">'
        '<div id="sel-placeholder" class="placeholder-black"></div>'
        '</div></div>'
    )

    comp_html = ('<div class="comparison-box"><div class="comparison-images">'
                 + ''.join(comp_items) + '</div>'
                 + '<div class="selector-box">'
                 + '<label for="sel-input">Select rank:</label>'
                 + f'<input type="number" id="sel-input" min="1" max="{max(len(df) if df is not None else 0,1)}" step="1" placeholder="e.g. 7">'
                 + '<span class="hint">Enter a rank number (row in the table) to display that molecule.</span>'
                 + '</div></div>')

    # JS: store rank->base64 map and swap image on input
    # To keep HTML size reasonable, the map contains only ranks that aren't already
    # shown above (fixed ranks); selecting a fixed rank still works by duplicating in map.
    js_data = json.dumps(all_ref_b64)
    js_labels = json.dumps({str(r): _slot_label(df, r) for r in range(1, len(df) + 1)}) if df is not None else '{}'
    tif_name_js = json.dumps(os.path.basename(tif_path))
    selector_js = f"""
    <script>
    const RANK_IMAGES = {js_data};
    const RANK_LABELS = {js_labels};
    const TIF_NAME = {tif_name_js};
    const STORAGE_CLASS = 'hdom_eval_class_' + TIF_NAME;
    const STORAGE_STAR = 'hdom_eval_star_' + TIF_NAME;
    const STORAGE_CHOSEN = 'hdom_eval_chosen_' + TIF_NAME;
    document.addEventListener('DOMContentLoaded', function() {{
        // --- rank preview selector ---
        const input = document.getElementById('sel-input');
        const img = document.getElementById('sel-img');
        const ph = document.getElementById('sel-placeholder');
        const label = document.getElementById('sel-label');
        function update() {{
            const v = parseInt(input.value, 10);
            if (!isNaN(v) && RANK_IMAGES[v]) {{
                img.src = 'data:image/png;base64,' + RANK_IMAGES[v];
                img.style.display = 'block';
                ph.style.display = 'none';
                label.textContent = RANK_LABELS[String(v)] || ('rank' + v);
            }} else {{
                img.style.display = 'none';
                ph.style.display = 'block';
                label.textContent = 'rank ?';
            }}
        }}
        input.addEventListener('input', update);
        input.addEventListener('change', update);

        // --- class buttons ---
        function refreshClassButtons() {{
            const cur = localStorage.getItem(STORAGE_CLASS) || '';
            document.querySelectorAll('.class-btn').forEach(function(b) {{
                if (b.dataset.value === cur) b.classList.add('active');
                else b.classList.remove('active');
            }});
        }}
        document.querySelectorAll('.class-btn').forEach(function(btn) {{
            btn.addEventListener('click', function() {{
                const v = btn.dataset.value;
                const cur = localStorage.getItem(STORAGE_CLASS) || '';
                if (cur === v) localStorage.removeItem(STORAGE_CLASS);
                else localStorage.setItem(STORAGE_CLASS, v);
                refreshClassButtons();
            }});
        }});
        refreshClassButtons();

        // --- star toggle ---
        const starBtn = document.getElementById('star-btn');
        function refreshStar() {{
            const v = localStorage.getItem(STORAGE_STAR) === '1';
            if (v) starBtn.classList.add('active'); else starBtn.classList.remove('active');
            starBtn.textContent = v ? '\u2605' : '\u2606';
        }}
        starBtn.addEventListener('click', function() {{
            const cur = localStorage.getItem(STORAGE_STAR) === '1';
            if (cur) localStorage.removeItem(STORAGE_STAR);
            else localStorage.setItem(STORAGE_STAR, '1');
            refreshStar();
        }});
        refreshStar();

        // --- matching-table row click = choose correct rank (override rank1) ---
        const rows = document.querySelectorAll('table.matching-table tbody tr');
        const chosenInfo = document.getElementById('chosen-info');
        function refreshChosen() {{
            const ch = parseInt(localStorage.getItem(STORAGE_CHOSEN), 10);
            rows.forEach(function(tr, i) {{
                if ((i + 1) === ch) tr.classList.add('chosen-row');
                else tr.classList.remove('chosen-row');
            }});
            if (!isNaN(ch)) {{
                var _lbl = (typeof RANK_LABELS !== 'undefined' && RANK_LABELS[String(ch)]) ? RANK_LABELS[String(ch)] : ('rank' + ch);
                chosenInfo.textContent = '(override \u2192 ' + _lbl + ')';
            }} else chosenInfo.textContent = '';
        }}
        rows.forEach(function(tr, i) {{
            tr.querySelectorAll('a').forEach(function(a) {{
                a.addEventListener('click', function(e) {{ e.stopPropagation(); }});
            }});
            tr.addEventListener('click', function() {{
                const rank = i + 1;
                const cur = parseInt(localStorage.getItem(STORAGE_CHOSEN), 10);
                if (cur === rank) localStorage.removeItem(STORAGE_CHOSEN);
                else localStorage.setItem(STORAGE_CHOSEN, String(rank));
                refreshChosen();
            }});
        }});
        refreshChosen();
    }});
    </script>
    """

    # Toolbar (class selector + star toggle + override indicator)
    classes = ['S', 'A', 'B', 'M', 'C', 'D', 'F']
    class_buttons = ''.join(
        f'<button type="button" class="class-btn" data-value="{c}">{c}</button>' for c in classes
    )
    toolbar_html = (
        '<div class="eval-toolbar">'
        '<span class="group-label">Class:</span>'
        + class_buttons +
        '<span class="group-label" style="margin-left:10px;">Mark:</span>'
        '<button type="button" id="star-btn" class="star-btn">\u2606</button>'
        '<span id="chosen-info" class="chosen-info"></span>'
        '<span class="hint" style="color:#888;font-size:0.85em;margin-left:auto;">'
        'Click a row in Matching Results to override rank1.</span>'
        '</div>'
    )

    # Matching table (same as original)
    tif_abs = os.path.abspath(tif_path)
    folder_name = os.path.basename(os.path.dirname(tif_abs))
    file_name = os.path.basename(tif_abs)
    table_html = get_table_html(df, folder_name=folder_name, file_name=file_name,
                               ref_map_path=ref_map_folder)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>hDOM Evaluate - {os.path.basename(tif_path)}</title>
    <style>{get_html_style()}{_extra_style()}</style>
</head>
<body>
    <div class="container">
        <h1>hDOM Evaluation Report</h1>
        <p><strong>File:</strong> {os.path.basename(tif_path)} ({molecule_length_bp / 1000:.1f} kbp)</p>
        {toolbar_html}
        {comp_html}
        <h2>Matching Results</h2>
        <div class="table-container">{table_html}</div>
        <h2>Information</h2>
        <div class="info-box">{info_content}</div>
    </div>
    {selector_js}
    <footer>Generated by hDOM_evaluate.py</footer>
</body>
</html>"""

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created: {html_path}")
    return summary_info


def generate_eval_index_html(summaries, output_path,
                             html_suffix='_eval.html',
                             page_title='hDOM Evaluation Summary Index',
                             generator='hDOM_evaluate.py'):
    """Summary index page for per-molecule reports (reuses index style/sort).
    Adds Star + Class columns and supports per-molecule rank override via localStorage.

    html_suffix: suffix of the per-molecule HTML files (e.g. '_eval.html', '_view.html').
    """
    rows_html = ""
    rank_data = {}
    for idx, s in enumerate(summaries, start=1):
        tif_name = s['tif_name']
        rank_data[tif_name] = s.get('ranks', [])
        ref_map = s['ref_map']
        aligned_bp = s['aligned_bp']
        length_bp = s['length_bp']
        log_p = s.get('score', s['log_p'])
        direction = s['direction']
        aligned_bp_fmt = f"{aligned_bp:,}" if isinstance(aligned_bp, (int, float)) else aligned_bp
        length_bp_fmt = f"{length_bp:,}" if isinstance(length_bp, (int, float)) else length_bp
        log_p_fmt = f"{log_p:.2f}" if isinstance(log_p, (int, float)) else log_p
        ref_map_low = str(ref_map).lower()
        rows_html += f"""
        <tr data-tif="{tif_name}">
            <td class="star-cell text-center" data-value="0">&#9734;</td>
            <td class="row-num text-right">{idx}</td>
            <td data-value="{tif_name.lower()}"><a href="{s['html_link']}">{tif_name}</a></td>
            <td class="ref-cell" data-value="{ref_map_low}">{ref_map}</td>
            <td class="text-right aligned-cell" data-value="{aligned_bp}">{aligned_bp_fmt}</td>
            <td class="text-right length-cell" data-value="{length_bp}">{length_bp_fmt}</td>
            <td class="text-right logp-cell" data-value="{log_p}">{log_p_fmt}</td>
            <td class="text-center dir-cell" data-value="{direction}">{direction}</td>
            <td class="text-center class-cell" data-value="">&nbsp;</td>
        </tr>"""

    # RANK_DATA is kept as a separate file (index_rank_data.json) and fetched
    # lazily so the HTML itself stays small and the browser can render the
    # table (and class/star cells) without waiting for the full JSON to parse.
    rank_data_path = os.path.join(os.path.dirname(os.path.abspath(output_path)),
                                  'index_rank_data.json')
    with open(rank_data_path, 'w', encoding='utf-8') as _f:
        json.dump(rank_data, _f, separators=(',', ':'))
    rank_data_js = '{}'  # placeholder; populated via fetch in the browser
    html_suffix_js = json.dumps(html_suffix)

    # Chromosome lengths (from first available ref_map_folder)
    chr_lengths = {}
    for s in summaries:
        rmf = s.get('ref_map_folder')
        if rmf and os.path.isdir(rmf):
            chr_lengths = get_chromosome_lengths(rmf)
            if chr_lengths:
                break
    sorted_chrs = sorted(chr_lengths.keys(), key=chromosome_sort_key) if chr_lengths else []
    max_chr_len = max(chr_lengths.values()) if chr_lengths else 0
    chr_lengths_js = json.dumps(chr_lengths)
    sorted_chrs_js = json.dumps(sorted_chrs)

    # Pre-build empty chr rows (filled by JS)
    chr_rows_html = ""
    for chr_name in sorted_chrs:
        chr_len = chr_lengths[chr_name]
        bar_width_pct = (chr_len / max_chr_len) * 100 if max_chr_len else 0
        chr_len_mbp = chr_len / 1_000_000
        chr_rows_html += f"""
        <div class="chr-row">
            <div class="chr-label">{chr_name}</div>
            <div class="chr-count" data-chr="{chr_name}"></div>
            <div class="chr-bar-wrapper">
                <div class="chr-bar" data-chr="{chr_name}" style="width:{bar_width_pct:.1f}%"></div>
            </div>
            <div class="chr-length">{chr_len_mbp:.1f} Mbp</div>
        </div>"""

    extra_css = """
    body { max-width: 1400px !important; }
    .summary-table { table-layout: auto; }
    .summary-table td:nth-child(3), .summary-table th:nth-child(3) {
        max-width: 300px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
    }
    .summary-table td:nth-child(3) a { display: inline-block; max-width: 100%;
        overflow: hidden; text-overflow: ellipsis; vertical-align: middle; }
    tr.starred { background-color: #ffd400 !important; }
    tr.starred td { color: #000 !important; }
    .star-cell { cursor: pointer; font-size: 1.1em; user-select: none; }
    .class-cell { cursor: pointer; font-weight: bold; user-select: none; min-width: 28px; }
    .class-cell[data-value="S"],
    .class-cell[data-value="A"],
    .class-cell[data-value="B"],
    .class-cell[data-value="M"],
    .class-cell[data-value="C"] { color: #00b050; }
    .class-cell[data-value="D"] { color: #ffc000; }
    .class-cell[data-value="F"] { color: #e74c3c; }
    td.overridden { color: #ff0000 !important; font-weight: bold; }
    .chr-cover-section { margin-top: 30px; padding: 20px; background-color: #fff;
                         border-radius: 8px; border: 1px solid #ddd; color: #333; }
    .chr-cover-section h2 { color: #2c3e50; border-bottom: 2px solid #3498db;
                            padding-bottom: 8px; margin-top: 0; }
    .chr-cover-stats { margin: 10px 0 20px 0; color: #555; font-size: 0.9em; }
    .chr-row { display: flex; align-items: center; margin: 3px 0; }
    .chr-label { width: 80px; min-width: 80px; font-size: 0.85em; font-weight: bold;
                 text-align: right; padding-right: 5px; white-space: nowrap; color: #333; }
    .chr-count { width: 30px; min-width: 30px; font-size: 0.8em; color: #666;
                 text-align: right; padding-right: 8px; }
    .chr-bar-wrapper { flex: 1; }
    .chr-bar { height: 16px; background-color: #e0e0e0; border-radius: 3px;
               position: relative; overflow: hidden; border: 1px solid #ccc; }
    .cov-segment { position: absolute; top: 0; height: 100%; display: block; }
    .cov-red { background-color: #e74c3c; }
    .cov-green { background-color: #27ae60; }
    .cov-black { background-color: #000000; }
    .cov-segment:hover { opacity: 0.5; z-index: 1; }
    .chr-length { font-size: 0.75em; color: #888; padding-left: 8px;
                  white-space: nowrap; width: 70px; }
    .export-btn { display: inline-block; padding: 8px 18px; background-color: #2980b9;
                  color: #fff; border: none; border-radius: 4px; cursor: pointer;
                  font-size: 0.95em; font-weight: bold; margin-right: 10px; }
    .export-btn:hover { background-color: #3498db; }
    """

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>{page_title}</title>
    <style>{get_html_style()}{extra_css}</style>
    <script src="https://cdn.sheetjs.com/xlsx-0.20.3/package/dist/xlsx.full.min.js"></script>
    <script>
    // RANK_DATA is lazy-loaded from index_rank_data.json (see fetch below).
    let RANK_DATA = {rank_data_js};
    let RANK_DATA_READY = false;
    const CHR_LENGTHS = {chr_lengths_js};
    const SORTED_CHRS = {sorted_chrs_js};
    const HTML_SUFFIX = {html_suffix_js};
    const CLASSES = ['', 'A', 'B', 'M', 'C', 'D', 'F'];

    function getEffectiveRank(tif) {{
        const raw = localStorage.getItem('hdom_eval_chosen_' + tif);
        const r = parseInt(raw, 10);
        const ranks = RANK_DATA[tif] || [];
        if (!isNaN(r) && r >= 1 && r <= ranks.length) return ranks[r - 1];
        return ranks[0] || null;
    }}

    function rebuildCoverage() {{
        if (!SORTED_CHRS.length) return;
        // collect per-chr coverage from effective rank of each tif
        const byChr = {{}};
        SORTED_CHRS.forEach(c => {{ byChr[c] = []; }});
        Object.keys(RANK_DATA).forEach(function(tif) {{
            const d = getEffectiveRank(tif);
            if (!d || !d.ref_name || byChr[d.ref_name] === undefined) return;
            const start = d.aligned_bp;
            const end = start + (d.cov_length_bp || 0);
            byChr[d.ref_name].push({{
                start_bp: start, end_bp: end,
                tif: tif, log_p: d.log_p, html_link: tif.replace(/\\.tiff?$/i, HTML_SUFFIX)
            }});
        }});
        // clear then fill bars
        document.querySelectorAll('.chr-bar').forEach(function(bar) {{
            bar.innerHTML = '';
            const chr = bar.getAttribute('data-chr');
            const chrLen = CHR_LENGTHS[chr];
            const segs = byChr[chr] || [];
            document.querySelector('.chr-count[data-chr="' + chr + '"]').textContent = segs.length ? segs.length : '';
            segs.forEach(function(cov) {{
                let startPct = (cov.start_bp / chrLen) * 100;
                let widthPct = ((cov.end_bp - cov.start_bp) / chrLen) * 100;
                if (startPct < 0) {{ widthPct += startPct; startPct = 0; }}
                if (startPct + widthPct > 100) widthPct = 100 - startPct;
                if (widthPct <= 0) return;
                const molClass = localStorage.getItem('hdom_eval_class_' + cov.tif) || '';
                let cls = 'cov-black';
                if (molClass === 'S' || molClass === 'A' || molClass === 'B' || molClass === 'M' || molClass === 'C') cls = 'cov-green';
                else if (molClass === 'D' || molClass === 'F') cls = 'cov-red';
                const a = document.createElement('a');
                a.href = cov.html_link;  // suffix-aware, set below
                a.className = 'cov-segment ' + cls;
                a.style.left = startPct.toFixed(2) + '%';
                a.style.width = widthPct.toFixed(2) + '%';
                a.title = cov.tif + ': ' + cov.start_bp.toLocaleString() + '-' + cov.end_bp.toLocaleString() + ' bp, log(P)=' + (cov.log_p != null ? cov.log_p.toFixed(1) : 'N/A');
                bar.appendChild(a);
            }});
        }});
        // stats
        const total = Object.values(byChr).reduce((n, a) => n + a.length, 0);
        const covered = Object.values(byChr).filter(a => a.length > 0).length;
        const totalBp = Object.values(byChr).reduce((n, a) => n + a.reduce((m, c) => m + (c.end_bp - c.start_bp), 0), 0);
        const genomeBp = Object.values(CHR_LENGTHS).reduce((n, v) => n + v, 0);
        const statsEl = document.getElementById('chr-cover-stats');
        if (statsEl && genomeBp > 0) {{
            statsEl.textContent = total + ' molecules mapped to ' + covered + ' chromosomes, ' +
                (totalBp / 1e9).toFixed(2) + ' Gb / ' + (genomeBp / 1e9).toFixed(2) + ' Gb (' +
                (totalBp / genomeBp).toFixed(2) + 'x)';
        }}
    }}

    function applyOverride(tr) {{
        const tif = tr.getAttribute('data-tif');
        const key = 'hdom_eval_chosen_' + tif;
        const raw = localStorage.getItem(key);
        const ranks = RANK_DATA[tif] || [];
        const refCell = tr.querySelector('.ref-cell');
        const alignedCell = tr.querySelector('.aligned-cell');
        const lengthCell = tr.querySelector('.length-cell');
        const logpCell = tr.querySelector('.logp-cell');
        const dirCell = tr.querySelector('.dir-cell');
        const cells = [refCell, alignedCell, lengthCell, logpCell, dirCell];
        // Reset first (restore original rank-1 values)
        const origRef = ranks[0] ? ranks[0].ref_name : refCell.textContent;
        const origAligned = ranks[0] ? ranks[0].aligned_bp : alignedCell.getAttribute('data-value');
        const origLength = (ranks[0] && ranks[0].cov_length_bp != null) ? ranks[0].cov_length_bp : lengthCell.getAttribute('data-value');
        const origLogp = (ranks[0] && ranks[0].score != null) ? ranks[0].score : logpCell.getAttribute('data-value');
        const origDir = ranks[0] ? ranks[0].dir : dirCell.textContent;
        cells.forEach(c => {{ c.classList.remove('overridden'); c.removeAttribute('data-rank'); }});
        refCell.textContent = origRef;
        refCell.setAttribute('data-value', String(origRef).toLowerCase());
        alignedCell.textContent = Number(origAligned).toLocaleString();
        alignedCell.setAttribute('data-value', origAligned);
        lengthCell.textContent = Number(origLength).toLocaleString();
        lengthCell.setAttribute('data-value', origLength);
        logpCell.textContent = (typeof origLogp === 'number') ? origLogp.toFixed(2) : origLogp;
        logpCell.setAttribute('data-value', origLogp);
        dirCell.textContent = origDir;
        dirCell.setAttribute('data-value', origDir);
        const r = parseInt(raw, 10);
        if (!isNaN(r) && r >= 1 && r <= ranks.length && r !== 1) {{
            const d = ranks[r - 1];
            refCell.textContent = d.ref_name;
            refCell.setAttribute('data-value', String(d.ref_name).toLowerCase());
            alignedCell.textContent = Number(d.aligned_bp).toLocaleString();
            alignedCell.setAttribute('data-value', d.aligned_bp);
            lengthCell.textContent = Number(d.cov_length_bp).toLocaleString();
            lengthCell.setAttribute('data-value', d.cov_length_bp);
            logpCell.textContent = (d.score != null) ? Number(d.score).toFixed(2) : 'N/A';
            logpCell.setAttribute('data-value', d.score != null ? d.score : '');
            dirCell.textContent = d.dir;
            dirCell.setAttribute('data-value', d.dir);
            cells.forEach(c => {{ c.classList.add('overridden'); c.setAttribute('data-rank', r); }});
        }}
    }}

    function applyStar(tr) {{
        const tif = tr.getAttribute('data-tif');
        const key = 'hdom_eval_star_' + tif;
        const starred = localStorage.getItem(key) === '1';
        const cell = tr.querySelector('.star-cell');
        if (starred) {{
            tr.classList.add('starred');
            cell.innerHTML = '&#9733;';
            cell.setAttribute('data-value', '1');
        }} else {{
            tr.classList.remove('starred');
            cell.innerHTML = '&#9734;';
            cell.setAttribute('data-value', '0');
        }}
    }}

    function applyClass(tr) {{
        const tif = tr.getAttribute('data-tif');
        const key = 'hdom_eval_class_' + tif;
        const v = localStorage.getItem(key) || '';
        const cell = tr.querySelector('.class-cell');
        cell.textContent = v || '\u00A0';
        cell.setAttribute('data-value', v);
        if (v) cell.classList.add('has-class'); else cell.classList.remove('has-class');
    }}

    document.addEventListener('DOMContentLoaded', function() {{
        const rows = document.querySelectorAll('#summaryTable tbody tr');
        // Fast path: class + star rely only on localStorage, so we can render
        // them immediately without waiting for the RANK_DATA fetch.
        rows.forEach(function(tr) {{
            applyStar(tr);
            applyClass(tr);
            tr.querySelector('.star-cell').addEventListener('click', function(e) {{
                e.stopPropagation();
                const tif = tr.getAttribute('data-tif');
                const key = 'hdom_eval_star_' + tif;
                if (localStorage.getItem(key) === '1') localStorage.removeItem(key);
                else localStorage.setItem(key, '1');
                applyStar(tr);
            }});
            tr.querySelector('.class-cell').addEventListener('click', function(e) {{
                e.stopPropagation();
                const tif = tr.getAttribute('data-tif');
                const key = 'hdom_eval_class_' + tif;
                const cur = localStorage.getItem(key) || '';
                const idx = CLASSES.indexOf(cur);
                const next = CLASSES[(idx + 1) % CLASSES.length];
                if (next === '') localStorage.removeItem(key);
                else localStorage.setItem(key, next);
                applyClass(tr);
                if (RANK_DATA_READY) rebuildCoverage();
            }});
        }});
        // Lazy-load RANK_DATA; once it arrives, apply overrides + coverage.
        fetch('index_rank_data.json').then(function(r) {{ return r.json(); }}).then(function(data) {{
            RANK_DATA = data;
            RANK_DATA_READY = true;
            rows.forEach(applyOverride);
            rebuildCoverage();
        }}).catch(function(e) {{ console.error('RANK_DATA load failed:', e); }});
    }});

    function exportClassifications() {{
        const result = {{}};
        for (let i = 0; i < localStorage.length; i++) {{
            const key = localStorage.key(i);
            if (key.startsWith('hdom_eval_class_') ||
                key.startsWith('hdom_eval_chosen_') ||
                key.startsWith('hdom_eval_star_')) {{
                result[key] = localStorage.getItem(key);
            }}
        }}
        const blob = new Blob([JSON.stringify(result, null, 2)], {{type: 'application/json'}});
        const a = document.createElement('a');
        a.href = URL.createObjectURL(blob);
        a.download = 'hDOM_Classifications.json';
        a.click();
    }}

    function exportXlsx() {{
        const table = document.getElementById('summaryTable');
        const rows = table.querySelectorAll('tbody tr');
        const data = [];
        rows.forEach(function(tr) {{
            const tif = tr.getAttribute('data-tif');
            const starred = tr.querySelector('.star-cell').getAttribute('data-value') === '1' ? 'Y' : '';
            const refMap = tr.querySelector('.ref-cell').textContent;
            const alignedBp = Number(tr.querySelector('.aligned-cell').getAttribute('data-value'));
            const lengthBp = Number(tr.querySelector('.length-cell').getAttribute('data-value'));
            const logP = Number(tr.querySelector('.logp-cell').getAttribute('data-value'));
            const dir = tr.querySelector('.dir-cell').textContent;
            const cls = tr.querySelector('.class-cell').getAttribute('data-value') || '';
            const chosenRaw = localStorage.getItem('hdom_eval_chosen_' + tif);
            const chosenRank = chosenRaw ? parseInt(chosenRaw, 10) : 1;
            data.push({{
                'TIF File': tif,
                'Chosen Rank': chosenRank,
                'Reference Map': refMap,
                'Aligned Position (bp)': alignedBp,
                'Molecule Length (bp)': lengthBp,
                'score': logP,
                'Dir': dir,
                'Class': cls,
                'Star': starred
            }});
        }});
        const ws = XLSX.utils.json_to_sheet(data);
        const wb = XLSX.utils.book_new();
        XLSX.utils.book_append_sheet(wb, ws, 'Summary');
        XLSX.writeFile(wb, 'hDOM_Classifications.xlsx');
    }}

    function sortTable(n) {{
        const table = document.getElementById("summaryTable");
        const tbody = table.tBodies[0];
        const rows = Array.from(tbody.rows);
        const curCol = tbody.getAttribute('data-sort-col');
        const curDir = tbody.getAttribute('data-sort-dir') || 'desc';
        const dir = (curCol === String(n) && curDir === 'asc') ? 'desc' : 'asc';
        const mult = dir === 'asc' ? 1 : -1;
        const CLASS_ORDER = {{'S':0,'A':1,'B':2,'M':3,'C':4,'D':5,'F':6,'':7}};
        // Decorate once: read each cell's sort key a single time instead of on
        // every comparison. With thousands of rows this removes the O(n log n)
        // getAttribute() calls that made sorting slow.
        const coll = new Intl.Collator(undefined, {{sensitivity: 'base'}});
        const deco = rows.map(function(r) {{
            const raw = r.cells[n].getAttribute('data-value') || '';
            return {{ r: r, raw: raw, lower: raw.toLowerCase(),
                     num: parseFloat(raw),
                     cls: (raw in CLASS_ORDER) ? CLASS_ORDER[raw] : 99 }};
        }});
        deco.sort(function(a, b) {{
            if (n === 8) return (a.cls - b.cls) * mult;
            if (!isNaN(a.num) && !isNaN(b.num) && a.raw !== '' && b.raw !== '') return (a.num - b.num) * mult;
            return coll.compare(a.lower, b.lower) * mult;
        }});
        const frag = document.createDocumentFragment();
        deco.forEach(function(d) {{ frag.appendChild(d.r); }});
        tbody.appendChild(frag);
        tbody.setAttribute('data-sort-col', String(n));
        tbody.setAttribute('data-sort-dir', dir);
        for (let i = 0; i < deco.length; i++) {{
            const c = deco[i].r.cells[1];
            if (c) c.textContent = String(i + 1);
        }}
    }}
    </script>
</head>
<body>
    <div class="container">
        <h1>{page_title}</h1>
        <p style="color:#888;font-size:0.9em;">
            Click <b>&#9734;</b> to star a row (yellow). Click the Class cell to cycle A&rarr;B&rarr;M&rarr;C&rarr;D&rarr;E&rarr;F.
            Overridden values are shown in red.
        </p>
        <div style="margin: 10px 0;">
            <button class="export-btn" onclick="exportXlsx()">Export (.xlsx)</button>
            <button class="export-btn" onclick="exportClassifications()">Export Class</button>
        </div>
        <table class="summary-table" id="summaryTable">
            <thead>
                <tr>
                    <th class="text-center" onclick="sortTable(0)">&#9733;</th>
                    <th class="text-right">#</th>
                    <th onclick="sortTable(2)">TIF File</th>
                    <th onclick="sortTable(3)">Reference Map</th>
                    <th class="text-right" onclick="sortTable(4)">Aligned Position (bp)</th>
                    <th class="text-right" onclick="sortTable(5)">Molecule Length (bp)</th>
                    <th class="text-right" onclick="sortTable(6)">score</th>
                    <th class="text-center" onclick="sortTable(7)">Dir</th>
                    <th class="text-center" onclick="sortTable(8)">Class</th>
                </tr>
            </thead>
            <tbody>
                {rows_html}
            </tbody>
        </table>
        <div class="chr-cover-section">
            <h2>Chromosome Coverage (override-aware)</h2>
            <div id="chr-cover-stats" class="chr-cover-stats"></div>
            {chr_rows_html}
        </div>
    </div>
    <footer>Generated by {generator}</footer>
</body>
</html>"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created Summary Index: {output_path}")


def _collect_tifs_from_input(path):
    """If path is a dir, return its .tif files. If it's an index.html, return tifs
    in the same directory. If it's a .tif, return [path]."""
    if os.path.isdir(path):
        return sorted(os.path.join(path, f) for f in os.listdir(path)
                      if f.lower().endswith(('.tif', '.tiff')))
    if path.lower().endswith('.html'):
        d = os.path.dirname(os.path.abspath(path)) or '.'
        return sorted(os.path.join(d, f) for f in os.listdir(d)
                      if f.lower().endswith(('.tif', '.tiff')))
    return [path]


def main():
    parser = argparse.ArgumentParser(
        description="Generate hDOM evaluation HTML reports (with rank 1,2,11,12 + selectable slot).")
    parser.add_argument("inputs", nargs='+',
                        help="One or more .tif files, a directory, or a generated index.html")
    parser.add_argument("--ref-map", default=None,
                        help="Path to ref_map folder (overrides automatic detection)")
    args = parser.parse_args()

    all_tif_files = []
    for p in args.inputs:
        all_tif_files.extend(_collect_tifs_from_input(p))

    # dedupe while preserving order
    seen = set()
    all_tif_files = [x for x in all_tif_files if not (x in seen or seen.add(x))]

    if not all_tif_files:
        print("No .tif files found to process.")
        sys.exit(1)

    # Resolve ref_map folder
    ref_map_override = None
    if args.ref_map:
        ref_map_override = os.path.abspath(args.ref_map)
        if not os.path.isdir(ref_map_override):
            print(f"Warning: --ref-map path does not exist: {ref_map_override}")

    summaries = []
    for tif in all_tif_files:
        s = tif2eval_html(tif, ref_map_override=ref_map_override)
        if s:
            summaries.append(s)

    if summaries:
        first_dir = os.path.dirname(os.path.abspath(all_tif_files[0]))
        index_path = os.path.join(first_dir, "index.html")
        keep_path = os.path.join(first_dir, "index_keep.html")
        if os.path.exists(index_path) and not os.path.exists(keep_path):
            os.rename(index_path, keep_path)
            print(f"Preserved original index.html as: {keep_path}")
        generate_eval_index_html(summaries, index_path)
        print(f"\nEvaluation reports generated.")
        print(f"Index page: {index_path}")


if __name__ == "__main__":
    main()
