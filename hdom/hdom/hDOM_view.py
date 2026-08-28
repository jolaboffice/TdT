#!/usr/bin/env python3
"""Per-molecule alignment view: the molecule over its candidate reference windows, as HTML.

One page per molecule, written next to the molecule as {base}_view.html:
  1. Comparison panel order: green profile, rank11, rank12, data, rank1, rank2, rank?
  2. A green intensity profile of the data TIF above rank11, with peak count
     and labels-per-100kb annotation.
  3. Matching Results table split into two sections:
        top    -> ranks 1, 2, 11, 12
        bottom -> ranks 3-10, 13-20

Usage:
    cd /path/to/result_folder
    hDOM_view.py            # current directory
    hDOM_view.py /path/to   # explicit folder
"""

import os
import sys
import re
import argparse
import json
import pandas as pd
import numpy as np

# Self-contained "trash" (discard) widget baked into every view. Pressing 't'
# (layout/IME-independent via e.code) or clicking the corner button toggles a
# discard mark; state is loaded from / saved to this folder's trash.json through
# save_trash.php. Relative paths resolve correctly for both the raw .html and the
# index.php?file= served page. __BASEJSON__ is replaced with the molecule base.
TRASH_WIDGET = r"""
<script>(function(){
  var BASE = __BASEJSON__;
  var trashed = false;
  var bar=document.createElement('div');
  bar.style.cssText='position:fixed;top:0;left:0;right:0;z-index:99999;text-align:center;font-weight:bold;padding:7px;font-size:15px;font-family:sans-serif;display:none;cursor:pointer;';
  var btn=document.createElement('button');
  btn.style.cssText='position:fixed;bottom:16px;right:16px;z-index:99999;padding:9px 14px;border-radius:7px;border:2px solid #c0392b;font-weight:bold;font-size:14px;cursor:pointer;font-family:sans-serif;box-shadow:0 2px 6px rgba(0,0,0,.25);';
  function render(){
    if(trashed){bar.style.display='block';bar.style.background='#c0392b';bar.style.color='#fff';bar.textContent='🗑 TRASH — marked for discard (press t to undo)';document.documentElement.style.outline='5px solid #c0392b';document.documentElement.style.outlineOffset='-5px';btn.style.background='#c0392b';btn.style.color='#fff';btn.textContent='🗑 Trashed (t)';}
    else{bar.style.display='none';document.documentElement.style.outline='none';btn.style.background='#fff';btn.style.color='#c0392b';btn.textContent='🗑 Trash (t)';}
  }
  function save(){var fd=new FormData();fd.append('base',BASE);fd.append('trash',trashed?'1':'0');fetch('save_trash.php',{method:'POST',body:fd}).catch(function(e){});}
  function toggle(){trashed=!trashed;render();save();}
  bar.addEventListener('click',toggle);btn.addEventListener('click',function(e){e.preventDefault();toggle();});
  document.addEventListener('keydown',function(e){if(e.code==='KeyT'||e.key==='t'||e.key==='T'){var x=e.target&&e.target.tagName?e.target.tagName.toLowerCase():'';if(x==='input'||x==='textarea')return;e.preventDefault();toggle();}},false);
  function attach(){document.body.appendChild(bar);document.body.appendChild(btn);render();}
  if(document.body)attach();else document.addEventListener('DOMContentLoaded',attach);
  fetch('trash.json',{cache:'no-store'}).then(function(r){return r.json();}).then(function(d){trashed=!!(d&&d[BASE]);render();}).catch(function(e){});
})();</script>
"""

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)
# Y-variant: resolve the original (production) view modules from python/hDOM,
# while refined_green_peaks is imported from this dir (hDOM_Y) below.
sys.path.append(_SCRIPT_DIR)

from hDOM_tif2html import (
    BPP,
    process_query_image,
    process_ref_match,
    get_html_style,
    find_ref_map_folder,
    read_ref_map_from_info,
    image_to_base64,
)
from PIL import Image, ImageDraw
from DOM_lib import tif2RGprofile, tif2RGimage
from hDOM_evaluate import generate_eval_index_html
from hDOM_mal_Y import refined_green_peaks  # maligner-consistent green peaks

# --- Render the green channel at HALF height (upper half) ONLY in the REF-MAP panels
#     (l1/l2, d1/d2, m1/m2, p1/p2), so the red backbone shows through below. The DATA panel
#     (the real molecule) keeps FULL green. process_ref_match returns the ref-segment PNG as
#     base64; we zero its green channel in the bottom half. Disable with env GREEN_HALF=0. ---
if os.environ.get('GREEN_HALF', '1') == '1':
    import base64 as _b64_gh
    import io as _io_gh
    _orig_process_ref_match = process_ref_match
    def _process_ref_match_greenhalf(*_a, **_k):
        _b = _orig_process_ref_match(*_a, **_k)
        if _b:
            try:
                _im = Image.open(_io_gh.BytesIO(_b64_gh.b64decode(_b))).convert('RGB')
                _arr = np.array(_im); _h = _arr.shape[0]
                _arr[_h // 2:, :, 1] = 0   # zero green channel in the bottom half rows
                _b = image_to_base64(Image.fromarray(_arr))
            except Exception:
                pass
        return _b
    process_ref_match = _process_ref_match_greenhalf

# Terminal labels dropped from scoring: removed_peaks.json holds their pixel positions, drawn as white dots.
_REMOVED_CACHE = {}
def _load_removed_peaks(folder):
    if folder in _REMOVED_CACHE:
        return _REMOVED_CACHE[folder]
    import json as _j
    d = {}
    for _p in (os.path.join(folder, 'removed_peaks.json'),):
        if os.path.exists(_p):
            try:
                d = _j.load(open(_p, encoding='utf-8'))
            except Exception:
                d = {}
            if d:
                break
    _REMOVED_CACHE[folder] = d
    return d

# --- homolog PAF cov/id for H_local / H_mDOM ref-map labels ---------------------
# Reads the minimap2 haplotype PAF (chr{num}.paf: query=PATERNAL, target=MATERNAL);
# cols 10/11 = residue matches / alignment-block length.  For a homolog region we
# report cov (fraction of the region spanned by alignment blocks; <100% => SV/gap)
# and id (overlap-weighted per-base identity, matches/blocklen).
_CHAINS = os.environ.get('HDOM_CHAINS', os.path.join(os.path.dirname(_SCRIPT_DIR), 'chains'))
_PAF_CACHE = {}
def _load_paf_blocks(chrnum):
    if chrnum in _PAF_CACHE:
        return _PAF_CACHE[chrnum]
    path = f'{_CHAINS}/chr{chrnum}.paf'
    blocks = []
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            for line in f:
                c = line.rstrip('\n').split('\t')
                if len(c) < 11:
                    continue
                try:
                    qs, qe = int(c[2]), int(c[3])
                    ts, te = int(c[7]), int(c[8])
                    matches, blen = int(c[9]), int(c[10])
                except ValueError:
                    continue
                blocks.append((qs, qe, ts, te, matches, blen))
    _PAF_CACHE[chrnum] = blocks
    return blocks

def _paf_cov_id(ref_name, start_bp, length_bp):
    """cov (int %) and id (float %) of region [start, start+length] on haplotype
    ref_name (e.g. 'chr07_m') vs its homolog. Returns (None, None) if no PAF/overlap."""
    m = re.match(r'chr(\w+?)_([pm])$', str(ref_name))
    if not m or not length_bp or length_bp <= 0:
        return None, None
    chrnum, hap = m.group(1), m.group(2)
    blocks = _load_paf_blocks(chrnum)
    if not blocks:
        return None, None
    s, e = int(start_bp), int(start_bp) + int(length_bp)
    covered = 0; idnum = 0.0; idden = 0.0
    for qs, qe, ts, te, matches, blen in blocks:
        bs, be = (qs, qe) if hap == 'p' else (ts, te)
        ov = min(e, be) - max(s, bs)
        if ov <= 0:
            continue
        covered += ov
        span = be - bs
        if span > 0:
            w = ov / span
            idnum += matches * w
            idden += blen * w
    if idden <= 0:
        return None, None
    cov = min(covered / length_bp, 1.0) * 100.0
    ident = idnum / idden * 100.0
    return int(round(cov)), ident

def _h_label(base_label, ref_name, reg):
    """Two-line H ref-map label: base_label (H_local/H_mDOM) + cov/id sub-line."""
    cov, ident = _paf_cov_id(ref_name, reg.get('aligned_bp') if reg else None,
                             reg.get('cov_length_bp') if reg else None) if reg else (None, None)
    if cov is None:
        return base_label
    id_str = '100' if ident >= 99.95 else f'{ident:.1f}'
    return (f'{base_label}<br><span style="font-size:.78em;color:#888;font-weight:normal">'
            f'cov{cov}% id{id_str}%</span>')

# Layout order around "data": rank11, rank12, data, rank1, rank2
TOP_RANKS = [11, 12]
BOTTOM_RANKS = [1, 2]
TABLE_TOP_RANKS = [1, 2, 11, 12]
TABLE_BOTTOM_RANKS = list(range(3, 11)) + list(range(13, 21))


def _row_to_match_info(row):
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


def _render_rank_item(label, b64, rtype=None):
    if rtype is None:
        rtype = label.rstrip('0123456789')  # local/mDOM/mal/hDOM  (for per-type toggle)
    if b64 is None:
        return (f'<div class="comparison-item" data-rtype="{rtype}"><div class="image-label">{label}</div>'
                f'<div class="image-container"><div class="placeholder-black"></div></div></div>')
    return (f'<div class="comparison-item" data-rtype="{rtype}"><div class="image-label">{label}</div>'
            f'<div class="image-container">'
            f'<img src="data:image/png;base64,{b64}"></div></div>')


def _build_green_profile_data(tif_path):
    """Return (data_g_raw_list, n_peaks, labels_per_100kb, peaks_px).

    Curve = green profile from tif2RGprofile (same as before). But n_peaks and
    peaks_px are the REFINED maligner labels (refined_green_peaks) so the
    displayed count and the chart markers match what maligner actually uses.
    """
    try:
        _r, g, _raw_peaks = tif2RGprofile(tif_path)
    except Exception as e:
        print(f"Error computing green profile for {tif_path}: {e}")
        return [], 0, 0.0, []

    try:
        peaks = refined_green_peaks(tif_path)
    except Exception as e:
        print(f"Error computing refined peaks for {tif_path}: {e}")
        peaks = _raw_peaks

    w = len(g)
    molecule_length_bp = int(w * BPP)
    n_peaks = int(len(peaks))
    labels_per_100kb = (n_peaks * 100000.0 / molecule_length_bp) if molecule_length_bp > 0 else 0.0

    g_list = g.tolist() if hasattr(g, 'tolist') else list(g)
    peaks_px = [int(p) for p in peaks]
    return g_list, n_peaks, labels_per_100kb, peaks_px


def _format_cell(df_fmt, cols):
    """Apply the same cell formatting as hDOM_tif2html.get_table_html."""
    if 'aligned_bp' in cols:
        df_fmt['aligned_bp'] = df_fmt['aligned_bp'].apply(
            lambda x: f"{int(x):,}" if pd.notna(x) else "")
    if 'shift_bp' in cols:
        df_fmt['shift_bp'] = df_fmt['shift_bp'].apply(
            lambda x: f"{int(x):,}" if pd.notna(x) else "")
    for c in ['scale', 'cc_r', 'cc_g', 'log(P)']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(
                lambda x: f"{x:.2f}" if pd.notna(x) and isinstance(x, (int, float)) else x)
    for c in ['cc_rg', 'cc_rg2']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(
                lambda x: f"{x:.1e}" if pd.notna(x) and isinstance(x, (int, float)) else x)
    for c in ['misc', 'dir']:
        if c in cols:
            df_fmt[c] = df_fmt[c].apply(
                lambda x: f"{int(x)}" if pd.notna(x) else "")
    return df_fmt


def _build_table_section(df_full, sel_indices, folder_name, file_name, ref_map_path):
    """Render a matching-results table for the given 1-based sel_indices
    (each referring to a row in df_full). Rows whose sel index is out of
    range are skipped.
    """
    valid = [s for s in sel_indices if 1 <= s <= len(df_full)]
    if not valid:
        return ""

    sub = df_full.iloc[[s - 1 for s in valid]].reset_index(drop=True)
    df_fmt = sub.copy()
    cols = df_fmt.columns

    df_fmt = _format_cell(df_fmt, cols)

    # Rank hyperlinks using the original sel index, not the new row position
    raw_ranks = [str(df_full.iloc[s - 1].get('rank', s)) for s in valid]
    if 'rank' in cols:
        new_ranks = []
        for i, sel_idx in enumerate(valid):
            rank_val = raw_ranks[i]
            display = rank_val if 'poly N' in rank_val else str(sel_idx)
            ref_param = f'&ref_map={ref_map_path}' if ref_map_path else ''
            new_ranks.append(
                f'<a href="../view.php?folder={folder_name}&file={file_name}'
                f'&sel={sel_idx}{ref_param}">{display}</a>')
        df_fmt['rank'] = new_ranks

    html = df_fmt.to_html(index=False, classes='matching-table',
                          border=0, escape=False)

    parts = html.split('<tbody>')
    if len(parts) == 2:
        header, body = parts[0] + '<tbody>', parts[1]
        rows = re.findall(r'<tr>.*?</tr>', body, re.DOTALL)
        new_rows = []
        for i, row in enumerate(rows):
            sel_idx = valid[i]
            rank_raw = raw_ranks[i]
            cls = []
            if 'method' in cols:
                src_row = df_full.iloc[sel_idx - 1]
                m_val = (str(src_row['method']).lower()
                         if pd.notna(src_row['method']) else '')
                if 'poly N' in rank_raw:
                    cls.append('method-polyn')
                elif 'pcc' in m_val:
                    cls.append('method-pcc')
                elif 'mal' in m_val:
                    cls.append('method-mal')
                try:
                    if int(rank_raw) in (1, 2):
                        cls.append('rank-top')
                except (TypeError, ValueError):
                    pass

            attrs = f' data-sel="{sel_idx}"'
            if cls:
                attrs += f' class="{" ".join(cls)}"'
            row = row.replace('<tr>', f'<tr{attrs}>', 1)
            new_rows.append(row)

        footer_match = re.search(r'(</tbody>.*?</table>)', body, re.DOTALL)
        footer = footer_match.group(1) if footer_match else '</tbody></table>'
        html = header + "".join(new_rows) + footer

    return html


def _build_split_matching_tables(df, folder_name, file_name, ref_map_path):
    if df is None or df.empty:
        return "<p>No matching results available.</p>"
    # Single table of all rows; index.php (hdom_table_from_csv) re-renders at serve
    # time with the F1 column + mDOM/hDOM/maligner Top + Other-ranks ordering.
    all_idx = list(range(1, len(df) + 1))
    html = _build_table_section(df, all_idx, folder_name, file_name, ref_map_path)
    out = []
    if html:
        out.append('<h3 class="table-section-title">Matching Results</h3>')
        out.append(f'<div class="table-container">{html}</div>')
    if not out:
        return "<p>No matching results available.</p>"
    return "".join(out)


def _nav_style():
    return """
    .nav-bar { display: flex; align-items: center; justify-content: space-between;
               margin-bottom: 20px; gap: 24px; }
    .nav-btn { display: inline-flex; align-items: center; justify-content: center;
               width: 48px; height: 48px; border: 2px solid #3498db; border-radius: 8px;
               background-color: #fff; color: #3498db; font-size: 1.5em; font-weight: bold;
               cursor: pointer; text-decoration: none; transition: all 0.15s; user-select: none;
               flex-shrink: 0; }
    .nav-btn:hover { background-color: #3498db; color: #fff; }
    .nav-btn.disabled { border-color: #ccc; color: #ccc; pointer-events: none; }
    .nav-counter { font-size: 1.1em; color: #555; font-weight: bold; }
    .title-row { display: flex; align-items: baseline; justify-content: space-between; gap: 16px; }
    .title-row h1 { margin-bottom: 0; white-space: nowrap; border-bottom: none; padding-bottom: 0; }
    .title-counter { font-size: 1.3em; color: #555; font-weight: bold; white-space: nowrap; }
    .jump-num { cursor: pointer; color: #3498db; text-decoration: underline dotted; }
    .jump-num:hover { color: #1d6fa5; }
    .jump-input { width: 5em; font-size: 1em; padding: 2px 4px; border: 1px solid #3498db;
                  border-radius: 3px; font-weight: bold; }
    .progress-bar-container { width: 100%; height: 4px; background-color: #ddd; margin-bottom: 16px; }
    .progress-bar-fill { height: 100%; background-color: #3498db; transition: width 0.3s; }
    .rank-separator { height: 1px; background-color: #444; margin: 4px 0; }
    .data-separator { height: 2px; background-color: #888; margin: 4px 0; }
    .profile-section { background: #fff; padding: 15px 20px; border-radius: 8px;
                       border: 1px solid #ddd; margin: 20px 0 10px 0; max-width: 100%;
                       overflow: hidden; box-sizing: border-box; }
    .profile-section h3 { margin: 0 0 8px 0; font-size: 1.0em; color: #2c3e50;
                          font-weight: bold; }
    .profile-section h3 .profile-stats { font-weight: normal; color: #666;
                                         font-size: 0.88em; margin-left: 12px; }
    .profile-chart-container { position: relative; height: 200px; width: 100%;
                               border: 1px solid #eee; padding: 10px; border-radius: 4px;
                               box-sizing: border-box; min-width: 0; }
    .comparison-item.profile-item .image-container { background-color: #fff; }
    .green-overlay { position: absolute; top: 0; left: 50%; transform: translateX(-50%);
        max-width: 100%; height: auto; display: none !important; pointer-events: none; }
    body.green-on .green-overlay { display: block !important; }
    .profile-chart-box { position: relative; height: 110px; width: 100%; }
    .selector-box { margin: 10px 0 0 0; padding: 10px; background-color: #222;
                    border-radius: 4px; color: #ddd; font-family: 'Segoe UI', sans-serif; }
    .selector-box label { margin-right: 8px; font-weight: bold; }
    .selector-box input { width: 60px; padding: 4px 6px; font-size: 0.95em;
                          border: 1px solid #555; background-color: #111; color: #eee;
                          border-radius: 3px; }
    .selector-box .hint { margin-left: 10px; color: #888; font-size: 0.85em; }
    .placeholder-black { width: 100%; height: 30px; background-color: #000; }
    .index-link { font-size: 0.95em; color: #3498db; text-decoration: none; }
    .index-link:hover { text-decoration: underline; }
    .table-section-title { color: #2c3e50; margin-top: 22px; margin-bottom: 6px;
                           font-size: 1.05em; }
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
    .class-btn.active[data-value="A"] { background-color: #00b050; border-color: #00b050; color: #fff; }
    .class-btn.active[data-value="B"],
    .class-btn.active[data-value="C"],
    .class-btn.active[data-value="D"],
    .class-btn.active[data-value="M"] { background-color: #ffc000; border-color: #ffc000; color: #fff; }
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


CLASS_LIST = ['S', 'A', 'B', 'C', 'D', 'M', 'F']


def _build_class_toolbar():
    buttons = ''.join(
        f'<button type="button" class="class-btn" data-value="{c}">{c}</button>'
        for c in CLASS_LIST
    )
    return (
        '<div class="eval-toolbar">'
        '<span class="group-label">Class:</span>'
        + buttons +
        '<span class="group-label" style="margin-left:10px;">Mark:</span>'
        '<button type="button" id="star-btn" class="star-btn">\u2606</button>'
        '<button type="button" id="green-toggle-btn" style="margin-left:12px;padding:4px 10px;'
        'cursor:pointer;border-radius:4px;border:1px solid #27ae60;font-weight:bold;'
        'font-size:0.85em;">\u26aa Green: OFF</button>'
        '<span id="chosen-info" class="chosen-info"></span>'
        '<span class="hint" style="color:#888;font-size:0.85em;margin-left:auto;">'
        'Click a row in Matching Results to override rank1.</span>'
        '</div>'
    )


def tif2view_html(tif_path, mol_index, total_mols, prev_html, next_html,
                  ref_map_override=None, view_filenames=None, label_dots=False):
    """Generate a _view.html for a single tif with prev/next navigation."""
    ext = os.path.splitext(tif_path)[1].lower()
    if ext not in ('.tif', '.tiff'):
        print(f"Error: Input must be a .tif/.tiff file. Got: {tif_path}")
        return None

    base_path = os.path.splitext(tif_path)[0]
    info_path = f"{base_path}_info.txt"
    csv_path = f"{base_path}_matching_results.csv"
    html_path = f"{base_path}_view.html"

    # Info
    info_content = f"Information file not found: {info_path}"
    if os.path.exists(info_path):
        with open(info_path, 'r', encoding='utf-8') as f:
            info_content = f.read()

    # Query image — alignment runs on the binarized (_samevec) tif, but the
    # "data" panel should show the REAL molecule. Prefer an original tif if one
    # sits beside the samevec tif (_orig.tif / .tif); fall back to tif_path.
    # label_dots mode (qlalf_acro only): data panel = REAL molecule (original tif)
    # with maligner-label green dots drawn ABOVE the rectangle. Default OFF keeps
    # the original behavior (data = the given tif, no dots) for all other folders.
    display_tif = tif_path
    if label_dots:
        _fold = os.path.dirname(os.path.abspath(tif_path))
        _db = base_path[:-len('_samevec')] if base_path.endswith('_samevec') else base_path  # full-path stem
        _stem = os.path.basename(_db)                                                        # basename stem
        # Prefer a LOCAL sibling _orig.tif (samevec-included name) — for straightened
        # molecules this is the linearized real molecule. Only if absent fall back to
        # the samevec-stripped lookup / /tmp/term_orig (acro/qlalf real-molecule panel).
        for _cand in (base_path + '_orig.tif', base_path + '_orig.tiff',
                      _db + '_orig.tif', _db + '_orig.tiff', _db + '.tif', _db + '.tiff'):
            if os.path.exists(_cand) and os.path.abspath(_cand) != os.path.abspath(tif_path):
                display_tif = _cand
                break
    query_img, query_b64 = process_query_image(display_tif)
    if not query_img:
        return None
    overlay_b64 = None   # synthetic green dots live on a separate transparent overlay so they can be toggled;
                         # they are never baked into the data image
    if label_dots:
        try:
            # Green dots = peaks kept in the scored tif, not the raw peaks of the display tif
            _pk = [int(round(10 + p)) for p in refined_green_peaks(tif_path)]
            _w, _h = query_img.size
            _r = max(2, round(_h * 0.116)); _dr = _r * 2; _mg = 2 * _dr + 4   # green=diamond 2x
            # Markers are drawn on a supersampled RGBA strip, then LANCZOS-downscaled.
            _ss = 4
            _strip = Image.new('RGBA', (_w * _ss, _mg * _ss), (0, 0, 0, 0))
            _sd = ImageDraw.Draw(_strip)
            _cyS = _dr * _ss; _R = _dr * _ss
            for _cx in _pk:
                _X = _cx * _ss
                _sd.polygon([(_X, _cyS-_R), (_X+_R, _cyS), (_X, _cyS+_R), (_X-_R, _cyS)],
                            fill=(0, 255, 0, 255), outline=(0, 140, 0, 255))
            # White dots = positions of the dropped terminal labels
            try:
                _rmv = _load_removed_peaks(_fold)
                _RW = max(3, _r * 2) * _ss
                for _cx in _rmv.get(_stem, []):
                    _X = int(round(10 + _cx)) * _ss
                    _sd.ellipse([_X-_RW, _cyS-_RW, _X+_RW, _cyS+_RW],
                                fill=(255, 255, 255, 255), outline=(30, 30, 30, 255), width=max(2, _ss))
            except Exception:
                pass
            _strip = _strip.resize((_w, _mg), Image.LANCZOS)
            # Clean data panel: black strip + molecule, with no markers baked in
            _cv = Image.new('RGB', (_w, _h + _mg), (0, 0, 0))
            _cv.paste(query_img, (0, _mg))
            # Overlay: transparent PNG the size of the data panel, markers only in the strip
            _ov = Image.new('RGBA', (_w, _h + _mg), (0, 0, 0, 0))
            _ov.paste(_strip, (0, 0), _strip)
            query_img, query_b64 = _cv, image_to_base64(_cv)
            overlay_b64 = image_to_base64(_ov)
        except Exception as _e:
            print(f"[label-dots] {os.path.basename(display_tif)}: {_e}")
    query_final_w = query_img.size[0]
    molecule_length_bp = int((query_final_w - 20) * BPP)

    # Green intensity profile data for the Chart.js plot above the rank stack
    data_g_raw, n_peaks, labels_per_100kb, peaks_px = _build_green_profile_data(tif_path)

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
        'labels_per_100kb': labels_per_100kb,
    }

    tif_dir = os.path.dirname(os.path.abspath(tif_path))
    ref_map_folder = (ref_map_override if ref_map_override
                      else read_ref_map_from_info(tif_path)
                      or find_ref_map_folder(tif_dir))

    try:
        _, raw_w, _ = tif2RGimage(tif_path)
    except Exception:
        raw_w = max(query_final_w - 20, 1)

    # Pre-render ref images
    all_ref_b64 = {}
    ranks_data = []
    if df is not None and not df.empty:
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
            r_aligned = row.get('aligned_bp')
            if pd.isna(r_aligned):
                r_aligned = row.get('shift_bp', 0)
            r_logp = row.get('log(P)')
            r_scale = row.get('scale', 1.0)
            if pd.isna(r_scale):
                r_scale = 1.0
            r_scale = float(r_scale)
            if not np.isfinite(r_scale):   # e.g. mal scale can be inf on degenerate fragment counts
                r_scale = 0.0
            _cov = raw_w * BPP * r_scale
            ranks_data.append({
                'ref_name': ('' if pd.isna(ref_name) else str(ref_name)) or 'N/A',
                'aligned_bp': int(r_aligned) if not pd.isna(r_aligned) else 0,
                'log_p': float(r_logp) if not pd.isna(r_logp) else None,
                'dir': '+' if row.get('dir', 1) == 1 else '-',
                'scale': r_scale,
                'cov_length_bp': int(_cov) if np.isfinite(_cov) else 0,
            })
            if pd.isna(ref_name) or not ref_name:
                continue
            b64 = process_ref_match(_row_to_match_info(row), ref_map_folder,
                                    query_final_w)
            if b64:
                all_ref_b64[i + 1] = b64
    summary_info['ranks'] = ranks_data
    summary_info['ref_map_folder'] = ref_map_folder

    # Build comparison block: profile, rank11, rank12, data, rank1, rank2, rank?
    comp_items = []

    # Green profile chart: shown by default (original behavior); omitted only in
    # label_dots mode (qlalf_acro), where green labels are shown as dots instead.
    if not label_dots:
        comp_items.append(
            '<div class="comparison-item profile-item">'
            '<div class="image-label">profile</div>'
            '<div class="image-container">'
            '<div class="profile-chart-box"><canvas id="chartDataG"></canvas></div>'
            '</div></div>'
        )
        comp_items.append('<div class="rank-separator"></div>')

    # Select PCC top-2 and maligner top-2 by METHOD (robust to a variable number
    # of PCC rows). The old code assumed exactly 10 PCC rows so maligner sat at
    # df positions 11/12; when PCC has fewer rows the maligner best hit shifts up
    # and would otherwise drop out of the comparison panel.
    def _pos_by_method(method, n=2):
        if df is None or df.empty or 'method' not in df.columns:
            return []
        out = []
        for i in range(len(df)):
            r = df.iloc[i]
            if str(r.get('method', '')).strip().lower() != method.lower():
                continue
            if 'poly N' in str(r.get('rank', '')):  # poly-N rows never fill a rank slot
                continue
            out.append(i + 1)
            if len(out) >= n:
                break
        return out
    local_pos = _pos_by_method('local')
    mdom_pos = _pos_by_method('mdom') or _pos_by_method('mal') or list(TOP_RANKS)
    pcc_pos = _pos_by_method('pcc') or list(BOTTOM_RANKS)

    # Panel order: data / mDOM1,2 / mal1,2 / hDOM1,2  (gray separator per group)
    # data first (inline so we can tag the img with an id)
    _ov_html = (f'<img class="green-overlay" src="data:image/png;base64,{overlay_b64}">'
                if overlay_b64 else '')
    comp_items.append(
        '<div class="comparison-item"><div class="image-label">data</div>'
        '<div class="image-container" style="position:relative">'
        f'<img id="data-img" src="data:image/png;base64,{query_b64}">{_ov_html}</div></div>'
    )
    comp_items.append('<div class="data-separator"></div>')

    # local rank1, rank2 (local-align method; only present in the 4-method comparison CSVs)
    for j, r in enumerate(local_pos, 1):
        comp_items.append(_render_rank_item(f'local{j}', all_ref_b64.get(r)))
    if local_pos:
        comp_items.append('<div class="data-separator"></div>')

    # mDOM rank1, rank2 (re-anchored; falls back to mal for CSVs that carry no mdom rows).
    # That fallback must not fire when the CSV already carries its own mal rows: the maligner
    # panels are rendered again below, so it would draw the identical two panels twice.
    _mdom_rows = _pos_by_method('mdom')
    if _mdom_rows or not _pos_by_method('mal'):
        _top_label = 'mDOM' if _mdom_rows else 'mal'
        for j, r in enumerate(mdom_pos, 1):
            comp_items.append(_render_rank_item(f'{_top_label}{j}', all_ref_b64.get(r)))
        comp_items.append('<div class="data-separator"></div>')

    # maligner rank1, rank2
    mal_pos = _pos_by_method('mal')
    for j, r in enumerate(mal_pos, 1):
        comp_items.append(_render_rank_item(f'mal{j}', all_ref_b64.get(r)))
    comp_items.append('<div class="data-separator"></div>')

    # hDOM = original PCC method -> label hDOM1, hDOM2
    for j, r in enumerate(pcc_pos, 1):
        comp_items.append(_render_rank_item(f'hDOM{j}', all_ref_b64.get(r)))

    # H_local panel (homolog of local rank1) -> tied to the 'local' ref_map toggle
    hloc_pos = _pos_by_method('h_local', n=1)
    if hloc_pos:
        comp_items.append('<div class="data-separator"></div>')
        for p in hloc_pos:
            _hn = str(df.iloc[p - 1].get('ref_name', '')) if df is not None else ''
            _reg = ranks_data[p - 1] if p - 1 < len(ranks_data) else None
            comp_items.append(_render_rank_item(_h_label('H_local', _hn, _reg), all_ref_b64.get(p), rtype='local'))
    # H_mDOM panel (homolog of mDOM rank1) -> tied to the 'mDOM' ref_map toggle
    hmd_pos = _pos_by_method('h_mdom', n=1)
    if hmd_pos:
        comp_items.append('<div class="data-separator"></div>')
        for p in hmd_pos:
            _hn = str(df.iloc[p - 1].get('ref_name', '')) if df is not None else ''
            _reg = ranks_data[p - 1] if p - 1 < len(ranks_data) else None
            comp_items.append(_render_rank_item(_h_label('H_mDOM', _hn, _reg), all_ref_b64.get(p), rtype='mDOM'))
    # a single 'homolog' row, where the method is not split by channel -> treat as H_local
    hom_pos = _pos_by_method('homolog', n=1)
    if hom_pos:
        comp_items.append('<div class="data-separator"></div>')
        for p in hom_pos:
            _hn = str(df.iloc[p - 1].get('ref_name', '')) if df is not None else ''
            _reg = ranks_data[p - 1] if p - 1 < len(ranks_data) else None
            comp_items.append(_render_rank_item(_h_label('H_local', _hn, _reg), all_ref_b64.get(p), rtype='local'))

    # Separator before selectable slot
    comp_items.append('<div class="rank-separator"></div>')

    # Selectable slot
    comp_items.append(
        '<div class="comparison-item" id="sel-item">'
        '<div class="image-label"><span id="sel-label">hDOM ?</span></div>'
        '<div class="image-container">'
        '<img id="sel-img" style="display:none;">'
        '<div id="sel-placeholder" class="placeholder-black"></div>'
        '</div></div>'
    )

    max_rank = max(len(df) if df is not None else 0, 1)
    # --- selectable-rank labels: d{n}=mDOM, m{n}=maligner, p{n}=hDOM(pcc) -> sel_idx ---
    _pref = {'local': 'l', 'mdom': 'd', 'mal': 'm', 'pcc': 'p'}
    _label_map = {}
    _lab_order = {'l': [], 'd': [], 'm': [], 'p': []}
    if df is not None and not df.empty and 'method' in df.columns:
        for _pos in range(len(df)):
            _m = str(df.iloc[_pos].get('method', '')).strip().lower()
            _rk = str(df.iloc[_pos].get('rank', '')).strip()
            if 'poly' in _rk.lower():
                continue
            _p = next((v for k, v in _pref.items() if k in _m), None)
            if _p is None:
                continue
            try:
                _n = int(float(_rk))
            except (TypeError, ValueError):
                continue
            _label_map[f'{_p}{_n}'] = _pos + 1
            _lab_order[_p].append((_n, f'{_p}{_n}'))
    _ordered_labels = [l for p in ('l', 'd', 'm', 'p') for _, l in sorted(_lab_order[p])]
    datalist_html = ''.join(f'<option value="{l}">' for l in _ordered_labels)
    comp_html = ('<div class="comparison-box"><div class="comparison-images">'
                 + ''.join(comp_items) + '</div>'
                 + '<div class="selector-box">'
                 + '<label for="sel-input">Select rank:</label>'
                 + '<input type="text" id="sel-input" list="sel-rank-list" autocomplete="off"'
                 + ' style="width:7em;" placeholder="e.g. d1, m3, p2">'
                 + f'<datalist id="sel-rank-list">{datalist_html}</datalist>'
                 + '<span class="hint">Enter d1.., m1.., p1.. to display.</span>'
                 + '</div></div>')

    # Green intensity profile section (placed above the comparison box /
    # above rank11), rendered with Chart.js to match view.php exactly.
    profile_data_js = json.dumps({
        'data_g': data_g_raw,
        'bpp': BPP,
        'molecule_length_bp': len(data_g_raw) * BPP,
        'peaks_px': peaks_px,
    })
    profile_html = (
        '<div class="profile-section">'
        f'<h3>Data Green<span class="profile-stats">peaks: {n_peaks} &nbsp;|&nbsp; '
        f'{labels_per_100kb:.2f} labels / 100 kb</span></h3>'
        '</div>'
    )

    # Molecule intensity profile for the click-to-view popup (label_dots/qlalf_acro
    # only; left empty elsewhere so the popup stays inert).
    _mol_prof = {'r': [], 'g': [], 'bpp': BPP}
    if label_dots:
        try:
            _pr, _pg, _ = tif2RGprofile(display_tif)
            _mol_prof = {'r': [round(float(v), 1) for v in _pr],
                         'g': [round(float(v), 1) for v in _pg], 'bpp': BPP}
        except Exception as _e:
            print(f"[mol-profile] {os.path.basename(display_tif)}: {_e}")

    # Selector + profile JS
    js_data = json.dumps(all_ref_b64)
    label_map_js = json.dumps(_label_map)
    tif_name_js = json.dumps(os.path.basename(tif_path))
    view_files_js = json.dumps(view_filenames or [])
    mol_prof_js = json.dumps(_mol_prof)
    selector_js = f"""
    <script>
    const RANK_IMAGES = {js_data};
    const LABEL_MAP = {label_map_js};
    const pData = {profile_data_js};
    const MOL_PROFILE = {mol_prof_js};
    const TIF_NAME = {tif_name_js};
    const VIEW_FILES = {view_files_js};
    let _molProfChart = null;
    function showMolProf() {{
        const ov = document.getElementById('mol-prof-overlay');
        if (!ov) return;
        if (!MOL_PROFILE || !MOL_PROFILE.r || !MOL_PROFILE.r.length) return;
        ov.style.display = 'block';
        if (!_molProfChart && typeof Chart !== 'undefined') {{
            const bpp = MOL_PROFILE.bpp || 200;
            const r = MOL_PROFILE.r || [], g = MOL_PROFILE.g || [];
            const xs = r.map(function(_v, i) {{ return (i * bpp) / 1000; }});
            _molProfChart = new Chart(
                document.getElementById('mol-prof-chart').getContext('2d'), {{
                type: 'line',
                data: {{ labels: xs, datasets: [
                    {{ label: 'red (DNA)', data: r, borderColor: 'rgba(255,70,70,1)',
                       borderWidth: 1.2, pointRadius: 0, fill: false }},
                    {{ label: 'green (label)', data: g, borderColor: 'rgba(0,220,0,1)',
                       borderWidth: 1.2, pointRadius: 0, fill: false }}
                ] }},
                options: {{ animation: false, maintainAspectRatio: false,
                    scales: {{ x: {{ title: {{ display: true, text: 'kb', color: '#bbb' }},
                                     ticks: {{ color: '#bbb', maxTicksLimit: 12 }} }},
                               y: {{ ticks: {{ color: '#bbb' }} }} }},
                    plugins: {{ legend: {{ labels: {{ color: '#ddd' }} }} }} }}
            }});
        }}
    }}
    const STORAGE_CLASS = 'hdom_eval_class_' + TIF_NAME;
    const STORAGE_STAR = 'hdom_eval_star_' + TIF_NAME;
    const STORAGE_CHOSEN = 'hdom_eval_chosen_' + TIF_NAME;
    document.addEventListener('DOMContentLoaded', function() {{
        const input = document.getElementById('sel-input');
        const img = document.getElementById('sel-img');
        const ph = document.getElementById('sel-placeholder');
        const label = document.getElementById('sel-label');
        function update() {{
            const key = (input.value || '').trim().toLowerCase();
            const v = LABEL_MAP[key];
            if (v && RANK_IMAGES[v]) {{
                img.src = 'data:image/png;base64,' + RANK_IMAGES[v];
                img.style.display = 'block';
                ph.style.display = 'none';
                label.textContent = key;
            }} else {{
                img.style.display = 'none';
                ph.style.display = 'block';
                label.textContent = 'hDOM ?';
            }}
        }}
        input.addEventListener('input', update);
        input.addEventListener('change', update);

        // --- class buttons (shared storage key with hDOM_evaluate.py) ---
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

        // --- green artifact-dot toggle (default off, body.green-on, persisted across navigation via localStorage) ---
        const greenBtn = document.getElementById('green-toggle-btn');
        function refreshGreen() {{
            const on = localStorage.getItem('greenMarkers') === '1';
            if (on) document.body.classList.add('green-on');
            else document.body.classList.remove('green-on');
            if (greenBtn) {{
                greenBtn.textContent = on ? '\U0001F7E2 Green: ON' : '⚪ Green: OFF';
                greenBtn.style.background = on ? '#27ae60' : '#eee';
                greenBtn.style.color = on ? '#fff' : '#444';
            }}
        }}
        if (greenBtn) greenBtn.addEventListener('click', function() {{
            const cur = localStorage.getItem('greenMarkers') === '1';
            if (cur) localStorage.removeItem('greenMarkers');
            else localStorage.setItem('greenMarkers', '1');
            refreshGreen();
        }});
        refreshGreen();

        // --- matching-table row click = choose correct rank (override rank1) ---
        // The view uses split tables; use data-sel attribute to get original rank.
        const rows = document.querySelectorAll('table.matching-table tbody tr[data-sel]');
        const chosenInfo = document.getElementById('chosen-info');
        function refreshChosen() {{
            const ch = parseInt(localStorage.getItem(STORAGE_CHOSEN), 10);
            rows.forEach(function(tr) {{
                const sel = parseInt(tr.getAttribute('data-sel'), 10);
                if (sel === ch) tr.classList.add('chosen-row');
                else tr.classList.remove('chosen-row');
            }});
            if (!isNaN(ch)) chosenInfo.textContent = '(override \u2192 rank ' + ch + ')';
            else chosenInfo.textContent = '';
        }}
        rows.forEach(function(tr) {{
            tr.querySelectorAll('a').forEach(function(a) {{
                a.addEventListener('click', function(e) {{ e.stopPropagation(); }});
            }});
            // rank name -> show this molecule's intensity profile (no rank override)
            tr.querySelectorAll('.rank-prof').forEach(function(sp) {{
                sp.addEventListener('click', function(e) {{ e.stopPropagation(); showMolProf(); }});
            }});
            tr.addEventListener('click', function() {{
                const rank = parseInt(tr.getAttribute('data-sel'), 10);
                if (isNaN(rank)) return;
                const cur = parseInt(localStorage.getItem(STORAGE_CHOSEN), 10);
                if (cur === rank) localStorage.removeItem(STORAGE_CHOSEN);
                else localStorage.setItem(STORAGE_CHOSEN, String(rank));
                refreshChosen();
            }});
        }});
        refreshChosen();

        // --- clickable jump-to-molecule counter ---
        function jumpTo(n) {{
            if (!Number.isInteger(n) || n < 1 || n > VIEW_FILES.length) return;
            const target = VIEW_FILES[n - 1];
            if (target) window.location.href = target;
        }}
        document.querySelectorAll('.jump-num').forEach(function(span) {{
            span.addEventListener('click', function() {{
                const parent = span.parentNode;
                const inp = document.createElement('input');
                inp.type = 'number';
                inp.min = 1;
                inp.max = VIEW_FILES.length;
                inp.value = span.textContent.trim();
                inp.className = 'jump-input';
                parent.replaceChild(inp, span);
                inp.focus();
                inp.select();
                let done = false;
                function restore() {{
                    if (done) return; done = true;
                    parent.replaceChild(span, inp);
                }}
                inp.addEventListener('keydown', function(e) {{
                    if (e.key === 'Enter') {{
                        const n = parseInt(inp.value, 10);
                        if (!isNaN(n)) {{ jumpTo(n); }} else {{ restore(); }}
                    }} else if (e.key === 'Escape') {{
                        restore();
                    }}
                }});
                inp.addEventListener('blur', restore);
            }});
        }});

        // ---- Data Green intensity profile, x-axis in real molecule kb ----
        const bpp = pData.bpp;
        const xMaxKb = pData.molecule_length_bp / 1000;

        // Raw green channel, scaled to half height (matches view.php)
        const data_g = (pData.data_g || []).map(v => v / 2);
        // raw molecule pixel width (without the 10 px tif padding on each side)
        const rawW = (pData.data_g || []).length;
        const paddedW = rawW + 20;

        function toKbSeries(arr) {{
            return arr.map((y, idx) => ({{ x: (idx * bpp) / 1000, y: y }}));
        }}

        const profileOptions = {{
            responsive: true,
            maintainAspectRatio: false,
            interaction: {{ intersect: false, mode: 'index' }},
            layout: {{ padding: {{ left: 0, right: 0, top: 4, bottom: 4 }} }},
            scales: {{
                x: {{
                    display: true,
                    type: 'linear',
                    min: 0,
                    max: xMaxKb,
                    grid: {{ display: true }},
                    ticks: {{
                        stepSize: 20,
                        callback: function(value) {{ return Math.round(value); }}
                    }}
                }},
                y: {{ beginAtZero: true, ticks: {{ display: false }} }}
            }},
            plugins: {{ legend: {{ display: false }} }},
            animation: false
        }};

        const gCanvas = document.getElementById('chartDataG');
        let gChart = null;
        if (gCanvas && typeof Chart !== 'undefined') {{
            gChart = new Chart(gCanvas.getContext('2d'), {{
                type: 'line',
                data: {{
                    datasets: [{{
                        label: 'Data Green',
                        data: toKbSeries(data_g),
                        borderColor: 'rgba(0, 255, 0, 1)',
                        borderWidth: 2,
                        fill: false,
                        backgroundColor: 'rgba(0, 255, 0, 0.1)',
                        pointRadius: 0
                    }}, {{
                        label: 'maligner labels',
                        data: (pData.peaks_px || []).map(function(i) {{ return {{ x: (i * bpp) / 1000, y: data_g[i] }}; }}),
                        showLine: false,
                        pointRadius: 4,
                        pointBackgroundColor: 'rgba(255, 45, 45, 1)',
                        pointBorderColor: 'rgba(255, 45, 45, 1)'
                    }}]
                }},
                options: profileOptions
            }});
        }}

        // Align the chart's plotting area with the rendered data image so
        // that each label's peak sits directly above the nick in the tif.
        // The data image uses `max-width: 100%; height: auto; margin: 0 auto`,
        // so short molecules are displayed at natural size (centered) and
        // long ones are shrunk to the container width. We measure the image
        // at runtime and translate that into Chart.js layout padding.
        function alignProfileToDataImg() {{
            if (!gChart) return;
            const dataImg = document.getElementById('data-img');
            const chartBox = gCanvas && gCanvas.parentNode;
            if (!dataImg || !chartBox) return;
            if (!dataImg.complete || dataImg.naturalWidth === 0) return;
            const imgRect = dataImg.getBoundingClientRect();
            const boxRect = chartBox.getBoundingClientRect();
            if (imgRect.width <= 0 || boxRect.width <= 0) return;
            const scale = imgRect.width / paddedW;       // CSS px per tif px
            const contentLeft = imgRect.left + 10 * scale;
            const contentRight = imgRect.right - 10 * scale;
            const padLeft = Math.max(0, contentLeft - boxRect.left);
            const padRight = Math.max(0, boxRect.right - contentRight);
            // guard: only update when padding actually changed (>0.5px) so a scrollbar
            // toggle -> resize -> re-align feedback loop can't jitter the page.
            if (Math.abs((gChart.options.layout.padding.left||0) - padLeft) > 0.5 ||
                Math.abs((gChart.options.layout.padding.right||0) - padRight) > 0.5) {{
                gChart.options.layout.padding.left = padLeft;
                gChart.options.layout.padding.right = padRight;
                gChart.update('none');
            }}
        }}

        const _dataImgEl = document.getElementById('data-img');
        if (_dataImgEl) {{
            if (_dataImgEl.complete) {{
                requestAnimationFrame(alignProfileToDataImg);
            }} else {{
                _dataImgEl.addEventListener('load', alignProfileToDataImg);
            }}
        }}
        window.addEventListener('resize', function() {{
            requestAnimationFrame(alignProfileToDataImg);
        }});
    }});

    // Keyboard navigation: left/right arrow keys
    document.addEventListener('keydown', function(e) {{
        if (e.target.tagName === 'INPUT') return;  // don't hijack input fields
        if (e.key === 'ArrowLeft') {{
            const prev = document.getElementById('nav-prev');
            if (prev && !prev.classList.contains('disabled')) prev.click();
        }} else if (e.key === 'ArrowRight') {{
            const next = document.getElementById('nav-next');
            if (next && !next.classList.contains('disabled')) next.click();
        }}
    }});
    </script>
    """

    # Navigation buttons
    prev_cls = '' if prev_html else ' disabled'
    next_cls = '' if next_html else ' disabled'
    prev_href = f'href="{prev_html}"' if prev_html else 'href="#"'
    next_href = f'href="{next_html}"' if next_html else 'href="#"'

    progress_pct = (mol_index / total_mols * 100) if total_mols > 0 else 0

    nav_html = (
        '<div class="title-row">'
        '<h1>hDOM Analysis Report</h1>'
        f'<span class="title-counter"><span class="jump-num" title="Click to jump to a molecule">{mol_index}</span> / {total_mols}</span>'
        '</div>'
        f'<div class="progress-bar-container"><div class="progress-bar-fill" style="width:{progress_pct:.1f}%"></div></div>'
        '<div class="nav-bar">'
        f'<a id="nav-prev" class="nav-btn{prev_cls}" {prev_href} title="Previous molecule (Left arrow)">&#9664;</a>'
        f'<span class="nav-counter"><span class="jump-num" title="Click to jump to a molecule">{mol_index}</span> / {total_mols}</span>'
        f'<a id="nav-next" class="nav-btn{next_cls}" {next_href} title="Next molecule (Right arrow)">&#9654;</a>'
        '</div>'
    )

    # Matching tables (split into two sections)
    tif_abs = os.path.abspath(tif_path)
    folder_name = os.path.basename(os.path.dirname(tif_abs))
    file_name = os.path.basename(tif_abs)
    tables_html = _build_split_matching_tables(
        df, folder_name=folder_name, file_name=file_name,
        ref_map_path=ref_map_folder)

    # self-contained trash widget: 't' key / corner button toggles a discard mark,
    # persisted in this folder's trash.json via save_trash.php. Works whether the
    # page is opened raw or served through index.php (relative fetch paths).
    trash_widget = TRASH_WIDGET.replace(
        '__BASEJSON__', json.dumps(os.path.splitext(os.path.basename(tif_path))[0]))

    green_profile_html = '' if label_dots else profile_html
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>hDOM View - {os.path.basename(tif_path)}</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>{get_html_style()}{_nav_style()}</style>
</head>
<body>
    <!-- The backslashes in the regex below are doubled: this is a Python string, and a lone
         backslash before a dot is not a Python escape, which warns at compile time. -->
    <script data-vredirect>(function(){{var p=location.pathname;if(/_view\\.html$/.test(p)){{location.replace('index.php?file='+p.split('/').pop().replace(/_view\\.html$/,''));}}}})();</script>
    <div class="container">
        {nav_html}
        <p><strong>File:</strong> {os.path.basename(tif_path)} ({molecule_length_bp / 1000:.1f} kbp)
           &nbsp;&nbsp;<a class="index-link" href="index.php">&#9776; Index</a></p>
        {_build_class_toolbar()}
        {green_profile_html}
        {comp_html}
        <h2>Matching Results</h2>
        {tables_html}
        <h2>Information</h2>
        <div class="info-box">{info_content}</div>
        {nav_html}
    </div>
    <div id="mol-prof-overlay" style="display:none;position:fixed;inset:0;background:rgba(0,0,0,0.6);z-index:9999;" onclick="if(event.target===this)this.style.display='none'">
      <div style="position:absolute;left:50%;top:50%;transform:translate(-50%,-50%);background:#111;padding:14px 16px;border-radius:8px;width:82%;max-width:1050px;box-shadow:0 4px 24px rgba(0,0,0,.6);">
        <div style="color:#eee;font-size:14px;margin-bottom:8px;font-family:sans-serif;">Molecule intensity profile
          <span style="float:right;cursor:pointer;color:#aaa;font-size:18px;" onclick="document.getElementById('mol-prof-overlay').style.display='none'">&#10005;</span></div>
        <div style="background:#000;padding:8px;border-radius:4px;height:260px;"><canvas id="mol-prof-chart"></canvas></div>
      </div>
    </div>
    {selector_js}
    <footer>Generated by hDOM_view.py</footer>
    {trash_widget}
</body>
</html>"""

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Created: {html_path}")
    return summary_info


def main():
    parser = argparse.ArgumentParser(
        description="Generate hDOM view HTML reports with prev/next navigation.")
    parser.add_argument('data_folder', nargs='?', default='.',
                        help='Folder with hDOM results (default: current directory)')
    parser.add_argument('--ref-map', default=None,
                        help='Path to ref_map folder (overrides automatic detection)')
    parser.add_argument('--label-dots', action='store_true',
                        help='Data panel = original tif with green-label dots above it, '
                             'green profile chart omitted (qlalf_acro style). Default off.')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip molecules whose _view.html already exists (incremental).')
    parser.add_argument('--shard', default=None,
                        help='k/N: process only molecules with index %% N == k (parallel sharding).')
    args = parser.parse_args()
    _shk, _shn = (int(x) for x in args.shard.split('/')) if args.shard else (0, 1)

    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        print(f"Error: Directory not found: {data_folder}")
        sys.exit(1)

    tif_files = sorted([
        os.path.join(data_folder, f)
        for f in os.listdir(data_folder)
        if f.lower().endswith(('.tif', '.tiff'))
        and not f.lower().endswith(('_orig.tif', '_orig.tiff'))
    ])

    if not tif_files:
        print(f"No .tif files found in {data_folder}")
        sys.exit(1)

    total = len(tif_files)
    print(f"Found {total} .tif files in {data_folder}")
    print("Generating view HTML reports (profile + ranks 11,12,data,1,2)...")
    print()

    ref_map_override = None
    if args.ref_map:
        ref_map_override = os.path.abspath(args.ref_map)

    view_filenames = []
    for tif in tif_files:
        base = os.path.splitext(os.path.basename(tif))[0]
        view_filenames.append(f"{base}_view.html")

    summaries = []
    for i, tif_file in enumerate(tif_files):
        mol_index = i + 1
        prev_html = view_filenames[i - 1] if i > 0 else None
        next_html = view_filenames[i + 1] if i < total - 1 else None

        if _shn > 1 and (i % _shn) != _shk:
            continue                                   # not this shard
        if args.skip_existing and os.path.exists(os.path.join(data_folder, view_filenames[i])):
            continue                                   # incremental: view already made

        summary = tif2view_html(
            tif_file, mol_index, total,
            prev_html, next_html,
            ref_map_override=ref_map_override,
            view_filenames=view_filenames,
            label_dots=args.label_dots,
        )
        if summary:
            summaries.append(summary)

    if not summaries:
        print("No reports generated.")
        sys.exit(1)

    index_path = os.path.join(data_folder, "index.html")
    keep_path = os.path.join(data_folder, "index_keep.html")
    if os.path.exists(index_path) and not os.path.exists(keep_path):
        os.rename(index_path, keep_path)
        print(f"Preserved original index.html as: {keep_path}")
    generate_eval_index_html(
        summaries, index_path,
        html_suffix='_view.html',
        page_title='hDOM View Summary Index',
        generator='hDOM_view.py',
    )

    print(f"\nDone! {len(summaries)} view reports generated.")
    print(f"Index page: {index_path}")


if __name__ == '__main__':
    main()
