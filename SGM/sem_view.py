#!/usr/bin/env python3
"""Build the SGM view -- SEM (QD) molecules in the standard hDOM align-page frame.

The SEM project already ships, per molecule, everything a viewer needs:
    {QD}.maps              query fragment sizes (bp)  -- what maligner aligned
    {QD}_mal.txt           maligner DP output, ranked candidates
    {QD}_points.json       the traced polyline in the stitched SEM image
    {QD}_adjusted.tif      the stitched, contrast-adjusted SEM image
    {QD}_rank{n}_strip.png rank 1..6 figures: Query / Ref(smoothed) / Ref(raw) label strips

so nothing has to be re-aligned. This script only re-frames it:

  data panel   the molecule itself -- the SEM image resampled along the traced polyline's arc
               length, so the panel's x axis is the measured bp axis
  labels       the Query band lifted out of the rank strip figure
  rank1..6     the Ref(smoothed) band lifted out of each rank strip figure
  table        {QD}_matching_results.csv in the hDOM schema, built from {QD}_mal.txt --
               index.php re-renders the matching table from it, so the page inherits the
               standard table, navigation, class buttons and marking UI for free

The rank strip figures are a fixed matplotlib layout (2400x600, bands at rows 177/222/267,
x 285..2255), verified identical across all 53 molecules, so the bands crop by constant.

SEM carries no AT signal, so there is no red channel and no cc_r: `mal` is the only method,
and the reference shown is the SEM project's own smoothed reference (rm100), not hDOM's
ref_map. The two are the same assembly and enzyme (HG002v1.1 / Nb.BssSI; chromosome lengths
agree to 0.01%), they differ only in how the map was derived, so the names map straight over.

Usage: python3 sem_view.py --data <data dir> --out <folder> [--template-index <index.php>] [--only QD025 ...]
"""
import argparse
import base64
import csv
import glob
import io
import json
import os
import re
import sys

import numpy as np
import tifffile as tiff
from PIL import Image
from scipy.ndimage import map_coordinates

# The style and toolbar are borrowed from the hDOM viewer, so an SGM page sits in the same frame
# as an alignment page. This reads the hdom package in this repository (hDOM/hdom); HDOM_DIR points
# elsewhere.
_HDOM = os.environ.get('HDOM_DIR') or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'hDOM', 'hdom')
sys.path.insert(0, _HDOM)
from hDOM_tif2html import get_html_style                      # noqa: E402
import hDOM_view as V                                         # noqa: E402


def toolbar():
    """Class / mark toolbar without the green overlay toggle (SGM has no green overlay)."""
    return re.sub(r'<button[^>]*id="green-toggle-btn".*?</button>', '', V._build_class_toolbar())

Image.MAX_IMAGE_PIXELS = None

# Every path is an argument, with no default. --data is the directory holding the <QD###>/ folders,
# --out the folder the view is written to, --template-index the index.php to copy in.
SRC = OUT = TEMPLATE_INDEX = None

# rank strip figure geometry: the bands sit in these row windows, but their exact extent is
# detected per figure -- the width varies (2400..7200 px) and the axes frame is anti-aliased,
# so a fixed crop drags in the grey border lines.
BAND = {'query': (170, 215), 'ref_sm': (215, 260), 'ref_raw': (260, 305)}
# Common display width for every panel. The source rank figures are 2400-3780 px wide, so a
# narrower W merges label ticks that sit a few pixels apart and the molecule then shows fewer
# labels than its .maps has (QD065: 18 -> 7). Keep this at or above the widest source figure.
W = int(os.environ.get('SEM_VIEW_W', 6000))
HALF = 26                        # data panel half-height in source pixels around the trace
N_PANEL = 2                      # ref panels drawn as fixed rows; the rest stay in the selector


def b64(img):
    buf = io.BytesIO()
    img.save(buf, format='PNG')
    return base64.b64encode(buf.getvalue()).decode()


# ---------------------------------------------------------------- data panel
def sem_strip(folder, qd, flip=False):
    """The SEM image resampled along the traced polyline's arc length.

    The measured length -- and therefore the bp axis every other panel uses -- is the arc
    length of this polyline, so sampling along it is what puts the photo on the same axis.
    A plain x-range crop would not do: 18 of the 53 traces double back in x, and one wanders
    1805 px in y.
    """
    j = json.load(open(f'{folder}/{qd}_points.json'))
    # measure points are the label positions; guide points only steer the tracer and must not
    # be joined into the path (13 of the 53 molecules carry them)
    types = j.get('point_types') or ['measure'] * len(j['points'])
    pts = np.array([p for p, t in zip(j['points'], types) if t == 'measure'], float)
    raw = f'{folder}/{qd}.tif'
    img = tiff.imread(raw if os.path.exists(raw) else f'{folder}/{qd}_adjusted.tif')
    if img.ndim == 3:
        img = img[..., 0] if img.shape[2] <= 4 else img[0]
    d = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0], np.cumsum(d)])
    g = np.linspace(0, s[-1], W)
    x = np.interp(g, s, pts[:, 0])
    y = np.interp(g, s, pts[:, 1])
    dx = np.gradient(x)
    dy = np.gradient(y)
    n = np.hypot(dx, dy)
    n[n == 0] = 1
    nx, ny = -dy / n, dx / n                       # unit normal
    k = np.arange(-HALF, HALF + 1)
    rr = y[None, :] + k[:, None] * ny[None, :]
    cc = x[None, :] + k[:, None] * nx[None, :]
    out = map_coordinates(img.astype(float), [rr.ravel(), cc.ravel()],
                          order=1, mode='nearest').reshape(len(k), W)
    lo, hi = np.percentile(out, 0.5), np.percentile(out, 99.5)
    out = np.clip((out - lo) / max(hi - lo, 1e-9), 0, 1) * 255
    if flip:
        out = out[:, ::-1]
    return Image.fromarray(out.astype(np.uint8)).convert('RGB')


def band_x(arr):
    """The band's drawn extent in x.

    The pure-black run gives the axis area, but the outermost label ticks are drawn a couple of
    pixels beyond it, and clipping to the black run alone drops them (QD025 loses its first and
    last of 14). Widen the range to cover any green tick as well.
    """
    c = (arr.max(axis=2) < 25).mean(axis=0)
    x = np.where(c > 0.5)[0]
    lo, hi = int(x.min()) + 1, int(x.max()) - 1
    g = ((arr[:, :, 1].astype(int) - arr[:, :, 0].astype(int)) > 60).any(axis=0)
    gx = np.where(g)[0]
    if len(gx):
        lo = min(lo, int(gx.min())); hi = max(hi, int(gx.max()))
    return lo, hi



def query_ticks(arr):
    """centres of the green label ticks, normalised to the band width."""
    x0, x1 = band_x(arr)
    g = ((arr[:, :, 1].astype(int) - arr[:, :, 0].astype(int)) > 60).any(axis=0)
    gx = np.where(g)[0]
    if not len(gx):
        return np.array([])
    cen, s0, prev = [], gx[0], gx[0]
    for v in gx[1:]:
        if v > prev + 1:
            cen.append((s0 + prev) / 2)
            s0 = v
        prev = v
    cen.append((s0 + prev) / 2)
    return (np.array(cen) - x0) / max(x1 - x0, 1)


def raw_band(png, which):
    """The band's pure-black rows only: the rows either side of it are the grey axes frame."""
    r0, r1 = BAND[which]
    arr = np.array(Image.open(png).convert('RGB'))
    win = arr[r0:r1]
    rows = np.where((win.max(axis=2) < 25).mean(axis=1) > 0.5)[0]
    return win[rows.min():rows.max() + 1]


def to_axis(arr, a=1.0, b=0.0, mirror=False):
    """Band resampled onto a common W-px axis; u_source = a*u + b, mirrored if asked.

    Columns that fall outside the band are painted black, not clamped: clamping samples the
    figure's left margin, which is where the white "Query"/"Ref" row label sits.
    """
    x0, x1 = band_x(arr)
    u = np.linspace(0, 1, W)
    us = a * u + b
    inside = (us >= 0) & (us <= 1)
    src = 1.0 - us if mirror else us
    idx = np.rint(np.clip(x0 + src * (x1 - x0), x0, x1)).astype(int)
    out = arr[:, idx].copy()
    out[:, ~inside] = 0
    return Image.fromarray(out)


def fit_axis(q_ref, q_n):
    """Least-squares u_n = a*u_ref + b from the same molecule's ticks in two figures.

    Each rank figure is drawn at its own alignment scale, so the molecule's own label track sits
    at different relative positions from figure to figure. Matching those ticks puts every rank's
    reference band back on rank 1's axis, which is the axis the Query row shows.
    """
    if len(q_ref) != len(q_n) or len(q_ref) < 2:
        return 1.0, 0.0
    a, b = np.polyfit(q_ref, q_n, 1)
    return float(a), float(b)


# ---------------------------------------------------------------- table
def norm_ref(name):
    """SEM 'Chr3_m' -> hDOM-style 'chr03_m' (same assembly, same enzyme)."""
    m = re.match(r'chr(\w+?)_(\w)$', name.strip(), re.I)
    if not m:
        return name
    c = m.group(1)
    return f'chr{c.zfill(2) if c.isdigit() else c.upper()}_{m.group(2).lower()}'


COLS = ['rank', 'ref_name', 'method', 'aligned_bp', 'dir', 'scale', 'cc_r', 'local_cc_r',
        'cc_g', 'cc_rg', 'cc_rg2', 'log(P)', 'misc', 'F1', 'm_score',
        'mal_matched', 'mal_FP', 'mal_FN', 'mal_sizing']


def mal_rows(folder, qd, limit):
    """{QD}_mal.txt -> hDOM-schema rows. query_misses = FP, ref_misses = FN.

    limit is the number of rank strip figures that exist, so the table never lists a candidate
    with no panel behind it -- _mal.txt carries 10, the figures stop at 6.
    """
    with open(f'{folder}/{qd}_mal.txt') as f:
        rd = list(csv.DictReader(f, delimiter='\t'))
    rows = []
    for i, r in enumerate(rd[:limit], 1):
        try:
            qs, qe = float(r['query_start_bp']), float(r['query_end_bp'])
            rs, re_ = float(r['ref_start_bp']), float(r['ref_end_bp'])
            matched = int(r['num_matched_chunks'])
            fp, fn = int(r['query_misses']), int(r['ref_misses'])
        except (KeyError, ValueError):
            continue
        # hDOM convention: molecule_bp * scale = ref span covered (hDOM_view: raw_w*BPP*scale)
        scale = (re_ - rs) / (qe - qs) if qe != qs else 1.0
        f1 = 2 * matched / (2 * matched + fp + fn) if (2 * matched + fp + fn) else 0.0
        rows.append([i, norm_ref(r['ref_map']), 'mal', int(rs),
                     1 if r['is_forward'].upper().startswith('F') else -1, round(scale, 2),
                     '', '', '', '', '', '', '', round(f1, 2), r.get('m_score', ''),
                     matched, fp, fn, r.get('sizing_score_rescaled', '')])
    return rows


# ---------------------------------------------------------------- page
def build_page(qd, panels, n_rank, mol_index, total, all_mols):
    items = [
        '<div class="comparison-item"><div class="image-label">SEM</div>'
        f'<div class="image-container"><img id="data-img" src="data:image/png;base64,{panels["sem"]}"></div></div>',
        '<div class="comparison-item"><div class="image-label">Query</div>'
        f'<div class="image-container"><img src="data:image/png;base64,{panels["query"]}"></div></div>',
        '<div class="data-separator"></div>',
    ]
    for i in range(1, min(n_rank, N_PANEL) + 1):
        items.append(
            f'<div class="comparison-item" data-rtype="mal"><div class="image-label">rank{i}</div>'
            f'<div class="image-container"><img src="data:image/png;base64,{panels["ref"][i]}"></div></div>')
    items.append('<div class="rank-separator"></div>')
    items.append(
        '<div class="comparison-item" id="sel-item"><div class="image-label">'
        '<span id="sel-label">rank ?</span></div><div class="image-container">'
        '<img id="sel-img" style="display:none;"><div id="sel-placeholder" class="placeholder-black"></div>'
        '</div></div>')

    opts = ''.join(f'<option value="m{i}">' for i in range(1, n_rank + 1))
    rank_imgs = ', '.join(f'"{i}": "{panels["ref"][i]}"' for i in range(1, n_rank + 1))
    label_map = ', '.join(f'"m{i}": {i}' for i in range(1, n_rank + 1))
    pct = 100.0 * mol_index / max(total, 1)
    counter = (f'<span class="jump-num" title="Click to jump to a molecule">{mol_index}</span>'
               f' / {total}')
    order = json.dumps(all_mols)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>SGM View - {qd}</title>
    <style>{get_html_style()}{V._nav_style()}</style>
</head>
<body>
    <div class="container">
        <div class="title-row"><h1>SGM Analysis Report</h1>
        <span class="title-counter">{counter}</span></div>
        <div class="progress-bar-container"><div class="progress-bar-fill" style="width:{pct:.1f}%"></div></div>
        <div class="nav-bar">
            <a id="nav-prev" class="nav-btn" href="#" title="Previous molecule (Left arrow)">&#9664;</a>
            <span class="nav-counter">{counter}</span>
            <a id="nav-next" class="nav-btn" href="#" title="Next molecule (Right arrow)">&#9654;</a>
        </div>
        <p><strong>File:</strong> {qd} ({panels['kb']:.1f} kbp, {panels['nlab']} labels)
           &nbsp;&nbsp;<a class="index-link" href="index.php">&#9776; Index</a></p>
        {toolbar()}
        <div class="comparison-box"><div class="comparison-images">{''.join(items)}</div>
        <div class="selector-box"><label for="sel-input">Select rank:</label>
        <input type="text" id="sel-input" list="sel-rank-list" autocomplete="off" style="width:7em;" placeholder="e.g. m1">
        <datalist id="sel-rank-list">{opts}</datalist>
        <span class="hint">Enter m1.. to display.</span></div></div>
        <h2>Matching Results</h2>
        <div class="table-container"><table class="dataframe matching-table"><thead><tr>
        <th>rank</th><th>ref_name</th><th>method</th><th>aligned_bp</th><th>dir</th><th>scale</th>
        </tr></thead><tbody></tbody></table></div>
    </div>
<script>
    const RANK_IMAGES = {{{rank_imgs}}};
    const LABEL_MAP = {{{label_map}}};
    const MOL_ORDER = {order};
    (function() {{
      const inp = document.getElementById('sel-input');
      const img = document.getElementById('sel-img');
      const ph = document.getElementById('sel-placeholder');
      const lab = document.getElementById('sel-label');
      function show() {{
        const k = LABEL_MAP[inp.value.trim().toLowerCase()];
        if (k && RANK_IMAGES[k]) {{
          img.src = 'data:image/png;base64,' + RANK_IMAGES[k];
          img.style.display = ''; ph.style.display = 'none';
          lab.textContent = 'rank ' + k;
        }} else {{ img.style.display = 'none'; ph.style.display = ''; lab.textContent = 'rank ?'; }}
      }}
      inp.addEventListener('change', show); inp.addEventListener('input', show);

      // n / N -> jump to the n-th molecule of this folder
      document.querySelectorAll('.jump-num').forEach(function(el) {{
        el.addEventListener('click', function() {{
          const v = prompt('Jump to molecule (1-' + MOL_ORDER.length + ')', el.textContent.trim());
          const n = parseInt(v, 10);
          if (n >= 1 && n <= MOL_ORDER.length) {{
            location.href = 'index.php?file=' + encodeURIComponent(MOL_ORDER[n - 1]);
          }}
        }});
      }});

      document.addEventListener('keydown', function(e) {{
        const t = (e.target.tagName || '').toLowerCase();
        if (t === 'input' || t === 'textarea') return;
        if (e.key === 'ArrowLeft') {{ const a = document.getElementById('nav-prev'); if (a && a.href) location.href = a.href; }}
        if (e.key === 'ArrowRight') {{ const a = document.getElementById('nav-next'); if (a && a.href) location.href = a.href; }}
      }});
    }})();
</script>
</body>
</html>
"""


# ---------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser(description='SEM molecules in the standard align-page frame')
    ap.add_argument('--data', required=True, help='directory holding the <QD###>/ folders')
    ap.add_argument('--out', required=True, help='directory to write the view into')
    ap.add_argument('--template-index', default='',
                    help='index.php to copy next to the views (default: none)')
    ap.add_argument('--only', nargs='*', help='subset of QD names')
    a = ap.parse_args()
    global SRC, OUT, TEMPLATE_INDEX
    SRC, OUT, TEMPLATE_INDEX = a.data, a.out, a.template_index
    os.makedirs(a.out, exist_ok=True)

    mols = a.only or sorted(
        os.path.basename(d) for d in glob.glob(f'{SRC}/QD*')
        if os.path.isdir(d) and glob.glob(f'{d}/{os.path.basename(d)}_rank1_strip.png'))
    summary, lengths = [], {}

    for i, qd in enumerate(mols, 1):
        folder = f'{SRC}/{qd}'
        try:
            strips = sorted(glob.glob(f'{folder}/{qd}_rank*_strip.png'),
                            key=lambda p: int(re.search(r'rank(\d+)_strip', p).group(1)))
            rows = mal_rows(folder, qd, len(strips))
            # the strip figure draws the query in the reference's orientation, so a reverse
            # alignment (dir = -1) mirrors it against the traced order; verified on all 53.
            flip = bool(rows) and rows[0][4] == -1
            q1 = raw_band(strips[0], 'query')
            qt1 = query_ticks(q1)
            ref = {}
            d1 = rows[0][4] if rows else 1
            for n, sp in enumerate(strips, 1):
                # each rank figure is drawn in ITS OWN reference orientation, so a rank whose
                # dir differs from rank 1's arrives mirrored against the axis Query is shown on
                mirror = bool(rows) and n <= len(rows) and rows[n - 1][4] != d1
                qn = query_ticks(raw_band(sp, 'query'))
                if mirror:
                    qn = 1.0 - qn[::-1]
                fa, fb = fit_axis(qt1, qn)
                ref[n] = b64(to_axis(raw_band(sp, 'ref_sm'), fa, fb, mirror))
            panels = {'sem': b64(sem_strip(folder, qd, flip)),
                      'query': b64(to_axis(q1)), 'ref': ref}
            total_bp, frags = 0, 0
            with open(f'{folder}/{qd}.maps') as f:
                fl = f.readline().split('\t')
                total_bp, frags = int(fl[1]), int(fl[2])
            panels['kb'] = total_bp / 1000.0
            panels['nlab'] = max(frags - 1, 0)

            with open(f'{a.out}/{qd}_matching_results.csv', 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(COLS)
                w.writerows(rows)

            open(f'{a.out}/{qd}_view_ytb.html', 'w').write(
                build_page(qd, panels, len(strips), i, len(mols), mols))

            top = rows[0] if rows else ['', ''] * 3
            summary.append([qd, top[14] if rows else '', '', '', top[1] if rows else '', '', '-', ''])
            lengths[qd] = round(total_bp / 1000.0)
            print(f'[{i:3d}/{len(mols)}] {qd:10s} ranks={len(strips)} '
                  f'{total_bp/1000:7.1f} kb  labels={panels["nlab"]:3d}  mal rows={len(rows)}', flush=True)
        except Exception as e:
            print(f'[{i:3d}/{len(mols)}] {qd:10s} FAIL {type(e).__name__}: {e}', flush=True)

    with open(f'{a.out}/_summary.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['molecule', 'score', 'local1', 'mDOM1', 'mal1', 'pcc1', 'H', 'mdom_class'])
        w.writerows(summary)
    json.dump(lengths, open(f'{a.out}/_lengths.json', 'w'))
    if not os.path.exists(f'{a.out}/index.php'):
        import shutil
        if TEMPLATE_INDEX:
            shutil.copy2(TEMPLATE_INDEX, f'{a.out}/index.php')
    print(f'\n{len(summary)}/{len(mols)} -> {a.out}')


if __name__ == '__main__':
    main()
