#!/usr/bin/env python3
"""Unified straighten+align pipeline (v5 'full_snap' logic, parameterized).

For each {mol}_orig.tif (2-frame RAW molecule, frame0=red DNA / frame1=green nick) in --input:
  detect ML labels -> backbone -> sweep(arc-length straighten) -> snap(labels to swept green
  peaks +-4px) -> widthprop(green binarize, height 846)  ==> {mol}.tif (align) + {mol}_orig.tif
  (display: straightened RAW real green)
  -> hDOM_calc_Q candidates -> local_align_pipeline.rerank (local_cc_r x F1 -> null log P)
  -> hDOM_view view  -> inject data(raw)+Linearized before/after panels (default on)
  -> deploy {mol}.{tif,_orig.tif,_matching_results.csv,_mlpeaks.json,_view.html} to --output

This is the single consolidated tool: replaces batch_sweep_full_snap/xy/cr (they differed ONLY
in the input/output paths). The data(raw) panel uses the INPUT {mol}_orig.tif as the raw source.

Usage:
  python3 straighten_align.py --input <raw_folder> --output <dest> [--shard k/N] [--no-before-after]
  # parallel: run N shards, e.g.  --shard 0/4, --shard 1/4, ...
"""
import os
os.environ['OPT_COMBINED'] = '1'
os.environ.setdefault('MATCH_WORKERS', '1'); os.environ.setdefault('MALIGNER_WORKERS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys, csv, glob, json, shutil, subprocess, argparse, warnings, numpy as np, tifffile
warnings.filterwarnings('ignore')
_PKG = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, _PKG)
import local_align_pipeline as P
import Label_Detect_ML_based as ML
from scipy.ndimage import map_coordinates, median_filter
from scipy.signal import peak_widths, find_peaks as _find_peaks
from hDOM_tif2html_Y import process_query_image

CALC = os.path.join(_PKG, 'candidates', 'hDOM_calc_Q.py')
VIEW = os.path.join(_PKG, 'hDOM_view.py')
C_AMP, FWHM2SIG, BPP = P.C_AMP, P.FWHM2SIG, P.BPP
TRIM_PX = P.TRIM_PX          # terminal 5 kb / BPP -- same window hDOM scores with

ap = argparse.ArgumentParser(description='Unified straighten+align (v5 snap logic, parameterized)')
ap.add_argument('--input', required=True, help='folder with {mol}_orig.tif (raw 2-frame molecules)')
ap.add_argument('--output', required=True, help='destination folder (deploy straightened+aligned+view)')
ap.add_argument('--shard', default='0/1', help='k/N shard for parallelism (default 0/1 = all)')
ap.add_argument('--no-before-after', action='store_true', help='skip data(raw)+Linearized 2-panel injection')
ap.add_argument('--logp', action='store_true', help='compute K=5000 random-null log(P) (h_score) for local rows + {base}_vsrandom.txt sidecar')
args = ap.parse_args()
B, V4 = args.input, args.output
shk, shn = (int(x) for x in args.shard.split('/'))
BEFORE_AFTER = not args.no_before_after
LOGP = args.logp
os.makedirs(V4, exist_ok=True)

CACHE = {}
def _key(f):
    """Molecule name from a tif path, with either `_orig.tif` or `.tif` stripped.

    The straightened image and its raw display copy share one name, and the ML label cache is keyed
    on that name, so both spellings have to collapse to the same key.
    """
    n = os.path.basename(f)
    for s in ('_orig.tif', '.tif'):
        if n.endswith(s): n = n[:-len(s)]; break
    return n
P.refined_green_peaks = lambda f: np.asarray(CACHE[_key(f)], dtype=int)

def backbone(red, win=4, reacq=20, maxstep=None, min_run=10):
    """Trace the DNA backbone: for each image column, the intensity centroid of the red strand.

    Walking outward in both directions from the brightest column, each step takes the centroid
    within +-win rows of where the strand sat in the previous column, which keeps the trace on the
    molecule instead of on a neighbouring speck. A column whose window holds nothing above
    threshold is retried once at +-reacq rows: a sharp terminal bend or a locally faint stretch
    would otherwise end the trace early and truncate the molecule.

    The threshold is set from the strand's own brightness -- the median of the upper 60% of column
    maxima -- rather than a fixed level, so dim and bright molecules are traced the same way.

    With maxstep set (the BBFIX pass, see chosen_backbone), a column whose centroid sits more than
    maxstep rows from the previous column is dropped instead of accepted: a DNA backbone in these
    crops rises 0.02-0.04 px per column, so a 2 px step means the tracer left the strand for a
    nearby speck and would otherwise drag the straightened image into a false hump. maxstep=None
    reproduces the original trace exactly.

    Returns one y per column, NaN where the strand was never acquired.
    """
    H, W = red.shape; r = np.clip(red - np.median(red), 0, None); cp = r.max(0)
    pos = cp[cp > 0]; la = np.median(pos[pos > np.percentile(pos, 40)]) if len(pos) else cp.max()
    thr = max(la * 0.15, r.max() * 0.02); x0 = int(np.argmax(cp))
    def cw(x, yp, wn=win):
        """Intensity centroid within +-wn rows of the previous column's y, or None if nothing is there."""
        lo, hi = max(0, int(round(yp)) - wn), min(H, int(round(yp)) + wn + 1); w = r[lo:hi, x]
        if w.sum() <= 0 or w.max() < thr: return None
        return lo + (np.arange(hi - lo) * w).sum() / w.sum()
    y = np.full(W, np.nan); y0 = cw(x0, int(np.argmax(r[:, x0]))) or float(np.argmax(r[:, x0]))
    for o in (range(x0, W), range(x0, -1, -1)):
        yp = y0
        for x in o:
            c = cw(x, yp)
            # curve/dim-gap re-acquire: if the +-win window lost the strand, retry with a wider
            # (+-reacq) window (thr-gated, bounded jump). Molecule crops place DNA near the image
            # edge, so a sharp terminal bend or a faint stretch shouldn't truncate the backbone.
            if c is None and reacq > win: c = cw(x, yp, wn=reacq)
            if maxstep is not None and c is not None and abs(c - yp) > maxstep: c = None
            if c is not None: y[x] = c; yp = c
    if maxstep is None: return y
    # A single survivor inside a rejected run would anchor the interpolation and pull the trace off
    # the strand, so drop points far from a local median. Wide window first, so a narrow excursion
    # cannot drag the median toward itself.
    for wm in (101, 61, 31):
        v = np.where(~np.isnan(y))[0]
        if len(v) < 10: break
        mf = median_filter(np.interp(np.arange(W), v, y[v]), size=wm, mode='nearest')
        bad = (~np.isnan(y)) & (np.abs(y - mf) > maxstep)
        if bad.any(): y[bad] = np.nan
    # Short surviving islands between rejected runs are re-acquisition residue.
    ok = ~np.isnan(y); i = 0
    while i < W:
        if not ok[i]: i += 1; continue
        j = i
        while j < W and ok[j]: j += 1
        gl = 0; k = i - 1
        while k >= 0 and not ok[k]: gl += 1; k -= 1
        gr = 0; k = j
        while k < W and not ok[k]: gr += 1; k += 1
        if (j - i) < min_run and (gl >= j - i or gr >= j - i) and (gl + gr) > 0: y[i:j] = np.nan
        i = j
    return y


# ---------- straightening QC: repair traces that produce a false hump ----------
# The tracer can step off the strand onto a nearby speck where the red signal is locally faint, and
# arc-length resampling then bakes that excursion into the straightened image as a "speed bump" the
# molecule never had. BBFIX (backbone(..., maxstep=...)) removes those steps, but applying it to
# every molecule would perturb ~75% of them for no visible gain, so it is used only where the defect
# is actually present. The test is made on the straightening this run produced, not on any stored
# result: straighten once, measure how far the DNA band wanders vertically in that image (ydev), and
# how far the BBFIX trace disagrees with the one used (d99). Their product is large only when a real
# excursion has been baked in. A pipeline re-run therefore judges its own output -- a first pass with
# different defects legitimately selects a different set of molecules.
BBFIX_MAXSTEP = float(os.environ.get('BBFIX_MAXSTEP', '2.0'))    # max backbone move per column (px)
BBFIX_MINRUN  = int(os.environ.get('BBFIX_MINRUN', '10'))        # shortest surviving island kept
BBFIX_THRESH  = float(os.environ.get('BBFIX_THRESH', '40'))      # d99 * ydev above which BBFIX is used
BBFIX_ON      = os.environ.get('BBFIX', '1') not in ('0', 'no', 'off')

def _span_fill(y):
    """Traced span plus the trace linearly interpolated across its internal gaps."""
    v = np.where(~np.isnan(y))[0]
    if len(v) < 20: return None, None
    xs = np.arange(v.min(), v.max() + 1)
    return xs, np.interp(xs, v, y[v])

def _trace_gap(y0, y1):
    """99th-percentile vertical disagreement between two traces over their common span (px)."""
    x0, f0 = _span_fill(y0); x1, f1 = _span_fill(y1)
    if f0 is None or f1 is None: return 0.0
    lo, hi = max(x0[0], x1[0]), min(x0[-1], x1[-1])
    if hi - lo < 50: return 0.0
    g = np.arange(lo, hi + 1)
    return float(np.percentile(np.abs(np.interp(g, x0, f0) - np.interp(g, x1, f1)), 99))

def _band_dev(swept_red):
    """How far the DNA band wanders from its own smooth trend in a straightened image (px, p99).

    A correct straightening leaves the strand on a level band; a trace excursion tilts the sampling
    normals there, so the band's intensity centroid rides up or down over that stretch. Measured
    against a 151-column median so genuine slow drift does not count.
    """
    r = np.clip(swept_red.astype(float) - np.median(swept_red), 0, None); tot = r.sum(0)
    if not (tot > 0).any(): return 0.0
    keep = tot > np.percentile(tot[tot > 0], 10); rows = np.arange(r.shape[0])
    c = np.full(r.shape[1], np.nan); m = tot > 0
    c[m] = (r[:, m] * rows[:, None]).sum(0) / tot[m]; c[~keep] = np.nan
    v = np.where(~np.isnan(c))[0]
    if len(v) < 60: return 0.0
    gx = np.arange(v.min(), v.max() + 1); cv = np.interp(gx, v, c[v])
    cs = np.array([cv[max(0, i - 4):i + 5].mean() for i in range(len(cv))])
    return float(np.percentile(np.abs(cs - median_filter(cs, size=151, mode='nearest')), 99))

SWEEP_SM = int(os.environ.get('SWEEP_SM', '6'))   # backbone smoothing half-window (columns)
def sweep(a, ml, sm=None):
    """Straighten the molecule by resampling along the backbone's arc length.

    The traced backbone is smoothed over +-sm columns and its arc length s(x) accumulated; the
    output grid is uniform in s, so a curved molecule comes out straight and, at the same time,
    comes out LONGER in pixels than its horizontal extent. That length change is the point: the
    curved image under-counts the molecule, and the scale the aligner then fits absorbs the error.
    Each output column is sampled along the local normal, so the strand's cross-section is carried
    over rather than smeared.

    Both frames are warped with the SAME coordinates, so red and green stay co-registered.

    Label positions come along by the same map -- x -> s(x) -- but note that np.interp clamps
    anything outside the traced span [x0, x1] onto an end pixel; those terminal labels are what
    TRIM_PX later drops.

    Returns (2-frame straightened image, label positions in the new grid, (x0, x1) traced span),
    or None if the trace was too short to be usable.
    """
    if sm is None: sm = SWEEP_SM
    red = a[0].astype(float)
    r0 = _warp(a, backbone(red), ml, sm)
    if r0 is None or not BBFIX_ON: return r0
    yb = backbone(red, maxstep=BBFIX_MAXSTEP, min_run=BBFIX_MINRUN)
    d = _trace_gap(backbone(red), yb)
    if d <= 0: return r0
    if d * _band_dev(r0[0][0]) < BBFIX_THRESH: return r0
    return _warp(a, yb, ml, sm) or r0

def _warp(a, y, ml, sm):
    """Resample both frames along the arc length of trace y (the body of sweep, given a trace)."""
    H, Wd = a[0].shape; v = np.where(~np.isnan(y))[0]
    if len(v) < 10: return None
    x0, x1 = v.min(), v.max(); xs = np.arange(x0, x1 + 1); yv = y[xs]; m = ~np.isnan(yv); yv = np.interp(xs, xs[m], yv[m])
    ys = np.array([yv[max(0, i - sm):i + sm + 1].mean() for i in range(len(xs))]); dy = np.gradient(ys)
    s = np.concatenate([[0], np.cumsum(np.sqrt(1 + dy**2))[:-1]]); Stot = s[-1]
    Wp = int(np.ceil(Stot)) + 1; sg = np.linspace(0, Stot, Wp); xf = np.interp(sg, s, xs); yf = np.interp(xf, xs, ys)
    dyf = np.interp(xf, xs, dy); Lf = np.sqrt(1 + dyf**2); nr, nc = 1.0 / Lf, -dyf / Lf; Bb = (H - 1) // 2; ks = np.arange(-Bb, Bb + 1)
    RR = yf[None, :] + ks[:, None] * nr[None, :]; CC = xf[None, :] + ks[:, None] * nc[None, :]
    out = np.zeros((2, len(ks), Wp), np.float32)
    for f in (0, 1): out[f] = map_coordinates(a[f].astype(float), [RR.ravel(), CC.ravel()], order=1, mode='constant', cval=0).reshape(len(ks), Wp)
    return out.astype(np.uint16), [float(np.interp(x, xs, s) / Stot * (Wp - 1)) for x in ml], (x0, x1)

def chosen_backbone(red, sm=None):
    """The trace sweep() would actually use for this molecule (stock, or BBFIX if the QC test fires).

    Kept in step with sweep() so the before/after panel crops the raw image over the same span that
    was straightened; a mismatch would misreport the unfolding as a length change.
    """
    if sm is None: sm = SWEEP_SM
    y = backbone(red)
    if not BBFIX_ON: return y
    a = np.stack([red, np.zeros_like(red)])
    r0 = _warp(a, y, [], sm)
    if r0 is None: return y
    yb = backbone(red, maxstep=BBFIX_MAXSTEP, min_run=BBFIX_MINRUN)
    d = _trace_gap(y, yb)
    if d <= 0 or d * _band_dev(r0[0][0]) < BBFIX_THRESH: return y
    return yb

def snap(oc, sweptgreen):
    """Move each warped label onto the nearest real green peak in the straightened image (+-4 px).

    Warping carries a label by interpolating its old coordinate, which can leave it a pixel or two
    off the peak it belongs to. Snapping puts it back on the intensity maximum. A label with no peak
    within +-4 px keeps its warped position rather than being discarded, so the label COUNT is
    preserved -- dropping one here would silently change F1 and the maligner query map.
    """
    g = np.max(sweptgreen, 0).astype(float)
    pk, _ = _find_peaks(g, height=g.max() * 0.25, distance=3)
    out = []
    for c in oc:
        cand = [p for p in pk if abs(p - c) <= 4]
        out.append(float(min(cand, key=lambda p: abs(p - c))) if cand else c)
    return out

def measure(g2, xs):
    """FWHM of each label in the straightened green frame, floored at 1.5 px.

    Fed to widthprop so the synthetic label keeps the width actually measured on this molecule.
    The floor, and the 3.0 px fallback when the fit fails, keep a degenerate peak from producing a
    delta function the correlation could not see.
    """
    g = np.max(g2, 0).astype(float); o = []
    for x in xs:
        try: w = float(peak_widths(g, [int(x)], rel_height=0.5)[0][0])
        except Exception: w = 3.0
        o.append(max(w, 1.5))
    return o

def widthprop(sr, oc, fw):
    """Rebuild the green frame as Gaussians: measured width, one fixed height (C_AMP).

    Green enters the pcc score through a cosine that does not subtract the mean, so a bright label
    counts for more than a faint one at the same position -- brightness, which carries no positional
    information, ends up steering the ranking. Redrawing every label at the same amplitude removes
    that, while the per-label width measured on the real image is kept, since width does reflect how
    well the site is localised.

    Only the green frame is synthetic; red stays the measured image.
    """
    Wp = sr.shape[1]; H = sr.shape[0]; g = np.zeros((H, Wp), np.float32); cols = np.arange(Wp); rows = np.arange(H); Bb = H // 2
    for c, f in zip(oc, fw):
        sig = max(f, 1e-3) / FWHM2SIG
        g = np.maximum(g, C_AMP * np.exp(-((cols - c)**2) / (2 * sig**2))[None, :] * np.exp(-((rows - Bb)**2) / (2 * 2.0**2))[:, None])
    return g

# ---------- [1] straighten ----------
mols = sorted(os.path.basename(f)[:-len('_orig.tif')] for f in glob.glob(f'{B}/*_orig.tif'))
mols = [k for k in mols if not k.endswith(('_L', '_R'))]
mine = [k for i, k in enumerate(mols) if i % shn == shk]
ST = f'{V4}/.sa_{shk}'; shutil.rmtree(ST, ignore_errors=True); os.makedirs(ST, exist_ok=True)
if not os.path.exists(f'{ST}/ref_map'): os.symlink(P.REFMAP, f'{ST}/ref_map')
rawmol = {}; staged = []; removed = {}; kept = {}
for k in mine:
  # One unusable molecule must not take the rest of the shard with it: a shard is a slice of
  # the input list, so an exception here would drop every molecule that happens to share it.
  try:
    src = f'{B}/{k}_orig.tif'
    if not os.path.exists(src): continue
    a = tifffile.imread(src).astype(float)
    if a.ndim != 3 or a.shape[0] < 2 or min(a.shape[1:]) < 1:
        print(f'  unusable image {k} (shape {a.shape}; need 2 frames with non-zero size)',
              flush=True); continue
    # HDOM_LABEL_EDIT: a curated label list for this molecule replaces ML detection. Positions are
    # RAW-frame pixels (the same frame detect_labels returns), so straightening maps them onward
    # exactly as it maps detected ones -- a hand-checked molecule differs from an automatic one only
    # in where the labels came from.
    ml = None
    _ed = os.environ.get('HDOM_LABEL_EDIT', '')
    if _ed:
        _ef = os.path.join(_ed, f'{k}_labeledit.json')
        if os.path.exists(_ef):
            try:
                _j = json.load(open(_ef, encoding='utf-8'))
                ml = sorted(int(round(float(x))) for x in _j.get('final', []))
                print(f'  [labeledit] {k}: {len(ml)} labels', flush=True)
            except Exception as e:
                print(f'  labeledit read fail {k}: {e}', flush=True); ml = None
    if ml is None:
        ml = sorted(int(x) for x in ML.detect_labels(src))
    if not ml: print(f'  no ML {k}', flush=True); continue
    r = sweep(a, ml)
    if r is None: print(f'  sweep fail {k}', flush=True); continue
    st, oc, ext = r; oc = snap(oc, st[1]); fw = measure(st[1], [int(round(c)) for c in oc])
    wp = np.stack([st[0], widthprop(st[0], oc, fw)]).astype(np.uint16)
    tifffile.imwrite(f'{ST}/{k}.tif', wp)              # align (widthprop green)
    tifffile.imwrite(f'{ST}/{k}_orig.tif', st)         # display: straightened RAW (real green)
    CACHE[k] = sorted(int(round(c)) for c in oc)
    json.dump(CACHE[k], open(f'{ST}/{k}_mlpeaks.json', 'w', encoding='utf-8'))
    # Labels inside the terminal 5 kb. Candidate generation (pcc's cc_r*cc_g**2 and the maligner
    # run inside hDOM_calc_Q) uses the FULL label set above; these are dropped only for the
    # second, reported maligner pass in [2b] and for the white-circle marking in the view.
    # Under 2 survivors there is no alignable map, so that molecule gets no trimmed pass at all.
    w2 = wp.shape[-1]
    # The reported maligner pass uses the same terminal rule as F1. If the two rules diverge, the
    # mal row and the local row in the view are computed from different label sets, and the same
    # molecule ends up described two different ways.
    _tl, _tr = P.terminal_cuts(np.asarray(CACHE[k], float), w2, TRIM_PX)
    kept[k] = [p for p in CACHE[k] if _tl <= p <= w2 - 1 - _tr]
    removed[k] = sorted(set(CACHE[k]) - set(kept[k])) if len(kept[k]) >= 2 else []
    rawmol[k] = a; staged.append(k)
  except Exception as e:
    print(f'  straighten fail {k}: {type(e).__name__}: {e}', flush=True)
print(f'[1] shard{shk} straightened+staged {len(staged)}/{len(mine)}', flush=True)

# ---------- [2] candidates + rerank ----------
# hDOM_calc_Q scans every .tif in the folder it is given, and stage [1] wrote two per molecule:
# {mol}.tif, whose green is redrawn at one fixed amplitude, and {mol}_orig.tif, which keeps the
# measured green for display. Only the first is meant to be aligned -- the fixed amplitude exists
# precisely because label brightness carries no positional information and distorts the score. So
# aligning the display copy searches the whole genome with the very signal that was removed, and
# its output was being deleted immediately afterwards. Holding those files aside for the duration
# of the call halves the candidate-generation work and changes nothing that is kept.
_hold = f'{ST}/_display'
os.makedirs(_hold, exist_ok=True)
for _f in glob.glob(f'{ST}/*_orig.tif'):
    shutil.move(_f, os.path.join(_hold, os.path.basename(_f)))
subprocess.run(['python3', CALC, ST, '--ref-map', P.REFMAP, '--bpp', str(BPP)], check=False)
for _f in glob.glob(f'{_hold}/*_orig.tif'):
    shutil.move(_f, os.path.join(ST, os.path.basename(_f)))
shutil.rmtree(_hold, ignore_errors=True)
# Positions the reported maligner pass dropped; the viewer marks them with a white circle.
def _stem(k):
    """Drop the `_samevec` suffix so the sidecar key matches the deployed molecule name."""
    return k[:-len('_samevec')] if k.endswith('_samevec') else k
json.dump({_stem(k): v for k, v in removed.items() if v}, open(f'{ST}/removed_peaks.json', 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)

# ---------- [2b] reported maligner: same tif, terminal-5 kb labels dropped ----------
# hDOM_calc_Q has already run maligner on the FULL label set, and those hits stay in the
# candidate pool that the local re-rank scores together with pcc's. The maligner row REPORTED
# for the molecule is a second run with the terminal labels dropped: straightening maps labels
# through np.interp over the traced backbone [x0,x1], so a label outside it clamps onto pixel 0
# or w2-1 and becomes a spurious query fragment maligner is forced to explain. Both runs are
# kept on disk (_maligner_dp.txt = candidates, _maligner_dp_trimmed.txt = reported).
sys.path.insert(0, os.path.dirname(CALC))
import hDOM_mal_Q as MAL
MALBIN = os.environ.get('MALIGNER_DP') or shutil.which('maligner_dp') or 'maligner_dp'
GENOME = os.path.join(P.REFMAP, 'genome.maps')

def trim_mal(k):
    """maligner_dp on the trimmed label set -> top-10 by total_rescaled_score (ascending=better)."""
    if len(kept.get(k, [])) < 2: return None      # not an alignable map; molecule keeps the full run
    wd = f'{ST}/_tm_{k[:24]}'
    shutil.rmtree(wd, ignore_errors=True); os.makedirs(wd)
    try:
        os.symlink(P.REFMAP, f'{wd}/ref_map')
        shutil.copy2(f'{ST}/{k}.tif', f'{wd}/{k}.tif')
        json.dump(kept[k], open(f'{wd}/{k}_mlpeaks.json', 'w', encoding='utf-8'))
        q = MAL.create_maps_from_tif(f'{wd}/{k}.tif', BPP)
        r = subprocess.run([MALBIN, '--max-alignments', '3000',
                            '--max-alignments-per-reference', '3', '--query-max-misses', '3',
                            '--ref-max-misses', '6', '--min-query-rescaling', '0.82',
                            '--max-query-rescaling', '1.18', q, GENOME],
                           capture_output=True, text=True, timeout=3600)
        ls = [l for l in r.stdout.split('\n') if l.strip()]
        if len(ls) < 2: return None
        h = {c: i for i, c in enumerate(ls[0].split('\t'))}
        al = sorted((l.split('\t') for l in ls[1:]),
                    key=lambda f: float(f[h['total_rescaled_score']]))
        out = f'{ST}/{k}_maligner_dp_trimmed.txt'
        with open(out, 'w', encoding='utf-8') as fh:
            fh.write(ls[0] + '\n')
            for f in al[:10]: fh.write('\t'.join(f) + '\n')
        return out
    except Exception as e:
        print(f'  trim-mal fail {k}: {e}', flush=True); return None
    finally:
        shutil.rmtree(wd, ignore_errors=True)
ok = []
import figure3_view as F3V   # Figure3 CSV (native mal cols from re-run-on-straightened maligner) + viewer
if LOGP:
    import logp_null as LP; LP.load_refs()
    print(f'[logp] refs loaded ({len(LP.ALL)} chr), K={LP.K} null/molecule', flush=True)
for k in staged:
    cf = f'{ST}/{k}_matching_results.csv'
    if not os.path.exists(cf): continue
    try: rows = P.rerank(f'{ST}/{k}.tif', cf)
    except Exception as e: print(f'  rerank fail {k}: {e}', flush=True); continue
    if rows is None: continue
    with open(cf, 'w', newline='', encoding='utf-8') as f: w = csv.writer(f); w.writerow(P.COLS); w.writerows(rows)
    # mal rows come from the trimmed pass; molecules with <2 surviving labels fall back to the
    # full-label run so they still get a maligner row rather than a blank one.
    F3V.to_figure3_csv(cf, trim_mal(k) or f'{ST}/{k}_maligner_dp.txt')   # 14->19 col, native mal
    if LOGP:
        # rerank has just rewritten cf, so its log(P) holds null_logp; force restores the K=5000 null.
        print(f'  [logp] {k[:40]} -> {LP.compute_logp(k, ST, force=True)}', flush=True)
    ok.append(k)
print(f'[2] shard{shk} reranked {len(ok)}', flush=True)

# ---------- [3] view + before/after + deploy ----------
subprocess.run(['python3', VIEW, ST, '--ref-map', P.REFMAP, '--label-dots'], check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
ANCH = '<div class="comparison-images"><div class="comparison-item"><div class="image-label">data</div>'
dep = 0
for k in ok:
    vf = f'{ST}/{k}_view.html'
    if BEFORE_AFTER and os.path.exists(vf):
        h = open(vf, encoding='utf-8').read()
        if ANCH in h:
            rs = rawmol[k]; y = chosen_backbone(rs[0].astype(float)); v = np.where(~np.isnan(y))[0]
            x0, x1 = (int(v.min()), int(v.max())) if len(v) else (0, rs.shape[-1] - 1)
            tp = f'{ST}/_raw_{k[:12]}.tif'; tifffile.imwrite(tp, rs[:, :, x0:x1 + 1].astype(np.uint16)); _, rb = process_query_image(tp); os.remove(tp)
            # length-honest display: raw(curved bbox width = x1-x0) and Linearized(arc-length = straightened
            # tif width) are drawn at a COMMON px-scale (raw scaled to raw_w/lin_w), center-aligned, so the
            # unfolding (straightened longer than the curved bounding box) is VISIBLE instead of both being
            # squeezed to equal container width by max-width:100%.
            raw_w = x1 - x0 + 1
            try: lin_w = tifffile.imread(f'{ST}/{k}.tif').shape[-1]
            except Exception: lin_w = raw_w
            pct = min(100.0, 100.0 * raw_w / lin_w) if lin_w else 100.0
            pan = ('<div class="comparison-item"><div class="image-label" style="color:#e67e22">data(raw)</div>'
                   '<div class="image-container" style="position:relative"><img class="rawcurved" src="data:image/png;base64,' + rb + '" style="image-rendering:pixelated"></div></div>')
            h = h.replace(ANCH, '<div class="comparison-images">' + pan + '<div class="comparison-item"><div class="image-label" style="color:#27ae60">Linearized</div>', 1)
            css = '<style id="lenfix">img.rawcurved{max-width:%.1f%% !important;}</style>' % pct
            h = h.replace('</head>', css + '</head>', 1) if '</head>' in h else css + h
            open(vf, 'w', encoding='utf-8').write(h)
    # _maligner_dp.txt is produced inside the staging dir, which is removed below, so it is
    # copied out with the rest: it holds the maligner metrics reported for the molecule.
    for suf in ('_matching_results.csv', '.tif', '_orig.tif', '_view.html', '_mlpeaks.json',
                '_vsrandom.txt', '_maligner_dp.txt', '_maligner_dp_trimmed.txt'):
        s = f'{ST}/{k}{suf}'
        if os.path.exists(s): shutil.copy2(s, f'{V4}/{k}{suf}')
    dep += 1
import figure3_view as F3V   # Figure3-style viewer + _lengths.json (length/lin. panel labels, shard-safe merge)
# removed_peaks.json is folder-scoped, so merge instead of overwrite (other shards write it too).
# Clearing the entry when nothing was trimmed keeps a re-run from leaving stale white circles.
_rp = {}
if os.path.exists(f'{V4}/removed_peaks.json'):
    try: _rp = json.load(open(f'{V4}/removed_peaks.json', encoding='utf-8'))
    except Exception: _rp = {}
for k in ok:
    if removed.get(k): _rp[_stem(k)] = removed[k]
    else: _rp.pop(_stem(k), None)
json.dump(_rp, open(f'{V4}/removed_peaks.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
F3V.add_lengths(V4, {k: F3V.mol_len_kb(f'{V4}/{k}.tif') for k in ok})
F3V.add_summary(V4, {k: F3V.summary_row(k, f'{V4}/{k}_matching_results.csv') for k in ok})   # list/nav
F3V.deploy_viewer(V4)
shutil.rmtree(ST, ignore_errors=True)
print(f'[3] shard{shk} DONE deployed {dep}  -> {V4}', flush=True)
