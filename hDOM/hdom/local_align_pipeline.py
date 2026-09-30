#!/usr/bin/env python3
"""Integrated local-align pipeline: raw molecule TIFF -> 4-method comparison result (+ HTML view later).

Purpose: compare LOCAL align against the other tools (mDOM / maligner / pcc) on the same molecule.

Stages:
  [1] green binarization = widthprop (height fixed 846, per-peak width = measured original FWHM)
  [2] candidates          = hDOM_calc_Q  (pcc cc_rg2 top-10  +  maligner top-10)     [global cc_r only finds loci]
  [3] rank (4 methods, each ranked independently over the SAME candidate pool):
        local : LOCAL cc_r (10kb tiles / 5kb / +-3kb) x F1     -> rows l1,l2,...
        mDOM  : global cc_r x F1^2   (pcc candidates)          -> rows d1,d2,...
        mal   : maligner m_score (kept from Q)                 -> rows m1,m2,...
        pcc   : cosine cc_r x cc_g^2 (kept from Q)             -> rows p1,p2,...
  [4] homolog H row = homolog of LOCAL rank1 (l1) if l1 not sex chr AND l1/l2 not a homolog pair (C_term style)
  [5] score column = per-method log(P) : local/mDOM/H -> per-molecule null ; pcc -> cosine null (Q) ; mal -> m_score (Q)
  [6] view = (built next)  -> CSV + origin.json + _summary.csv for now

Original hDOM_calc_Q.py is NEVER modified (invoked as subprocess; ref_map symlink lets its maligner find genome.maps).

Usage:
  python3 local_align_pipeline.py <raw_tif | folder> <out_folder> [--tag NAME] [--no-view]
"""
import os, sys, csv, glob, shutil, subprocess, argparse, json
import numpy as np
import tifffile as tiff
from scipy import signal, stats

CODE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODE)
import DOM_lib_Y as lib
import mdom_reanchor as R
from hDOM_mal_Y import refined_green_peaks
from paf_liftover import PafLiftover
import pandas as pd
import warnings; warnings.filterwarnings('ignore')

# calc_Q candidates with --ml-labels: the maligner map builder (candidates/hDOM_mal_Q.create_maps_from_tif)
# reads an authoritative {base}_mlpeaks.json sidecar and overrides its fixed-threshold peaks with the EXACT
# integer ML label positions; pcc reads the widthprop green channel (ML peaks baked as Gaussians). Production
# hDOM_Q lacks that sidecar override, so ML mode MUST use CALC_Q_ML (and write the STAGE sidecar).
CALC_Q_PROD = os.path.join(CODE, 'candidates', 'hDOM_calc_Q.py')
CALC_Q_ML   = os.path.join(CODE, 'candidates', 'hDOM_calc_Q.py')  # hDOM_mal_Q reads {base}_mlpeaks.json
CALC_Q   = CALC_Q_PROD           # default; --ml-labels switches to CALC_Q_ML
MLPATH   = CODE  # Label_Detect_ML_based lives here
REFMAP   = os.environ.get('HDOM_REFMAP', os.path.join(os.path.dirname(CODE), 'ref_map'))
CHAINS   = os.environ.get('HDOM_CHAINS', os.path.join(os.path.dirname(CODE), 'chains'))
BPP      = int(os.environ.get('HDOM_BPP', 200))   # bases per reference pixel; must match the ref_map
# SCALE_MIN/MAX below stay ratios, not bpp: `scale` is how far the molecule is stretched relative
# to BPP, so the same 0.905-1.10 window is the right physical range whatever BPP is set to.
C_AMP    = 846
FWHM2SIG = 1.6651
WIN, STEP, MAXS = 50, 25, 15
TRIM_PX = int(os.environ.get('HDOM_TRIM_PX', 5000 // BPP))   # terminal bases dropped from scoring
                             # (both ends, R & G). 25 px = 5 kb; HDOM_TRIM_PX=0 keeps terminal labels.
SCALE_MIN = float(os.environ.get('LOCAL_SCALE_MIN', 181/200))   # env override; default = pcc bpp 181
SCALE_MAX = float(os.environ.get('LOCAL_SCALE_MAX', 220/200))   # env override; default = pcc bpp 220
# note: maligner already searches scale 0.82-1.18, so lowering LOCAL_SCALE_MIN unlocks mal's
# existing low-scale candidates for over-stretched long molecules (no calc_Q re-run needed).
SCALE_STEP = float(os.environ.get('LOCAL_SCALE_STEP', 0.025))   # scale re-optimization grid step
# Two-stage scale search, the way the pcc scan does it: sweep the coarse grid, then refine around
# the winner at LOCAL_SCALE_FINE. Cost is what dominates the re-rank -- opt_ps calls local_ccr once
# per scale point -- so a flat 0.005 grid is 40 points and 101 s per molecule, while coarse+refine
# is 13 points. 0 (the default) keeps the single coarse pass.
SCALE_FINE = float(os.environ.get('LOCAL_SCALE_FINE', 0))
# How many coarse points to refine around. 1 follows only the coarse winner, which misses the
# molecules whose scale-vs-score curve has more than one peak: on a 50-molecule check, refining
# the winner alone reproduced the flat-0.005 answer for 48/50, and the two it missed had their
# true optimum in a different coarse cell (scale 1.030 vs 0.910).
SCALE_TOPK = int(os.environ.get('LOCAL_SCALE_TOPK', 1))
# 0.025 is 5x coarser than the pcc scan, which sweeps 1 bpp = 0.005. Half a grid cell (0.0125)
# displaces the far end of a 1,000 px molecule by 12 px, which is why every tile needs a shift:
# the median |shift| a true locus asks for is 6 px, and refining the grid to 0.005 halves it.
_SCALE_ANCHOR = 181 / 200                  # grid anchored here (=default SCALE_MIN 0.905) so that
                                           # lowering SCALE_MIN only ADDS low scales without shifting
                                           # existing grid points (avoids discretization ranking flips).
def _scale_grid(smin=None, smax=None, step=None):
    """The scale values the re-rank is allowed to try, anchored so the grid is reproducible.

    Steps of SCALE_STEP away from _SCALE_ANCHOR in both directions, clipped to [smin, smax]. The
    anchor matters: without it the grid would depend on where the range happens to start, and two
    runs with different bounds would not share candidate scales.
    """
    import math
    smin = SCALE_MIN if smin is None else smin
    smax = SCALE_MAX if smax is None else smax
    step = SCALE_STEP if step is None else step
    k0 = int(math.ceil((smin - _SCALE_ANCHOR) / step - 1e-9))
    k1 = int(math.floor((smax - _SCALE_ANCHOR) / step + 1e-9))
    return [round(_SCALE_ANCHOR + k * step, 6) for k in range(k0, k1 + 1)]
OFFSET_PX, OFF_STEP = 25, 5   # position re-optimization: search +-5 kb in 1 kb steps
HOM_TOL_BP = 30000
# opt_ps picks the scale by local_cc_r*F1 rather than F1 alone. Selecting on F1 alone lets a
# scale that matches more labels win even where the red profile disagrees, which costs rank-1
# accuracy; set OPT_COMBINED=0 for the F1-only behaviour.
OPT_COMBINED = os.environ.get('OPT_COMBINED', '1') == '1'
# Re-centre the offset window on the shift a scale change implies (see _centers).
# SCALE_ANCHORED=0 restores the pre-2026-09 behaviour (window always at the candidate).
SCALE_ANCHORED = os.environ.get('SCALE_ANCHORED', '1') != '0'
# MAL_KEEP_SCALE=1 (env): for maligner-derived candidates, keep maligner's own scale value
# (re-optimize OFFSET only, not scale) instead of grid-searching scale in opt_ps. Default off.
MAL_KEEP_SCALE = os.environ.get('MAL_KEEP_SCALE', '0') == '1'
# poly-N exclusion: a candidate whose mapped reference span overlaps an assembly N-gap
# (>=100bp N-run in the haplotype FASTA, e.g. acrocentric rDNA/satellite) has no real DNA or
# labels there, so any alignment landing on it is spurious and is dropped from EVERY method's
# rank (local/mDOM/mal/pcc). Uses the ported hDOM_kts overlaps_n. Set POLYN_EXCLUDE=0 to disable.
try:
    import polyn_filter as _PN
    _NREG = _PN.load_n_regions()
except Exception:
    _PN = None; _NREG = {}
POLYN_EXCLUDE = os.environ.get('POLYN_EXCLUDE', '1') != '0'
def _is_polyn(ref, ab, w_px, scale):
    """True if the span this candidate would occupy runs into an assembly N-gap.

    Such a candidate is dropped from every method's ranking: the reference has no sequence there,
    so neither backbone nor labels can be matched. Disabled by POLYN_EXCLUDE.
    """
    return bool(POLYN_EXCLUDE and _PN is not None
                and _PN.overlaps_n(str(ref), int(ab), int(w_px), float(scale), _NREG, bpp=BPP))
# shift_gap and w_coh are the shift roughness the score multiplies in (shift_pen), carried in the
# CSV so it shows the penalty the ranking applied. local_ccr_j returns both, so carrying them
# costs nothing.
COLS = ['rank','ref_name','method','aligned_bp','dir','scale','cc_r','local_cc_r',
        'cc_g','cc_rg','cc_rg2','log(P)','misc','F1','shift_gap','w_coh']

# ---------- ML green-label detection (optional: --ml-labels) ----------
# Default = classic refined_green_peaks. --ml-labels swaps in the adaptive RandomForest
# detector (Label_Detect_ML_based.detect_labels, a drop-in) AND the ML calc_Q (candidates/,
# which reads a per-molecule _mlpeaks.json sidecar). Each molecule is detected ONCE on its
# raw tif (real green) and reused for the widthprop bake + rerank (cached by base name).
ML_MODE = False
_ml_detect = None          # bound to Label_Detect_ML_based.detect_labels when --ml-labels
_ML_CACHE = {}
def _ml_key(f):
    """Molecule name from a tif path, with the _orig/.tif suffixes removed.

    Used as the cache key so a molecule is detected once and reused for both the widthprop bake
    and the re-rank, rather than being re-detected on each.
    """
    n = os.path.basename(f)
    for s in ('_orig.tif', '_orig.tiff', '.tif', '.tiff'):
        if n.endswith(s): return n[:-len(s)]
    return n
def _green_peaks(f):
    """Green nick-label x-positions: ML (cached per molecule) if --ml-labels, else classic."""
    if ML_MODE:
        k = _ml_key(f)
        if k not in _ML_CACHE:
            _ML_CACHE[k] = sorted(int(round(x)) for x in _ml_detect(f))
        return np.asarray(_ML_CACHE[k], dtype=int)
    return refined_green_peaks(f)

# ---------- [1] widthprop ----------
def build_widthprop_tif(src_tif, out_tif):
    """Rewrite the green channel as one Gaussian per detected label, all at the same amplitude.

    Each label keeps the width it was measured at -- the Gaussian sigma comes from that peak's FWHM
    -- but every peak is given the same height C_AMP. Brightness therefore stops carrying weight in
    any score computed on this image, while position and width survive. Red is copied through
    untouched, and the synthesised green is written as a band across the middle rows.
    """
    d = tiff.imread(src_tif); red = d[0].astype(np.uint16); H, W = d[1].shape
    g = np.max(d[1], axis=0).astype(float)
    pk = _green_peaks(src_tif)
    x = np.arange(W); prof = np.zeros(W)
    if len(pk):
        fwhm = signal.peak_widths(g, pk, rel_height=0.5)[0]
        for p, fw in zip(pk, fwhm):
            sig = max(fw, 1e-3) / FWHM2SIG
            prof += C_AMP * np.exp(-((x - p) / sig) ** 2)
    prof = np.clip(prof, 0, 65535)
    newg = np.zeros((H, W), dtype=np.uint16); newg[H//2-2:H//2+3, :] = prof.astype(np.uint16)
    tiff.imwrite(out_tif, np.stack([red, newg]))

# ---------- red profiles ----------
_refprof = {}
def ref_rg(ref):
    """cached (red, green) reference profile — avoids re-loading big chromosome tifs per candidate."""
    if ref not in _refprof:
        p = f'{REFMAP}/{ref}.tif'
        if os.path.exists(p):
            r, g, _ = lib.tif2RGprofile(p)
            _refprof[ref] = (np.asarray(r, float), np.asarray(g, float))
        else:
            _refprof[ref] = (None, None)
    return _refprof[ref]
def ref_red(ref):
    """Red profile of a reference chromosome, from the same cache as ref_rg."""
    return ref_rg(ref)[0]
_refpk = {}
def ref_peaks(ref):
    """cached (ref_peaks, ref_len) — avoids re-reading the chromosome tif per candidate."""
    if ref not in _refpk:
        _refpk[ref] = R.ref_green_peaks(f'{REFMAP}/{ref}.tif')
    return _refpk[ref]


# 'outermost' keeps every terminal label but the outermost one at each end. A fixed 5 kb window
# can also take a second terminal label that is a real interior site; HDOM_TRIM_RULE=window
# selects that fixed window instead.
TRIM_RULE = os.environ.get('HDOM_TRIM_RULE', 'outermost')
def terminal_cuts(dp, w2, trim):
    """How much to cut at each end, in molecule pixels.

    'window' (the original) cuts the whole terminal `trim`. 'outermost' cuts just past the single
    outermost label at each end, so a second terminal label -- which the fixed window also removed,
    though it is a real site -- keeps its say. A clamped end is exempt and keeps the full cut: two
    labels on one pixel, or two on the image edge, are the arc-length resampling artefact the trim
    exists for, not evidence.

    The exemption is itself bounded. Labels within a pixel of the outermost one are that same pile
    -- one label, however many columns the resampling stacked onto it -- and a label beyond the
    pile is an interior site that was never terminal. Where such a label exists, the cut stops just
    past the pile and keeps it; where none does, the full window is cut as before, so the scored
    reference window is unchanged for every end that has nothing to save.

    Returns (cut_left, cut_right); both default to `trim`, so callers are unchanged under 'window'.
    """
    if TRIM_RULE != 'outermost' or trim <= 0:
        return trim, trim
    L = [float(x) for x in dp if x < trim]
    Rr = [float(x) for x in dp if x > w2 - 1 - trim]
    def keep_full(side, edge):
        if len(side) < 2: return True
        q = [int(round(v)) for v in side]
        if len(q) != len(set(q)): return True                     # duplicate pixels = clamped
        return sum(1 for v in q if abs(v - edge) <= 1) >= 2
    def past_pile(side, left):
        """Cut past the outermost pile, or None when the pile is the whole window (nothing saved)."""
        q = [int(round(v)) for v in side]
        o = min(q) if left else max(q)
        if all(abs(v - o) <= 1 for v in q): return None           # no interior label to rescue
        e = max(v for v in q if abs(v - o) <= 1) if left else min(v for v in q if abs(v - o) <= 1)
        return min(trim, (e + 1) if left else (w2 - 1 - e + 1))
    def cut(side, left, edge):
        if not keep_full(side, edge):
            return (int(min(side)) + 1) if left else (int(w2 - 1 - max(side)) + 1)
        p = past_pile(side, left) if side else None
        return trim if p is None else p
    return cut(L, True, 0), cut(Rr, False, w2 - 1)

def _mapped_ref(ref_prof, MOL, aligned_bp, scale, direction):
    """The reference resampled onto the molecule's pixel grid at a given placement.

    Returns an array of length MOL: reference pixel aligned_bp/BPP + scale*i for molecule pixel i,
    reversed when direction is -1, linearly interpolated. Note the direction of the mapping -- the
    REFERENCE is brought to the molecule, so molecule coordinates stay fixed and every candidate is
    comparable on the same axis. The pcc scan does the opposite, resampling the molecule instead.
    """
    sh = aligned_bp / BPP; p = np.arange(MOL)
    xs = sh + scale * (MOL - 1 - p) if direction == -1 else sh + scale * p
    # fast linear interp at integer ref positions (avoids allocating np.arange(len(ref)) each call)
    N = len(ref_prof)
    xi = np.clip(xs, 0.0, N - 1.0)
    lo = xi.astype(np.intp); frac = xi - lo
    hi = np.minimum(lo + 1, N - 1)
    return ref_prof[lo] * (1.0 - frac) + ref_prof[hi] * frac

def local_ccr_j(mol_red, ref_prof, aligned_bp, scale, direction):
    """tile-wise red Pearson with local +-3kb shift: (mean of per-tile best, shift roughness).

    The roughness is the mean |difference| between the shifts chosen by neighbouring tiles, in px.
    It costs nothing: the shifts are the argmax the loop already takes. A true locus needs a shift
    that drifts slowly (residual stretch), so its neighbouring tiles agree; a wrong locus finds its
    per-tile optimum in an unrelated place each time, so its shifts jump. Measured as an ADJACENT
    difference, not a spread: a locus whose shift ramps across the whole search range is smooth.
    """
    MOL = len(mol_red); mf = _mapped_ref(ref_prof, MOL, aligned_bp, scale, direction)
    # green-only trimming: RED uses the full molecule (backbone at ends is real signal);
    # only the green label-set (F1) drops the terminal 5 kb. So no red trim here.
    vals = []; shifts = []
    for st in range(0, MOL - WIN, STEP):
        a = mol_red[st:st+WIN]
        if np.std(a) < 1e-6: continue
        best = -9.0; bs = None
        for s in range(-MAXS, MAXS+1, 3):
            l = st + s
            if l < 0 or l + WIN > len(mf): continue
            bb = mf[l:l+WIN]
            if np.std(bb) < 1e-6: continue
            c = np.corrcoef(a, bb)[0, 1]
            if c > best: best = c; bs = s
        if best > -9: vals.append(best); shifts.append(bs)
    if not vals: return np.nan, np.nan
    jump = float(np.abs(np.diff(shifts)).mean()) if len(shifts) > 1 else 0.0
    return float(np.mean(vals)), jump

def local_ccr(mol_red, ref_prof, aligned_bp, scale, direction):
    """tile-wise red Pearson with local +-3kb shift; mean of per-tile best."""
    return local_ccr_j(mol_red, ref_prof, aligned_bp, scale, direction)[0]

# SHIFT_PEN=1 multiplies the ranking score by (1 - roughness / full search width), so a candidate
# only keeps the correlation its local shifts EARNED COHERENTLY. Reported cc_rL is untouched: this
# changes the order of candidates, not the number printed in the table. Off by default.
SHIFT_PEN = os.environ.get('SHIFT_PEN', '1') == '1'    # production default; SHIFT_PEN=0 restores the plain product
def shift_pen(jump):
    # env is read per call so a driver can score both arms in one process
    if os.environ.get('SHIFT_PEN', '1' if SHIFT_PEN else '0') != '1' or jump != jump: return 1.0
    return max(0.0, 1.0 - float(jump) / (2.0 * MAXS))

def global_ccr(mol_red, ref_prof, aligned_bp, scale, direction):
    """whole-molecule red Pearson at the alignment (single corrcoef)."""
    MOL = len(mol_red); mf = _mapped_ref(ref_prof, MOL, aligned_bp, scale, direction)
    if np.std(mol_red) < 1e-6 or np.std(mf) < 1e-6: return np.nan
    return float(np.corrcoef(mol_red, mf)[0, 1])

def global_cc_rg(mol_red, mol_green, ref, aligned_bp, scale, direction):
    """cosine cc_r (red) + cc_g (green dcc) at an alignment, matching the pcc-scan method.
    Used to fill cc_g/cc_rg/cc_rg2 for maligner-only candidates (which lack cosine green)."""
    mr, mg = ref_rg(ref)
    if mr is None: return np.nan, np.nan
    cr = lib.normalize(mr); cg = lib.normalize(mg)
    yr = lib.normalize(mol_red); yg = lib.normalize(mol_green)
    n = len(yr); nn = int(scale * n)
    if nn <= 1: return np.nan, np.nan
    yyr = signal.resample(yr, nn); yyg = signal.resample(yg, nn)
    tr = yyr if direction == 1 else yyr[::-1]; tg = yyg if direction == 1 else yyg[::-1]
    sh = int(round(aligned_bp / BPP))
    if sh < 0 or sh + nn > len(cr): return np.nan, np.nan
    return float(lib.corrcoef(cr[sh:sh+nn], tr)), float(lib.dcc(cg[sh:sh+nn], tg))

# ---------- helpers ----------
def _f(v):
    """float(v), or NaN if it will not convert. CSV fields are often blank."""
    try: return float(v)
    except Exception: return float('nan')
def is_sex(ref):
    """True for chrX/chrY references."""
    return ref.startswith('chrX') or ref.startswith('chrY')
_lifters = {}
def lifter(chrnum):
    """Cached PafLiftover for one chromosome number, or None if that PAF is absent.

    Used to decide whether two candidate loci are the two haplotypes of the same place.
    """
    if chrnum not in _lifters:
        p = f'{CHAINS}/chr{chrnum}.paf'
        _lifters[chrnum] = PafLiftover(p) if os.path.exists(p) else None
    return _lifters[chrnum]

def _cbase(r):
    """Chromosome name without the _m/_p haplotype suffix."""
    r = str(r); return r[:-2] if r.endswith(('_m', '_p')) else r

def paf_lift_point(lo, pos, src, tol_edge=20000):
    """Robust PAF-block-interpolation lift (survives point-lift CIGAR gaps). src='p'(query,PATERNAL)
    or 'm'(target,MATERNAL). Prefers a block that CONTAINS pos; else the nearest edge-matching block.
    Returns {'pos':int, 'strand':+1/-1} or None."""
    if lo is None or not getattr(lo, 'blocks', None): return None
    best = None
    for blk in lo.blocks:
        qs, qe, ts, te, st = blk.qstart, blk.qend, blk.tstart, blk.tend, blk.strand
        strand = 1 if st == '+' else -1
        if src == 'p':                                   # query -> target
            if not (qs - tol_edge <= pos <= qe + tol_edge): continue
            frac = (min(max(pos, qs), qe) - qs) / max(qe - qs, 1)
            mapped = ts + frac * (te - ts) if st == '+' else te - frac * (te - ts)
            inside = qs <= pos <= qe
        else:                                            # target -> query
            if not (ts - tol_edge <= pos <= te + tol_edge): continue
            frac = (min(max(pos, ts), te) - ts) / max(te - ts, 1)
            mapped = qs + frac * (qe - qs) if st == '+' else qe - frac * (qe - qs)
            inside = ts <= pos <= te
        if inside: return {'pos': int(mapped), 'strand': strand}
        if best is None: best = {'pos': int(mapped), 'strand': strand}
    return best

def paf_homolog(l1, l2, tol_edge=20000, tol_map=50000):
    """PAF-block homolog test (robust to gaps): l1,l2 opposite haplotypes of same chr are homolog iff
    both fall in corresponding query/target ranges of one PAF block and l1 block-maps near l2."""
    r1, a1 = l1; r2, a2 = l2
    if not (r1 and r2 and str(r1)[:3] == 'chr' and str(r2)[:3] == 'chr' and _cbase(r1) == _cbase(r2)
            and str(r1)[-2:] in ('_m', '_p') and str(r2)[-2:] in ('_m', '_p') and str(r1)[-2:] != str(r2)[-2:]):
        return False
    try: a1, a2 = int(float(a1)), int(float(a2))
    except Exception: return False
    lo = lifter(str(r1)[3:-2])
    if lo is None or not getattr(lo, 'blocks', None): return False
    qpos, tpos = (a1, a2) if str(r1)[-1] == 'p' else (a2, a1)   # query=paternal, target=maternal
    for blk in lo.blocks:
        qs, qe, ts, te, st = blk.qstart, blk.qend, blk.tstart, blk.tend, blk.strand
        if not (qs - tol_edge <= qpos <= qe + tol_edge): continue
        if not (ts - tol_edge <= tpos <= te + tol_edge): continue
        frac = (min(max(qpos, qs), qe) - qs) / max(qe - qs, 1)
        mapped_t = ts + frac * (te - ts) if st == '+' else te - frac * (te - ts)
        if abs(mapped_t - tpos) <= tol_map: return True
    return False
def null_logp(scores):
    """log10 P from a normal fitted to `scores`, the score being converted included.

    Two caveats, both deliberate and both worth knowing before the column is read as a probability.

    The population is the candidate list, so the score being converted is INSIDE the sample its own
    mean and sd come from. With three to six candidates per molecule that is not a null in the usual
    sense; it reports how far rank 1 stands out from the others, rewritten on a log scale.

    The conversion fits a normal by moment matching and reads its upper tail (base 10), as in
    logp_null.py -- but here without even a 5000-draw null behind the two estimates. The tail is not
    itself fitted, so this is not "tail extrapolation" in the extreme-value sense. Use the column to
    order candidates, not as a calibrated p-value. The proper per-molecule null is `--logp`
    (logp_null.py), which replaces this for local rows; see METHODS.md section 6.

    sd falls back to 1.0 when every candidate scores the same (logp_null.py uses 1e-9 for its own
    guard; the two differ, which only matters in that degenerate case).
    """
    mu = float(np.mean(scores)); sd = float(np.std(scores)) or 1.0
    return lambda s: f'{float(stats.norm.logsf((s-mu)/sd))/np.log(10):.2f}'

_SCALES_MOD = _scale_grid()

def homolog_H(method, ranked, mol_red, mol_green, dp, w2, lpfn, hom_tol=HOM_TOL_BP):
    """Return a 14-cell H row (method 'H_local'/'H_mDOM') for `ranked` (list of dicts with
    ref,ab,di,sc and the method's score fields), or None.
    Add H iff rank1 is not a sex chr AND rank2 is not rank1's homolog.
    Prefer an existing ranked candidate at the homolog locus (already scored by the method,
    incl. local re-optimization); else lift rank1 -> homolog and evaluate by the method.
    (mDOM has no re-opt, so an existing candidate row is especially preferred there.)"""
    if not ranked: return None
    a = ranked[0]; ref1 = str(a['ref'])
    if is_sex(ref1) or ref1[:3] != 'chr' or ref1[-2:] not in ('_p', '_m'): return None
    chrnum = ref1[3:-2]; hap = ref1[-1]; hom = f'chr{chrnum}_{"m" if hap == "p" else "p"}'
    lo = lifter(chrnum)
    if lo is None or not os.path.exists(f'{REFMAP}/{hom}.tif'): return None
    res = paf_lift_point(lo, int(a['ab']), hap)   # PAF-block lift (robust to point-lift gaps)
    if not res: return None
    hpos = int(res['pos']); hdir = int(a['di']) * int(res['strand'])
    # already a homolog pair at rank2 (PAF-block homolog, same criterion as classify)? then no H needed
    if len(ranked) >= 2 and paf_homolog((ref1, a['ab']), (str(ranked[1]['ref']), ranked[1]['ab'])):
        return None
    mlab = 'H_local' if method == 'local' else 'H_mDOM'
    # (1) prefer an existing candidate at the homolog locus (rank 2..N)
    for c in ranked[1:]:
        if str(c['ref']) == hom and abs(int(c['ab']) - hpos) <= hom_tol:
            if method == 'local':
                return ['H', hom, mlab, c['ab'], c['di'], c['sc'], c.get('gcc', ''), c.get('lc', ''),
                        c.get('cc_g', ''), c.get('cc_rg', ''), c.get('cc_rg2', ''), c.get('logp', ''), c.get('FP', ''), c.get('f1', '')]
            return ['H', hom, mlab, c['ab'], c['di'], c['sc'], c.get('cc_r', ''), '',
                    c.get('cc_g', ''), c.get('cc_rg', ''), c.get('cc_rg2', ''), c.get('logp', ''), c.get('FP', ''), c.get('f1', '')]
    # (2) lift + evaluate at the homolog position
    rp_red = ref_red(hom)
    if rp_red is None: return None
    rp, rl = ref_peaks(hom)
    def f1_at(ab, di, sc, trim):
        M, FP, FN, f1 = R.label_f1(dp, w2, rp.astype(float), rl, ab, di, sc, trim_px=trim)
        return (0.0 if (f1 is None or np.isnan(f1)) else float(f1)), (int(FP) if not np.isnan(FP) else 0)
    if method == 'local':
        best = (-1.0, hpos, float(a['sc']), 0)              # re-optimize offset+scale (green trim), like local
        for scv in _SCALES_MOD:
            for do in range(-OFFSET_PX, OFFSET_PX + 1, OFF_STEP):
                ab2 = hpos + do * BPP
                if ab2 < 0: continue
                f, FP = f1_at(ab2, hdir, scv, TRIM_PX)
                if f > best[0]: best = (f, ab2, scv, FP)
        f1, ab, sc, FP = best
        if f1 < 0: return None
        try:
            lc, _jmp = local_ccr_j(mol_red, rp_red, ab, sc, hdir)
        except Exception:
            return None
        if np.isnan(lc): return None
        lc = max(lc, 0.0); score = lc * f1 * shift_pen(_jmp)
        gcc = global_ccr(mol_red, rp_red, ab, sc, hdir)
        _, cg2 = global_cc_rg(mol_red, mol_green, hom, ab, sc, hdir)
        ok2 = (gcc == gcc and cg2 == cg2)
        return ['H', hom, mlab, ab, hdir, round(sc, 2), round(gcc, 2) if gcc == gcc else '', round(lc, 2),
                round(cg2, 2) if cg2 == cg2 else '', round(gcc * cg2, 3) if ok2 else '',
                round(max(gcc, 0.0) * cg2 * cg2, 3) if ok2 else '', lpfn(score), FP, round(f1, 2),
                '' if _jmp != _jmp else f'{_jmp:.2f}', f'{shift_pen(_jmp):.3f}']
    # mDOM (C1): no re-opt, evaluate at lifted position with rank1's scale
    sc = float(a['sc']); ab = hpos
    f1, FP = f1_at(ab, hdir, sc, 0)
    gcc = global_ccr(mol_red, rp_red, ab, sc, hdir)
    if gcc != gcc: gcc = 0.0
    score = max(gcc, 0.0) * f1 * f1
    _, cg2 = global_cc_rg(mol_red, mol_green, hom, ab, sc, hdir)
    ok2 = (cg2 == cg2)
    return ['H', hom, mlab, ab, hdir, round(sc, 2), round(gcc, 2), '',
            round(cg2, 2) if ok2 else '', round(gcc * cg2, 3) if ok2 else '',
            round(max(gcc, 0.0) * cg2 * cg2, 3) if ok2 else '', lpfn(score), FP, round(f1, 2), '', '']

# ---------- [3]-[5] rerank one molecule ----------
def rerank(tif, csvp):
    df = pd.read_csv(csvp)
    if 'method' not in df.columns: return None
    pcc = df[df['method'] == 'pcc']; mal = df[df['method'] == 'mal']
    if pcc.empty and mal.empty: return None
    try:
        dp = _green_peaks(tif).astype(float)
        dr, dg, _ = lib.tif2RGprofile(tif); mol_red = np.asarray(dr, float); mol_green = np.asarray(dg, float); w2 = len(dg)
    except Exception:
        return None

    # long-molecule offset widening (env-gated experiment): maligner's GLOBAL position error grows
    # with molecule length; the default +-5 kb opt_ps offset window can't correct it, so the true
    # homolog locus is missed (lands on a wrong scale). For molecules >= LONG_MOL_KB, widen the
    # offset re-optimization window to LONG_OFFSET_PX. Default (env unset) = no change.
    _long_kb = float(os.environ.get('LONG_MOL_KB', '1e9'))
    off_px = int(os.environ.get('LONG_OFFSET_PX', str(OFFSET_PX))) if (w2 * BPP >= _long_kb * 1000) else OFFSET_PX

    def f1_gen(ref, ab, di, sc, trim):
        rp, rl = ref_peaks(ref)
        _tl, _tr = terminal_cuts(dp, w2, trim)
        M, FP, FN, f1 = R.label_f1(dp, w2, rp.astype(float), rl, ab, di, sc,
                                   trim_px=trim, trim_l=_tl, trim_r=_tr)
        return (0.0 if (f1 is None or np.isnan(f1)) else float(f1)), (int(FP) if not np.isnan(FP) else 0)
    def f1_at(ref, ab, di, sc):        # local uses green-only 5 kb trim
        return f1_gen(ref, ab, di, sc, TRIM_PX)

    _SCALES = _scale_grid()
    def _centers(ab0, sc0, scv):
        """Where to centre the offset window when the trial scale differs from the candidate's.

        The placement maps molecule pixel x to reference pixel  shift + scale*u,  u = x (dir +1) or
        w2-1-x (dir -1), with shift = aligned_bp/BPP. So the u=0 end sits at aligned_bp and the
        u=w2-1 end at aligned_bp + scale*(w2-1): position and scale are COUPLED. Changing the scale
        by ds moves the far end by ds*(w2-1) pixels, i.e. ds*(w2-1)*BPP bases. Searching the scale
        grid while holding the window at ab0 therefore only ever evaluates the trial scale at the
        wrong place -- for a 392 kb molecule a 0.074 scale step needs a 29 kb shift, six times the
        +-5 kb window, so the true locus is unreachable and a poor scale wins by default.

        Which end is anchored is not known (it depends on where the candidate's evidence actually
        sits), so both ends and the midpoint are tried: t = 0, 0.5, 1 of the implied shift. Cost is
        3x the offset loop, not a wider window, so a wrong candidate cannot wander further than
        before -- every centre stays on the line the scale change itself defines.
        """
        d = -(scv - sc0) * (w2 - 1) * BPP
        # a maligner row whose scale came out non-finite (degenerate query map: too few fragments)
        # would make d inf and int(round(...)) raise; such a candidate just keeps its own centre.
        if not np.isfinite(d): return [ab0]
        return [ab0, ab0 + 0.5 * d, ab0 + d] if SCALE_ANCHORED and abs(d) >= BPP else [ab0]
    def opt_ps(ref, ab0, di, trim, sc0=None):
        """re-optimize aligned_bp (+-5 kb per centre) and scale ([0.905,1.10]). Returns (bestF1, ab, scale, FP).
        OPT_COMBINED: pick the SCALE by local_cc_r*F1 (not F1 alone) so a marginally-higher-F1 but
        red-poor scale isn't chosen (fixes homolog scale mismatch). Offset stays F1-optimized per scale."""
        _sc0 = 1.0 if sc0 is None else float(sc0)
        if OPT_COMBINED:
            # exact: for each scale, F1-best offset, then local_cc_r there; pick max local_cc_r*F1
            rpr = ref_red(ref)
            def _eval(scv):
                bf1 = -1.0; bab = ab0; bfp = 0
                for _c in _centers(ab0, _sc0, scv):
                  for doff in range(-off_px, off_px + 1, OFF_STEP):
                    ab2 = int(round(_c)) + doff * BPP
                    if ab2 < 0: continue
                    try:
                        f12, FP2 = f1_gen(ref, ab2, di, scv, trim)
                    except Exception:
                        continue
                    if f12 > bf1: bf1, bab, bfp = f12, ab2, FP2
                if bf1 < 0 or rpr is None: return None
                try:
                    lc = local_ccr(mol_red, rpr, bab, scv, di)
                except Exception:
                    return None
                if np.isnan(lc): return None
                return (max(lc, 0.0) * bf1, bab, float(scv), bfp, bf1)
            best = (-1.0, ab0, 1.0, 0, -1.0)   # combscore, ab, sc, FP, f1
            seen = set(); coarse = []
            for scv in _SCALES:
                seen.add(round(float(scv), 6))
                r = _eval(scv)
                if not r: continue
                coarse.append(r)
                if r[0] > best[0]: best = r
            if SCALE_FINE > 0 and coarse:
                # refine one coarse cell either side of the TOPK best coarse points, so a curve with
                # more than one peak is followed in every candidate basin, not just the leading one
                coarse.sort(key=lambda r: -r[0])
                for r0 in coarse[:max(1, SCALE_TOPK)]:
                    lo = max(SCALE_MIN, r0[2] - SCALE_STEP)
                    hi = min(SCALE_MAX, r0[2] + SCALE_STEP)
                    v = lo
                    while v <= hi + 1e-9:
                        key = round(v, 6)
                        if key not in seen:
                            seen.add(key)
                            r = _eval(v)
                            if r and r[0] > best[0]: best = r
                        v += SCALE_FINE
            return (best[4], best[1], best[2], best[3])
        best = (-1.0, ab0, 1.0, 0)
        for scv in _SCALES:
          for _c in _centers(ab0, _sc0, scv):
            for doff in range(-off_px, off_px + 1, OFF_STEP):
                ab2 = int(round(_c)) + doff * BPP
                if ab2 < 0: continue
                try:
                    f12, FP2 = f1_gen(ref, ab2, di, scv, trim)
                except Exception:
                    continue
                if f12 > best[0]:
                    best = (f12, ab2, float(scv), FP2)
        return best
    def opt_offset(ref, ab0, di, sc, trim):
        """re-optimize aligned_bp (+-5 kb) only, keeping scale (for mDOM: pcc scale already limited)."""
        best = (-1.0, ab0, 0)
        for doff in range(-off_px, off_px + 1, OFF_STEP):
            ab2 = ab0 + doff * BPP
            if ab2 < 0: continue
            try:
                f12, FP2 = f1_gen(ref, ab2, di, sc, trim)
            except Exception:
                continue
            if f12 > best[0]:
                best = (f12, ab2, FP2)
        return best

    # ----- mDOM (C1): pcc top-10, NO re-optimization (PCC anchor kept), dedup by locus -----
    # F1 = label_f1 at the PCC alignment (no offset/scale re-opt, no trim); cc_r = PCC's cc_r.
    # mDOM = max(cc_r,0) * F1^2.  (Ablation: offset re-opt + cc_r recompute broke red/green
    # co-registration -> A_group homolog rate 98%->81%; C1 reproduces the original ~99%.)
    CLUSTER_BP = 1_000_000
    def _key(ref, ab): return (ref, int(ab) // CLUSTER_BP)
    MOL_LEN_BP = w2 * BPP   # molecule's genomic footprint; two same-ref candidates within one
                            # molecule-length overlap -> same locus; farther apart -> distinct (repeat).
    def _dedup_locus(items, merge_bp=MOL_LEN_BP):
        """distance-based dedup: drop a candidate within merge_bp (= molecule length) of a
        higher-scored same-ref one. Fixes the fixed-1Mb-block boundary artifact (near-duplicates
        straddling a block edge, e.g. 3,999,200 vs 4,001,000, survived as separate rank1/rank2)
        while keeping genuine repeat loci (>1 molecule-length apart) separate. items must be
        pre-sorted by score desc."""
        kept = []
        for d in items:
            if any(k['ref'] == d['ref'] and abs(int(k['ab']) - int(d['ab'])) < merge_bp for k in kept):
                continue
            kept.append(d)
        return kept
    mdom_by = {}
    for _, r in pcc.iterrows():
        ref = str(r['ref_name']); ab = int(r['aligned_bp']); di = int(r['dir']); sc = float(r['scale'])
        if ref_red(ref) is None: continue                # ref tif must exist
        f1, FP = f1_gen(ref, ab, di, sc, 0)              # F1 at PCC position, no re-opt, no trim
        cc_r = _f(r.get('cc_r'))                         # PCC's cc_r (co-registered red+green)
        if cc_r != cc_r: cc_r = 0.0
        score = max(cc_r, 0.0) * f1 * f1
        d = dict(ref=ref, ab=ab, di=di, sc=round(sc,2), cc_r=round(cc_r,2), f1=round(f1,2), FP=FP,
                 cc_g=r.get('cc_g',''), cc_rg=r.get('cc_rg',''), cc_rg2=r.get('cc_rg2',''), sc_m=score)
        k = _key(ref, ab)
        if k not in mdom_by or score > mdom_by[k]['sc_m']:
            mdom_by[k] = d              # dedup by locus, keep best mDOM score
    _mcand = [d for d in mdom_by.values() if not _is_polyn(d['ref'], d['ab'], w2, d['sc'])]
    mdom = _dedup_locus(sorted(_mcand, key=lambda d: -d['sc_m']))
    if mdom:
        lpm = null_logp([d['sc_m'] for d in mdom])
        for d in mdom: d['logp'] = lpm(d['sc_m'])

    # ----- local: pcc UNION mal (dedup by ref+LOCUS), local cc_r x F1 -----
    # A ref (e.g. chr05_p) can appear at multiple distinct loci (>1 Mb apart); each is a
    # separate candidate. Only near-duplicate positions (same locus, pcc vs mal) are merged.
    # (CLUSTER_BP / _key defined above with the mDOM loop; reused here.)
    pcc_by_key = {}
    for _, r in pcc.iterrows():
        pcc_by_key.setdefault(_key(str(r['ref_name']), r['aligned_bp']), r)
    cand = {}
    for _, r in pd.concat([pcc, mal]).iterrows():
        ref = str(r['ref_name']); rp_red = ref_red(ref)
        if rp_red is None: continue
        ab0 = int(r['aligned_bp']); di = int(r['dir'])
        # re-optimize offset (+-5 kb) AND scale ([0.905,1.10]) for max F1 (green trim); this also
        # re-fits extreme-scale maligner candidates into the valid range (replaces the scale filter).
        # MAL_KEEP_SCALE=1: maligner candidates keep maligner's scale (offset-only re-opt).
        if MAL_KEEP_SCALE and str(r.get('method','')).strip().lower() == 'mal':
            try: _msc = float(r['scale'])
            except Exception: _msc = 1.0
            f1, ab, FP = opt_offset(ref, ab0, di, _msc, TRIM_PX); sc = _msc
        else:
            try: _c_sc = float(r['scale'])
            except Exception: _c_sc = 1.0
            f1, ab, sc, FP = opt_ps(ref, ab0, di, TRIM_PX, _c_sc)
        if f1 < 0: continue
        try:
            lc, _jmp = local_ccr_j(mol_red, rp_red, ab, sc, di)
        except Exception:
            continue
        if np.isnan(lc): continue
        lc = max(lc, 0.0); score = lc * f1 * shift_pen(_jmp)
        key = _key(ref, ab)
        pr = pcc_by_key.get(key)
        if pr is not None:
            gcc = _f(pr.get('cc_r'))
            ccg = pr.get('cc_g',''); ccrg = pr.get('cc_rg',''); ccrg2 = pr.get('cc_rg2','')
        else:
            # maligner-only locus: robust red via interp; cosine green at this locus (selection context)
            gcc = global_ccr(mol_red, rp_red, ab, sc, di)
            _, cg2 = global_cc_rg(mol_red, mol_green, ref, ab, sc, di)
            ok2 = (gcc == gcc and cg2 == cg2)
            ccg = round(cg2, 2) if cg2 == cg2 else ''
            ccrg = round(gcc * cg2, 3) if ok2 else ''
            ccrg2 = round(max(gcc, 0.0) * cg2 * cg2, 3) if ok2 else ''
        if key not in cand or score > cand[key]['sc_l']:
            cand[key] = dict(ref=ref, ab=ab, di=di, sc=round(sc,2), lc=round(lc,2), f1=round(f1,2), FP=FP,
                             gcc=round(gcc,2) if gcc==gcc else '', cc_g=ccg, cc_rg=ccrg, cc_rg2=ccrg2, sc_l=score,
                             jmp=_jmp, wj=shift_pen(_jmp))
    _lcand = [d for d in cand.values() if not _is_polyn(d['ref'], d['ab'], w2, d['sc'])]
    local = _dedup_locus(sorted(_lcand, key=lambda d: -d['sc_l']))
    if not local: return None
    lpl = null_logp([d['sc_l'] for d in local])
    for d in local: d['logp'] = lpl(d['sc_l'])

    # ----- [4] homolog H rows: H_local (homolog of local rank1) and H_mDOM (of mDOM rank1) -----
    Hloc = homolog_H('local', local, mol_red, mol_green, dp, w2, lpl)
    Hmd = homolog_H('mDOM', mdom, mol_red, mol_green, dp, w2, lpm) if mdom else None

    # ----- assemble rows: local, H_local, mDOM, H_mDOM, mal, pcc -----
    out = []
    for i, d in enumerate(local, 1):
        out.append([i, d['ref'], 'local', d['ab'], d['di'], d['sc'], d['gcc'], d['lc'],
                    d['cc_g'], d['cc_rg'], d['cc_rg2'], d['logp'], d['FP'], d['f1'],
                    '' if d['jmp'] != d['jmp'] else f"{d['jmp']:.2f}", f"{d['wj']:.3f}"])
    if Hloc:
        out.append(Hloc)
    for i, d in enumerate(mdom, 1):
        out.append([i, d['ref'], 'mDOM', d['ab'], d['di'], d['sc'], d['cc_r'], '',
                    d['cc_g'], d['cc_rg'], d['cc_rg2'], d['logp'], d['FP'], d['f1'], '', ''])
    if Hmd:
        out.append(Hmd)
    _i = 0
    for _, r in mal.iterrows():
        if _is_polyn(r['ref_name'], int(r['aligned_bp']), w2, float(r['scale'])): continue
        _i += 1
        out.append([_i, r['ref_name'], 'mal', int(r['aligned_bp']), int(r['dir']), r['scale'],
                    '', '', '', '', '', r.get('log(P)',''), r.get('misc',''), '', '', ''])
    _i = 0
    for _, r in pcc.iterrows():
        if _is_polyn(r['ref_name'], int(r['aligned_bp']), w2, float(r['scale'])): continue
        _i += 1
        out.append([_i, r['ref_name'], 'pcc', int(r['aligned_bp']), int(r['dir']), r['scale'],
                    r.get('cc_r',''), '', r.get('cc_g',''), r.get('cc_rg',''), r.get('cc_rg2',''),
                    r.get('log(P)',''), r.get('misc',''), '', '', ''])
    return out

# ---------- driver ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input'); ap.add_argument('out')
    ap.add_argument('--tag', default='local'); ap.add_argument('--no-view', action='store_true')
    ap.add_argument('--origins', default=None, help='json {base: origin} to record per-molecule origin')
    ap.add_argument('--ml-labels', action='store_true',
                    help='use adaptive ML (RandomForest) green-label detection instead of classic peaks; '
                         'also switches candidates to the ML build (reads _mlpeaks.json sidecar)')
    ap.add_argument('--logp', action='store_true',
                    help='compute K=5000 random-null log(P) (h_score) for local rows + write '
                         '{base}_vsrandom.txt sidecar (same Null A as full_logp_run; slow: loads all ref chr)')
    a = ap.parse_args()

    global ML_MODE, _ml_detect, CALC_Q
    if a.ml_labels:
        sys.path.insert(0, MLPATH)
        import Label_Detect_ML_based as _ML
        _ml_detect = _ML.detect_labels
        # v3 opt-in: prefer an authoritative {base}_mlpeaks.json sidecar next to the INPUT tif
        # (labels detected on the ORIGINAL bent molecule + arc-length-mapped) over re-detecting
        # on the swept image. Keeps the real swept green while preserving the exact label set.
        if os.environ.get('USE_INPUT_MLPEAKS') == '1':
            _base_detect = _ML.detect_labels
            def _ml_detect(f, _bd=_base_detect):
                _sc = os.path.join(os.path.dirname(os.path.abspath(f)), _ml_key(f) + '_mlpeaks.json')
                if os.path.exists(_sc):
                    try:
                        return sorted(float(x) for x in json.load(open(_sc, encoding='utf-8')))
                    except Exception:
                        pass
                return _bd(f)
            print('[ml-labels] USE_INPUT_MLPEAKS=1: input sidecar labels preferred over re-detection')
        ML_MODE = True
        CALC_Q = CALC_Q_ML
        print('[ml-labels] adaptive ML green-label detection ON (ML candidate build)')

    raw = ([a.input] if a.input.endswith(('.tif', '.tiff'))
           else sorted(t for t in glob.glob(os.path.join(a.input, '*.tif')) if '_orig' not in t))
    if not raw: print('no input tif'); return
    OUT = os.path.abspath(a.out); os.makedirs(OUT, exist_ok=True)
    STAGE = OUT + '/.stage'
    if os.path.isdir(STAGE): shutil.rmtree(STAGE)
    os.makedirs(STAGE)
    # let hDOM_mal_Q find genome.maps: STAGE/ref_map -> production ref_map
    link = os.path.join(STAGE, 'ref_map')
    if not os.path.exists(link): os.symlink(REFMAP, link)

    origins = json.load(open(a.origins, encoding='utf-8')) if a.origins and os.path.exists(a.origins) else {}
    bases = []
    for t in raw:
        base = os.path.splitext(os.path.basename(t))[0]
        build_widthprop_tif(t, os.path.join(STAGE, base + '.tif'))
        shutil.copy2(t, os.path.join(OUT, base + '_orig.tif'))
        if ML_MODE:  # STAGE: candidates/hDOM_mal_Q reads it for EXACT ML maligner labels; OUT: metadata (views/QC)
            _pk = _ML_CACHE.get(base, [])
            json.dump(_pk, open(os.path.join(STAGE, base + '_mlpeaks.json'), 'w', encoding='utf-8'))
            json.dump(_pk, open(os.path.join(OUT, base + '_mlpeaks.json'), 'w', encoding='utf-8'))
        bases.append(base)
    print(f'[1] widthprop: {len(bases)} tif')

    subprocess.run(['python3', CALC_Q, STAGE, '--ref-map', REFMAP, '--bpp', str(BPP)], check=False)
    print('[2] hDOM_calc_Q done')

    origin_out = {}; ok = 0; summ = []
    import figure3_view as F3V   # Figure3 CSV (native mal cols) + viewer + _lengths.json
    if a.logp:
        import logp_null as LP; LP.load_refs()
        print(f'[logp] refs loaded ({len(LP.ALL)} chr, genome {LP.GENOME_BP/1e9:.2f}Gb), K={LP.K} null/molecule')
    for base in bases:
        tif = os.path.join(STAGE, base + '.tif'); csvp = os.path.join(STAGE, base + '_matching_results.csv')
        if not os.path.exists(csvp): print(f'  [skip] no candidates: {base}'); continue
        rows = rerank(tif, csvp)
        if rows is None: print(f'  [skip] rerank fail: {base}'); continue
        with open(os.path.join(OUT, base + '_matching_results.csv'), 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f); w.writerow(COLS); w.writerows(rows)
        shutil.copy2(tif, os.path.join(OUT, base + '.tif'))
        F3V.to_figure3_csv(os.path.join(OUT, base + '_matching_results.csv'),
                           os.path.join(STAGE, base + '_maligner_dp.txt'))  # 14->19 col, maligner-native mal
        if a.logp:
            # rerank has just rewritten the CSV, so its log(P) holds null_logp; force restores the K=5000 null.
            print(f'      [logp] {base[:40]:40s} -> {LP.compute_logp(base, OUT, force=True)}')
        org = origins.get(base, a.tag); origin_out[base] = org; ok += 1
        def r1(m): return next((r[1] for r in rows if r[2] == m), '-')
        def sc(m): return next((r[11] for r in rows if r[2] == m), '')
        score_val = sc('local')
        if a.logp:  # on-disk CSV now holds the 5000-null log_p for local rows
            try:
                for _r in csv.DictReader(open(os.path.join(OUT, base + '_matching_results.csv'), encoding='utf-8')):
                    if _r['method'] == 'local': score_val = _r['log(P)']; break
            except Exception: pass
        hasH = any(r[2] == 'homolog' for r in rows)
        summ.append((base, org, score_val, r1('local'), r1('mDOM'), r1('mal'), r1('pcc'), 'Y' if hasH else '-'))
        print(f'  [ok] {base[:40]:40s} [{org:8s}] local={r1("local")} mDOM={r1("mDOM")} mal={r1("mal")} pcc={r1("pcc")} H={"Y" if hasH else "-"}')
    json.dump(origin_out, open(os.path.join(OUT, 'origin.json'), 'w', encoding='utf-8'))   # origin kept here (not in _summary)
    F3V.add_summary(OUT, {b: F3V.summary_row(b, os.path.join(OUT, b + '_matching_results.csv'))
                          for b in origin_out})   # Figure3-format _summary.csv (list/nav/score)
    shutil.rmtree(STAGE, ignore_errors=True)
    print(f'[3-5] reranked: {ok}/{len(bases)}')

    # [6] A_group-design view: comparison panels via hDOM_view --label-dots (adds local1,2)
    #     + extended index.php (4-method table: local/mDOM/mal/pcc + local_cc_r + H)
    import figure3_view as F3V   # Figure3-style viewer + _lengths.json (length/lin. panel labels)
    F3V.add_lengths(OUT, {b: F3V.mol_len_kb(os.path.join(OUT, b + '.tif')) for b in bases})
    if not a.no_view:
        subprocess.run(['python3', os.path.join(CODE, 'hDOM_view.py'),
                        OUT, '--ref-map', REFMAP, '--label-dots'], check=False)
        vi = os.path.join(OUT, 'index.html')
        if os.path.exists(vi): shutil.move(vi, os.path.join(OUT, '_viewgen_index.html'))
        F3V.deploy_viewer(OUT)   # 2-section local+mal viewer (was local_align_index.php)
    print(f'\n[done] {ok}/{len(bases)} -> {OUT}')

if __name__ == '__main__':
    main()
