#!/usr/bin/env python3
"""sv_la.py — single-molecule SV confidence (local-align scoring line; local_ccr x F1).

Shared scoring engine for both SV tests:
 - Inversion (pure-red verdict): pcc rank-1 anchor, hybrid localization
   (hybrid_localize.hybrid_scan_A) over (center, width, gamma), segment inverted in place.
   present iff red_d > 0 and p_T <= ALPHA (T = max-z vs the matched null-of-T, null_bank_T).
   Green (galign) drives discovery/localization only and is excluded from the verdict.
 - Translocation (junction IUT): implemented in confidence.py / null_cross.py; sv_la.py
   supplies the shared placement + scoring.

Reuses local_align_pipeline (LA), mdom_reanchor (R), DOM_lib_Y (lib) by import.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
config.ensure_sys_path()
import csv, json
import numpy as np

try:
    import local_align_pipeline as LA      # local_ccr, global_ccr, ref_red, ref_peaks, BPP, WIN...
    import mdom_reanchor as R              # label_f1, ref_green_peaks, TOL_PX
    import DOM_lib_Y as lib                # tif2RGprofile
except ImportError as e:
    raise SystemExit(f"cannot import shared analysis modules ({e}); "
                     "set SV_LA_CODE / SV_HOMOLOG_CODE / SV_ML_CODE (see README)")

# green nick detection: authoritative ML detect_labels, pipeline refined_green_peaks fallback
_ML_PATH = config.ML_CODE
_ml_detect = None; _MLC = {}
def mol_green_peaks(tif):
    """green nick px. Priority: {base}_mlpeaks.json sidecar -> ML detect_labels (cached)
       -> refined_green_peaks fallback."""
    global _ml_detect
    if tif in _MLC: return _MLC[tif]
    base = tif[:-4] if tif.endswith('.tif') else (tif[:-5] if tif.endswith('.tiff') else tif)
    sc = base + '_mlpeaks.json'
    if os.path.exists(sc):
        try:
            pk = np.sort(np.asarray([float(x) for x in json.load(open(sc))], float))
            _MLC[tif] = pk; return pk
        except Exception:
            pass
    try:
        if _ml_detect is None:
            if _ML_PATH not in sys.path: sys.path.insert(0, _ML_PATH)
            from Label_Detect_ML_based import detect_labels as _dl
            _ml_detect = _dl
        pk = np.sort(np.asarray([float(x) for x in _ml_detect(tif)], float))
        if pk.size >= 1:
            _MLC[tif] = pk; return pk
    except Exception:
        pass
    pk = np.sort(np.asarray(LA._green_peaks(tif), float)); _MLC[tif] = pk; return pk

# ---- constants (inherited from local_align_pipeline) ----
BPP      = LA.BPP                # 200 bp/px
WIN      = LA.WIN               # 50 px = 10 kb tile (local_ccr)
TRIM_PX  = LA.TRIM_PX          # 25 px = 5 kb end trim (whole-molecule F1)
ROOT     = config.ROOT

# ---- breakpoint scan / registration / verdict parameters ----
ENDPOINT_R   = 10000     # left/right endpoint independent scan radius (+/-bp)
STEP_FINE    = 200       # fine breakpoint step (bp) = 1 px
STEP_COARSE  = 2000      # coarse breakpoint step (bp)
MIN_SEG_WIDTH = 8000     # min segment width (bp); endpoints may not cross
REG_TRANS_R  = 6000      # registration search radius (+/-bp); region_score R_trans default
ALPHA        = 0.05
POOL_DIR     = config.POOL_DIR
POOL_CAP     = 3000      # null-pool placement cap


# ============================ placement ============================
class Anchor:
    __slots__ = ('stem','chrom','aligned_bp','direction','scale','W','mol_red',
                 'gpk','ref_red','ref_pk','ref_len','cc_r','logp')

def place(stem, group='M_group', base_dir=None):
    """Use the pcc rank-1 alignment as anchor (fixes the genomic locus). The score
       is computed separately as local_ccr x F1."""
    base = f'{base_dir}/{stem}' if base_dir else f'{ROOT}/{group}/{stem}'
    with open(base + '_matching_results.csv') as f:
        pcc = [r for r in csv.DictReader(f) if r['method'] == 'pcc']
    if not pcc:
        raise ValueError(f'{stem}: no pcc candidate')
    r1 = pcc[0]
    P = Anchor()
    P.stem = stem
    P.chrom = r1['ref_name']
    P.aligned_bp = int(float(r1['aligned_bp']))
    P.direction = int(float(r1['dir']))
    P.scale = float(r1['scale'])
    P.cc_r = _f(r1.get('cc_r')); P.logp = _f(r1.get('log(P)'))
    tif = base + '.tif'
    dr, dg, _ = lib.tif2RGprofile(tif)
    P.mol_red = np.asarray(dr, float)
    P.W = len(P.mol_red)
    P.gpk = mol_green_peaks(tif)   # data green nick px (ML)
    rr = LA.ref_red(P.chrom)
    if rr is None:
        raise ValueError(f'{stem}: ref tif missing for {P.chrom}')
    P.ref_red = rr
    P.ref_pk, P.ref_len = LA.ref_peaks(P.chrom)
    return P


def _f(v):
    try: return float(v)
    except Exception: return float('nan')

def make_anchor(tif, chrom, aligned_bp, direction, scale):
    """Build an Anchor for an arbitrary tif (e.g. a translocation half) under an explicit
       alignment (chrom, ab, dir, scale)."""
    P = Anchor(); P.stem = os.path.basename(tif); P.chrom = chrom
    P.aligned_bp = int(aligned_bp); P.direction = int(direction); P.scale = float(scale)
    dr, dg, _ = lib.tif2RGprofile(tif)
    P.mol_red = np.asarray(dr, float); P.W = len(P.mol_red)
    P.gpk = mol_green_peaks(tif)
    P.ref_red = LA.ref_red(chrom)
    if P.ref_red is None: raise ValueError(f'ref missing {chrom}')
    P.ref_pk, P.ref_len = LA.ref_peaks(chrom)
    P.cc_r = float('nan'); P.logp = float('nan')
    return P

def _aligned_eff(P, slide_bp):
    """Absorb the svmark local slide (slide_bp) into the anchor aligned_bp
       (both channels share one frame)."""
    if not slide_bp: return P.aligned_bp
    return P.aligned_bp + (P.scale * slide_bp if P.direction == 1 else -P.scale * slide_bp)

def _mol_bp_to_ref_px(P, mol_bp, ab):
    p = np.asarray(mol_bp, float) / BPP
    if P.direction == 1:
        return ab / BPP + P.scale * p
    return ab / BPP + P.scale * (P.W - 1 - p)


# ============================ region score (global_ccr x F1, in-place inversion) ============================
# red discriminator = global_ccr (local_ccr's +/-3kb tile is orientation-robust).
# Segment inverted in place at offset 0: fwd/inv share the same position (pure orientation contrast).
REG_SIGMA_PX = 4.0   # green Gaussian sigma (px)

def _score_at(content, labels, ab_o, scale, P, w_reg):
    gcc = LA.global_ccr(content, P.ref_red, ab_o, scale, 1)
    gcc = 0.0 if gcc != gcc else float(gcc)
    if len(labels) == 0:
        return (0.0, gcc, 0.0)
    _, _, _, f1 = R.label_f1(labels, w_reg, P.ref_pk, P.ref_len, ab_o, 1, scale, trim_px=0)
    f1 = 0.0 if (f1 is None or f1 != f1) else float(f1)
    return (max(gcc, 0.0) * f1, gcc, f1)

def _green_align(labels, ab_reg, scale, P, w_reg):
    """Region labels' overlap with ref nicks at offset 0: sum over labels of the nearest-ref
       Gaussian (0..nlab). Measured in place (no offset slide)."""
    if len(labels) == 0:
        return 0.0
    lo = ab_reg / BPP
    qref = lo + scale * np.asarray(labels, float)              # labels (region-local px) -> ref px
    rp = P.ref_pk[(P.ref_pk >= lo - 10) & (P.ref_pk <= lo + scale * w_reg + 10)]
    if len(rp) == 0:
        return 0.0
    # = sqrt(2) (literal, bit-reproducible)
    sig = REG_SIGMA_PX * 1.41421356
    d = rp[None, :] - qref[:, None]
    return float(np.exp(-(d / sig) ** 2).max(axis=1).sum())     # sum of each label's nearest-ref overlap

def region_score(P, l_bp, r_bp, slide_bp=0, R_trans=REG_TRANS_R):
    """Invert the segment [l_bp,r_bp] in place (offset 0) and score forward / inverted at the
       same position: red flipped, green labels mirrored about the center. galign_inv/fwd =
       inverted/forward labels' agreement with the reference."""
    ab = _aligned_eff(P, slide_bp)
    p_l = min(l_bp, r_bp) / BPP; p_r = max(l_bp, r_bp) / BPP
    i0 = max(int(round(p_l)), 0); i1 = min(int(round(p_r)), P.W)
    w_reg = i1 - i0
    if w_reg < WIN:                       # cannot form a tile (<10kb)
        return None
    r_a = _mol_bp_to_ref_px(P, i0 * BPP, ab); r_b = _mol_bp_to_ref_px(P, i1 * BPP, ab)
    rlo_px = float(min(r_a, r_b)); rhi_px = float(max(r_a, r_b))
    ab_reg = rlo_px * BPP
    seg = P.mol_red[i0:i1]
    R_reg = seg if P.direction == 1 else seg[::-1]      # ref-increasing order
    dp = P.gpk[(P.gpk >= i0) & (P.gpk < i1)]
    loc = (dp - i0) if P.direction == 1 else ((i1 - 1) - dp)
    nlab = int(len(loc))
    labels_inv = (w_reg - 1) - loc                     # mirror about the segment center (in place)
    ab_o = ab_reg                                      # in place: offset 0
    s_inv, gcc_inv, f1_inv = _score_at(R_reg[::-1], labels_inv, ab_o, P.scale, P, w_reg)
    s_fwd, gcc_fwd, f1_fwd = _score_at(R_reg,       loc,        ab_o, P.scale, P, w_reg)   # same position, orientation only
    galign_inv = _green_align(labels_inv, ab_reg, P.scale, P, w_reg)   # inverted labels vs reference (search objective)
    galign_fwd = _green_align(loc,        ab_reg, P.scale, P, w_reg)
    return dict(score_fwd=max(s_fwd, 0.0), score_inv=max(s_inv, 0.0),
                gcc_fwd=gcc_fwd, gcc_inv=gcc_inv, f1_fwd=f1_fwd, f1_inv=f1_inv,
                galign_inv=galign_inv, galign_fwd=galign_fwd,
                off=0.0,                               # in place
                nlab=nlab, i0=i0, i1=i1, rlo_px=rlo_px, rhi_px=rhi_px, w_reg=w_reg)


# ============================ whole-molecule score ============================
def reg_shift_molpx(P, off_bp):
    """Convert a registration o* (off_bp, canonical ref-increasing frame) into a mol-px
       shift (accounting for dir / scale).
       dir+1: mol px prop to ref px -> +o*/scale. dir-1: mol px prop to -ref px -> -o*/scale."""
    d = (off_bp / BPP) / P.scale
    return int(round(d if P.direction == 1 else -d))

def _rev_shift_1d(arr, i0, i1, d):
    """1D arr: reverse [i0,i1], then shift the reversed block by d px (placed at
       [i0+d,i1+d]; the rest of arr is unchanged)."""
    out = np.asarray(arr, float).copy(); seg = out[i0:i1][::-1]; W = len(out); d = int(round(d))
    a0 = max(0, i0 + d); a1 = min(W, i1 + d)
    if a1 > a0: out[a0:a1] = seg[(a0 - (i0 + d)):(a1 - (i0 + d))]
    return out

def whole_score(P, l_bp=None, r_bp=None, slide_bp=0, off_bp=0.0):
    """Whole-molecule local_ccr x F1. If l_bp/r_bp are given, score the 'inverted
       hypothesis' with that segment inverted (+ registration o* shift).
       off_bp = the shared offset found by region_score (aligning labels to ref);
       applied to whole too so scoring and visualization agree.
       slide_bp is accepted for call-signature compatibility but not applied."""
    mol = P.mol_red
    gpk = P.gpk
    if l_bp is not None and r_bp is not None:
        p_l = min(l_bp, r_bp) / BPP; p_r = max(l_bp, r_bp) / BPP
        i0 = max(int(round(p_l)), 0); i1 = min(int(round(p_r)), P.W)
        if i1 - i0 >= 2:
            d = reg_shift_molpx(P, off_bp)
            mol = _rev_shift_1d(P.mol_red, i0, i1, d)         # inversion + o* shift
            m = (gpk >= i0) & (gpk < i1)
            gpk = gpk.copy(); gpk[m] = np.clip((i0 + i1 - 1) - gpk[m] + d, 0, P.W - 1)   # mirror + o* shift
            gpk = np.sort(gpk)
    ccr = LA.local_ccr(mol, P.ref_red, P.aligned_bp, P.scale, P.direction)
    _, _, _, f1 = R.label_f1(gpk, P.W, P.ref_pk, P.ref_len, P.aligned_bp, P.direction, P.scale, trim_px=TRIM_PX)
    ccr = 0.0 if ccr != ccr else max(float(ccr), 0.0)
    f1 = 0.0 if (f1 is None or f1 != f1) else float(f1)
    return dict(score=ccr * f1, ccr=ccr, f1=f1)


# ============================ inversion: independent breakpoint scan ============================
def _nlab_in(P, l_bp, r_bp):
    """Label count in the segment, same basis as region_score
       ([round(l/BPP), round(r/BPP)) px, i1 exclusive)."""
    i0 = max(int(round(min(l_bp, r_bp) / BPP)), 0); i1 = min(int(round(max(l_bp, r_bp) / BPP)), P.W)
    return int(((P.gpk >= i0) & (P.gpk < i1)).sum())

def _neighbor_bounds(P, a_bp, b_bp, R):
    """Label-preserving bounds: constrain the scan so an endpoint does not cross a
       neighbouring label. Returns (lmin, lmax, rmin, rmax, nlab0)."""
    qbp = P.gpk * BPP
    # same mask as _nlab_in (rounded-px half-open) so boundary rounding agrees
    i0 = max(int(round(min(a_bp, b_bp) / BPP)), 0); i1 = min(int(round(max(a_bp, b_bp) / BPP)), P.W)
    inb = (P.gpk >= i0) & (P.gpk < i1)
    seg = qbp[inb]; nlab0 = int(inb.sum())
    L = b_bp - a_bp; Reff = min(R, L)
    Wbp = P.W * BPP
    if nlab0 > 0:
        first, last = float(seg[0]), float(seg[-1])
        prevL = float(qbp[qbp < first].max()) if np.any(qbp < first) else 0.0
        nextL = float(qbp[qbp > last].min()) if np.any(qbp > last) else Wbp
        lmin, lmax = max(a_bp - Reff, prevL + 1), min(a_bp + Reff, first - 1)
        rmin, rmax = max(b_bp - Reff, last + 1),  min(b_bp + Reff, nextL - 1)
    else:
        lmin, lmax = max(a_bp - Reff, 0.0), a_bp + Reff
        rmin, rmax = b_bp - Reff, min(b_bp + Reff, Wbp)
    return lmin, lmax, rmin, rmax, nlab0

def inversion_scan(P, a_bp, b_bp, slide_bp=0, R=ENDPOINT_R,   # scan radius +/-R (same for obs and null)
                   coarse=STEP_COARSE, fine=STEP_FINE):
    """Scan breakpoints [a,b] within +/-R, inverting each candidate in place, and pick the one
       whose inverted labels best match the reference (max galign_inv). R<=0 keeps the marked
       segment (no scan). Returns region_score dict + l_bp/r_bp + dreg (= inv - fwd)."""
    if R <= 0:                                        # keep the marked segment (no scan)
        rs = region_score(P, a_bp, b_bp, slide_bp)
        if rs is None: return None
        return dict(l_bp=float(min(a_bp, b_bp)), r_bp=float(max(a_bp, b_bp)),
                    dreg=rs['score_inv'] - rs['score_fwd'], **rs)
    lmin, lmax, rmin, rmax, nlab0 = _neighbor_bounds(P, a_bp, b_bp, R)

    def grid(lo, hi, step):
        if hi <= lo: return np.array([np.clip((lo + hi) / 2, 0.0, P.W * BPP)])
        return np.arange(lo, hi + 1e-6, step)

    def scan(lg, rg):
        best = None
        for l in lg:
            for r in rg:
                if r - l < MIN_SEG_WIDTH: continue
                if _nlab_in(P, l, r) != nlab0: continue   # label preservation
                rs = region_score(P, float(l), float(r), slide_bp)
                if rs is None: continue
                obj = rs['galign_inv']                   # objective: inverted labels vs reference
                if best is None or obj > best[0]:
                    best = (obj, float(l), float(r), rs)
        return best

    b = scan(grid(lmin, lmax, coarse), grid(rmin, rmax, coarse))
    if b is None:                                     # coarse failed -> keep the mark
        rs = region_score(P, a_bp, b_bp, slide_bp)
        if rs is None: return None
        b = (rs['galign_inv'], float(a_bp), float(b_bp), rs)
    _, lc, rc, _ = b
    b2 = scan(grid(max(lmin, lc - coarse), min(lmax, lc + coarse), fine),
              grid(max(rmin, rc - coarse), min(rmax, rc + coarse), fine))
    _, l_bp, r_bp, rs = (b2 or b)
    dreg = rs['score_inv'] - rs['score_fwd']
    out = dict(l_bp=l_bp, r_bp=r_bp, dreg=dreg, **rs)
    return out


# ============================ null (cross-molecule) ============================
_POOL = None
POOL_CACHE = config.POOL_CACHE
def _pack(P):   # exclude ref_red/ref_pk (reloaded from chrom)
    return dict(stem=P.stem, chrom=P.chrom, aligned_bp=P.aligned_bp, direction=P.direction,
                scale=P.scale, W=P.W, mol_red=P.mol_red, gpk=P.gpk)
def _unpack(d):
    P = Anchor()
    for k, v in d.items(): setattr(P, k, v)
    P.ref_red = LA.ref_red(P.chrom); P.ref_pk, P.ref_len = LA.ref_peaks(P.chrom)
    P.cc_r = float('nan'); P.logp = float('nan')
    return P
def build_pool(cap=POOL_CAP, verbose=True):
    global _POOL
    if _POOL is not None: return _POOL
    import glob, pickle
    if os.path.exists(POOL_CACHE):
        try:
            packs = pickle.load(open(POOL_CACHE, 'rb'))
            if len(packs) >= min(cap, 50):
                _POOL = [_unpack(d) for d in packs]
                if verbose: print(f'[pool] loaded {len(_POOL)} from cache', flush=True)
                return _POOL
        except Exception: pass
    fs = sorted(glob.glob(POOL_DIR + '/*_matching_results.csv'))
    stems = [os.path.basename(f)[:-len('_matching_results.csv')] for f in fs]
    if len(stems) > cap:
        idx = np.unique(np.linspace(0, len(stems) - 1, cap).round().astype(int))
        stems = [stems[i] for i in idx]
    Ps = []
    for s in stems:
        try: Ps.append(place(s, base_dir=POOL_DIR))
        except Exception: pass
    _POOL = Ps
    try: pickle.dump([_pack(P) for P in Ps], open(POOL_CACHE, 'wb'))
    except Exception: pass
    if verbose: print(f'[pool] placed {len(Ps)} molecules from {POOL_DIR}', flush=True)
    return Ps

# ============================ diagnostic data (for visualization) ============================
def _mapped_ref_arr(ref_prof, n, ab_px, scale):
    """Sample ref_prof at n points over [ab_px, ab_px+scale*n) (dir=+1). For plotting."""
    x = ab_px + scale * np.arange(n)
    N = len(ref_prof); xi = np.clip(x, 0, N - 1)
    lo = xi.astype(int); frac = xi - lo; hi = np.minimum(lo + 1, N - 1)
    return ref_prof[lo] * (1 - frac) + ref_prof[hi] * frac

def diag_data(P, obs, null_T, T_obs, slide_bp=0, green_d=None, red_d=None, whole_d=None):
    """At the chosen breakpoint: segment red (mol fwd / mol inv / ref) + labels + null
       distributions (region and whole). In place (offset 0): fwd/inv share the same
       position, so there is a single ref."""
    i0, i1 = obs['i0'], obs['i1']; w_reg = i1 - i0
    ab_reg = obs['rlo_px'] * BPP; off = obs['off']
    seg = P.mol_red[i0:i1]
    R_reg = seg if P.direction == 1 else seg[::-1]
    ref_o = _mapped_ref_arr(P.ref_red, w_reg, (ab_reg + off) / BPP, P.scale)
    def _n(a):
        a = np.asarray(a, float); mn, mx = a.min(), a.max()
        return ((a - mn) / (mx - mn)) if mx > mn else a * 0
    dp = P.gpk[(P.gpk >= i0) & (P.gpk < i1)]
    loc = (dp - i0) if P.direction == 1 else ((i1 - 1) - dp)   # ref-increasing local px
    rp = P.ref_pk[(P.ref_pk >= (ab_reg + off) / BPP - 3) & (P.ref_pk <= (ab_reg + off) / BPP + P.scale * w_reg + 3)] - (ab_reg + off) / BPP
    return dict(
        w_reg=int(w_reg), scale=float(P.scale), off_bp=float(off),
        dir=int(P.direction), i0=int(i0), i1=int(i1),
        mol_fwd=[round(float(v), 4) for v in _n(R_reg)],
        mol_inv=[round(float(v), 4) for v in _n(R_reg[::-1])],
        ref=[round(float(v), 4) for v in _n(ref_o)],
        data_labels_px=[round(float(x), 1) for x in loc],
        data_labels_inv_px=[round(float(w_reg - 1 - x), 1) for x in loc],
        ref_labels_px=[round(float(x), 1) for x in rp],
        green_d=(None if green_d is None else round(float(green_d), 4)),
        red_d=(None if red_d is None else round(float(red_d), 4)),
        whole_d=(None if whole_d is None else round(float(whole_d), 4)),
        T=(None if T_obs is None else round(float(T_obs), 4)),
        null_T=[round(float(x), 4) for x in null_T])

# ============================ verdict ============================
def analyze_inversion(stem, a_bp, b_bp, group='M_group', slide_bp=0, R=ENDPOINT_R, do_null=True, dump_dir=None, base_dir=None):   # R = +/-10kb breakpoint search
    P = place(stem, group, base_dir=base_dir)
    # hybrid localization over (center, width, gamma). Lazy import avoids a circular import;
    # hybrid_scan_A falls back to inversion_scan when the region has no labels.
    import hybrid_localize as _HL
    obs = _HL.hybrid_scan_A(P, a_bp, b_bp, slide_bp=slide_bp, R=R)
    if obs is None:
        return dict(stem=stem, verdict='n/a (segment too short)')
    wf = whole_score(P, slide_bp=slide_bp)                                    # forward whole
    wi = whole_score(P, obs['l_bp'], obs['r_bp'], slide_bp)                   # segment inverted in place
    green_d = obs['galign_inv'] - obs['galign_fwd']  # green Delta = galign (reversed vs forward)
    red_d   = obs['gcc_inv'] - obs['gcc_fwd']        # red (backbone) Delta
    dwhole = wi['score'] - wf['score']
    res = dict(stem=stem, chrom=P.chrom, nlab=obs['nlab'],
               l_bp=int(obs['l_bp']), r_bp=int(obs['r_bp']), w_kb=round((obs['r_bp']-obs['l_bp'])/1000, 1),
               gamma=round(float(obs.get('gamma', 1.0)), 3),   # local stretch (1.0 = global scale)
               region_fwd=round(obs['score_fwd'], 4), region_inv=round(obs['score_inv'], 4),
               gcc_fwd=round(obs['gcc_fwd'], 3), gcc_inv=round(obs['gcc_inv'], 3),
               f1_fwd=round(obs['f1_fwd'], 3), f1_inv=round(obs['f1_inv'], 3),
               galign_fwd=round(obs['galign_fwd'], 3), galign_inv=round(obs['galign_inv'], 3),
               green_d=round(green_d, 4), red_d=round(red_d, 4), dreg=round(obs['dreg'], 4),
               whole_fwd=round(wf['score'], 4), whole_inv=round(wi['score'], 4),
               dwhole=round(dwhole, 4))
    if not do_null:
        res['verdict'] = '(null not run)'; return res
    # verdict = matched null-of-T: T (hybrid max-z) tested against null pool windows run through
    # the full hybrid, keyed by marked width. Green does not enter the verdict (green_d shown only).
    null_T = _HL.null_bank_T(b_bp - a_bp)
    n = len(null_T)
    T_obs = obs.get('z')
    p_T = (1 + int(np.sum(null_T >= T_obs))) / (1 + n) if (n and T_obs is not None) else float('nan')
    p_conf = p_T
    res.update(n_null=n, T=round(float(T_obs), 4) if T_obs is not None else None,
               p_T=round(p_T, 4) if n else None, p_conf=round(p_conf, 4) if n else None)
    if not (n >= 5):
        verdict = 'n/a (insufficient null)'
    elif p_conf <= ALPHA and red_d > 0:                  # present: p_T <= ALPHA and red_d > 0
        verdict = 'present'
    else:
        verdict = 'rejected'
    res['verdict'] = verdict
    res['a_bp'] = int(a_bp); res['b_bp'] = int(b_bp); res['slide_bp'] = int(slide_bp)
    if dump_dir:
        try:
            dd = diag_data(P, obs, null_T, T_obs, slide_bp,
                           green_d=green_d, red_d=red_d, whole_d=dwhole)
            json.dump({**{k: v for k, v in res.items()}, 'sv': 'inversion', 'diag': dd},
                      open(os.path.join(dump_dir, f'{stem}_svla.json'), 'w'))
        except Exception as e:
            print(f'  [svla dump failed] {stem}: {e}')
    return res
