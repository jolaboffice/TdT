"""Hybrid inversion breakpoint localization + matched null-of-T (joint over center, width, gamma).

For each local stretch gamma (grid 1.0 +- 0.15), within the +/-10kb (ENDPOINT_R) label-preserving band:
green re-localizes the center (max galign_inv) and red refines the symmetric half-width (max
length-normalized score z). Statistic T = max-z over (center, width, gamma), tested against a matched
null-of-T (null pool windows run through the identical hybrid, keyed by marked width).
"""
import json
import os
import numpy as np
import sv_la as SL
import config

BPP = SL.BPP


# --- length -> null(red_d) calibration curve (for width_z) ---
_CALIB = None


def load_calib(path=None):
    """Load and smooth the precomputed length->null(red_d) quantile curve."""
    global _CALIB
    if _CALIB is None:
        d = json.load(open(path or config.CALIB_LENGTHNULL))
        Ls = np.array(sorted(int(k) for k in d), float)
        med = np.array([d[str(int(L))]['med'] for L in Ls])
        p95 = np.array([d[str(int(L))]['p95'] for L in Ls])
        # median ~flat vs length: light moving-average smoothing
        med = _movavg(med, 5)
        # upper tail must be non-increasing in length: isotonic (PAVA) fit
        p95 = _isotonic_decreasing(p95)
        _CALIB = (Ls, med, p95)
    return _CALIB


def _movavg(a, w):
    a = np.asarray(a, float); k = w // 2; out = a.copy()
    for i in range(len(a)):
        lo = max(0, i - k); hi = min(len(a), i + k + 1); out[i] = a[lo:hi].mean()
    return out


def _isotonic_decreasing(y):
    """Pool-adjacent-violators fit for a non-increasing sequence (equal weights)."""
    y = np.asarray(y, float).copy()
    val = list(y); wt = [1.0] * len(y)
    i = 0
    while i < len(val) - 1:
        if val[i] < val[i + 1] - 1e-12:
            nv = (val[i] * wt[i] + val[i + 1] * wt[i + 1]) / (wt[i] + wt[i + 1])
            val[i] = nv; wt[i] += wt[i + 1]; del val[i + 1]; del wt[i + 1]
            if i > 0: i -= 1
        else:
            i += 1
    out = []
    for v, w in zip(val, wt): out.extend([v] * int(round(w)))
    return np.asarray(out[:len(y)], float)


def width_z(red_d, L_bp):
    """Length-fair score: (red_d - null_median(L)) / (null_p95(L) - null_median(L))."""
    Ls, med, p95 = load_calib()
    Lc = float(np.clip(L_bp, Ls[0], Ls[-1]))
    m = float(np.interp(Lc, Ls, med)); s = float(np.interp(Lc, Ls, p95)) - m
    return (red_d - m) / max(s, 1e-3)


# ---- local-stretch (gamma) axis: per-window correction beyond the global scale ----
GAMMAS = np.round(1.0 + np.arange(-0.15, 0.15 + 1e-9, 0.025), 4)   # 13 values, incl. 1.0 (baseline nested)


def _score_window_gamma(P, l_bp, r_bp, gamma, slide_bp=0):
    """region_score at [l_bp,r_bp] under local stretch gamma (center-anchored; eff_scale =
    P.scale*gamma, symmetric about the center; == sv_la.region_score at gamma=1.0)."""
    ab = SL._aligned_eff(P, slide_bp)
    p_l = min(l_bp, r_bp) / BPP; p_r = max(l_bp, r_bp) / BPP
    i0 = max(int(round(p_l)), 0); i1 = min(int(round(p_r)), P.W); w_reg = i1 - i0
    if w_reg < SL.WIN:
        return None
    r_a = SL._mol_bp_to_ref_px(P, i0 * BPP, ab); r_b = SL._mol_bp_to_ref_px(P, i1 * BPP, ab)
    rlo = float(min(r_a, r_b)); rhi = float(max(r_a, r_b)); rc = 0.5 * (rlo + rhi)
    eff = P.scale * gamma
    ab_o = (rc - eff * (w_reg / 2.0)) * BPP                     # center-anchored, symmetric stretch
    seg = P.mol_red[i0:i1]
    R_reg = seg if P.direction == 1 else seg[::-1]
    dp = P.gpk[(P.gpk >= i0) & (P.gpk < i1)]
    loc = (dp - i0) if P.direction == 1 else ((i1 - 1) - dp)
    nlab = int(len(loc)); labels_inv = (w_reg - 1) - loc
    s_inv, gcc_inv, f1_inv = SL._score_at(R_reg[::-1], labels_inv, ab_o, eff, P, w_reg)
    s_fwd, gcc_fwd, f1_fwd = SL._score_at(R_reg,       loc,        ab_o, eff, P, w_reg)
    galign_inv = SL._green_align(labels_inv, ab_o, eff, P, w_reg)
    galign_fwd = SL._green_align(loc,        ab_o, eff, P, w_reg)
    return dict(score_fwd=max(s_fwd, 0.0), score_inv=max(s_inv, 0.0),
                gcc_fwd=gcc_fwd, gcc_inv=gcc_inv, f1_fwd=f1_fwd, f1_inv=f1_inv,
                galign_inv=galign_inv, galign_fwd=galign_fwd, off=0.0,
                nlab=nlab, i0=i0, i1=i1, rlo_px=rlo, rhi_px=rhi, w_reg=w_reg, gamma=float(gamma))


def hybrid_scan_A(P, a_bp, b_bp, slide_bp=0, R=SL.ENDPOINT_R, gammas=GAMMAS):
    """Joint hybrid inversion localization. For each local stretch gamma, green re-localizes the center
    (max galign_inv at the marked width) and red refines the symmetric width, both scanned 1px over the
    +/-10kb (ENDPOINT_R) label-preserving band (endpoints within lmin/lmax, rmin/rmax from
    _neighbor_bounds). T = max width_z over (center, width, gamma). Green diagnostics read at gamma=1,
    red (gcc/red_d/T) at gamma*."""
    lmin, lmax, rmin, rmax, nlab0 = SL._neighbor_bounds(P, a_bp, b_bp, R)
    if nlab0 < 1:                                     # no labels -> green cannot anchor; fall back
        o = SL.inversion_scan(P, a_bp, b_bp, slide_bp=slide_bp, R=R)
        return None if o is None else {**o, 'cstar': (a_bp + b_bp) / 2.0, 'z': None, 'gamma': 1.0, 'hybrid': False}
    W0 = float(b_bp - a_bp); c0 = (a_bp + b_bp) / 2.0
    c_lo = max(lmin + W0 / 2, rmin - W0 / 2); c_hi = min(lmax + W0 / 2, rmax - W0 / 2)

    def _green_center(g):
        """green center at stretch g: max galign_inv over the +-10kb, label-preserving band."""
        bc = (None, c0)
        c = c_lo
        while c <= c_hi + 1e-6:
            if SL._nlab_in(P, c - W0 / 2, c + W0 / 2) == nlab0:
                rs = _score_window_gamma(P, float(c - W0 / 2), float(c + W0 / 2), g, slide_bp)
                if rs is not None and (bc[0] is None or rs['galign_inv'] > bc[0]):
                    bc = (rs['galign_inv'], float(c))
            c += BPP
        return bc[1]

    best = (None, 1.0, c0, W0 / 2, None)     # (z, gamma, cstar, h, rs)
    for g in gammas:
        cstar = _green_center(g)
        h_lo = max(SL.MIN_SEG_WIDTH / 2.0, cstar - lmax, rmin - cstar)
        h_hi = min(cstar - lmin, rmax - cstar)
        h = h_lo
        while h <= h_hi + 1e-6:
            l, r = cstar - h, cstar + h
            if SL._nlab_in(P, l, r) == nlab0:
                rs = _score_window_gamma(P, float(l), float(r), g, slide_bp)
                if rs is not None:
                    z = width_z(rs['gcc_inv'] - rs['gcc_fwd'], 2 * h)
                    if best[0] is None or z > best[0]:
                        best = (z, float(g), cstar, float(h), rs)
            h += BPP
    zmax, gstar, cstar, hstar, rs = best
    if rs is None:
        return None
    l_bp, r_bp = cstar - hstar, cstar + hstar
    out = {**rs}
    # green diagnostics (f1/galign) read at gamma=1; red (gcc/red_d/T) stays at gamma*
    if abs(gstar - 1.0) > 1e-9:
        g1 = _score_window_gamma(P, float(l_bp), float(r_bp), 1.0, slide_bp)
        if g1 is not None:
            for kk in ('f1_fwd', 'f1_inv', 'galign_fwd', 'galign_inv'):
                out[kk] = g1[kk]
    out.update(l_bp=l_bp, r_bp=r_bp, dreg=rs['score_inv'] - rs['score_fwd'],
               cstar=cstar, gamma=gstar, z=zmax, hybrid=True)
    return out


# --- matched null of the hybrid statistic T = max-z ---
_TBANK = {}


def null_bank_T(W0, R=SL.ENDPOINT_R, cap=1000):
    """Matched null of T. Each pool window of length W0 (>=1 label) is run through the full joint
    hybrid (center, width, gamma; max-z). Keyed by marked width W0; cached on disk
    (config.TNULL_CACHE_DIR)."""
    k = int(round(W0 / 2000.0))
    if k in _TBANK:
        return _TBANK[k]
    cf = os.path.join(config.TNULL_CACHE_DIR, 'tnull_%d.npy' % k)
    if os.path.exists(cf):
        out = np.load(cf); _TBANK[k] = out; return out
    cand = []
    for P in SL.build_pool(verbose=False):
        Wb = P.W * BPP; c = 0.0
        while c + W0 <= Wb:
            if SL._nlab_in(P, c, c + W0) >= 1:
                cand.append((P, c, c + W0))
            c += W0
    rng = np.random.RandomState(k * 131 + 7)
    if len(cand) > cap:
        cand = [cand[i] for i in rng.choice(len(cand), cap, replace=False)]
    Ts = []
    for P, ca, cb in cand:
        try:
            h = hybrid_scan_A(P, ca, cb, R=R)
        except Exception:
            continue
        if h is not None and h.get('z') is not None:
            Ts.append(float(h['z']))
    out = np.asarray(Ts, float)
    try:
        os.makedirs(config.TNULL_CACHE_DIR, exist_ok=True); np.save(cf, out)
    except Exception:
        pass
    _TBANK[k] = out
    return out
