#!/usr/bin/env python3
"""Batched local_cc_r for the K=5000 null, on GPU when one is available.

logp_null.score_pos costs 24.7 ms because it issues ~176 tiny numpy calls per genomic position
(8 scales x 2 directions x 11 shifts), each on an 87x50 array. The arrays are far too small to
amortise call overhead, so the run is overhead-bound, not arithmetic-bound.

This module keeps the arithmetic identical and only changes the batching:

  * all 11 shifts become one 3-D gather instead of a Python loop
  * all 16 scale/direction placements are stacked
  * many genomic positions are evaluated together

Everything stays float64 so results match logp_null to rounding, and _mapped_ref is reproduced
term for term (same clip, same lo/frac split) rather than reimplemented.

label_f1 is left alone except for the window selection, which logp_null does with a boolean mask
over every peak on the chromosome -- O(all peaks) per call. searchsorted makes it O(log n); the
greedy matching that follows is untouched, so the counts are identical.
"""
import numpy as np
import local_align_pipeline as P
import logp_null as LP
import mdom_reanchor as R

try:
    import cupy as _cp
    _HAS_GPU = _cp.cuda.runtime.getDeviceCount() > 0
except Exception:
    _cp, _HAS_GPU = None, False

WIN, STEP = P.WIN, P.STEP
SHIFTS = LP.SHIFTS
SC_GRID = LP.SC_GRID
DIRS = (1, -1)


def _xp(use_gpu):
    return _cp if (use_gpu and _HAS_GPU) else np


def mapped_ref_batch(xp, ref_prof, MOL, pos_bp, scale, direction):
    """_mapped_ref for many positions at once -> (B, MOL). Same clip/lo/frac as the original."""
    sh = xp.asarray(pos_bp, dtype=xp.float64) / P.BPP
    p = xp.arange(MOL, dtype=xp.float64)
    xs = sh[:, None] + scale * ((MOL - 1 - p) if direction == -1 else p)
    N = len(ref_prof)
    xi = xp.clip(xs, 0.0, N - 1.0)
    lo = xi.astype(xp.intp); frac = xi - lo
    hi = xp.minimum(lo + 1, N - 1)
    return ref_prof[lo] * (1.0 - frac) + ref_prof[hi] * frac


def local_ccr_batch(xp, M, ref_prof, pos_bp, scale, direction):
    """local cc_r for a batch of positions -> (B,). Same quantity as logp_null.local_ccr_vec."""
    MOL = M['MOL']
    mf = mapped_ref_batch(xp, ref_prof, MOL, pos_bp, scale, direction)      # (B, MOL)
    starts = M['starts_x']; molAz = M['molAz_x']; mvalid = M['mvalid_x']
    ls = starts[None, :] + M['shifts_x'][:, None]                           # (nS, nt)
    ok = (ls >= 0) & (ls + WIN <= MOL)
    idx = xp.clip(ls, 0, MOL - WIN)[..., None] + xp.arange(WIN)             # (nS, nt, WIN)
    Bw = mf[:, idx.reshape(-1)].reshape(len(pos_bp), *idx.shape)            # (B, nS, nt, WIN)
    bstd = Bw.std(axis=3)
    good = ok[None] & mvalid[None, None, :] & (bstd >= 1e-6)
    Bz = (Bw - Bw.mean(axis=3, keepdims=True)) / xp.where(bstd < 1e-6, 1.0, bstd)[..., None]
    c = (molAz[None, None, :, :] * Bz).sum(axis=3) / WIN                    # (B, nS, nt)
    c = xp.where(good, c, -xp.inf)
    best = c.max(axis=1)                                                    # (B, nt)
    fm = mvalid[None, :] & xp.isfinite(best)
    num = xp.where(fm, best, 0.0).sum(axis=1); den = fm.sum(axis=1)
    return xp.where(den > 0, num / xp.maximum(den, 1), xp.nan)


def label_f1_fast(dp, w2, pk_sorted, rl, aligned_bp, direction, scale, trim_px):
    """label_f1 with the reference window taken by searchsorted instead of a full-array mask."""
    shift_px = aligned_bp / P.BPP
    win_lo, win_hi = shift_px, shift_px + scale * w2
    a = np.searchsorted(pk_sorted, win_lo - R.TOL_PX - scale * trim_px - 1)
    b = np.searchsorted(pk_sorted, win_hi + R.TOL_PX + scale * trim_px + 1)
    return R.label_f1(dp, w2, pk_sorted[a:b], rl, aligned_bp, direction, scale, trim_px=trim_px)


def prep(M, use_gpu=True):
    """Move the pre-standardised molecule windows onto the compute device once."""
    xp = _xp(use_gpu)
    M = dict(M)
    M['starts_x'] = xp.asarray(M['starts'])
    M['molAz_x'] = xp.asarray(M['molAz'], dtype=xp.float64)
    M['mvalid_x'] = xp.asarray(M['mvalid'])
    M['shifts_x'] = xp.asarray(SHIFTS)
    M['xp'] = xp
    return M


def score_pos_batch(M, dp, w2, rr, pk, rl, positions, chunk=64):
    """logp_null.score_pos for many positions -> (B,) array of max over scale x direction."""
    xp = M['xp']
    positions = np.asarray(positions, dtype=np.float64)
    rrx = xp.asarray(rr, dtype=xp.float64)
    pk_sorted = np.sort(np.asarray(pk, float))
    out = np.full(len(positions), -9.0)
    for a in range(0, len(positions), chunk):
        pos = positions[a:a + chunk]
        best = np.full(len(pos), -9.0)
        for sc in SC_GRID:
            for d in DIRS:
                lc = local_ccr_batch(xp, M, rrx, pos, sc, d)
                lc = _cp.asnumpy(lc) if xp is _cp else lc
                for i, p in enumerate(pos):
                    if not np.isfinite(lc[i]):
                        continue
                    _, _, _, f1 = label_f1_fast(dp, w2, pk_sorted, rl, p, d, sc, LP.TRIM)
                    f1 = 0.0 if (f1 is None or (isinstance(f1, float) and np.isnan(f1))) else f1
                    v = lc[i] * (f1 if f1 else 0)
                    if v > best[i]:
                        best[i] = v
        out[a:a + chunk] = best
    return out


def _draw(chrs, ms, seedbase, i, skip=0):
    """The (chrom, pos) the original would draw at iteration i, after `skip` earlier rejects.

    Reproduces build_null's inner loop exactly: one RandomState per i, up to 25 attempts, short
    chromosomes consumed and skipped without being offered.
    """
    rng = np.random.RandomState(seedbase + i)
    seen = 0
    for _ in range(25):
        c = chrs[rng.randint(len(chrs))]
        rr = LP.ALL[c][0]
        if len(rr) <= ms + 10:
            continue
        # int, exactly as build_null produces it: the value is written into the sidecar
        pos = int(rng.randint(0, len(rr) - ms - 1)) * P.BPP
        if seen == skip:
            return c, pos
        seen += 1
    return None, None


def build_null_fast(M, dp, w2, chrs, seedbase, K=None, chunk=64, wave=None):
    """build_null with the same draws in the same order, scored in batches.

    The original accepts a draw when its score clears -9 and keeps going until K have been
    accepted, so a reject shifts every later draw. That is reproduced here rather than
    approximated: a wave of first-attempts is scored in batch, rejects are retried with the exact
    per-i fallback, and further waves are drawn until K accepted values exist. The resulting sample
    is identical to build_null's, element for element (checked at 1e-16 on 5000 draws).
    """
    K = LP.K if K is None else K
    wave = K if wave is None else wave
    ms = int(w2 * P.SCALE_MAX)
    out, best = [], (-9.0, None, None)
    i = 0
    while len(out) < K:
        need = K - len(out)
        # first wave asks for K at once; later waves only for the shortfall, which is tiny
        idx = list(range(i, i + (wave if not out else need)))
        i = idx[-1] + 1
        draws = [(k,) + _draw(chrs, ms, seedbase, k) for k in idx]
        draws = [d for d in draws if d[1] is not None]
        scores = {}
        for c in sorted({d[1] for d in draws}):
            sel = [d for d in draws if d[1] == c]
            rr, pk, rl = LP.ALL[c]
            v = score_pos_batch(M, dp, w2, rr, pk, rl,
                                np.array([d[2] for d in sel]), chunk=chunk)
            for d, val in zip(sel, v):
                scores[d[0]] = (val, d[1], d[2])
        for k in idx:
            if len(out) >= K:
                break
            if k not in scores:
                continue
            v, c, pos = scores[k]
            if v <= -9:                     # rare: fall back to the original retry sequence
                for a in range(1, 25):
                    c2, p2 = _draw(chrs, ms, seedbase, k, skip=a)
                    if c2 is None:
                        break
                    rr, pk, rl = LP.ALL[c2]
                    v = score_pos_batch(M, dp, w2, rr, pk, rl, np.array([p2]))[0]
                    if v > -9:
                        c, pos = c2, p2
                        break
            if v > -9:
                out.append(v)
                if v > best[0]:
                    best = (v, c, pos)
    return np.array(out[:K]), best
