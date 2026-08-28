#!/usr/bin/env python3
"""Green-channel label-overlap statistic (galign) for translocation halves."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sv_la as S


def _green_overlap(P):
    """galign: continuous Gaussian nearest-reference-nick label overlap
       (direction-aware; sigma = REG_SIGMA_PX*sqrt(2))."""
    g = np.asarray(P.gpk, float)
    if len(g) == 0: return 0.0
    lo = P.aligned_bp / S.BPP
    qref = lo + P.scale * (g if P.direction == 1 else (P.W - 1 - g))   # label px -> ref px (direction-aware)
    rp = P.ref_pk
    if rp is None or len(rp) == 0: return 0.0
    sig = S.REG_SIGMA_PX * 1.41421356
    d = np.asarray(rp, float)[None, :] - qref[:, None]
    return float(np.exp(-(d / sig) ** 2).max(axis=1).sum())            # sum over labels of the nearest-ref overlap
