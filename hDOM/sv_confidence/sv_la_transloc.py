#!/usr/bin/env python3
"""Single-locus continuous-placement geometry for translocation half-segments."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sv_la as S

BPP = S.BPP


def _cont_aligned_bp(ab_S, dir_S, scale_S, W, p0, wh):
    """aligned_bp of the molecule-px [p0, p0+wh) sub-segment under the single-locus
       continuous placement (direction-aware)."""
    if dir_S == 1:
        return ab_S + scale_S * p0 * BPP
    return ab_S + scale_S * (W - p0 - wh) * BPP
