#!/usr/bin/env python3
"""The label-agreement score F1, and the mDOM prototype that first used it.

Two things live here. `label_f1` is the one the whole pipeline depends on:
local_align_pipeline imports it as `R.label_f1` and it supplies the green half of the deployed
objective `local cc_r x F1`. `reanchor_molecule` is the earlier prototype that motivated it, kept
because the mDOM rows in the output CSV are produced this way.

The prototype in one line: take the pcc alignment as given -- red and green co-registered, red
localising well -- and replace the green score at that same placement, keeping its (shift, scale,
direction) untouched:

    F1   = 2M / (2M + FP + FN)      M matched, FP data-only peaks, FN reference-only peaks
    mDOM = max(cc_r, 0) * F1**2

What that buys over the cosine is a symmetric penalty. dcc does not subtract the mean, so it
rewards a bright label that happens to land near a reference site and charges nothing for a label
the reference lacks, or for a site the molecule failed to label. F1 charges both at the same rate
and ignores brightness entirely, since only positions enter.

Red and green stay on ONE alignment, as in the legacy pcc scan; only the green term changes.
"""
import os, sys
_REFMAP = os.environ.get('HDOM_REFMAP', os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ref_map'))
import numpy as np
import pandas as pd
from scipy import signal
import DOM_lib_Y as lib
from hDOM_mal_Y import refined_green_peaks
import polyn_filter as _pn

_N_REGIONS = _pn.load_n_regions()  # {chr: [(start_bp,end_bp)]}, cached

BPP = 200
TOL_PX = 12          # label match tolerance in reference pixels (~2.4 kb @200bpp)
REF_PROM_FRAC = 0.10 # ref green peak prominence as fraction of max
TOPK = 6             # PCC candidates to rescore

_ref_cache = {}


def ref_green_peaks(ref_file):
    """Cached (peak positions, profile length) for a reference chromosome's green channel.

    Peaks are called by prominence at REF_PROM_FRAC of the profile maximum, not by a fixed height,
    so a chromosome with a few very strong sites does not lose its ordinary ones. Cached because a
    genome-wide run asks for the same chromosome once per candidate.
    """
    if ref_file in _ref_cache:
        return _ref_cache[ref_file]
    _, g, _ = lib.tif2RGprofile(ref_file)
    g = np.asarray(g, dtype=float)
    peaks, _ = signal.find_peaks(g, prominence=REF_PROM_FRAC * (g.max() if g.max() > 0 else 1))
    _ref_cache[ref_file] = (peaks, len(g))
    _ref_cache[ref_file] = (peaks, len(g))
    return _ref_cache[ref_file]


def label_f1(data_peaks_px, w2, ref_peaks_px, ref_len, aligned_bp, direction, scale, trim_px=0):
    """Label agreement at one placement: returns (matched, FP, FN, F1).

    The molecule's label positions are mapped into reference pixels at the given alignment --
    shift + scale*x, or shift + scale*(w2-1-x) when the direction is -1 -- and then matched
    greedily to the reference sites, nearest first, each reference site usable once, accepting a
    pair only within TOL_PX (12 reference pixels, about 2.4 kb at 200 bp/pixel). What is left over
    on each side is the error: FP is a molecule label with no reference partner, FN a reference
    site the molecule did not label. F1 = 2M / (2M + FP + FN) weights the two the same.

    Only positions matter. Peak height never enters, which is the point of using this instead of
    the cosine: a label counts once whether it is faint or saturated.

    trim_px > 0 drops the terminal trim_px pixels at both ends before scoring, because arc-length
    resampling clamps any label outside the traced backbone onto pixel 0 or w2-1 and those are not
    real. The cut is applied to both sides: molecule labels outside [trim_px, w2-1-trim_px] go, and
    so do reference sites that map back outside it, or the reference would be charged FN for sites
    the molecule was never allowed to match. If trimming would empty the molecule's label set the
    trim is abandoned for that molecule -- the rule exists to drop unreliable terminal labels, not
    to leave a molecule with nothing and rank it on red alone.

    Returns (0, 0, 0, nan) when both sides are empty; the caller treats nan as no green evidence.
    """
    shift_px = aligned_bp / BPP
    s = scale
    # map data peaks (data-pixel) -> reference-pixel frame
    if direction == 1:
        mapped = shift_px + s * data_peaks_px
    else:
        mapped = shift_px + s * (w2 - 1 - data_peaks_px)
    win_lo, win_hi = shift_px, shift_px + s * w2
    _trim_ok = False
    if trim_px and trim_px > 0 and (w2 - 2 * trim_px) > 0:
        lo_d, hi_d = trim_px, w2 - 1 - trim_px
        mp = mapped[(data_peaks_px >= lo_d) & (data_peaks_px <= hi_d)]
        ref_dpx = (w2 - 1 - (ref_peaks_px - shift_px) / s) if direction == -1 else (ref_peaks_px - shift_px) / s
        rp = ref_peaks_px[(ref_dpx >= lo_d) & (ref_dpx <= hi_d)]
        # A molecule whose every label sits inside the terminal window would be scored on an
        # empty label set (F1 = 0, log P = nan) and end up ranked on red alone. Keep its labels
        # instead: the trim exists to drop unreliable terminal labels, not to erase a molecule.
        _trim_ok = not (len(mp) == 0 and len(data_peaks_px) > 0)
    if not _trim_ok:
        # reference peaks inside the aligned window (with tolerance slack)
        rp = ref_peaks_px[(ref_peaks_px >= win_lo - TOL_PX) & (ref_peaks_px <= win_hi + TOL_PX)]
        mp = mapped[(mapped >= win_lo - TOL_PX) & (mapped <= win_hi + TOL_PX)]
    if len(rp) == 0 and len(mp) == 0:
        return 0, 0, 0, np.nan
    # greedy nearest matching within tolerance
    used_ref = np.zeros(len(rp), dtype=bool)
    M = 0
    for x in mp:
        if len(rp) == 0:
            break
        d = np.abs(rp - x)
        d[used_ref] = np.inf
        j = int(np.argmin(d)) if len(d) else -1
        if j >= 0 and d[j] <= TOL_PX:
            used_ref[j] = True
            M += 1
    FP = len(mp) - M          # data peaks with no ref partner (extra labels)
    FN = len(rp) - M          # ref peaks with no data partner (missed labels)
    denom = 2 * M + FP + FN
    f1 = (2 * M) / denom if denom > 0 else np.nan
    return M, FP, FN, f1


def reanchor_molecule(base, folder, dedup_by_ref=False):
    """Re-score this molecule's top pcc candidates with F1 and return them ranked by mDOM.

    Reads the molecule's existing CSV, takes the first TOPK pcc rows, and recomputes green at each
    one's alignment. The pcc cosine values are carried through into their own columns so the CSV
    keeps both scores for the same placement rather than overwriting one with the other. Candidates
    whose span meets an assembly N-gap are flagged, not dropped, at this stage.
    """
    csv = os.path.join(folder, base + '_matching_results.csv')
    tif = os.path.join(folder, base + '.tif')
    df = pd.read_csv(csv)
    pcc = df[df['method'] == 'pcc'].head(TOPK).copy()
    if pcc.empty:
        return None

    data_peaks = refined_green_peaks(tif).astype(float)
    _, data_g, _ = lib.tif2RGprofile(tif)
    w2 = len(data_g)

    rows = []
    for _, r in pcc.iterrows():
        ref_name = r['ref_name']
        ref_file = os.path.join(_REFMAP, f'{ref_name}.tif')
        if not os.path.exists(ref_file):
            continue
        ref_peaks, ref_len = ref_green_peaks(ref_file)
        M, FP, FN, f1 = label_f1(data_peaks, w2, ref_peaks.astype(float), ref_len,
                                 float(r['aligned_bp']), int(r['dir']), float(r['scale']))
        cc_r = float(r['cc_r'])
        mdom = max(cc_r, 0.0) * (f1 ** 2) if not np.isnan(f1) else 0.0   # dDOM = max(cc_r,0)*F1^2
        is_polyn = _pn.overlaps_n(ref_name, float(r['aligned_bp']), w2, float(r['scale']), _N_REGIONS)
        rows.append({'ref_name': ref_name, 'aligned_bp': int(r['aligned_bp']), 'dir': int(r['dir']),
                     'scale': float(r['scale']),
                     'cc_r': round(cc_r, 2), 'M': M, 'FP': FP, 'FN': FN,
                     'F1': round(f1, 3) if not np.isnan(f1) else np.nan,
                     'mDOM': round(mdom, 3), 'is_polyn': bool(is_polyn),
                     # carry the PCC (cosine) green so the CSV keeps real cc_g/cc_rg/cc_rg2
                     # for this same alignment; F1/mDOM go in their own columns.
                     'cos_ccg': r.get('cc_g', np.nan), 'cos_ccrg': r.get('cc_rg', np.nan),
                     'cos_ccrg2': r.get('cc_rg2', np.nan),
                     'old_ccrg2': round(float(r['cc_rg2']), 3)})
    out = pd.DataFrame(rows).sort_values('mDOM', ascending=False).reset_index(drop=True)
    if dedup_by_ref:
        out = out.drop_duplicates(subset='ref_name', keep='first').reset_index(drop=True)
    return out


if __name__ == '__main__':
    FAILS = [
        '0202_3_4th_26_27_3',
        '0315_5_1_merged_t033_t034_mol003_121kb',
        '0317_HG002_NbBssSI_TdT_4_merged_t026_t027_mol001_267kb',
        '0329_5_merged_t023_t024_mol001_270kb',
        'NbBssSIHG002_1uM_SA_EGFP_3_merged_t024_t025_mol001_286kb',
    ]
    folder = os.environ.get('MDOM_DEMO_FOLDER', '.')
    for base in (sys.argv[1:] or FAILS):
        out = reanchor_molecule(base, folder)
        print(f'\n=== {base[:48]} ===')
        if out is None:
            print('  no pcc rows'); continue
        print(out.to_string(index=False))
        top = out.iloc[0]
        r2 = out.iloc[1]['ref_name'] if len(out) > 1 else None
        homol = (str(top['ref_name'])[:-2] == str(r2)[:-2]) if r2 else False
        print(f'  -> re-anchored mDOM rank1={top.ref_name}, rank2={r2}  homologous_pair={homol}')
