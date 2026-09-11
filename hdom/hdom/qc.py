#!/usr/bin/env python3
"""Per-molecule image quality checks, written beside the alignment as `{base}_qc.json`.

Neither check looks at the alignment. They ask whether the IMAGE can support one, which is a
different question and worth asking first: a molecule whose green labels sit barely above their
own background will still produce a rank 1, and nothing downstream will say the labels it was
ranked on are mostly noise.

  green_snr    median(height of detected green peaks) / median(green profile)
               A ratio, not (signal - background) / noise. Fails below GREEN_SNR_MIN.

  red_medroll  backbone width measured per column, then the standard deviation of that width
               over a rolling +-30 column window, then the median over columns. Large when the
               stained backbone thickens and thins along its length. Fails above RED_MEDROLL_MAX.
               Returns None when fewer than MIN_COLS columns are bright enough to trace, which is
               itself a failure: the molecule is too short or too broken to measure.

The thresholds are empirical. They were set by looking at molecules that had already been judged
by eye, not derived, so treat them as a boundary drawn around known-good data rather than a
calibrated cut. What they are good at is catching data that does not resemble that set at all.
On the 3,079 classified molecules of this study they remove 8 (0.3%), i.e. they change nothing;
on a set imaged with an unoptimised green fluorophore they removed 37 of 39 (94.9%), with green
SNR clustered at 1.2-1.8 against a median of 6.6 here. That contrast is the point of shipping
them: they tell a new user whether this pipeline is being applied to the kind of image it was
built for.

Nothing is filtered by default. `hDOM_pipeline.py` writes the sidecar and warns; passing
`--qc-filter` is what actually drops a molecule, and the alignment viewer shows a banner when a
sidecar reports a failure.
"""
import os, json
import numpy as np
import tifffile as tiff
from scipy.ndimage import median_filter

GREEN_SNR_MIN = 2.5      # green_snr below this fails
RED_MEDROLL_MAX = 0.33   # red_medroll above this fails
MIN_COLS = 70            # columns needed before red_medroll is meaningful


def green_snr(src_tif, peaks=None):
    """median peak height / median background, on the green column-max profile.

    `peaks` lets the caller pass the label set the pipeline actually used (ML detection); without
    it the fixed-threshold peaks are recomputed, which is what the original QC used.
    """
    g = np.max(tiff.imread(src_tif)[1], axis=0).astype(float)
    bg = max(float(np.median(g)), 1.0)
    if peaks is None:
        from hDOM_mal_Y import refined_green_peaks
        peaks = refined_green_peaks(src_tif)
    peaks = np.asarray([p for p in peaks if 0 <= p < len(g)], dtype=int)
    return float(np.median(g[peaks]) / bg) if peaks.size else 0.0


def red_medroll(src_tif):
    """Rolling variability of the backbone's width. None when the molecule cannot be traced."""
    red = tiff.imread(src_tif).astype(float)[0]
    H, W = red.shape
    bg = np.percentile(red, 20)
    # Centre line: intensity-weighted centroid within +-4 px of each column's brightest row,
    # then median-filtered, so a single bright speck cannot drag the trace off the backbone.
    peak = np.argmax(red, 0)
    offs = np.arange(-4, 5)
    ridx = np.clip(peak[None, :] + offs[:, None], 0, H - 1)
    ci = np.tile(np.arange(W), (9, 1))
    vals = np.clip(red[ridx, ci] - bg, 0, None)
    s = vals.sum(0)
    # A column entirely at or below background sums to zero. Divide only where it is safe and fall
    # back to that column's brightest row, rather than letting the 0/0 raise a warning per molecule.
    cent = np.divide((ridx * vals).sum(0), s, out=peak.astype(float), where=s > 0)
    sm = median_filter(cent, 15, mode='nearest')
    # Width per column, over the brighter 60% of columns only.
    bright = red.max(0) > np.percentile(red.max(0), 40)
    sp = []
    for c in np.where(bright)[0]:
        t = int(round(sm[c]))
        lo, hi = max(0, t - 7), min(H, t + 8)
        w = np.clip(red[lo:hi, c] - bg, 0, None)
        if w.sum() < 1e-6:
            continue
        pos = np.arange(lo, hi)
        m = (pos * w).sum() / w.sum()
        sp.append(np.sqrt(((pos - m) ** 2 * w).sum() / w.sum()))
    if len(sp) < MIN_COLS:
        return None
    sp = np.array(sp)
    n, h = len(sp), 30
    rs = np.empty(n)
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        rs[i] = sp[a:b].std()
    return float(np.median(rs))


def evaluate(src_tif, peaks=None):
    """Both checks plus their verdicts, as the dict written to `{base}_qc.json`."""
    snr = green_snr(src_tif, peaks=peaks)
    med = red_medroll(src_tif)
    fail_green = snr < GREEN_SNR_MIN
    fail_red = (med is None) or (med > RED_MEDROLL_MAX)
    reasons = []
    if fail_green:
        reasons.append(f'green SNR {snr:.2f} is below {GREEN_SNR_MIN} (labels are close to their '
                       f'own background, so detected label positions are unreliable)')
    if med is None:
        reasons.append(f'the red backbone could not be traced over {MIN_COLS} columns '
                       f'(molecule too short or too broken to measure)')
    elif fail_red:
        reasons.append(f'red medRoll {med:.3f} is above {RED_MEDROLL_MAX} (backbone width varies '
                       f'along the molecule, so its length scale is uneven)')
    return {
        'green_snr': round(snr, 3),
        'red_medroll': None if med is None else round(med, 4),
        'green_snr_min': GREEN_SNR_MIN,
        'red_medroll_max': RED_MEDROLL_MAX,
        'fail_green': bool(fail_green),
        'fail_red': bool(fail_red),
        'pass': not (fail_green or fail_red),
        'reasons': reasons,
    }


def write_sidecar(src_tif, out_dir=None, base=None, peaks=None):
    """Run both checks and write `{base}_qc.json`. Returns the dict."""
    q = evaluate(src_tif, peaks=peaks)
    d = out_dir or os.path.dirname(os.path.abspath(src_tif))
    b = base or os.path.basename(src_tif).rsplit('.', 1)[0].removesuffix('_orig')
    p = os.path.join(d, f'{b}_qc.json')
    tmp = p + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(q, fh, indent=1)
    os.replace(tmp, p)
    return q


if __name__ == '__main__':
    import sys, glob
    args = sys.argv[1:]
    if not args:
        print(__doc__.strip().split('\n')[0])
        print('usage: python3 qc.py MOL_orig.tif [...]   |   python3 qc.py FOLDER')
        raise SystemExit(2)
    tifs = []
    for a in args:
        tifs += sorted(glob.glob(os.path.join(a, '*_orig.tif'))) if os.path.isdir(a) else [a]
    print(f'{"molecule":52s} {"SNR":>6s} {"medRoll":>8s}  verdict')
    for t in tifs:
        q = write_sidecar(t)
        med = 'untraceable' if q['red_medroll'] is None else f'{q["red_medroll"]:8.3f}'
        tag = 'pass' if q['pass'] else 'FAIL: ' + ', '.join(
            x for x, f in (('green', q['fail_green']), ('red', q['fail_red'])) if f)
        print(f'{os.path.basename(t)[:52]:52s} {q["green_snr"]:6.2f} {med:>11s}  {tag}')
