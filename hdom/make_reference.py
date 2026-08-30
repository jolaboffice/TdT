#!/usr/bin/env python3
"""Build an hDOM reference map from chromosome FASTA.

The pipeline needs a reference folder holding one `{name}.tif` per chromosome plus a
`genome.maps`. This builds both from FASTA:

    python3 make_reference.py --fasta FASTA_DIR --out REF_MAP_DIR --bpp 200

FASTA_DIR    one `{name}.fasta` per chromosome (single record each). The stem becomes the
             reference name used throughout the pipeline.
REF_MAP_DIR  written: `{name}.tif` (RGB uint8, 10 identical rows) and `genome.maps`.

Recipe, at `bpp` bases per pixel and `n = ceil(len(seq)/bpp)` pixels:

  R channel   histogram of A/T positions -> low-pass -> rescaled to 5..255
              This is the AT-richness profile the red backbone image is matched against.
  G channel   histogram of the nicking-enzyme recognition sites (motif + reverse complement)
              -> low-pass -> rescaled to 5..255. The label track.
  genome.maps the same G channel reduced to fragment lengths, which is the only form
              maligner_dp reads: `name<TAB>length_bp<TAB>n_fragments<TAB>frag1<TAB>frag2...`

`--motif` defaults to Nb.BssSI (CTCGTG). Give the recognition sequence of whatever nicking
enzyme was used; the reverse complement is searched too.

Two low-pass choices are provided because they are not interchangeable -- a reference built
with one will not reproduce alignments made against the other:

  --smooth gaussian   gaussian_filter1d, sigma 1.0 on both channels (default)
  --smooth butter     butter(4, 0.4) filtfilt

Check with --validate whenever the recipe or its parameters change: two references that differ
even slightly disagree about which homologue a molecule came from.

Use `--validate EXISTING_REF_DIR` to byte-compare the output against a reference you already
have, which is the only reliable way to confirm a recipe matches one used before.
"""
import argparse, glob, os, sys
import numpy as np
from PIL import Image
from scipy import signal, ndimage

Image.MAX_IMAGE_PIXELS = None
HEIGHT = 10


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgt', 'TGCATGCA'))[::-1].upper()


def read_sequence(path):
    """FASTA sequence as upper-case bytes, newlines and the header stripped."""
    with open(path, 'rb') as f:
        raw = f.read()
    if raw[:1] == b'>':
        raw = raw[raw.find(b'\n') + 1:]
    return raw.replace(b'\n', b'').replace(b'\r', b'').upper()


def site_positions(seq, motif):
    """0-based positions of the motif and its reverse complement, deduped and sorted."""
    found = set()
    for pat in {motif.upper().encode(), revcomp(motif).encode()}:
        i = seq.find(pat)
        while i >= 0:
            found.add(i)
            i = seq.find(pat, i + 1)
    return np.array(sorted(found), dtype=np.int64)


def _binned(indicator_or_pos, n, bpp, is_positions):
    if is_positions:
        p = indicator_or_pos
        return np.bincount(np.minimum(p // bpp, n - 1), minlength=n).astype(np.float64)
    a = indicator_or_pos
    pad = n * bpp - a.size
    if pad > 0:
        a = np.concatenate([a, np.zeros(pad, a.dtype)])
    return a[:n * bpp].reshape(n, bpp).sum(axis=1).astype(np.float64)


def smooth(v, how, sigma):
    if how == 'gaussian':
        return ndimage.gaussian_filter1d(v, sigma)
    b, a = signal.butter(4, 0.4)
    return np.clip(signal.filtfilt(b, a, v), 0, None)


def scale_5_255(v):
    lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        return np.full(v.shape, 5.0)
    return (v - lo) / (hi - lo) * 250.0 + 5.0


def build_channels(seq, bpp, motif, how, sigma, gsigma):
    """R and G uint8 profiles.

    The G channel is binarised to site presence, smoothed, and then square-rooted before
    scaling. The square root is what makes overlapping sites combine sub-linearly: two sites in
    adjacent pixels reach 1.27x one isolated site rather than the 1.78x plain summation gives.
    It also widens the apparent kernel -- sqrt of a sigma-1 gaussian is a sigma-sqrt(2) one --
    so an isolated site prints the 9-pixel footprint [7,21,63,127,162,127,63,21,7].
    """
    edges = np.arange(0, len(seq) + bpp, bpp)      # trailing edge included, so n = ceil(L/bpp)
    n = len(edges) - 1
    if n < 2:
        raise ValueError('sequence too short for this bpp')
    arr = np.frombuffer(seq, dtype=np.uint8)
    at = np.flatnonzero((arr == ord('A')) | (arr == ord('T')))
    r_hist, _ = np.histogram(at, bins=edges)
    g_hist, _ = np.histogram(site_positions(seq, motif), bins=edges)

    r_sm = smooth(r_hist.astype(np.float64), how, sigma)
    g_sm = np.sqrt(smooth((g_hist > 0).astype(np.float64), how, gsigma))
    r8 = np.clip(scale_5_255(r_sm), 0, 255).astype(np.uint8)
    g8 = np.clip(scale_5_255(g_sm), 0, 255).astype(np.uint8)
    return r8, g8


def save_tif(path, r8, g8):
    w = r8.size
    row = np.zeros((w, 3), dtype=np.uint8)
    row[:, 0] = r8
    row[:, 1] = g8
    Image.frombuffer('RGB', (w, HEIGHT),
                     np.repeat(row[None], HEIGHT, axis=0).tobytes(), 'raw', 'RGB', 0, 0).save(path)


def maps_row(name, g8, bpp):
    """One genome.maps line: the G channel as inter-site fragment lengths.

    Mirrors what the pipeline's own peak reader does, so the fragments maligner_dp aligns
    against are the same sites the pcc stage correlates against.
    """
    g = g8.astype(np.float64).copy()
    # Green floor: drop anything below 10% of the way from the mean to the maximum, so find_peaks
    # does not return the smoothing shoulders around a real site. The 0.1 is inherited, not derived
    # -- it has no measurement behind it and was never tuned. It is repeated verbatim in DOM_lib.py,
    # DOM_lib_Y.py, DOM_lib_Q.py and hDOM_calc_Q.py (three times each), and changing it here without
    # changing them would put genome.maps out of step with the profiles the pcc stage correlates
    # against. See METHODS.md section 7.
    g[g < g.mean() + (g.max() - g.mean()) * 0.1] = 0
    peaks, _ = signal.find_peaks(g)
    total = g8.size * bpp
    sites = np.concatenate(([0], peaks.astype(np.int64) * bpp, [total]))
    frags = np.concatenate((np.diff(np.concatenate(([0], sites))), [total - sites[-1]]))
    return '\t'.join([name, str(total), str(len(frags))] + [str(int(x)) for x in frags]), len(peaks)


def main():
    ap = argparse.ArgumentParser(description='Build an hDOM reference map from FASTA')
    ap.add_argument('--fasta', required=True, help='folder of {name}.fasta')
    ap.add_argument('--out', required=True, help='destination reference folder')
    ap.add_argument('--bpp', type=int, default=200, help='bases per pixel (default 200)')
    ap.add_argument('--motif', default='CTCGTG',
                    help='nicking enzyme recognition site (default CTCGTG = Nb.BssSI)')
    ap.add_argument('--smooth', choices=('gaussian', 'butter'), default='gaussian')
    ap.add_argument('--sigma', type=float, default=1.0,
                    help='gaussian sigma for the R (AT) channel, pixels (default 1.0)')
    ap.add_argument('--gsigma', type=float, default=1.0,
                    help='gaussian sigma for the G (label) channel, pixels (default 1.0)')
    ap.add_argument('--only', nargs='*', help='subset of names to build')
    ap.add_argument('--validate', help='existing reference folder to byte-compare against')
    a = ap.parse_args()

    out = os.path.expanduser(a.out)
    os.makedirs(out, exist_ok=True)
    names = a.only or sorted(os.path.basename(p)[:-len('.fasta')]
                             for p in glob.glob(os.path.join(os.path.expanduser(a.fasta), '*.fasta')))
    if not names:
        sys.exit(f'no *.fasta in {a.fasta}')

    rows, ok = [], True
    for name in names:
        fa = os.path.join(os.path.expanduser(a.fasta), name + '.fasta')
        if not os.path.exists(fa):
            print(f'  !! missing {fa}', file=sys.stderr)
            continue
        seq = read_sequence(fa)
        r8, g8 = build_channels(seq, a.bpp, a.motif, a.smooth, a.sigma, a.gsigma)
        tif = os.path.join(out, name + '.tif')
        save_tif(tif, r8, g8)
        row, npk = maps_row(name, g8, a.bpp)   # bare name: what the aligners key on
        rows.append(row)
        print(f'{name}: {len(seq):,} bp -> {r8.size:,} px  sites={npk:,}', flush=True)
        if a.validate:
            ref = os.path.join(os.path.expanduser(a.validate), name + '.tif')
            same = os.path.exists(ref) and open(ref, 'rb').read() == open(tif, 'rb').read()
            print(f'   VALIDATE tif_identical={same}')
            ok = ok and same
        del seq, r8, g8

    # A partial build writes genome.maps.partial, so it cannot replace a complete genome.maps
    # with one covering only the chromosomes just built.
    fn = 'genome.maps' if not a.only else 'genome.maps.partial'
    with open(os.path.join(out, fn), 'w', encoding='utf-8') as f:
        f.write('\n'.join(rows) + '\n')
    print(f'\n{len(rows)} references -> {out}/{fn}')
    if a.validate:
        print('ALL IDENTICAL' if ok else 'MISMATCH')


if __name__ == '__main__':
    main()
