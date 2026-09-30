#!/usr/bin/env python
"""The pcc scan: whole-molecule two-channel correlation against every reference, on the GPU.

This is the conventional half of candidate generation. A molecule is compared with each reference
at every shift, both directions, and a range of scales, and the ten best loci are handed on. What
it optimises is

    cc_rg2 = cc_r * cc_g**2

with cc_r the Pearson correlation of the red backbone profile and cc_g the mean-free cosine of the
green label profile (both from DOM_lib_Q). Green enters squared, so a locus with a good backbone
but poor label agreement is pushed down hard.

Two details of that ranking are worth knowing because they are not obvious from the formula:

  * candidates are re-ranked on cc_rg2 recomputed from the returned cc_r and cc_g, not on the raw
    FFT score, so a peak that scores well while its green has collapsed to ~0 cannot win;
  * a cc_g of exactly zero is floored to 0.001 rather than left at zero, so the product does not
    annihilate an otherwise usable candidate.

Scale search is two-stage, over bases-per-pixel rather than a ratio: a coarse pass every
--coarse-step bpp across [--start-bpp, --end-bpp], then a 1 bpp sweep around every coarse point
within --fine-threshold of the best. Both range ends are always probed, because a molecule whose
true scale sits at the edge would otherwise only be reached if a fine pass happened to be
triggered next to it.

All shifts are evaluated at once by FFT: correlating with a window of ones gives the running sums
that the Pearson denominator needs, so one transform covers every offset. CuPy is required.

Profiles are cached per file with an LRU bound (DATA_CACHE_LIMIT / MAP_CACHE_LIMIT) so that
scanning many molecules against many references does not re-read the same TIFFs or exhaust memory.
"""
import sys
import time
import os
import argparse
import numpy as np
try:
    import cupy as cp
except Exception as exc:
    print(f"Error: CuPy is required for hDOM_pcc_ft.py but is not available: {exc}")
    sys.exit(1)
from collections import OrderedDict
from scipy import signal
import DOM_lib_Q as lib


def parse_args():
    """Command line. --start-bpp/--end-bpp bound the scale search and can also come from the
    PCC_START_BPP / PCC_END_BPP environment variables.
    """
    parser = argparse.ArgumentParser(description='PCC calculation (FFT Accelerated)')
    parser.add_argument('--map', required=True, help='Reference map tif OR folder of tif maps')
    parser.add_argument('--data-folder', required=True, help='Folder containing tif data files')
    parser.add_argument('--bpp', type=int, default=200)
    parser.add_argument('--start-bpp', type=int, default=int(os.environ.get('PCC_START_BPP','181')))
    parser.add_argument('--end-bpp', type=int, default=int(os.environ.get('PCC_END_BPP','220')))
    parser.add_argument('--coarse-step', type=int, default=5)
    parser.add_argument('--fine-threshold', type=float, default=0.01)
    parser.add_argument('--circular', action='store_true',
                        help='Treat reference genomes as circular (default: linear)')
    return parser.parse_args()


# Cache limits to avoid memory blow-up when many data TIFs are used.
# Set to None for unlimited, or a small integer to cap memory usage.
DATA_CACHE_LIMIT = 20
MAP_CACHE_LIMIT = None


def _cache_get(cache, key):
    """Fetch from an LRU dict and mark the entry as most recently used. None if absent or disabled."""
    if cache is None:
        return None
    if key in cache:
        cache.move_to_end(key)
        return cache[key]
    return None


def _cache_set(cache, key, value, limit):
    """Insert into an LRU dict, evicting the oldest entries once `limit` is exceeded."""
    if cache is None:
        return
    cache[key] = value
    cache.move_to_end(key)
    if limit is not None:
        while len(cache) > limit:
            cache.popitem(last=False)


def _next_pow2(n):
    """Smallest power of two >= n, for FFT lengths."""
    return 1 << (n - 1).bit_length()


def _fft_valid_corr(x, y, num_shifts):
    """Return dot products for valid shifts using FFT (length = num_shifts)."""
    nfft = _next_pow2(len(x) + len(y) - 1)
    corr = cp.fft.ifft(cp.fft.fft(x, nfft) * cp.fft.fft(y[::-1], nfft)).real
    start = len(y) - 1
    return corr[start:start + num_shifts]


def _fft_window_sums(x, n, num_shifts):
    ones = cp.ones(n, dtype=x.dtype)
    return _fft_valid_corr(x, ones, num_shifts)


def _compute_scores_ft(xx_r, xx_g, y_r, y_g):
    """Compute scores for all shifts using FFT-based correlations."""
    L = len(xx_r)
    n = len(y_r)
    num_shifts = L - n
    if num_shifts <= 0:
        return cp.array([]), cp.array([]), cp.array([])

    # Red: Pearson correlation using sums (mean/variance), not min-max.
    xy_r = _fft_valid_corr(xx_r, y_r, num_shifts)
    sum_x_r = _fft_window_sums(xx_r, n, num_shifts)
    sum_x2_r = _fft_window_sums(xx_r * xx_r, n, num_shifts)
    sum_y_r = cp.sum(y_r)
    sum_y2_r = cp.sum(y_r * y_r)

    den_r = (n * sum_x2_r - sum_x_r * sum_x_r) * (n * sum_y2_r - sum_y_r * sum_y_r)
    den_r = cp.sqrt(cp.maximum(den_r, 0.0))
    cc_r = cp.zeros_like(xy_r)
    valid_r = den_r != 0
    cc_r[valid_r] = (n * xy_r[valid_r] - sum_x_r[valid_r] * sum_y_r) / den_r[valid_r]

    # Green: DCC (cosine similarity)
    xy_g = _fft_valid_corr(xx_g, y_g, num_shifts)
    sum_x2_g = _fft_window_sums(xx_g * xx_g, n, num_shifts)
    sum_y2_g = cp.sum(y_g * y_g)
    den_g = cp.sqrt(cp.maximum(sum_x2_g * sum_y2_g, 0.0))
    dcc_g = cp.zeros_like(xy_g)
    valid_g = den_g != 0
    dcc_g[valid_g] = xy_g[valid_g] / den_g[valid_g]

    # Zero correlations are replaced so they cannot dominate the product score
    dcc_g[dcc_g == 0.0] = 0.001

    score = cc_r * dcc_g * dcc_g
    return score, cc_r, dcc_g


def pcc_ft(xx_r, xx_g, y_r, y_g):
    """Best placement of one molecule against one reference, at one scale.

    Scores every shift in both directions and returns (shift, direction, score, cc_r, cc_g) for the
    winner. Forward and reverse scores are interleaved into one array so a single argmax picks both
    the offset and the orientation. cc_r and cc_g are then recomputed at that one placement with the
    library functions, so the reported pair matches what the rest of the pipeline would measure
    rather than the FFT approximation. A cc_g of zero is floored to 0.001.
    """
    L = len(xx_r)
    n = len(y_r)
    num_shifts = L - n
    if num_shifts <= 0:
        return 0, 1, 0.0, 0.0, 0.0

    # Move to GPU
    xx_r_gpu = cp.asarray(xx_r)
    xx_g_gpu = cp.asarray(xx_g)
    y_r_gpu = cp.asarray(y_r)
    y_g_gpu = cp.asarray(y_g)
    y_r_rev_gpu = cp.asarray(y_r[::-1].copy())
    y_g_rev_gpu = cp.asarray(y_g[::-1].copy())

    # Forward scores
    score_fwd, _, _ = _compute_scores_ft(xx_r_gpu, xx_g_gpu, y_r_gpu, y_g_gpu)
    # Reverse scores
    score_rev, _, _ = _compute_scores_ft(xx_r_gpu, xx_g_gpu, y_r_rev_gpu, y_g_rev_gpu)

    # Interleave forward/reverse like hDOM_pcc.py
    cc = cp.empty(num_shifts * 2, dtype=cp.float64)
    cc[0::2] = score_fwd
    cc[1::2] = score_rev

    best_idx = int(cp.argmax(cc).get())
    best_score = float(cc[best_idx].get()) if cc.size > 0 else 0.0

    shift = best_idx // 2
    direction = 1 if best_idx % 2 == 0 else -1

    # Compute final result for the best match only using lib functions
    target_y_r = y_r if direction == 1 else y_r[::-1]
    target_y_g = y_g if direction == 1 else y_g[::-1]
    max_r = lib.corrcoef(xx_r[shift:n + shift], target_y_r)
    max_g = lib.dcc(xx_g[shift:n + shift], target_y_g)
    if max_g == 0.0:
        max_g = 0.001

    return shift, direction, best_score, max_r, max_g


def scan_pcc(xx_r, xx_g, y_r, y_g, log, args):
    """Search scale for the best placement, coarse then fine, and return it.

    Coarse pass: every --coarse-step bpp over [--start-bpp, --end-bpp], plus both endpoints.
    Fine pass: 1 bpp steps around each coarse point whose cc_rg2 is within --fine-threshold of the
    best, which keeps the fine work near the plausible scales without committing to the single
    coarse winner. Every trial is written to `log` as a CSV row.

    Returns (shift, direction, scale, score, cc_r, cc_g) for the best cc_rg2 found.
    """
    n = len(y_r)
    best_final = {"score": -1.0}

    # Coarse & Fine search logic remains same, but uses pcc_ft
    coarse_scales = np.arange(args.start_bpp, args.end_bpp + 1, args.coarse_step)
    # Make sure the start/end boundaries are always probed (the true scale can sit
    # at the edge of the range; otherwise it is only reachable by a lucky fine trigger).
    for edge in (args.start_bpp, args.end_bpp):
        if edge not in coarse_scales:
            coarse_scales = np.append(coarse_scales, edge)
    coarse_scales = np.unique(coarse_scales)
    coarse_results = []

    # Rank candidates by cc_rg2 = cc_r * cc_g**2 (recomputed) rather than the raw
    # FFT score, so a high-score-but-green-degenerate peak (cc_g ~ 0) is not chosen.
    def _ccrg2(cr, cg):
        return cr * cg * cg

    for bpp in coarse_scales:
        s = bpp / args.bpp
        yy_r = signal.resample(y_r, int(s * n))
        yy_g = signal.resample(y_g, int(s * n))
        res = pcc_ft(xx_r, xx_g, yy_r, yy_g)
        shift, direct, score, cr, cg = res
        if cg == 0.0:
            cg = 0.001
        print(f"{s:.3f},{bpp},{shift},{direct},{score:.3f},{cr:.3f},{cg:.3f}", file=log)
        coarse_results.append((bpp, _ccrg2(cr, cg), res, s))

    coarse_results.sort(key=lambda x: x[1], reverse=True)
    best_coarse_metric = coarse_results[0][1]

    # Fine search
    fine_scales = set()
    for bpp, metric, res, s in coarse_results:
        if best_coarse_metric - metric <= args.fine_threshold:
            for val in range(bpp - args.coarse_step, bpp + args.coarse_step + 1):
                if args.start_bpp <= val <= args.end_bpp:
                    fine_scales.add(val)

    best_metric = -1.0
    final_res = None

    for bpp in sorted(list(fine_scales)):
        s = bpp / args.bpp
        yy_r = signal.resample(y_r, int(s * n))
        yy_g = signal.resample(y_g, int(s * n))
        res = pcc_ft(xx_r, xx_g, yy_r, yy_g)
        shift, direct, score, cr, cg = res
        if cg == 0.0:
            cg = 0.001
        print(f"{s:.3f},{bpp},{shift},{direct},{score:.3f},{cr:.3f},{cg:.3f}", file=log)
        metric = _ccrg2(cr, cg)
        if metric > best_metric:
            best_metric = metric
            final_res = (shift, direct, s, score, cr, cg)

    return final_res


def process_file(RGmap, data_file, args, ref_name=None, map_cache=None, data_cache=None):
    """Scan one molecule against one reference and write the per-scale trace to a .log beside it.

    Both profiles are min-max normalised before scanning. Reference and molecule profiles are held
    in the supplied LRU caches, since a full run pairs every molecule with every chromosome.
    """
    # Generate reference-specific log filename if ref_name is provided
    data_dir = os.path.dirname(data_file)
    data_name = os.path.splitext(os.path.basename(data_file))[0]
    if ref_name:
        data_log = os.path.join(data_dir, f"{data_name}_{ref_name}_pcc.log")
    else:
        data_log = data_file.replace(".tif", "_pcc.log")
    with open(data_log, 'w', encoding='utf-8') as log:
        print("scale,bpp,shift,direction,score,cc_r,cc_g", file=log)
        cached_map = _cache_get(map_cache, RGmap)
        if cached_map is not None:
            map_r, map_g = cached_map
        else:
            map_r, map_g, _ = lib.tif2RGprofile(RGmap)
            _cache_set(map_cache, RGmap, (map_r, map_g), MAP_CACHE_LIMIT)

        cached_data = _cache_get(data_cache, data_file)
        if cached_data is not None:
            data_r, data_g = cached_data
        else:
            data_r, data_g, _ = lib.tif2RGprofile(data_file)
            _cache_set(data_cache, data_file, (data_r, data_g), DATA_CACHE_LIMIT)

        clip1r = lib.normalize(map_r)
        clip1g = lib.normalize(map_g)
        orig_ref_len = len(clip1r)

        # Circular genome: concatenate reference with itself to handle boundary-spanning matches
        if args.circular:
            clip1r = np.concatenate([clip1r, clip1r])
            clip1g = np.concatenate([clip1g, clip1g])

        result = scan_pcc(clip1r, clip1g, lib.normalize(data_r), lib.normalize(data_g), log, args)

        # Circular genome: wrap shift back into original reference length
        if args.circular and result is not None:
            shift = result[0] % orig_ref_len
            result = (shift,) + result[1:]

        return [data_file] + list(result)


if __name__ == '__main__':
    args = parse_args()
    t_start = time.time()

    # Resolve reference maps (single tif or folder)
    if os.path.isdir(args.map):
        ref_files = lib.Read_TIFF_folder(args.map)
        if not ref_files:
            print(f"Error: No reference tif files found in {args.map}")
            sys.exit(1)
    else:
        ref_files = [args.map]

    # Resolve data files
    if os.path.isdir(args.data_folder):
        tif_files = lib.Read_TIFF_folder(args.data_folder)
        output_dir = args.data_folder
    else:
        tif_files = [args.data_folder]
        output_dir = os.path.dirname(args.data_folder)

    results = []
    map_cache = OrderedDict()
    data_cache = OrderedDict()
    total_tasks = len(ref_files) * len(tif_files)
    processed = 0
    last_progress_ts = time.time()
    print(f"Start hDOM_pcc_ft: {len(ref_files)} reference(s) x {len(tif_files)} data file(s) = {total_tasks} task(s)")
    for ref_file in ref_files:
        ref_name = os.path.splitext(os.path.basename(ref_file))[0]
        RGmap_path = lib.use_tif_RGmap(ref_file, quiet=True)
        for f in tif_files:
            res = process_file(
                RGmap_path,
                f,
                args,
                ref_name=ref_name,
                map_cache=map_cache,
                data_cache=data_cache
            )
            results.append(res)
            processed += 1
            now = time.time()
            if now - last_progress_ts >= 30:
                elapsed = now - t_start
                print(
                    f"Progress: {processed}/{total_tasks} "
                    f"({processed * 100.0 / max(total_tasks, 1):.1f}%) "
                    f"- {int(elapsed // 60)}m {int(elapsed % 60)}s"
                )
                last_progress_ts = now

    t_end = time.time() - t_start
    print(f"\nDone. Total Time: {int(t_end//3600)}h {int((t_end%3600)//60)}m {int(t_end%60)}s")
