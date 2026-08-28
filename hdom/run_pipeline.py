#!/usr/bin/env python3
"""hDOM single-molecule optical-map pipeline -- one entry point.

Give it a folder of molecule TIFFs and a reference map; it returns, for every molecule, the
ranked chromosomal alignment from both scorers (hDOM red+label, maligner label-only) plus a
self-contained HTML view.

    python3 run_pipeline.py --molecules RAW_DIR --ref-map REF_MAP_DIR --output OUT_DIR

RAW_DIR       2-frame 16-bit TIFFs, frame 0 = red (AT-dye DNA backbone), frame 1 = green (nick
              label).  Either `{name}.tif` or `{name}_orig.tif`; both are accepted.
REF_MAP_DIR   one `chr*.tif` per reference chromosome, same 2-channel convention, at --bpp
              bp/pixel (the in-silico reference maps).
OUT_DIR       per molecule: `_matching_results.csv` (ranked hits, both scorers),
              `_view.html` (alignment view), `.tif` (straightened), plus `index.php`.

Stages, all inside hdom/straighten_align.py:
  1  ML label detection -> backbone trace -> arc-length straightening -> label snap -> widthprop
  2  candidate generation (hdom/candidates/hDOM_calc_Q.py: pcc top-10 by cc_r*cc_g**2 over every
     chromosome, and maligner_dp top-10 over the label maps -- both on the full label set)
     -> local re-rank of the pooled candidates (local_cc_r x label F1 -> null log P)
  2b the maligner row reported per molecule: a second maligner_dp run with the terminal-5 kb
     labels dropped (see README, "The terminal 5 kb")
  3  HTML view + before/after panels + deploy

Requires `maligner_dp` on PATH (Mendelowitz et al., maligner) and python3 with
numpy/scipy/pandas/tifffile/Pillow/scikit-learn/joblib.
"""
import os, sys, glob, shutil, argparse, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, 'hdom')


def die(msg):
    """Stop with a message.

    Used for the pre-flight checks below. Every one of them guards a dependency whose absence would
    otherwise produce a run that LOOKS successful -- fewer candidates, or none -- rather than an
    error, so they are all checked before any molecule is touched.
    """
    sys.exit(f'error: {msg}')


def main():
    """Check the environment, normalise the input folder, then run straighten_align over shards.

    Nothing scientific happens here. The work is in hdom/straighten_align.py; this wrapper exists so
    that a folder of TIFFs and a reference map are all a user needs to supply, and so that a missing
    GPU, a missing maligner_dp or a missing genome.maps stops the run with an explanation instead of
    quietly producing a partial candidate pool.
    """
    ap = argparse.ArgumentParser(description='hDOM straighten + align pipeline')
    ap.add_argument('--molecules', required=True, help='folder of molecule TIFFs')
    ap.add_argument('--ref-map', required=True, help='folder of reference chr*.tif')
    ap.add_argument('--output', required=True, help='destination folder')
    ap.add_argument('--shards', type=int, default=1, help='parallel worker processes (default 1)')
    ap.add_argument('--bpp', type=int, default=200,
                    help='bases per reference pixel; must match how --ref-map was built (default 200)')
    ap.add_argument('--logp', action='store_true',
                    help='also compute the K=5000 random-null log(P) (h_score) per molecule')
    ap.add_argument('--chains', default='', help='optional minimap2 haplotype PAF dir (homolog panel)')
    ap.add_argument('--keep-staged', action='store_true', help='keep the normalised input folder')
    a = ap.parse_args()

    if not shutil.which('maligner_dp'):
        die('maligner_dp not found on PATH -- the label-only scorer cannot run.')
    refs = sorted(glob.glob(os.path.join(a.ref_map, 'chr*.tif')))
    if not refs:
        die(f'no chr*.tif in --ref-map {a.ref_map}')
    # maligner_dp reads the reference as fragment lengths, not pixels, so the folder needs
    # genome.maps beside the TIFFs. make_reference.py writes both.
    if not os.path.exists(os.path.join(a.ref_map, 'genome.maps')):
        die(f'{a.ref_map}/genome.maps is missing.\n'
            '       It holds one line per reference: name, length_bp, then the inter-label\n'
            '       fragment lengths. maligner_dp needs it; the chr*.tif alone are not enough.')
    # The pcc scanner runs on the GPU and exits when CuPy cannot load. hDOM_calc_Q does not
    # distinguish that from a scanner that found nothing, so the run would finish with maligner
    # candidates only and no error.
    try:
        import cupy  # noqa: F401
    except Exception as e:
        die(f'CuPy is unavailable ({e}).\n'
            '       The pcc stage needs an NVIDIA GPU; without it only maligner candidates\n'
            '       would be generated. Install the CuPy wheel matching your CUDA version.')
    if not os.path.exists(os.path.join(PKG, 'label_model.joblib')):
        die('hdom/label_model.joblib missing (ML label detector).')
    # hDOM_calc_Q runs these by path rather than importing them, and skips a candidate
    # generator whose script is absent, so verify all of them before starting.
    for _s in ('hDOM_calc_Q.py', 'hDOM_pcc_ft_Q.py', 'hDOM_mal_Q.py', 'DOM_lib_Q.py'):
        if not os.path.exists(os.path.join(PKG, 'candidates', _s)):
            die(f'hdom/candidates/{_s} missing -- candidate generation would be incomplete.')

    # straighten_align consumes `{mol}_orig.tif`; accept plain `{mol}.tif` by staging a copy so a
    # user does not have to rename their data to match an internal convention.
    src = sorted(p for p in glob.glob(os.path.join(a.molecules, '*.tif'))
                 + glob.glob(os.path.join(a.molecules, '*.tiff')))
    if not src:
        die(f'no TIFFs in --molecules {a.molecules}')
    # Rebuilt from scratch every run. straighten_align aligns whatever it finds here, so a
    # file left behind by an interrupted run would be aligned again alongside the real input.
    staged = os.path.join(a.output, '_input')
    shutil.rmtree(staged, ignore_errors=True)
    os.makedirs(staged)
    # Cropping can leave behind a file with a zero-size frame. Name it here rather than let it
    # surface as a numpy reduction error from inside the label detector.
    import tifffile
    n, skipped = 0, []
    for p in src:
        b = os.path.basename(p)
        try:
            shp = tifffile.TiffFile(p).series[0].shape
        except Exception as e:
            skipped.append((b, f'unreadable ({type(e).__name__})')); continue
        if len(shp) != 3 or shp[0] < 2 or min(shp[1:]) < 1:
            skipped.append((b, f'shape {tuple(shp)}, need 2 frames with non-zero size')); continue
        stem = b[:-len('_orig.tif')] if b.endswith('_orig.tif') else os.path.splitext(b)[0]
        shutil.copy2(p, os.path.join(staged, stem + '_orig.tif'))
        n += 1
    if skipped:
        print(f'skipping {len(skipped)} unusable file(s):', flush=True)
        for b, why in skipped:
            print(f'  {b}  --  {why}', flush=True)
    if not n:
        die('no usable molecules left after that.')

    env = dict(os.environ, HDOM_REFMAP=os.path.abspath(a.ref_map), HDOM_BPP=str(a.bpp))
    if a.chains:
        env['HDOM_CHAINS'] = os.path.abspath(a.chains)
    print(f'{n} molecules | {len(refs)} reference chromosomes | {a.bpp} bp/px | '
          f'{a.shards} shard(s)', flush=True)

    base = [sys.executable, os.path.join(PKG, 'straighten_align.py'),
            '--input', staged, '--output', os.path.abspath(a.output)]
    if a.logp:
        base.append('--logp')
    # Shards partition the molecule list (k/N) and write disjoint files, so they run concurrently;
    # the folder-level registries they share are merged under a lock.
    procs = [subprocess.Popen(base + ['--shard', f'{k}/{a.shards}'], env=env)
             for k in range(a.shards)]
    rc = [p.wait() for p in procs]

    done = len(glob.glob(os.path.join(a.output, '*_view.html')))
    if not a.keep_staged:
        shutil.rmtree(staged, ignore_errors=True)
    print(f'\ndone: {done}/{n} molecules -> {a.output}')
    if done < n:
        print(f'  {n - done} produced no result (no ML label detected, or backbone trace failed)')
    sys.exit(max(rc) if any(rc) else 0)


if __name__ == '__main__':
    main()
