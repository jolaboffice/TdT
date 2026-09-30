"""Translocation IUT null (cross baseline; global_ccr + galign), with junction search.
   A pool molecule gets a random locus pair (la=L, lb=R) and a random mark; the same +/-10kb
   junction search picks jx and cross_red/cross_green are computed there.
   Both halves must be >= MIN_FRAG_KB (=60kb), so only molecules with W >= 2*MIN_FRAG qualify;
   the mark is uniform on [MIN_FRAG, W-MIN_FRAG].
   Single sequential process only (parallel shards OOM).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import tempfile, json, shutil
import numpy as np, tifffile
os.environ.setdefault('SVLA_TIFDIR', config.SVLA_TIFDIR); os.environ.setdefault('SVLA_BASE', config.SVLA_BASE)
import sv_la as S, sv_la_transloc as T, run_transloc_fixed as T2, junction_search as JS
BPP = S.BPP; ca = T._cont_aligned_bp
RAD = int(JS.SEARCH_KB * 1000 / BPP)   # 50 px = +/-10 kb
MIN_FRAG_KB = 60                        # min half-fragment length floor


def wscore_g(tif, chrom, ab, di, sc):
    P = S.make_anchor(tif, chrom, int(round(ab)), int(di), float(sc))
    g = S.LA.global_ccr(P.mol_red, P.ref_red, int(round(ab)), float(sc), int(di))
    return (0.0 if g != g else float(g)), float(T2._green_overlap(P))


def trial(arr, Wmol, la, lb, m0px, wd, tag=''):
    """la, lb = (chrom, aligned_bp, dir, scale) whole-molecule pool loci. The junction is searched over m0px +/- RAD."""
    lo = max(S.TRIM_PX + 30, m0px - RAD); hi = min(Wmol - S.TRIM_PX - 30, m0px + RAD)
    if hi <= lo:
        return None
    # junction search: max sum of left/right global_ccr
    best = None
    for jx in range(lo, hi + 1):
        cL = JS.piece_gccr(arr[:, :, 0:jx], la[0], ca(la[1], la[2], la[3], Wmol, 0, jx), la[2], la[3])
        cR = JS.piece_gccr(arr[:, :, jx:Wmol], lb[0], ca(lb[1], lb[2], lb[3], Wmol, jx, Wmol - jx), lb[2], lb[3])
        if best is None or cL + cR > best[0]:
            best = (cL + cR, jx)
    jx = best[1]
    # at the best jx: own and alt scores (global_ccr + galign)
    fdL = os.path.join(wd, f'L{tag}.tif'); fdR = os.path.join(wd, f'R{tag}.tif')
    tifffile.imwrite(fdL, arr[:, :, 0:jx]); tifffile.imwrite(fdR, arr[:, :, jx:Wmol])
    goL, galoL = wscore_g(fdL, la[0], ca(la[1], la[2], la[3], Wmol, 0, jx), la[2], la[3])          # own_L @ la
    goR, galoR = wscore_g(fdR, lb[0], ca(lb[1], lb[2], lb[3], Wmol, jx, Wmol - jx), lb[2], lb[3])    # own_R @ lb
    gaL, galaL = wscore_g(fdL, lb[0], ca(lb[1], lb[2], lb[3], Wmol, 0, jx), lb[2], lb[3])            # alt_L @ lb
    gaR, galaR = wscore_g(fdR, la[0], ca(la[1], la[2], la[3], Wmol, jx, Wmol - jx), la[2], la[3])    # alt_R @ la
    os.remove(fdL); os.remove(fdR)
    return ((goL - gaL) + (goR - gaR), (galoL - galaL) + (galoR - galaR))


def calibrate(n, seed=3):
    pool = S.build_pool(verbose=True)
    rng = np.random.RandomState(seed)
    loci_pool = [(P.chrom, int(P.aligned_bp), int(P.direction), float(P.scale)) for P in pool]
    wd = tempfile.mkdtemp(prefix='ncrsrch_')
    cr, cg = [], []; done = 0
    order = list(range(len(pool))); rng.shuffle(order)
    MFP = int(MIN_FRAG_KB * 1000 / BPP)                          # min fragment floor in px
    for idx in order:
        if done >= n:
            break
        P = pool[idx]
        if P.W < 2 * MFP + 2:                                    # need W >= 2*MFP so both halves qualify
            continue
        tif = f'{S.POOL_DIR}/{P.stem}.tif'
        if not os.path.exists(tif):
            continue
        la = loci_pool[rng.randint(len(loci_pool))]; lb = loci_pool[rng.randint(len(loci_pool))]
        try:
            arr = tifffile.imread(tif); Wmol = arr.shape[-1]
            if Wmol < 2 * MFP + 2:
                continue
            m0px = int(rng.randint(MFP, Wmol - MFP + 1))         # mark uniform on [MFP, Wmol-MFP]
            r = trial(arr, Wmol, la, lb, m0px, wd, tag=f'_{done}')
            if r is None:
                continue
            cr.append(r[0]); cg.append(r[1]); done += 1
            if done % 25 == 0:
                print(f'  {done}/{n} trials...', flush=True)
        except Exception:
            continue
    shutil.rmtree(wd, ignore_errors=True)
    return [np.array(x, float) for x in (cr, cg)]


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    cr, cg = calibrate(n)
    json.dump({'cross_red_null': [round(float(x), 4) for x in cr],
               'cross_green_null': [round(float(x), 4) for x in cg],
               'n': len(cr), 'seed': 3, 'junction_search': True},
              open(os.path.join(config.NULL_DIR, 'transloc_iut_null_cross.json'), 'w'))

    def summ(name, a):
        print(f"{name:12}: median={np.median(a):+.3f} p90={np.percentile(a, 90):+.3f} "
              f"p95={np.percentile(a, 95):+.3f} p99={np.percentile(a, 99):+.3f} max={np.max(a):+.3f}")
    print(f"\n=== null n={len(cr)} (cross baseline, junction search included) ===")
    summ('cross_red', cr); summ('cross_green', cg)

    def emp_p(a, obs):
        return (1 + int((a >= obs).sum())) / (1 + len(a))
    N = len(cr)
    fp = sum(1 for i in range(N)
             if max(emp_p(cr, cr[i]), emp_p(cg, cg[i])) <= 0.05 and cr[i] > 0 and cg[i] > 0)
    print(f"self-FP (cross IUT p<=0.05 and both>0): {fp}/{N} = {fp / N:.2%}")
