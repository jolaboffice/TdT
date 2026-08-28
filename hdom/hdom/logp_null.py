"""Shared log_p (h_score) computation — K=5000 random-position null (Null A).

Extracted from full_logp_run.py so BOTH local_align_pipeline.py and straighten_align.py
can compute the same per-molecule log(P) inline (instead of a separate post-process).

Method (identical to the sidecar generator):
  - null = K random genome positions on OTHER chromosomes; each scored by score_pos
    (max over scale-grid x direction of  local_cc_r x label-F1).
  - for each local rank: obs=score_pos at its aligned_bp; z=(obs-mu)/sd;
    log_p = log10(sf(z))  (<=0, per-position);  log_E = log_p + log10(N_independent).
  - CSV "log(P)" for LOCAL rows := log_p (display).  mal/pcc rows are left untouched.
  - writes {base}_vsrandom.txt sidecar (keeps both log_p and log_E, plus null stats).

Usage:
  import logp_null as LP;  LP.load_refs()      # once (loads all ref chr profiles)
  LP.compute_logp(base, folder)                # per molecule (reads/writes folder/{base}_*)
"""
import os, csv, glob, json, io, hashlib
import numpy as np
from scipy import stats
import local_align_pipeline as P
import DOM_lib_Y as lib
import mdom_reanchor as R

TRIM = 25
K = 5000
WIN, STEP, MAXS = P.WIN, P.STEP, P.MAXS
SHIFTS = np.arange(-MAXS, MAXS + 1, 3)
SC_GRID = list(np.arange(P.SCALE_MIN, P.SCALE_MAX + 1e-9, P.SCALE_STEP))
REF_GLOB = os.path.join(P.REFMAP, 'chr*.tif')

ALL = {}
GENOME_BP = 0.0


def load_refs():
    """Load every reference chromosome profile once into ALL {chr:(red,peaks,len)}; sum GENOME_BP."""
    global GENOME_BP
    if ALL:
        return
    for f in glob.glob(REF_GLOB):
        c = os.path.basename(f)[:-4]
        rr = P.ref_red(c)
        if rr is None:
            continue
        GENOME_BP += len(rr) * P.BPP
        pk, rl = P.ref_peaks(c)
        if pk is not None:
            ALL[c] = (rr, np.asarray(pk, float), rl)


def prep_mol(mol_red):
    """Pre-standardise the molecule's windows once, since the null re-scores it 5000 times.

    Cuts the red profile into WIN-wide windows every STEP and z-scores each one, which is what a
    per-window Pearson does internally. Doing it here means the inner loop only has to standardise
    the reference side. Windows with no variation are flagged in `mvalid` and skipped later.
    """
    MOL = len(mol_red)
    starts = np.arange(0, MOL - WIN, STEP)
    molA = mol_red[starts[:, None] + np.arange(WIN)].astype(float)
    mstd = molA.std(1)
    mvalid = mstd >= 1e-6
    molAz = np.zeros_like(molA)
    molAz[mvalid] = (molA[mvalid] - molA[mvalid].mean(1, keepdims=True)) / mstd[mvalid, None]
    return dict(MOL=MOL, starts=starts, molAz=molAz, mvalid=mvalid, nt=len(starts))


def local_ccr_vec(M, ref_prof, aligned_bp, scale, direction):
    """local cc_r at one placement, vectorised over windows: the same quantity as
    local_align_pipeline.local_ccr, arranged so all windows are evaluated per shift at once.

    For each of the shifts in SHIFTS the whole stack of reference windows is gathered, z-scored and
    dotted against the pre-standardised molecule stack, keeping a running per-window maximum. The
    mean of those maxima is returned; NaN if no window survives.
    """
    mf = P._mapped_ref(ref_prof, M['MOL'], aligned_bp, scale, direction)
    L = len(mf)
    starts = M['starts']; molAz = M['molAz']; mvalid = M['mvalid']
    best = np.full(M['nt'], -np.inf)
    for s in SHIFTS:
        ls = starts + s
        ok = (ls >= 0) & (ls + WIN <= L)
        if not ok.any():
            continue
        oki = np.where(ok)[0]; idx = ls[oki]
        B = mf[idx[:, None] + np.arange(WIN)]; bstd = B.std(1)
        mv = mvalid[oki] & (bstd >= 1e-6)
        if not mv.any():
            continue
        Bs = B[mv]; Bz = (Bs - Bs.mean(1, keepdims=True)) / bstd[mv][:, None]
        c = (molAz[oki][mv] * Bz).sum(1) / WIN
        t = oki[mv]; best[t] = np.maximum(best[t], c)
    fm = mvalid & np.isfinite(best)
    return float(best[fm].mean()) if fm.any() else np.nan


def score_pos(M, dp, w2, rr, pk, rl, pos):
    """The pipeline's ranking score at one genomic position: max over scale and direction of
    local cc_r x label F1.

    This is the function the null is built from, so a decoy is scored exactly as a real candidate
    is -- including the free choice of scale and orientation, which a decoy also gets.
    """
    best = -9.0
    for sc in SC_GRID:
        for d in (1, -1):
            lc = local_ccr_vec(M, rr, pos, sc, d)
            if np.isnan(lc):
                continue
            _, _, _, f1 = R.label_f1(dp, w2, pk, rl, pos, d, sc, trim_px=TRIM)
            f1 = 0.0 if (f1 is None or (isinstance(f1, float) and np.isnan(f1))) else f1
            best = max(best, lc * (f1 if f1 else 0))
    return best


def build_null(M, dp, w2, chrs, seedbase):
    """Returns (null_scores, best_decoy). best_decoy=(score, chrom, pos_bp) = the single
    highest-scoring random position (= best 'wrong answer' among the sampled decoys), saved
    so a best-decoy log(P) variant can be derived later WITHOUT re-scoring/re-running pcc."""
    ms = int(w2 * P.SCALE_MAX)
    out = []; i = 0
    best = (-9.0, None, None)     # (score, chrom, pos_bp)
    while len(out) < K:
        rng = np.random.RandomState(seedbase + i); i += 1; got = False
        for _ in range(25):
            c = chrs[rng.randint(len(chrs))]
            rr, pk, rl = ALL[c]
            if len(rr) <= ms + 10:
                continue
            pos = int(rng.randint(0, len(rr) - ms - 1)) * P.BPP
            v = score_pos(M, dp, w2, rr, pk, rl, pos)
            if v > -9:
                out.append(v); got = True
                if v > best[0]:
                    best = (v, c, pos)
                break
        if not got and i > K * 3:
            break
    return np.array(out), best


def _build_null_auto(M, dp, w2, chrs, seedb):
    """build_null, routed through logp_fast when that module is importable.

    logp_fast issues the same arithmetic in batches instead of ~176 tiny array calls per genomic
    position, which is where this stage spends its time. It returns the same sample element for
    element (checked to 1e-16 on 5000 draws against build_null), so this is a speed switch and not
    a change of method; set HDOM_LOGP_FAST=0 to force the reference implementation, which is what
    that check is run against.

    Imported here rather than at module level because logp_fast imports this module.
    """
    if os.environ.get('HDOM_LOGP_FAST', '1') == '0':
        return build_null(M, dp, w2, chrs, seedbase=seedb)
    try:
        import logp_fast as LF
        return LF.build_null_fast(LF.prep(M), dp, w2, chrs, seedb)
    except Exception as e:
        print(f'  [logp] fast path unavailable ({type(e).__name__}: {e}); using build_null',
              flush=True)
        return build_null(M, dp, w2, chrs, seedbase=seedb)


def compute_logp(base, folder, dp=None, force=False):
    """K=5000 null log_p for LOCAL rows of folder/{base}_matching_results.csv.
    Writes log(P) (local rows only, = log_p) + folder/{base}_vsrandom.txt. Returns status."""
    if not ALL:
        load_refs()
    side = f'{folder}/{base}_vsrandom.txt'
    if not force and os.path.exists(side):
        try:
            if 'DONE' in open(side, encoding='utf-8').read():
                return 'skip'
        except Exception:
            pass
    tif = f'{folder}/{base}.tif'
    if not os.path.exists(tif):
        return 'notif'
    mr, mg, _ = lib.tif2RGprofile(tif)
    mol_red = np.asarray(mr, float); w2 = len(mol_red)
    if dp is None:
        mlp = f'{folder}/{base}_mlpeaks.json'
        dp = (np.array(sorted(float(x) for x in json.load(open(mlp, encoding='utf-8'))))
              if os.path.exists(mlp) else np.asarray(P._green_peaks(tif), float))
    else:
        dp = np.asarray(sorted(float(x) for x in dp), float)
    csvp = f'{folder}/{base}_matching_results.csv'
    rows = list(csv.DictReader(open(csvp, encoding='utf-8')))
    if not rows:
        return 'nocsv'
    hd = list(rows[0].keys())
    locs = [r for r in rows if r['method'] == 'local']
    if not locs:
        return 'nolocal'
    sc0 = float(locs[0]['scale'])
    truechr = locs[0]['ref_name'].rsplit('_', 1)[0]
    chrs = [c for c in ALL if c.rsplit('_', 1)[0] != truechr]
    # The null draws decoy windows from every chromosome except the one the molecule matched,
    # so a reference holding only that chromosome leaves no population to draw from.
    if not chrs:
        return f'nodecoy: ref_map has no chromosome other than {truechr} to draw a null from'
    M = prep_mol(mol_red)
    seedb = int(hashlib.md5(base.encode()).hexdigest()[:7], 16) % 1000000
    null, best_decoy = _build_null_auto(M, dp, w2, chrs, seedb)
    if len(null) < K * 0.5:
        return 'nullfail'
    mu, sd = null.mean(), (null.std() or 1e-9)
    mol_bp = w2 * sc0 * P.BPP
    N = GENOME_BP / mol_bp
    logN = np.log10(N)
    per = []
    for r in locs:
        ref = r['ref_name']
        if ref not in ALL:
            continue
        rr, pk, rl = ALL[ref]
        obs = score_pos(M, dp, w2, rr, pk, rl, int(r['aligned_bp']))
        z = (obs - mu) / sd
        lp = float(stats.norm.logsf(z) / np.log(10))
        lE = lp + logN
        r['log(P)'] = f'{lp:.2f}'   # display = log_p (<=0); mal/pcc rows untouched
        per.append((r['rank'], ref, obs, z, lp, lE))
    buf = io.StringIO(); w = csv.DictWriter(buf, fieldnames=hd); w.writeheader(); w.writerows(rows)
    tmp = csvp + '.tmp'; open(tmp, 'w', encoding='utf-8').write(buf.getvalue()); os.replace(tmp, csvp)
    txt = (f'mol={base}\nnull n={len(null)} mean={mu:.4f} sd={sd:.4f} max={null.max():.4f}\n'
           f'N_ind={N:.0f} logN={logN:.3f} mol_bp={mol_bp:.0f}\n'
           f'best_decoy chrom={best_decoy[1]} pos={best_decoy[2]} score={best_decoy[0]:.4f}\n')
    for rk, ref, obs, z, lp, lE in per:
        txt += f'rank{rk} {ref} obs={obs:.4f} z={z:.3f} log_p={lp:.3f} log_E={lE:.3f}\n'
    txt += 'DONE\n'
    tmp2 = side + '.tmp'; open(tmp2, 'w', encoding='utf-8').write(txt); os.replace(tmp2, side)
    return 'done'
