"""Junction search over the continuous molecule with fixed loci.
   substrate = FIG4/{stem}.tif (Wmol px); pieces left=[0:WL], right=[Wmol-WR:Wmol].
   Each piece's locus is fixed; ab_full gives the full-molecule anchor and _cont_aligned_bp
   the sub-segment aligned_bp. jx scans mark_px +/- SEARCH_KB; objective = per-side whole-piece
   global_ccr (sum of left/right).
   Usage: python3 junction_search.py <stem> <Lchrom> <Lbp> <Ldir> <Lscale> <Rchrom> <Rbp> <Rdir> <Rscale> <mark_bp>
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import tempfile
import tifffile
os.environ.setdefault('SVLA_TIFDIR', config.SVLA_TIFDIR); os.environ.setdefault('SVLA_BASE', config.SVLA_BASE)
import sv_la as S, sv_la_transloc as T
BPP = S.BPP
FIG4 = config.FIG4
SEARCH_KB = 10

def piece_gccr(seg, chrom, ab, di, sc):
    """seg=(2,H,w) piece array. Whole-piece global_ccr (single Pearson over the piece)."""
    fd = os.path.join(tempfile.gettempdir(), f'js_{os.getpid()}_{seg.shape[-1]}_{int(ab)%100000}.tif')
    tifffile.imwrite(fd, seg)
    P = S.make_anchor(fd, chrom, int(round(ab)), int(di), float(sc)); os.remove(fd)
    g = S.LA.global_ccr(P.mol_red, P.ref_red, int(round(ab)), float(sc), int(di))
    return 0.0 if g != g else float(g)

def ab_full(ab_piece, di, sc, mstart, Wp, Wmol):
    """Convert a piece's (mstart..mstart+Wp) aligned_bp to the full-molecule (Wmol) make_anchor anchor."""
    if di == 1:
        return ab_piece - sc * mstart * BPP
    return ab_piece + sc * (Wp + mstart - Wmol) * BPP

def piece_width(stem, side):
    for tag in [f'_{side}', f'_{side.capitalize()}']:
        p = f'{FIG4}/{stem}{tag}.tif'
        if os.path.exists(p):
            return tifffile.imread(p).shape[-1]
    raise FileNotFoundError(f'{stem} {side}')

def search(stem, L, R, mark_bp):
    arr = tifffile.imread(f'{FIG4}/{stem}.tif'); Wmol = arr.shape[-1]
    WL = piece_width(stem, 'left'); WR = piece_width(stem, 'right'); mstartR = Wmol - WR
    Lc, Lab, Ldi, Lsc = L; Rc, Rab, Rdi, Rsc = R
    abfL = ab_full(Lab, Ldi, Lsc, 0, WL, Wmol)
    abfR = ab_full(Rab, Rdi, Rsc, mstartR, WR, Wmol)
    mpx = int(round(mark_bp / BPP)); rad = int(SEARCH_KB * 1000 / BPP)
    lo = max(S.TRIM_PX + 30, mpx - rad); hi = min(Wmol - S.TRIM_PX - 30, mpx + rad)
    rows = []
    for jx in range(lo, hi + 1):
        labL = T._cont_aligned_bp(abfL, Ldi, Lsc, Wmol, 0, jx)
        cL = piece_gccr(arr[:, :, 0:jx], Lc, labL, Ldi, Lsc)
        labR = T._cont_aligned_bp(abfR, Rdi, Rsc, Wmol, jx, Wmol - jx)
        cR = piece_gccr(arr[:, :, jx:Wmol], Rc, labR, Rdi, Rsc)
        rows.append((jx, cL, cR, cL + cR))
    return Wmol, mpx, rows, (WL, WR, mstartR)

if __name__ == '__main__':
    a = sys.argv
    stem = a[1]; L = (a[2], float(a[3]), int(a[4]), float(a[5])); R = (a[6], float(a[7]), int(a[8]), float(a[9])); mark = float(a[10])
    Wmol, mpx, rows, geom = search(stem, L, R, mark)
    best = max(rows, key=lambda r: r[3])
    print(f"{stem}  Wmol={Wmol}px mark={mpx}px  geom(WL,WR,mstartR)={geom}")
    print(f"  L={L[0]}@{L[1]:.0f}(s{L[3]},d{L[2]})  R={R[0]}@{R[1]:.0f}(s{R[3]},d{R[2]})")
    print(f"{'jx':>4} {'kb':>6} {'L_gcc':>7} {'R_gcc':>7} {'sum':>7}")
    for jx, cL, cR, s in rows:
        if jx % 3 == 0 or jx == best[0]:
            print(f"{jx:4d} {jx*BPP/1000:6.1f} {cL:7.3f} {cR:7.3f} {s:7.3f}{'  *' if jx==best[0] else ''}")
    print(f"\nbest junction: jx={best[0]} ({best[0]*BPP/1000:.1f}kb)  L_gcc={best[1]:.3f} R_gcc={best[2]:.3f} sum={best[3]:.3f}  (vs mark {(best[0]-mpx)*BPP/1000:+.1f}kb)")
