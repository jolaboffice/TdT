"""Confidence signals at a fixed junction (cross baseline).
   cross_red  = sum_half(own global_ccr - alt global_ccr);
   cross_green = sum_half(own galign - alt galign).   alt = opposite flank's locus across the junction.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import tempfile
import tifffile
os.environ.setdefault('SVLA_TIFDIR', config.SVLA_TIFDIR); os.environ.setdefault('SVLA_BASE', config.SVLA_BASE)
import sv_la as S, sv_la_transloc as T, run_transloc_fixed as T2
import junction_search as JS
BPP = S.BPP; FIG4 = JS.FIG4; ca = T._cont_aligned_bp

def scores(seg, chrom, ab, di, sc):
    fd = os.path.join(tempfile.gettempdir(), f'cf_{os.getpid()}_{seg.shape[-1]}_{int(ab)%100000}.tif')
    tifffile.imwrite(fd, seg)
    P = S.make_anchor(fd, chrom, int(round(ab)), int(di), float(sc)); os.remove(fd)
    g = S.LA.global_ccr(P.mol_red, P.ref_red, int(round(ab)), float(sc), int(di))
    return (0.0 if g != g else float(g)), float(T2._green_overlap(P))

def confidence(stem, L, R, cont, jx_bp):
    arr = tifffile.imread(f'{FIG4}/{stem}.tif'); Wmol = arr.shape[-1]
    WL = JS.piece_width(stem, 'left'); WR = JS.piece_width(stem, 'right'); mstartR = Wmol - WR
    jx = int(round(jx_bp / BPP))
    Lc, Lab, Ldi, Lsc = L; Rc, Rab, Rdi, Rsc = R
    abfL = JS.ab_full(Lab, Ldi, Lsc, 0, WL, Wmol); abfR = JS.ab_full(Rab, Rdi, Rsc, mstartR, WR, Wmol)
    segL = arr[:, :, 0:jx]; segR = arr[:, :, jx:Wmol]
    oL, gL = scores(segL, Lc, ca(abfL, Ldi, Lsc, Wmol, 0, jx), Ldi, Lsc)          # own
    oR, gR = scores(segR, Rc, ca(abfR, Rdi, Rsc, Wmol, jx, Wmol - jx), Rdi, Rsc)
    aL, galaL = scores(segL, Rc, ca(abfR, Rdi, Rsc, Wmol, 0, jx), Rdi, Rsc)           # alt: L@R locus
    aR, galaR = scores(segR, Lc, ca(abfL, Ldi, Lsc, Wmol, jx, Wmol - jx), Ldi, Lsc)   # alt: R@L locus
    return dict(own_L=oL, own_R=oR, alt_L=aL, alt_R=aR,
                cross_red=round((oL - aL) + (oR - aR), 3), cross_green=round((gL - galaL) + (gR - galaR), 3),
                gal_own=(round(gL, 2), round(gR, 2)),
                alt_gal=(round(galaL, 2), round(galaR, 2)))

if __name__ == '__main__':
    A = 'A026'; O = '0317_HG002_NbBssSI_TdT_5_merged_t058_t059_mol001_240kb'
    Acont = ('chr01_m', 28031800, 1, 0.98)
    CASES = [
        (A, 'A026 intra maternal (m.m)', ('chr01_m', 108699600, 1, 0.96), ('chr01_m', 108465600, 1, 1.0), Acont, 59600),
        (A, 'A026 intra paternal (p.p)', ('chr01_p', 108782400, 1, 0.96), ('chr01_p', 108548200, 1, 1.0), Acont, 59600),
        (A, 'A026 inter  (p·m)',    ('chr01_p', 108782400, 1, 0.96), ('chr01_m', 108465600, 1, 1.0), Acont, 59600),
        (A, 'A026 inter  (m·p)',    ('chr01_m', 108699600, 1, 0.96), ('chr01_p', 108548200, 1, 1.0), Acont, 59600),
        (O, '0317 (chr10_m·chr04_p)', ('chr10_m', 33697400, -1, 1.02), ('chr04_p', 185156800, -1, 0.99),
         ('chr04_m', 184390200, -1, 0.98), 100870),
    ]
    print(f"{'case':26} {'own_L':>6} {'own_R':>6} {'cross_red':>9} {'cross_green':>11}")
    for stem, name, L, R, cont, jx in CASES:
        d = confidence(stem, L, R, cont, jx)
        print(f"{name:26} {d['own_L']:6.3f} {d['own_R']:6.3f} {d['cross_red']:+9.3f} {d['cross_green']:+11.3f}")
        print(f"{'':26} galign own={d['gal_own']} alt={d['alt_gal']}")
