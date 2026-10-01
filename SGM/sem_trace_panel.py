"""Redraw the SEM panel of an SGM view as the raw image with the trace and QD marks lifted off it.

`_result.tif` paints the traced backbone and the QD circles directly onto the DNA, where they sit
on top of the signal they describe and are hard to read at page scale. This script keeps the photo
untouched and lifts the annotation clear of it: the trace is drawn as a continuous line a fixed
distance above (or below) the molecule, and each QD becomes a vertical arrow pointing back at it.

Two things are worth knowing about the inputs.

The red trace in `_result.tif` is broken wherever a QD circle was painted over it, so the mask
cannot be used directly; the trace is rebuilt by interpolating its y per column. And QD positions
are taken from `<QD>.maps`, not from `_points.json`: the latter stores every click, including the
guide points that only steer the least-cost path, so it holds more points than there are labels
(QD010: 76 clicks, 43 labels). `.maps` is what the Query strip is drawn from, so using it keeps
the two in the same frame.

The Query/rank strips are drawn on rank1's alignment axis, so a reverse alignment (dir = -1)
mirrors them against the photo. The photo stays as it is and the three strips are flipped
together, which keeps Query and the reference maps corresponding to each other.
"""
import os, re, io, json, csv, base64, shutil, argparse
import numpy as np, tifffile
from scipy import ndimage as ndi
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

ap = argparse.ArgumentParser()
ap.add_argument('mol')
ap.add_argument('--data', required=True, help='directory holding the <QD###>/ folders')
ap.add_argument('--out',  required=True, help='the view folder sem_view.py wrote')
ap.add_argument('--base', default='_adjusted', help='suffix of the background tif')
ap.add_argument('--off', type=int, default=125)
ap.add_argument('--trace-r', type=int, default=6)
ap.add_argument('--tick-len', type=int, default=85)
ap.add_argument('--tick-th', type=int, default=5)
ap.add_argument('--head', type=int, default=30)
ap.add_argument('--head-deg', type=float, default=26)
ap.add_argument('--gap', type=int, default=28)
ap.add_argument('--margin', type=int, default=250)
# The page scales each photo to the column width, so a molecule twice as wide shows an arrow of
# the same pixel size at half the apparent size. Arrow dimensions are therefore scaled with the
# crop width, which fixes how large they look on screen. The offset that lifts the trace off the
# DNA is deliberately NOT scaled: it stays a fixed distance in original pixels, so a molecule that
# is reduced more shows the trace closer to the DNA.
ap.add_argument('--ref-width', type=int, default=10990)
ap.add_argument('--ref-scale', type=float, default=1.5)
ap.add_argument('--side', choices=('auto','up','down'), default='auto')
a = ap.parse_args()
k = a.mol; D = f'{a.data}/{k}'; F = a.out
IC = r'<div class="image-container"[^>]*>'          # the container may already carry attributes

raw = tifffile.imread(f'{D}/{k}{a.base}.tif')
if raw.ndim == 3: raw = raw[..., 0]
res = tifffile.imread(f'{D}/{k}_result.tif')
R, G, B = res[..., 0].astype(int), res[..., 1].astype(int), res[..., 2].astype(int)
colm = (np.abs(R-G) > 25) | (np.abs(G-B) > 25) | (np.abs(R-B) > 25)
red = colm & (R > G+40) & (R > B+40) & (G < 150)
fr = np.array([float(x) for x in open(f'{D}/{k}.maps').read().split('\t')[3:]])
H, W = raw.shape
ys, xs = np.nonzero(red); o = np.argsort(xs); xs, ys = xs[o], ys[o]
_u = np.unique(xs)
x0 = max(0, int(_u.min()) - a.margin); x1 = min(W-1, int(_u.max()) + a.margin); CW = x1 - x0 + 1
f = a.ref_scale * CW / a.ref_width
# f scales the arrow and the trace thickness only; a.off stays in original pixels (see above).
OFF = a.off; TR = max(1, int(round(a.trace_r*f)))
TL = int(round(a.tick_len*f)); HD = int(round(a.head*f))
TH = max(1, int(round(a.tick_th*f)) | 1); GAP = int(round(a.gap*f))
ux = np.unique(xs); uy = np.array([ys[xs == x].mean() for x in ux])
_xx = np.arange(int(ux.min()), int(ux.max())+1); _yy = np.interp(_xx, ux, uy)
# Does the molecule double back? Everything below models the trace as one y per column, which a
# molecule that folds on itself is not: its columns hold two branches, their mean follows neither,
# and an arc length indexed by x cannot carry bp onto the image. Count the columns that hold more
# than one run of trace, and take the measured mask itself where there are many.
_nmulti = 0
for _x in ux:
    _c = np.nonzero(red[:, _x])[0]
    if len(_c) > 1 and int((np.diff(_c) > 1).sum()):
        _nmulti += 1
FOLDED = _nmulti > 0.02*len(ux)
_my = np.nonzero(red.any(1))[0]
_ylo, _yhi = (float(_my.min()), float(_my.max())) if FOLDED else (float(_yy.min()), float(_yy.max()))
# Which side the annotation goes. A trace that runs near the top of the frame leaves no room
# above, and padding the canvas makes the whole panel smaller on screen; drop it below instead
# when there is room there. The arrow still points back at the trace either way.
NEED = OFF + GAP + TL + TR
if a.side == 'auto':
    SIDE = -1 if (_ylo >= NEED or (H - _yhi) < NEED) else +1
else:
    SIDE = -1 if a.side == 'up' else +1
if SIDE < 0:
    PAD = max(0, int(np.ceil(NEED - _ylo)) + 10); PADB = 0
else:
    PAD = 0; PADB = max(0, int(np.ceil(NEED - (H - _yhi))) + 10)
H = H + PAD + PADB
_b = raw
if PAD:  _b = np.vstack([np.zeros((PAD, W), raw.dtype), _b])
if PADB: _b = np.vstack([_b, np.zeros((PADB, W), raw.dtype)])
out = np.dstack([_b, _b, _b]).copy()
def put(yy, xx, rgb):
    yy = yy + PAD
    ok = (yy >= 0) & (yy < H) & (xx >= 0) & (xx < W); out[yy[ok], xx[ok]] = rgb
_disk = np.hypot(*np.mgrid[-TR:TR+1, -TR:TR+1]) <= TR
if FOLDED:
    # Lift the measured mask instead of a reconstructed line, so both branches survive. The QD
    # circles painted over the trace leave short gaps in it; a horizontal closing bridges those
    # without joining anything that is really apart.
    _m = ndi.binary_dilation(ndi.binary_closing(red, np.ones((1, 41), bool)), _disk)
    _fy, _fx = np.nonzero(_m)
    put(_fy + SIDE*OFF, _fx, (139, 0, 0))
else:
    for oy in range(-TR, TR+1):
        for ox in range(-TR, TR+1):
            if oy*oy + ox*ox > TR*TR: continue
            put(np.round(_yy+oy).astype(int) + SIDE*OFF, _xx+ox, (139, 0, 0))
# QD positions. `_points.json` holds every click with a type, so the 'measure' ones are exactly the
# labels and the 'guide' ones -- which only steer the least-cost path -- drop out; they are stored
# in click order, i.e. along the molecule, which is the order a folded trace needs. Falling back to
# bp laid on the arc length indexed by x places a label wherever that index is not monotonic, which
# on a molecule that doubles back is most of it.
QXY = None
try:
    _j = json.load(open(f'{D}/{k}_points.json', encoding='utf-8'))
    _mp = [p for p, t in zip(_j['points'], _j['point_types']) if t == 'measure']
    if len(_mp) == len(fr) + 1: QXY = np.asarray(_mp, float)
except Exception as _e:
    print(f'  [{k}] points read failed ({type(_e).__name__}), falling back to the bp projection',
          flush=True)
if QXY is not None:
    QX, QY = QXY[:, 0], QXY[:, 1]
else:
    _ds = np.hypot(np.diff(_xx), np.diff(_yy))
    _s = np.concatenate([[0.0], np.cumsum(_ds)])
    _cum = np.concatenate([[0.0], np.cumsum(fr)]); _cum /= _cum[-1]
    QX = np.interp(_cum*_s[-1], _s, _xx)      # x of each QD
    QY = np.interp(QX, ux, uy)
d = np.array([0.0, float(SIDE)]); th = np.deg2rad(a.head_deg)
def stroke(pts):
    for oy in range(-(TH//2), TH//2+1):
        for ox in range(-(TH//2), TH//2+1):
            put(np.round(pts[:, 1]+oy).astype(int), np.round(pts[:, 0]+ox).astype(int), (255, 255, 0))
for cx, cy in zip(QX, QY):
    p0 = np.array([cx, float(cy) + SIDE*OFF]) + GAP*d
    stroke(p0[None, :] + np.outer(np.linspace(0, 1, TL*3)*TL, d))
    rp = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]) @ d
    rm = np.array([[np.cos(-th), -np.sin(-th)], [np.sin(-th), np.cos(-th)]]) @ d
    tri = np.array([p0, p0 + HD*rp, p0 + HD*rm])
    y0i, y1i = int(np.floor(tri[:, 1].min())), int(np.ceil(tri[:, 1].max()))
    x0i, x1i = int(np.floor(tri[:, 0].min())), int(np.ceil(tri[:, 0].max()))
    gy2, gx2 = np.mgrid[y0i:y1i+1, x0i:x1i+1]
    v0 = tri[1]-tri[0]; v1 = tri[2]-tri[0]; den = v0[0]*v1[1] - v1[0]*v0[1]
    if abs(den) > 1e-9:
        u = ((gx2-tri[0, 0])*v1[1] - v1[0]*(gy2-tri[0, 1]))/den
        w = (v0[0]*(gy2-tri[0, 1]) - (gx2-tri[0, 0])*v0[1])/den
        ins = (u >= 0) & (w >= 0) & (u+w <= 1)
        put(gy2[ins].astype(int), gx2[ins].astype(int), (255, 255, 0))
png = f'{F}/{k}_sem_trace.png'
Image.fromarray(out[:, x0:x1+1]).save(png, 'PNG', optimize=True); os.chmod(png, 0o644)

# ---- strips: on a reverse alignment flip all three together to match the photo ----
v = f'{F}/{k}_view_ytb.html'
if not os.path.exists(v + '.bak_sem'): shutil.copy2(v, v + '.bak_sem')
s = open(v + '.bak_sem', encoding='utf-8', errors='ignore').read()      # always rebuild from the original
r1 = next(iter(csv.DictReader(open(f'{F}/{k}_matching_results.csv', encoding='utf-8'))), None)
flip = bool(r1) and str(r1.get('dir', '1')).strip().startswith('-')
strips = {}
for mo in re.finditer(r'<div class="image-label"[^>]*>([^<]*)</div>\s*(' + IC + r')(.*?)</div>', s, re.S):
    lb = mo.group(1).strip()
    b = re.search(r'src="data:image/png;base64,([A-Za-z0-9+/=]+)"', mo.group(3))
    if not b or lb not in ('Query', 'rank1', 'rank2'): continue
    im = Image.open(io.BytesIO(base64.b64decode(b.group(1))))
    if flip: im = im.transpose(Image.Transpose.FLIP_LEFT_RIGHT)
    fp = f'{F}/{k}_{lb}_strip.png'; im.save(fp, 'PNG', optimize=True); os.chmod(fp, 0o644)
    strips[lb] = (os.path.basename(fp), im.size[0], np.array(im.convert('RGB')).astype(int))
qim = strips['Query'][2]
g = ((qim[:, :, 1] - qim[:, :, 0]) > 60).any(0); gx = np.nonzero(g)[0]
runs = []; st = gx[0]; pv = gx[0]
for vv in gx[1:]:
    if vv > pv+1: runs.append((st+pv)/2); st = vv
    pv = vv
runs.append((st+pv)/2)
Wq = strips['Query'][1]
a0 = (QX[0]-x0)/CW; a1 = (QX[-1]-x0)/CW
b0 = runs[0]/Wq; b1 = runs[-1]/Wq
kk = (a1-a0)/(b1-b0); dd = a0 - kk*b0
# Resample the strips into the SEM crop's coordinate system rather than stretching them with CSS
# width/margin: an over-wide image with a negative margin changes how the flex row sizes itself,
# which moves the strips relative to the photo.
WOUT = 4000
for lb, (fn, wq, arr) in list(strips.items()):
    src = Image.open(f'{F}/{fn}'); A = np.array(src.convert('RGB'))
    j = np.arange(WOUT)
    sx = ((j/WOUT) - dd)/kk*wq                      # screen fraction -> source strip x
    si = np.clip(np.round(sx), 0, wq-1).astype(int)
    ok = (sx >= 0) & (sx <= wq-1)
    outA = np.zeros((A.shape[0], WOUT, 3), np.uint8)
    outA[:, ok] = A[:, si[ok]]
    fp = f'{F}/{k}_{lb}_strip.png'
    Image.fromarray(outA).save(fp, 'PNG', optimize=True); os.chmod(fp, 0o644)
_v = int(os.path.getmtime(png))        # cache-bust: the file names do not change between runs
mi = re.search(r'<img id="data-img"[^>]*>', s)
s = s[:mi.start()] + f'<img id="data-img" src="{k}_sem_trace.png?v={_v}" style="max-width:100%;height:auto">' + s[mi.end():]
def fix(mo):
    lb = mo.group(1).strip()
    if lb not in strips: return mo.group(0)
    return (f'<div class="image-label">{mo.group(1)}</div>'
            f'<div class="image-container"><img src="{strips[lb][0]}?v={_v}" '
            f'style="width:100%;height:auto"></div>')
s = re.sub(r'<div class="image-label"[^>]*>([^<]*)</div>\s*' + IC + r'(.*?)</div>', fix, s, flags=re.S)
open(v, 'w', encoding='utf-8').write(s)
print(f'{k}  {CW}x{H} {os.path.getsize(png)/1e6:.1f}MB  QD {len(QX)}  query ticks {len(runs)}  '
      f'f{f:.2f} {"below" if SIDE>0 else "above"} pad{PAD+PADB} '
      f'dir{"-1 flipped" if flip else "+1"}  width {kk*100:.3f}% margin-left {dd*100:.3f}%')
