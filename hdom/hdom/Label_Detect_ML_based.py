#!/usr/bin/env python3
"""
Label_Detect_ML_based.py — adaptive (per-image) green nick-label detection for hDOM molecules.

Replaces the fixed-threshold detector (hDOM_mal_Y.refined_green_peaks, which cuts at
green_mean + 0.10*(green_max-green_mean)) with a RandomForest that decides label/no-label per
candidate from image-relative features, so the effective threshold ADAPTS to each molecule.

Trained on user-curated labels (edited FP/FN, plus molecules confirmed as needing no correction
confirmations). Group (leave-molecule-out) CV: F1 ~0.98, precision ~0.98 vs fixed-threshold ~0.964.

Features (12), mostly image-relative so the model generalizes across molecules:
  h_norm         (green_height - mean) / (max - mean)     # where the fixed 0.10 cut lives
  h_over_med     green_height / median-of-signal
  prom_frac      peak prominence / green_max
  fwhm           peak width (px)
  contrast_abs   green_height - local background
  contrast_ratio green_height / local background
  dleft, dright  distance to nearest neighbour candidate
  density        # candidates within +-30 px
  pos_norm       x / width
  voff           chromatic vertical offset: green-row-center - RED-BACKBONE-row (backbone traced,
                 not per-column, so it is robust to noise/gaps)  ->  real labels sit at a constant
                 offset from the backbone; noise deviates
  voff_dev       |voff - molecule-median voff|   (the discriminative form)

USAGE
  # as a library (drop-in replacement for refined_green_peaks):
  from Label_Detect_ML_based import detect_labels
  peaks = detect_labels('mol_orig.tif')          # -> sorted list of x (px)

  # CLI:
  python Label_Detect_ML_based.py mol_orig.tif            # print detected peaks + vs fixed
  python Label_Detect_ML_based.py --retrain               # rebuild train_data + retrain from Label_ML
"""
import os, sys, json, glob, warnings, argparse
import numpy as np, tifffile, joblib
from scipy import signal
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(HERE, 'label_model.joblib')
TRAIN_CSV = os.path.join(HERE, 'label_train_data.csv')
LABEL_ML = os.environ.get('HDOM_LABEL_ML', '')     # curated label dir (--retrain only)
OPT = os.environ.get('HDOM_LABEL_SRC', '')         # source *_orig.tif  (--retrain only)
SKIP = {'0128_Re2_22_23_2_samevec'}   # too-clumped molecule excluded from training
TOL = 3.0                              # candidate<->label match tolerance (px)
DEFAULT_THRESH = 0.5
FEATS = ['h_norm', 'h_over_med', 'prom_frac', 'fwhm', 'contrast_abs', 'contrast_ratio',
         'dleft', 'dright', 'density', 'pos_norm', 'voff', 'voff_dev']

# ---------------------------------------------------------------- feature extraction (shared) ----
def _col_com(col):
    """background-subtracted center-of-mass (row index) of a 1-D column, or nan."""
    c = col - np.median(col); c[c < 0] = 0; s = c.sum()
    return float((np.arange(len(c)) * c).sum() / s) if s > 1e-6 else np.nan

def _backbone_trace(red, win=4):
    """Continuity-tracked red-DNA backbone y(x) (same routine used for straightening): windowed COM
    that follows the previous row, robust to per-column noise/gaps. Returns y for every x."""
    H, Wd = red.shape; r = np.clip(red - np.median(red), 0, None); cp = r.max(0)
    pos = cp[cp > 0]; la = np.median(pos[pos > np.percentile(pos, 40)]) if len(pos) else cp.max()
    thr = max(la * 0.15, r.max() * 0.02); x0 = int(np.argmax(cp))
    def cw(x, yp):
        """Intensity centroid within +-win rows of the previous column's y, or None if nothing is there."""
        lo, hi = max(0, int(round(yp)) - win), min(H, int(round(yp)) + win + 1); w = r[lo:hi, x]
        if w.sum() <= 0 or w.max() < thr: return None
        return lo + (np.arange(hi - lo) * w).sum() / w.sum()
    y = np.full(Wd, np.nan); y0 = cw(x0, int(np.argmax(r[:, x0]))) or float(np.argmax(r[:, x0]))
    for order in (range(x0, Wd), range(x0, -1, -1)):
        yp = y0
        for x in order:
            c = cw(x, yp)
            if c is not None: y[x] = c; yp = c
    v = np.isfinite(y)
    if v.any(): y = np.interp(np.arange(Wd), np.where(v)[0], y[v])
    return y

def candidate_features(orig_tif):
    """Every plausible green peak in the molecule, with the 12 features the classifier reads.

    The candidate cut is deliberately far below the fixed detector's: peaks are taken at 2% of the
    molecule's green range rather than 10%. That is not the decision -- it is the pool the decision
    is made over, and it has to be permissive enough to contain the faint real labels the fixed
    threshold was missing. The classifier then removes the noise the loose cut let in.

    Two of the features need the red frame. `voff` is the vertical distance from the green peak's
    centroid to the traced backbone, which is dominated by chromatic aberration and is therefore
    nearly CONSTANT for real labels on a given molecule; `voff_dev`, the deviation from the
    molecule's own median voff, is the form that discriminates, since the constant itself varies
    with the optics of the acquisition.

    Returns (xs, feats, medv): integer x positions, one feature dict each, and the molecule's median
    voff. Features that could not be computed come back NaN and are imputed downstream.
    """
    d = tifffile.imread(orig_tif); red2d = d[0].astype(float); green2d = d[1].astype(float)
    g = np.max(green2d, axis=0); W = len(g); mean_g = g.mean(); max_g = g.max(); rng = max(max_g - mean_g, 1e-6)
    med_pos = np.median(g[g > mean_g]) if (g > mean_g).any() else mean_g
    red_by = _backbone_trace(red2d)
    cand, props = signal.find_peaks(g, height=mean_g + 0.02 * rng, prominence=0.02 * max_g, width=1, distance=3)
    cand = list(cand); proms = list(props['prominences']); wids = list(props['widths'])
    o = np.argsort(cand); cand = list(np.array(cand)[o]); proms = list(np.array(proms)[o]); wids = list(np.array(wids)[o])
    carr = np.array(cand, float); xs = []; feats = []; voffs = []
    for i, x in enumerate(cand):
        xi = int(x); h = g[xi]
        lo, hi = max(0, xi - 20), min(W, xi + 21)
        win = np.concatenate([g[lo:max(lo, xi - 3)], g[min(hi, xi + 4):hi]]); bg = float(np.median(win)) if win.size else mean_g
        dl = (x - carr[i - 1]) if i > 0 else 999.0; dr = (carr[i + 1] - x) if i < len(cand) - 1 else 999.0
        dens = int(np.sum(np.abs(carr - x) <= 30)) - 1
        gy = _col_com(green2d[:, max(0, xi - 1):xi + 2].sum(axis=1))
        ry = np.nanmedian(red_by[max(0, xi - 3):xi + 4]); voff = (gy - ry) if (np.isfinite(gy) and np.isfinite(ry)) else np.nan
        voffs.append(voff)
        xs.append(xi)
        feats.append({'h_norm': (h - mean_g) / rng, 'h_over_med': h / (med_pos + 1e-6), 'prom_frac': proms[i] / max_g,
                      'fwhm': wids[i], 'contrast_abs': h - bg, 'contrast_ratio': h / (bg + 1e-6),
                      'dleft': min(dl, 999), 'dright': min(dr, 999), 'density': dens, 'pos_norm': xi / W, 'voff': voff})
    medv = np.nanmedian(voffs) if np.isfinite(voffs).any() else np.nan
    for f in feats: f['voff_dev'] = abs(f['voff'] - medv) if (np.isfinite(f['voff']) and np.isfinite(medv)) else np.nan
    return xs, feats, medv

def _feature_matrix(feats):
    """Feature dicts -> array in FEATS order, NaNs filled with that feature's median.

    A NaN reaches here whenever a feature was undefined for one candidate -- no backbone under it,
    say -- and scikit-learn will not accept one. Imputing with the column median keeps the candidate
    in the pool, which matters because dropping it would silently lose a label.
    """
    X = np.array([[(f.get(k, np.nan) if np.isfinite(f.get(k, np.nan)) else np.nan) for k in FEATS] for f in feats], float)
    for j in range(X.shape[1]):
        col = X[:, j]; m = np.nanmedian(col); X[np.isnan(col), j] = (m if np.isfinite(m) else 0.0)
    return X

# ---------------------------------------------------------------- detection (main use) ----
_MODEL = None
def _model():
    """The trained RandomForest, loaded once and cached.

    A genome-wide run calls detect_labels for every molecule; reloading the joblib each time would
    dominate the stage.
    """
    global _MODEL
    if _MODEL is None: _MODEL = joblib.load(MODEL_PATH)
    return _MODEL

def detect_labels(orig_tif, thresh=DEFAULT_THRESH, return_prob=False):
    """Detect this molecule's green labels. Returns sorted x positions in pixels.

    Drop-in replacement for the fixed-threshold hDOM_mal_Y.refined_green_peaks. The difference that
    matters is that the cut is per-candidate and image-relative, so a dim molecule is not stripped
    of its labels and a bright one is not filled with noise, where one global fraction of the peak
    height cannot do both.

    return_prob gives (x, probability) for EVERY candidate rather than the thresholded set -- used
    when reviewing a molecule by hand, since it shows what sat just under the cut.
    """
    M = _model(); clf = M['model']
    xs, feats, _ = candidate_features(orig_tif)
    if not xs: return [] if not return_prob else []
    prob = clf.predict_proba(_feature_matrix(feats))[:, 1]
    if return_prob: return sorted(zip(xs, [float(p) for p in prob]))
    return sorted(int(xs[i]) for i in range(len(xs)) if prob[i] >= thresh)

# ---------------------------------------------------------------- retraining (when labels grow) ----
def build_training_data(out_csv=TRAIN_CSV):
    """Rebuild the training table from the user-curated Label_ML dir: molecules whose labels were
    edited, plus those confirmed as needing no correction."""
    import csv
    def is_labeled(b):
        """True if a curator actually ruled on this molecule -- edited a peak, or marked it reviewed.

        A molecule confirmed as needing no correction is as much a training example as an edited
        one; a molecule nobody looked at is not, and must not enter as if its labels were checked.
        """
        ef = f'{LABEL_ML}/{b}_labeledit.json'
        if not os.path.exists(ef): return False
        j = json.load(open(ef, encoding='utf-8'));  return isinstance(j, dict) and (bool(j.get('reviewed')) or
            any(p.get('status') in ('rm', 'add') for p in j.get('peaks', [])))
    def finals(b):
        """The curator's final label positions for this molecule."""
        return sorted(float(x) for x in json.load(open(f'{LABEL_ML}/{b}_labeledit.json', encoding='utf-8')).get('final', []))
    rows = []
    for mp in sorted(glob.glob(f'{LABEL_ML}/*_meta.json')):
        b = os.path.basename(mp)[:-len('_meta.json')]
        if b in SKIP or not is_labeled(b): continue
        orig = f'{OPT}/{b}_orig.tif'
        if not os.path.exists(orig): continue
        xs, feats, _ = candidate_features(orig)
        fin = np.array(finals(b), float)
        # ensure every final label is represented (add missed FN as forced candidates)
        have = set()
        for xi, f in zip(xs, feats):
            lab = 1 if (fin.size and np.min(np.abs(fin - xi)) <= TOL) else 0
            rows.append([b, xi, lab] + [f[k] if f[k] == f[k] else '' for k in FEATS])
            have.add(xi)
        for fx in fin:
            if not any(abs(int(fx) - xi) <= TOL for xi in xs):
                xs2, feats2, _ = candidate_features(orig)   # leftover no-op: result discarded
                # A curated label with no candidate under it, even at the 2% cut. It is kept as a
                # positive so the target count stays right, but its features are left blank and get
                # median-imputed in train(), so it carries no usable signal -- it is rare by
                # construction (the loose cut is chosen to make it so).
                rows.append([b, int(fx), 1] + [''] * len(FEATS))
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh); w.writerow(['mol', 'x', 'label'] + FEATS); w.writerows(rows)
    pos = sum(1 for r in rows if r[2] == 1)
    print(f'train_data: {len(set(r[0] for r in rows))} mols, {len(rows)} candidates (pos {pos}) -> {out_csv}')
    return out_csv

def train(train_csv=TRAIN_CSV, out_model=MODEL_PATH):
    """Retrain the RandomForest and report leave-molecule-out CV beside the fixed threshold.

    The CV is grouped BY MOLECULE, not by candidate. Candidates from one molecule share its
    brightness, its background and its chromatic offset, so a random split would let the model see
    the answer for a molecule it is then tested on and report a score it will not reproduce on new
    data. The fixed-threshold row is scored on the same folds so the comparison is like for like.
    """
    import csv
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import GroupKFold, cross_val_predict
    from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score
    rows = list(csv.reader(open(train_csv, encoding='utf-8'))); hd = rows[0]; body = rows[1:]
    fi = [hd.index(c) for c in FEATS]
    def _f(v):
        """CSV cell -> float; a blank (a forced positive with no features) becomes NaN."""
        try: return float(v)
        except Exception: return np.nan
    X = np.array([[_f(r[i]) for i in fi] for r in body]); y = np.array([int(r[hd.index('label')]) for r in body])
    groups = np.array([r[hd.index('mol')] for r in body])
    for j in range(X.shape[1]):
        col = X[:, j]; m = np.nanmedian(col); X[np.isnan(col), j] = (m if np.isfinite(m) else 0.0)
    hn = X[:, FEATS.index('h_norm')]
    def rep(n, p):
        """Print accuracy/precision/recall/F1 for one set of predictions."""
        print(f'  {n:24} acc {accuracy_score(y,p):.3f}  P {precision_score(y,p):.3f}  R {recall_score(y,p):.3f}  F1 {f1_score(y,p):.3f}')
    print(f'samples {len(y)} pos {y.sum()} neg {len(y)-y.sum()} mols {len(set(groups))}')
    rep('fixed threshold (0.1)', (hn >= 0.1).astype(int))
    clf = RandomForestClassifier(n_estimators=300, min_samples_leaf=3, random_state=0)
    rep('RandomForest (groupCV)', cross_val_predict(clf, X, y, cv=GroupKFold(5), groups=groups))
    clf = RandomForestClassifier(n_estimators=400, min_samples_leaf=3, random_state=0).fit(X, y)
    joblib.dump({'model': clf, 'feats': FEATS}, out_model)
    print('feature importance:', {f: round(v, 3) for f, v in sorted(zip(FEATS, clf.feature_importances_), key=lambda t: -t[1])})
    print(f'model -> {out_model}')

# ---------------------------------------------------------------- CLI ----
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Adaptive ML green-label detection (hDOM).')
    ap.add_argument('tif', nargs='?', help='molecule *_orig.tif to detect labels on')
    ap.add_argument('--thresh', type=float, default=DEFAULT_THRESH)
    ap.add_argument('--retrain', action='store_true', help='rebuild training data from Label_ML and retrain')
    a = ap.parse_args()
    if a.retrain:
        build_training_data(); train()
    elif a.tif:
        sys.path.insert(0, HERE)
        try:
            from hDOM_mal_Y import refined_green_peaks
            fixed = sorted(int(x) for x in refined_green_peaks(a.tif))
        except Exception:
            fixed = []
        ml = detect_labels(a.tif, thresh=a.thresh)
        near = lambda x, s: any(abs(x - p) <= 3 for p in s)
        new = [x for x in ml if not near(x, fixed)]; drop = [x for x in fixed if not near(x, ml)]
        print(f'{os.path.basename(a.tif)}')
        print(f'  ML labels : {len(ml)}   fixed-threshold : {len(fixed)}')
        print(f'  ML new (FN recovered): {len(new)}  {new}')
        print(f'  ML dropped (FP)      : {len(drop)}  {drop}')
        print(f'  peaks: {ml}')
    else:
        ap.print_help()
