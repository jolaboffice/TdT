"""Deploy the Figure3-style viewer (2-section local+mal table, length/lin. panel labels,
rank1/2 + mal1/2 panels) into a pipeline output folder, and build the _lengths.json the
viewer needs. Used by local_align_pipeline.py (standard) and straighten_align.py (linearized).

The viewer template `figure3_index.php` is the generic form of zave/Figure3/index.php
(length relabel gated on _lengths.json existence, not on folder name).
mdom/pcc/homolog rows stay in the CSV (kept for later use); the viewer just shows local+mal.
"""
import os, json, shutil, fcntl

CODE = os.path.dirname(os.path.abspath(__file__))
VIEWER_SRC = os.path.join(CODE, 'figure3_index.php')


def deploy_viewer(folder):
    """Copy the Figure3-style viewer into folder/index.php."""
    try:
        shutil.copy2(VIEWER_SRC, os.path.join(folder, 'index.php'))
        return True
    except Exception:
        return False


def mol_len_kb(tif, bpp=200):
    """Molecule length in kb = round(tif_width_px * bpp / 1000).  (matches Figure3 _lengths.json)"""
    try:
        import tifffile
        a = tifffile.imread(tif)
        w = a.shape[-1]
        return int(round(w * bpp / 1000.0))
    except Exception:
        return None


# 19-col Figure3 CSV = 14 rerank cols + 5 maligner-native cols (mal rows only)
FIG3_COLS = ['rank', 'ref_name', 'method', 'aligned_bp', 'dir', 'scale', 'cc_r', 'local_cc_r',
             'cc_g', 'cc_rg', 'cc_rg2', 'log(P)', 'misc', 'F1',
             'm_score', 'mal_matched', 'mal_FP', 'mal_FN', 'mal_sizing']


def parse_mal_native(dp_txt, topn=10):
    """Parse a saved maligner_dp .txt (produced by calc_Q) -> top-N alignments by m_score
    with native fields. Returns list of dicts, [] if empty, None if missing/unparseable
    (query_misses=FP, ref_misses=FN). No maligner re-run."""
    if not dp_txt or not os.path.exists(dp_txt):
        return None
    lines = [l for l in open(dp_txt, encoding='utf-8').read().strip().split('\n') if l]
    if len(lines) < 2:
        return []
    h = {c: i for i, c in enumerate(lines[0].split('\t'))}
    need = ('ref_map', 'ref_start_bp', 'is_forward', 'query_scaling_factor', 'm_score',
            'num_matched_chunks', 'query_misses', 'ref_misses', 'sizing_score_rescaled')
    if not all(k in h for k in need):
        return None
    out = []
    for l in lines[1:]:
        f = l.split('\t')
        try:
            out.append(dict(
                ref=f[h['ref_map']].replace('.tif', ''),
                ab=int(round(float(f[h['ref_start_bp']]))),
                dir=(1 if f[h['is_forward']].strip().upper().startswith('F') else -1),
                scale=round(float(f[h['query_scaling_factor']]), 2),
                m=round(float(f[h['m_score']]), 2),
                matched=f[h['num_matched_chunks']].strip(),
                FP=f[h['query_misses']].strip(),
                FN=f[h['ref_misses']].strip(),
                sizing=round(float(f[h['sizing_score_rescaled']]), 2)))
        except Exception:
            continue
    out.sort(key=lambda x: x['m'])   # lower m_score = better maligner alignment
    return out[:topn]


def to_figure3_csv(csv_path, dp_txt, topn=10):
    """Rewrite a 14-col rerank CSV into the 19-col Figure3 format IN PLACE:
    non-mal rows (local/mDOM/pcc/homolog) are kept & padded with 5 blank cols;
    mal rows are replaced with maligner-native fields (m_score, matched, FP, FN, sizing)
    parsed from dp_txt, with cc_*/log(P) blank. mdom/pcc/homolog data is preserved."""
    import csv as _csv
    import io as _io
    try:
        rows = list(_csv.DictReader(open(csv_path, encoding='utf-8')))
    except Exception:
        return 'nocsv'
    if not rows:
        return 'nocsv'
    native = parse_mal_native(dp_txt, topn)
    out = []
    for r in rows:                       # keep every non-mal row, padded to 19 cols
        if (r.get('method') or '').strip() == 'mal':
            continue
        row = {c: r.get(c, '') for c in FIG3_COLS[:14]}
        for c in FIG3_COLS[14:]:
            row[c] = ''
        out.append(row)
    if native is None:                   # dp_txt missing/bad: keep original mal rows (padded)
        for r in rows:
            if (r.get('method') or '').strip() != 'mal':
                continue
            row = {c: r.get(c, '') for c in FIG3_COLS[:14]}
            for c in FIG3_COLS[14:]:
                row[c] = ''
            out.append(row)
    else:                                # native mal rows (log(P) blank, m_score+matched/FP/FN/sizing)
        for i, d in enumerate(native, 1):
            row = {c: '' for c in FIG3_COLS}
            row.update({'rank': str(i), 'ref_name': d['ref'], 'method': 'mal',
                        'aligned_bp': str(d['ab']), 'dir': str(d['dir']), 'scale': f"{d['scale']:.2f}",
                        'm_score': f"{d['m']:.2f}", 'mal_matched': d['matched'], 'mal_FP': d['FP'],
                        'mal_FN': d['FN'], 'mal_sizing': f"{d['sizing']:.2f}"})
            out.append(row)
    buf = _io.StringIO()
    w = _csv.DictWriter(buf, fieldnames=FIG3_COLS)
    w.writeheader()
    w.writerows(out)
    tmp = csv_path + '.tmp'
    open(tmp, 'w', encoding='utf-8').write(buf.getvalue())
    os.replace(tmp, csv_path)
    return 'done'


# _summary.csv drives the viewer's molecule LIST + nav order + score column
SUMMARY_COLS = ['molecule', 'score', 'local1', 'mDOM1', 'mal1', 'pcc1', 'H', 'mdom_class']


def summary_row(base, csv_path):
    """Build a Figure3 _summary.csv row [molecule,score,local1,mDOM1,mal1,pcc1,H,mdom_class]
    from a molecule's matching_results.csv. score = local rank1 log(P). Returns None on failure."""
    import csv as _csv
    try:
        rows = list(_csv.DictReader(open(csv_path, encoding='utf-8')))
    except Exception:
        return None
    if not rows:
        return None

    def first(m):
        for r in rows:
            if (r.get('method') or '').strip().lower() == m:
                return r
        return None

    l1 = first('local')
    score = (l1.get('log(P)', '') if l1 else '')

    def ref1(m):
        r = first(m)
        return (r.get('ref_name', '') if r else '-')

    hasH = any(('homolog' in (r.get('method') or '').lower())
               or (r.get('method') or '').lower().startswith('h_')
               or (str(r.get('rank', '')).strip().upper() == 'H') for r in rows)
    return [base, score, ref1('local'), ref1('mdom'), ref1('mal'), ref1('pcc'),
            'Y' if hasH else '-', '']


def add_summary(folder, mapping):
    """Merge {molecule: row_list} into folder/_summary.csv under an exclusive lock (shard-safe)."""
    import csv as _csv
    import io as _io
    mapping = {k: v for k, v in (mapping or {}).items() if v}
    if not mapping:
        return
    p = os.path.join(folder, '_summary.csv')
    fd = os.open(p, os.O_RDWR | os.O_CREAT, 0o664)
    try:
        with os.fdopen(fd, 'r+') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            txt = f.read()
            d = {}
            if txt.strip():
                for row in list(_csv.reader(_io.StringIO(txt)))[1:]:   # skip header
                    if row:
                        d[row[0]] = row
            d.update(mapping)
            buf = _io.StringIO()
            w = _csv.writer(buf)
            w.writerow(SUMMARY_COLS)
            w.writerows(d.values())
            f.seek(0)
            f.truncate()
            f.write(buf.getvalue())
            fcntl.flock(f, fcntl.LOCK_UN)
    except Exception:
        pass


def add_lengths(folder, mapping):
    """Merge {base: length_kb} into folder/_lengths.json under an exclusive lock (shard-safe)."""
    mapping = {k: v for k, v in (mapping or {}).items() if v is not None}
    if not mapping:
        return
    p = os.path.join(folder, '_lengths.json')
    fd = os.open(p, os.O_RDWR | os.O_CREAT, 0o664)
    try:
        with os.fdopen(fd, 'r+') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            txt = f.read()
            try:
                d = json.loads(txt) if txt.strip() else {}
            except Exception:
                d = {}
            d.update(mapping)
            f.seek(0)
            f.truncate()
            f.write(json.dumps(d, ensure_ascii=False, indent=1))
            fcntl.flock(f, fcntl.LOCK_UN)
    except Exception:
        pass
