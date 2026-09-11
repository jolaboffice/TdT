#!/usr/bin/env python3
"""Structural checks on the package. No reference genome, no molecules, seconds to run.

    python3 tests/test_release.py          # or: pytest tests/

These cover the ways an incomplete copy of the package still runs to completion: a candidate
generator whose script is missing is skipped rather than reported, and a module carrying an
absolute path works only on the machine it came from.

Whether a known molecule lands on its known locus needs the reference and the example
molecules; that is tests/test_example.py.
"""
import os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, 'hdom')
CAND = os.path.join(PKG, 'candidates')


def test_layout():
    """Every file the entry point checks for at startup is present."""
    # LICENSE is not in this list: the licence sits at the root of the toolkit, one level up,
    # because maligner is bundled there under the same terms and a second copy here would only
    # invite the two to drift apart.
    for rel in ('hDOM_pipeline.py', 'make_reference.py', 'README.md',
                'hdom/straighten_align.py', 'hdom/local_align_pipeline.py',
                'hdom/mdom_reanchor.py', 'hdom/hDOM_view.py',
                'hdom/Label_Detect_ML_based.py', 'hdom/label_model.joblib',
                'hdom/figure3_view.py', 'hdom/qc.py',
                # One viewer, not two. figure3_view.py copies this file into every output
                # folder, so a second copy here would ship runs a viewer older than the repo.
                'hdom/index.php'):
        assert os.path.exists(os.path.join(ROOT, rel)), f'missing {rel}'


def test_subprocess_scripts_present():
    """hDOM_calc_Q invokes these by filename, and skips a step whose script is absent.

    The names are read out of the source rather than hard-coded, so a candidate generator
    added to hDOM_calc_Q is checked for without editing this test.
    """
    src = open(os.path.join(CAND, 'hDOM_calc_Q.py'), encoding='utf-8').read()
    names = set(re.findall(r'["\']([A-Za-z_0-9]+\.py)["\']', src))
    names = {n for n in names if n.startswith(('hDOM_', 'DOM_'))}
    assert names, 'expected hDOM_calc_Q to name its subprocess scripts'
    for n in names:
        assert os.path.exists(os.path.join(CAND, n)), \
            f'hDOM_calc_Q.py invokes {n} but hdom/candidates/{n} is missing'


def test_no_absolute_paths():
    """A module carrying an absolute path works only on the machine it was written on."""
    bad = []
    for base, _, files in os.walk(ROOT):
        if '__pycache__' in base or os.sep + 'tests' in base:
            continue
        for f in files:
            if not f.endswith('.py'):
                continue
            p = os.path.join(base, f)
            for i, line in enumerate(open(p, encoding='utf-8', errors='ignore'), 1):
                if line.lstrip().startswith('#'):
                    continue
                if re.search(r'''["'](/home/|/var/www/|/tmp/)''', line):
                    bad.append(f'{os.path.relpath(p, ROOT)}:{i}')
    assert not bad, 'absolute paths: ' + ', '.join(bad)


def test_modules_import():
    """Import the pipeline modules the way straighten_align does, from the package dir."""
    code = ('import sys; sys.path.insert(0, %r);'
            'import local_align_pipeline as P, mdom_reanchor, figure3_view, polyn_filter;'
            'assert P.BPP > 0 and P.TRIM_PX >= 0' % PKG)
    r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_bpp_is_configurable():
    """--bpp has to reach the scoring constants; TRIM_PX is 5 kb expressed in pixels."""
    code = ('import sys; sys.path.insert(0, %r);'
            'import local_align_pipeline as P;'
            'print(P.BPP, P.TRIM_PX)' % PKG)
    for bpp, trim in ((200, 25), (250, 20)):
        env = dict(os.environ, HDOM_BPP=str(bpp))
        out = subprocess.run([sys.executable, '-c', code], capture_output=True,
                             text=True, env=env).stdout.split()
        assert out == [str(bpp), str(trim)], f'HDOM_BPP={bpp} gave {out}'


def test_terminal_trim_rule():
    """Terminal-5 kb trim: drop labels within TRIM_PX of either end, but never all of them."""
    trim, w2 = 25, 600
    for peaks, keep in (([2, 100, 300, 598], [100, 300]),      # both ends dropped
                        ([100, 300], [100, 300]),              # nothing to drop
                        ([2, 598], [2, 598])):                 # all terminal -> keep as-is
        got = [p for p in peaks if trim <= p <= w2 - 1 - trim]
        if len(got) < 2:
            got = list(peaks)
        assert got == keep, f'{peaks} -> {got}, expected {keep}'


def test_declared_dependencies_cover_imports():
    """Every third-party module imported at module level must be in requirements.txt.

    DOM_lib imports numba and matplotlib at the top of the file, so the package cannot even be
    imported without them however little of their surface it uses. On a machine that happens to
    have them installed the omission is invisible.
    """
    import re
    declared = {re.split(r'[<>=!\[]', l.strip())[0].lower()
                for l in open(os.path.join(ROOT, 'requirements.txt'), encoding='utf-8')
                if l.strip() and not l.startswith('#')}
    declared |= {'os', 'sys', 'json', 'csv', 'glob', 'shutil', 'subprocess', 'argparse', 're',
                 'warnings', 'math', 'time', 'io', 'base64', 'collections', 'itertools', 'random',
                 'datetime', 'hashlib', 'threading', 'concurrent', 'fcntl', 'functools', 'pickle',
                 'tempfile', 'traceback', 'copy', 'string', 'textwrap', 'statistics', 'bisect',
                 'pathlib', 'html', 'urllib', 'xml', 'typing', 'importlib', 'platform', 'ast',
                 'PIL', 'pillow', 'sklearn', 'scikit-learn', 'cv2', 'skimage'}
    alias = {'pil': 'pillow', 'sklearn': 'scikit-learn'}
    local = set()
    for base, _, files in os.walk(PKG):
        local |= {f[:-3] for f in files if f.endswith('.py')}
    missing = set()
    for base, _, files in os.walk(PKG):
        if '__pycache__' in base:
            continue
        for f in files:
            if not f.endswith('.py'):
                continue
            import ast as _ast
            for node in _ast.walk(_ast.parse(open(os.path.join(base, f), encoding='utf-8').read())):
                mods = ([a.name for a in node.names] if isinstance(node, _ast.Import)
                        else [node.module] if isinstance(node, _ast.ImportFrom) and node.module
                        and node.level == 0 else [])
                for m in mods:
                    top = m.split('.')[0]
                    if top in local or top in declared:
                        continue
                    if alias.get(top.lower(), top.lower()) in declared:
                        continue
                    missing.add(top)
    assert not missing, f'imported but not in requirements.txt: {sorted(missing)}'


def test_make_reference_roundtrip():
    """Build a reference from a synthetic FASTA and check genome.maps is well formed."""
    motif, bpp = 'CTCGTG', 200
    seq = ('ATATATATAT' * 20 + motif) * 300                    # a site every ~206 bp
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'fa'); os.makedirs(fa)
        open(os.path.join(fa, 'chrT.fasta'), 'w', encoding='utf-8').write('>chrT\n' + seq + '\n')
        out = os.path.join(d, 'ref')
        r = subprocess.run([sys.executable, os.path.join(ROOT, 'make_reference.py'),
                            '--fasta', fa, '--out', out, '--bpp', str(bpp), '--motif', motif],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, 'chrT.tif'))
        line = open(os.path.join(out, 'genome.maps'), encoding='utf-8').read().strip().split('\n')[0]
        f = line.split('\t')
        assert f[0] == 'chrT', 'the name field is the bare reference name, no extension'
        assert int(f[1]) == int(-(-len(seq) // bpp)) * bpp, 'length_bp must be n_px * bpp'
        assert int(f[2]) == len(f) - 3, 'declared fragment count must match the fields present'
        assert sum(int(x) for x in f[3:]) == int(f[1]), 'fragments must sum to the total length'


if __name__ == '__main__':
    fns = [(n, f) for n, f in sorted(globals().items()) if n.startswith('test_')]
    bad = 0
    for n, f in fns:
        try:
            f(); print(f'  PASS  {n}')
        except AssertionError as e:
            bad += 1; print(f'  FAIL  {n}: {e}')
        except Exception as e:
            bad += 1; print(f'  ERROR {n}: {type(e).__name__}: {e}')
    print(f'\n{len(fns) - bad}/{len(fns)} passed')
    sys.exit(1 if bad else 0)
