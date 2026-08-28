#!/usr/bin/env python3
"""End-to-end check: run the bundled example and compare against example/expected.json.

    python3 tests/test_example.py

Needs `maligner_dp` on PATH. Takes about a minute -- the example reference is two chromosomes,
not a whole genome. Unlike tests/test_release.py this exercises the scoring itself, so it
catches a change that runs cleanly but moves molecules to different loci.

The expected values come from running this pipeline against the bundled reference; all three
molecules land on the locus a genome-wide search assigns them.
"""
import json, os, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, 'example')


def main():
    """Run the bundled two-chromosome example and compare it with example/expected.json."""
    import shutil, subprocess, tempfile
    if not shutil.which('maligner_dp'):
        print('SKIP: maligner_dp not on PATH')
        return 0
    with tempfile.TemporaryDirectory() as out:
        r = subprocess.run([sys.executable, os.path.join(ROOT, 'run_pipeline.py'),
                            '--molecules', os.path.join(EX, 'molecules'),
                            '--ref-map', os.path.join(EX, 'reference'),
                            '--output', out, '--shards', '3'],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print('FAIL: pipeline exited', r.returncode)
            print(r.stdout[-2000:], r.stderr[-2000:])
            return 1
        # One checker for both examples, so the two-chromosome and the genome-wide runs are held
        # to the same criteria: chromosome, locus, F1, and the haplotype ranked second. h_score is
        # not among them -- it is a null-model p-value whose scale depends on how much genome the
        # null could draw from, so the same correct alignment scores about -1.4 here and about -9
        # against the full reference.
        chk = subprocess.run([sys.executable, os.path.join(ROOT, 'example_genome', 'check.py'),
                              out, os.path.join(EX, 'expected.json')],
                             capture_output=True, text=True)
        print(chk.stdout.rstrip())
        return chk.returncode


def test_example():
    """pytest entry point, so `pytest tests/` covers the example and not only the structural checks.

    Without this, nothing in this file is named test_*, so pytest collects it and runs nothing --
    the suite passes while the part that actually exercises the scoring never executes. Skipping
    rather than passing when maligner_dp is absent keeps that distinction visible: half the
    candidate pool is missing, so a green result would not mean what it appears to mean.

    pytest is imported inside the function because this file is also run as a plain script
    (`python tests/test_example.py`), including by CI, on machines that have no pytest.
    """
    import shutil
    import pytest
    if not shutil.which('maligner_dp'):
        pytest.skip('maligner_dp not on PATH: the label-only scorer cannot run')
    assert main() == 0


if __name__ == '__main__':
    sys.exit(main())
