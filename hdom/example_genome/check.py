#!/usr/bin/env python3
"""Compare a run against expected.json.

    python3 example_genome/check.py out_genome                       # genome-wide run
    python3 example_genome/check.py out example/expected.json        # two-chromosome run

What must reproduce is the placement: the chromosome and haplotype, the locus within a couple of
kilobases, the label-agreement F1, and the haplotype that comes second. h_score is not checked --
it is a null-model p-value whose scale depends on how much genome the null could draw from, so the
same correct alignment scores about -1.4 against two chromosomes and about -9 against all 47.
Holding it to a fixed value would fail runs that are right.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1] if len(sys.argv) > 1 else 'out_genome'
expp = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'expected.json')
exp = json.load(open(expp, encoding='utf-8'))
tol = exp.get('_tolerance', {'local1_bp': 2000, 'local1_F1': 0.03})
bad = 0
mols = [k for k in exp if not k.startswith('_')]
for mol in sorted(mols):
    e = exp[mol]
    p = os.path.join(out, mol + '_matching_results.csv')
    if not os.path.exists(p):
        print(f'  FAIL  {mol[:46]:<46} no result'); bad += 1; continue
    loc = [r for r in csv.DictReader(open(p, encoding='utf-8')) if r['method'] == 'local']
    if not loc:
        print(f'  FAIL  {mol[:46]:<46} no local rows'); bad += 1; continue
    ref, bp = loc[0]['ref_name'], int(float(loc[0]['aligned_bp']))
    f1 = float(loc[0]['F1']); r2 = loc[1]['ref_name'] if len(loc) > 1 else ''
    errs = []
    if ref != e['local1_ref']:
        errs.append(f"rank1 {ref} != {e['local1_ref']}")
    elif abs(bp - e['local1_bp']) > tol['local1_bp']:
        errs.append(f"rank1 at {bp:,} != {e['local1_bp']:,}")
    if abs(f1 - e['local1_F1']) > tol['local1_F1']:
        errs.append(f"F1 {f1} != {e['local1_F1']}")
    if r2 != e['local2_ref']:
        errs.append(f"rank2 {r2} != {e['local2_ref']}")
    if errs:
        bad += 1; print(f"  FAIL  {mol[:46]:<46} " + '; '.join(errs))
    else:
        print(f"  PASS  {mol[:46]:<46} {ref} {bp:,}  F1={f1}  rank2={r2}")
print(f'\n{len(mols) - bad}/{len(mols)} molecules reproduced')
sys.exit(1 if bad else 0)
