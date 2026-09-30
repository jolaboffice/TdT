#!/usr/bin/env python3
"""Poly-N region detection & rank-exclusion for the mDOM pipeline.

Ported from hDOM_kts (hDOM_calc_kts.py). A reference alignment whose mapped
molecule span overlaps an assembly N-gap (>=100bp run of 'N' in the haplotype
FASTA) is spurious and must be excluded from the rank. Works in bp coordinates
(aligned_bp and FASTA positions share the same reference bp system).
"""
import os
import re
import glob
import json

BPP = 200
MIN_N_BP = 100
FASTA_DIR = os.environ.get('HDOM_FASTA_DIR', '')   # per-chromosome FASTA; poly-N filter only
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'n_regions.json')


def find_n_regions(fasta_dir=FASTA_DIR, min_n_bp=MIN_N_BP):
    """{chr_name: [(start_bp, end_bp), ...]} for N-runs >= min_n_bp."""
    out = {}
    for fa in glob.glob(os.path.join(fasta_dir, '*.fasta')):
        chrom = os.path.splitext(os.path.basename(fa))[0]
        parts = []
        with open(fa, encoding='utf-8') as f:
            for line in f:
                if not line.startswith('>'):
                    parts.append(line.strip().upper())
        seq = ''.join(parts)
        regions = [(m.start(), m.end()) for m in re.finditer('N{%d,}' % min_n_bp, seq)]
        if regions:
            out[chrom] = regions
    return out


def load_n_regions(use_cache=True):
    if use_cache and os.path.exists(CACHE):
        with open(CACHE, encoding='utf-8') as f:
            return {k: [tuple(r) for r in v] for k, v in json.load(f).items()}
    nr = find_n_regions()
    with open(CACHE, 'w', encoding='utf-8') as f:
        json.dump(nr, f)
    return nr


def overlaps_n(ref_name, aligned_bp, mol_width_px, scale, n_regions, bpp=BPP):
    """True if molecule's mapped reference span overlaps any N-region of ref_name.
    Span is direction-independent (start..start+len), matching hDOM_kts."""
    regions = n_regions.get(ref_name)
    if not regions:
        return False
    mol_start = float(aligned_bp)
    mol_end = mol_start + float(mol_width_px) * float(scale) * bpp
    for n_start, n_end in regions:
        if max(mol_start, n_start) < min(mol_end, n_end):
            return True
    return False


if __name__ == '__main__':
    nr = load_n_regions(use_cache=False)
    tot = sum(len(v) for v in nr.values())
    bp = sum(e - s for v in nr.values() for s, e in v)
    print('chromosomes with N-regions: %d' % len(nr))
    print('total N-regions: %d, total %d bp' % (tot, bp))
    for k in sorted(nr):
        print('  %-10s %d region(s), %d bp' % (k, len(nr[k]), sum(e - s for s, e in nr[k])))
