# TdT: hDOM and SGM code

Code from a single study of single-molecule genome mapping, in two parts by imaging arm, with the
aligner both of them need.

| | what it does |
|---|---|
| [`hDOM/`](hDOM/) | Fluorescence (hDOM). Assigns single DNA molecules to chromosomal loci from two-channel optical maps. Candidates come from an FFT correlation scan and from maligner; they are re-ranked by a windowed backbone correlation multiplied by a label-set F1. |
| [`hDOM/sv_confidence/`](hDOM/sv_confidence/) | Decides whether an inversion or a translocation in one hDOM molecule is supported. Imports the scoring engine from `hDOM/hdom/` rather than reimplementing it, so both lines score with the same code. |
| [`SGM/`](SGM/) | SEM (SGM). Traces and measures DNA contours in SEM images, and turns them into `.maps` fragment profiles that maligner can align. |
| [`maligner/`](maligner/) | Not our work. Mendelowitz et al.'s ordered restriction map aligner, bundled so the pipelines run without hunting for a matching build. Shared by `hDOM/` and `SGM/`. See the licence note below. |

## Licence

This repository is GPL v3. It bundles maligner, which is GPL v3, and `SGM` imports its
Python library directly, so a permissive licence would not have been honest here.

## Citing

The dual-channel mapping this pipeline implements, including the `cc_rg2 = cc_r · cc_g²`
candidate score, is from:

> Jo K et al., *DOM: Dual Optical Mapping Integrating Sequence-Specific Barcodes and A/T Density Profiles*, Small Methods, 2026. doi:10.1002/smtd.70883

Please cite that paper alongside this software.

`maligner/` carries its own `LICENSE` and belongs to its authors:
Mendelowitz L, Schatz MC et al., *Maligner: a fast ordered restriction map aligner*,
Bioinformatics 32(7):1016-1022, 2016. https://github.com/LeeMendelowitz/maligner

## Installing

Python 3.9 or later. Each tool has its own `requirements.txt`; the mapping pipeline additionally
needs CuPy and an NVIDIA GPU, because its candidate scan runs there.

    pip install -r hDOM/requirements.txt                 # CPU-only; the GPU scan needs CuPy, see below
    pip install -r hDOM/sv_confidence/requirements.txt
    pip install -r SGM/requirements.txt

One environment holds all three. `sv_confidence` pins exact versions where the other two give
floors, and every pin satisfies both floors (numpy 1.23.5, scipy 1.10.1, matplotlib 3.9.4,
Pillow 9.3.0, tifffile 2021.7.2), so installing them together resolves to those versions rather
than failing. The pins are there because the structural-variant nulls were calibrated against
them; installing `sv_confidence` last keeps them.

maligner is already built for Linux under `maligner/build/bin`. Put it on `PATH`:

    export PATH="$PWD/maligner/build/bin:$PATH"

On other platforms, rebuild from `maligner/`'s own README. `SGM` finds the bundled Python
library on its own; `MALIGNER_PY` overrides it if you point at another build.

## Checking the install

    cd hDOM && python3 -m pytest tests/

Nine tests: eight structural checks, and one that runs the bundled example end to end and compares
it with the recorded result. About half a minute, no download needed.

## What the examples show

`hDOM/example/` searches two chromosome 20 maps; `hDOM/example_genome/` searches all 47. The
molecules are the same three, so running both answers one question: does the answer depend on
having narrowed the search first. It does not — chromosome, haplotype, locus, F1 and the runner-up
haplotype come out identical.

The 47 reference maps are not in this repository: 331 MB to download, about 900 MB unzipped.
They are deposited as their own Zenodo record, separate from this software, because they are data
derived from a public assembly rather than part of the program. The same record carries the SGM
maps `SGM` aligns against (`SGM_reference_maps_HG002v12.zip`):

**hDOM and SGM reference maps for HG002 v1.2 (Nb.BssSI)**
<https://doi.org/10.5281/zenodo.22929196>

This is the map `hDOM/example/` and `hDOM/example_genome/expected.json` are recorded against.
An earlier assembly version does not reproduce those values.

`hDOM/example_genome/README.md` shows how to unpack and use them.
