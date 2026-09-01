# hDOM-toolkit

Three tools from a single study of single-molecule genome mapping, kept together because two of
them depend on the third and all of them need the same aligner.

| | what it does |
|---|---|
| [`hdom/`](hdom/) | Assigns single DNA molecules to chromosomal loci from two-channel optical maps. Candidates come from an FFT correlation scan and from maligner; they are re-ranked by a windowed backbone correlation multiplied by a label-set F1. |
| [`sv_confidence/`](sv_confidence/) | Decides whether an inversion or a translocation in one molecule is supported. Imports the scoring engine from `hdom/` rather than reimplementing it, so both lines score with the same code. |
| [`sem_length/`](sem_length/) | Traces and measures DNA contours in SEM images, and turns them into `.maps` fragment profiles that maligner can align. |
| [`maligner/`](maligner/) | Not our work. Mendelowitz et al.'s ordered restriction map aligner, bundled so the pipelines run without hunting for a matching build. See the licence note below. |

## Licence

This repository is GPL v3. It bundles maligner, which is GPL v3, and `sem_length` imports its
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

    pip install -r hdom/requirements.txt        # add cupy-cuda12x or the wheel for your CUDA
    pip install -r sv_confidence/requirements.txt
    pip install -r sem_length/requirements.txt

One environment holds all three. `sv_confidence` pins exact versions where the other two give
floors, and every pin satisfies both floors (numpy 1.23.5, scipy 1.10.1, matplotlib 3.9.4,
Pillow 9.3.0, tifffile 2021.7.2), so installing them together resolves to those versions rather
than failing. The pins are there because the structural-variant nulls were calibrated against
them; installing `sv_confidence` last keeps them.

maligner is already built for Linux under `maligner/build/bin`. Put it on `PATH`:

    export PATH="$PWD/maligner/build/bin:$PATH"

On other platforms, rebuild from `maligner/`'s own README. `sem_length` finds the bundled Python
library on its own; `MALIGNER_PY` overrides it if you point at another build.

## Checking the install

    cd hdom && python3 -m pytest tests/

Nine tests: eight structural checks, and one that runs the bundled example end to end and compares
it with the recorded result. About half a minute, no download needed.

## What the examples show

`hdom/example/` searches two chromosome 20 maps; `hdom/example_genome/` searches all 47. The
molecules are the same three, so running both answers one question: does the answer depend on
having narrowed the search first. It does not — chromosome, haplotype, locus, F1 and the runner-up
haplotype come out identical.

The 47 reference maps are not in this repository: 330 MB to download, about 900 MB unzipped. They
are deposited as their own Zenodo record, separate from this software, because they are data
derived from a public assembly rather than part of the program:

**hDOM reference map for HG002 v1.1 (Nb.BssSI, 200 bp/pixel)**
<https://doi.org/10.5281/zenodo.22121923>

`hdom/example_genome/README.md` shows how to unpack and use them.
