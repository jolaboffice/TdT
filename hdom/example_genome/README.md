# Genome-wide example

The three molecules are named `A`, `B` and `C`. They are real HG002 molecules; the experiment
names they were collected under are kept in the `_provenance` block of `expected.json`.

The same three molecules as `example/`, searched against the whole HG002 reference — 47 chromosome
maps, both haplotypes — instead of the two chromosome 20 maps bundled with the package.

That is the point of running both. `example/` shows the pipeline works; this shows it works when
it is not told where to look. All three molecules land on the same locus either way, so the answer
does not depend on having narrowed the search first.

## The reference map is not in this repository

47 chromosome maps are about 860 MB unzipped, more than a repository should hold. Download
`hDOM_reference_map_HG002v12_1bp.zip` from the reference-map Zenodo record, which the top-level
README points to, unzip it, and point `--ref-map` at the folder it produces. The v1.1 map is a
different assembly and does not reproduce `expected.json`.

## Running it

    unzip hDOM_reference_map_HG002v12_1bp.zip -d ref_map
    python3 hDOM_pipeline.py \
        --molecules example_genome/molecules \
        --ref-map   ref_map \
        --output    out_genome \
        --shards 3 --logp
    python3 example_genome/check.py out_genome

Each molecule is scanned against every chromosome at every scale in the grid, so expect minutes
rather than the seconds the two-chromosome example takes.

## What reproduces, and what does not

`check.py` compares the chromosome and haplotype, the locus to within 2 kb, the F1, and the
haplotype ranked second. Those are identical between the two examples.

`h_score` is not compared. It is a null-model p-value: the null draws random positions from the
reference, so its scale depends on how much genome there is to draw from, and on whether `--logp`
was given at all. The same correct alignment scores about -1.4 against two chromosomes and about
-9 to -26 against all 47. The larger value is the more meaningful one -- it is the statement that
the molecule was located in ~6 Gb rather than in 130 Mb -- but it is not a number a test can hold
fixed across setups.
