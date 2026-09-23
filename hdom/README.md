# hDOM pipeline

Aligns single-molecule optical-map images to a reference genome using **both** channels of the
molecule: the AT-dye red backbone profile and the green nick labels. Give it molecule TIFFs and a
reference map; it returns the ranked chromosomal assignment per molecule plus an HTML view.

```
python3 hDOM_pipeline.py --molecules RAW_DIR --ref-map REF_MAP_DIR --output OUT_DIR
python3 hDOM_pipeline.py RAW_DIR          # short form: ref_map/ and out/ beside RAW_DIR
```

The short form fills in `--molecules` from the positional argument and defaults `--ref-map` to
`RAW_DIR/../ref_map` and `--output` to `RAW_DIR/../out`. Either form accepts the options below.

## Input

| | |
|---|---|
| `--molecules` | 2-frame 16-bit TIFFs. Frame 0 = red (AT-dye DNA backbone), frame 1 = green (nick label). `{name}.tif` or `{name}_orig.tif`. |
| `--ref-map` | One `chr*.tif` per reference chromosome, same 2-channel convention, at 200 bp/pixel — **plus `genome.maps`** (see below). |
| `--output` | Destination. Created if absent. |

`genome.maps` is the same references expressed as fragment lengths, which is the only form
`maligner_dp` reads: one tab-separated line per chromosome — name, length in bp, then the
inter-label fragment lengths. The `chr*.tif` alone are not sufficient, so it must sit in the
`--ref-map` folder next to the TIFFs.

## Building a reference

```
python3 make_reference.py --fasta FASTA_DIR --out REF_MAP_DIR --bpp 200 --motif CTCGTG
```

One `{name}.fasta` per chromosome in; `{name}.tif` and `genome.maps` out. `--motif` is the
nicking enzyme's recognition sequence (default Nb.BssSI); the reverse complement is searched
too. Per pixel the red channel is the A/T count and the green channel is site presence, each
smoothed with a sigma-1 gaussian, the green then square-rooted, and both rescaled to 5..255.

`--validate EXISTING_REF_DIR` byte-compares the result against a reference you already have.
Run it whenever the recipe or its parameters change: two references that differ only slightly
still disagree about which homologue a molecule came from, since homologous chromosomes are
the hardest pair to tell apart.

## Haplotype chains (optional, `--chains`)

The homologue panel in the view, and the homologous-pair criterion used to evaluate the ranking,
need to know when two candidate loci are the same place on the two haplotypes. That comes from a
per-chromosome minimap2 alignment of one haplotype against the other, with the CIGAR kept:

```
for n in 01 02 ... 22; do
    minimap2 -cx asm5 -t 8 chr${n}_m.fasta chr${n}_p.fasta > chains/chr${n}.paf
done
```

Target is maternal, query is paternal, so a lift with `src='p'` goes query to target and `src='m'`
goes the other way (`hdom/paf_liftover.py`). `-c` is required: the lift interpolates inside
alignment blocks and needs the `cg:Z` CIGAR, which minimap2 only emits with `-c`. The chains used
in the paper were built with minimap2 2.17-r941 and `asm5`, the preset for sequences a few percent
divergent, which is the right register for two haplotypes of one genome.

Pass the folder as `--chains chains/`. Without it the pipeline still runs; only the homologue
panel and the homologue-aware locus test are skipped.

## Example

`example/` holds three molecules and a two-chromosome reference (chr20 maternal and paternal),
built by `make_reference.py` from FASTA. It runs in about a minute:

```
python3 hDOM_pipeline.py --molecules example/molecules --ref-map example/reference \
                        --output /tmp/hdom_example --shards 3
```

All three land on the locus a full 47-chromosome run put them on. `example/expected.json`
records that, and `tests/test_example.py` checks it.

Note this example searches two chromosomes, not a genome, so it demonstrates the mechanics
rather than the specificity of a genome-wide search. `--logp` needs more than that: the null
draws decoy windows from every chromosome except the one the molecule matched, so on a chr20-
only reference it reports `nodecoy` and leaves `log(P)` as the re-rank value.

## Tests

```
python3 tests/test_release.py    # structure and packaging, seconds, no data needed
python3 tests/test_example.py    # end to end on example/, about a minute, needs maligner_dp
```

Options: `--shards N` (parallel workers), `--logp` (K=5000 random-null log P per molecule),
`--chains DIR` (minimap2 haplotype PAF, enables the homologue panel in the view).

## Output, per molecule

| file | contents |
|---|---|
| `_matching_results.csv` | Ranked hits from every scorer. `method` = `pcc` / `local` / `mal`. |
| `_view.html` | Alignment view: molecule over the reference window, red profile overlay, label toggle. |
| `.tif` | Straightened molecule, green re-baked as width-preserving Gaussians (the tif that was scored). |
| `_orig.tif` | Straightened molecule, real green channel (display). |
| `_mlpeaks.json` | ML label positions, px, on the straightened molecule. |
| `_maligner_dp.txt` | Raw `maligner_dp` output. |
| `index.php` | Browsable list of all molecules in the folder. |

Key `_matching_results.csv` columns: `local_cc_r` (windowed red correlation), `F1` (label
agreement), `log(P)` (**h_score**, null-model significance -- the headline hDOM score), and for
`method=mal` the maligner-native `Matched / FP / FN / F1 / sizing / m_score`.

## Requirements

- **An NVIDIA GPU with a working CUDA runtime.** The pcc stage correlates each molecule against
  every reference on the GPU through CuPy and exits if CuPy cannot load, which leaves the run
  with maligner candidates only.
- `maligner_dp` on `PATH` — the label-only scorer, from Mendelowitz et al., *maligner*
  (https://github.com/LeeMendelowitz/maligner). The results in the paper were produced with
  commit `f6a1f34` (`maligner_paper-33-gf6a1f34`), which reports itself as version 0.4. The
  version matters: `hdom/candidates/hDOM_mal_Q.py` widens several of its defaults
  (`--query-max-misses 3`, `--ref-max-misses 6`, rescaling 0.82–1.18,
  `--max-alignments-per-reference 3`), and a build that parses those flags differently would
  change which candidates reach the re-rank.
- Python 3 with `numpy scipy pandas tifffile Pillow scikit-learn joblib numba matplotlib cupy`
  (`pip install -r requirements.txt`; CuPy is listed separately in `requirements-gpu.txt`
  because the plain package builds from source, so install the wheel matching your CUDA
  version instead, e.g. `pip install cupy-cuda12x`)

## How it works

1. **Label detection** — a trained random forest (`hdom/label_model.joblib`) calls green peaks per
   molecule instead of one global threshold, so exposure differences between fields don't change
   the label count.
2. **Straightening** — the red backbone is traced, then the molecule is resampled along its own
   arc length. A curved molecule is measured at its true contour length, not its bounding box.
   Labels are carried through the same transform and snapped (±4 px) to the swept green peaks, so
   the label count is preserved rather than warped.
3. **Candidates** — `hdom/candidates/hDOM_calc_Q.py` produces two independent top-10 lists over
   the whole genome: pcc, ranked by `cc_r·cc_g²` (red profile correlation weighted by the green
   label correlation), and `maligner_dp` on the label map. These are the two scorers the paper
   compares; their union is the pool the re-rank scores.
4. **Re-rank** — every candidate is re-scored locally: `local_cc_r` (red correlation in sliding
   windows, tolerant of local stretch) × label `F1`, converted to a null-model `log(P)`.
   4b. **Reported maligner row** — a second `maligner_dp` run on the same molecule with the
   terminal-5 kb labels dropped (see below).
5. **View** — per-molecule HTML with the molecule drawn over its reference window.

### The terminal 5 kb

Straightening maps each label through `np.interp` over the traced backbone `[x0,x1]`, so a label
falling outside that span clamps onto pixel 0 or `w2-1`. The terminal 5 kb is exactly where that
artifact lands, and a clamped label becomes a query fragment maligner is forced to explain.

The window is applied at two specific points, and **not** during candidate generation:

| step | terminal 5 kb |
|---|---|
| pcc candidates, top 10 by `cc_r·cc_g²` | **kept** |
| maligner candidates, top 10 (inside `hDOM_calc_Q`) | **kept** |
| local re-rank `local_cc_r × F1` | **dropped** (`mdom_reanchor.label_f1`, `TRIM_PX`) |
| the maligner row reported for the molecule | **dropped** (second pass, stage 4b) |

Both scorers therefore search over the same candidate pool, built from every label, while the
scores that decide the ranking ignore the unreliable terminal ones. The reported maligner row
comes from a second `maligner_dp` run on the trimmed label set (`_maligner_dp_trimmed.txt`); the
candidate run is kept alongside it as `_maligner_dp.txt`.

A molecule left with fewer than two labels after trimming has no alignable map, so it gets no
trimmed pass and reports its full-label maligner row instead. Dropped positions are recorded in
`removed_peaks.json` and drawn as white circles over the green label diamonds in the view.

## Layout

```
hDOM_pipeline.py           entry point
make_reference.py          builds a reference map from FASTA
METHODS.md                 every stage, every constant, and what has been measured
hdom/
  straighten_align.py      stages 1-5, the pipeline proper
  local_align_pipeline.py  local re-rank (local_cc_r, F1, offset/scale optimisation)
  mdom_reanchor.py         label_f1 + the 5 kb trim rule
  logp_null.py             K=5000 random-null log(P)
  Label_Detect_ML_based.py ML label detector       (+ label_model.joblib)
  hDOM_view.py             per-molecule HTML view
  figure3_view.py          folder viewer + registries  (+ index.php)
  qc.py                    per-molecule image quality checks
  polyn_filter.py          poly-N reference masking (optional; set HDOM_FASTA_DIR)
  paf_liftover.py          haplotype PAF lift-over for the homologue panel
  candidates/              red-PCC + maligner candidate generation, vendored unmodified
```

Environment overrides. All optional, and the defaults are what the paper used, so a normal run
sets none of them.

| | |
|---|---|
| `HDOM_REFMAP`, `HDOM_CHAINS` | reference and chain folders; `hDOM_pipeline.py` sets both for you |
| `HDOM_FASTA_DIR` | enables poly-N reference masking (`polyn_filter.py`) |
| `HDOM_LABEL_ML` / `HDOM_LABEL_SRC` | only for `Label_Detect_ML_based.py --retrain` |
| `HDOM_LABEL_EDIT` | folder of `{molecule}_labeledit.json` curated label lists (see §1 of METHODS.md) |
| `HDOM_TRIM_PX` | terminal pixels dropped from scoring; default `5000/bpp` = 25 px. `0` keeps terminal labels |
| `BBFIX` | `0` disables the backbone-repair pass entirely (see §2 of METHODS.md) |
| `BBFIX_MAXSTEP`, `BBFIX_MINRUN`, `BBFIX_THRESH` | its three constants: `2.0` px, `10` columns, `40` |
| `SWEEP_SM` | backbone smoothing half-window in columns, default `6` |

## Reproducibility

Straightening is deterministic: re-running the pipeline on the same input reproduces the
straightened TIFF byte-for-byte. One caveat — arc length is summed in floating point and the
output width is `ceil(arc_length)+1`, so a molecule whose arc length lands very near an integer
can differ by one pixel between machines or BLAS thread counts. The pipeline pins
`OMP_NUM_THREADS=1` to keep this stable.

## Defaults

These are the settings the HG002 v1.2 run of all 3,164 molecules used. The last column gives the
switch back to the simpler behaviour each one replaces.

| Setting | Default | Simpler alternative |
|---|---|---|
| `genome.maps` site position | snapped to the true motif bp | `make_reference.py --grid binned` (`peak_index × bpp`, a 200 bp lattice) |
| `BBFIX_THRESH` | **10** | `BBFIX_THRESH=40` |
| `HDOM_TRIM_RULE` | **`outermost`** | `HDOM_TRIM_RULE=window`, a fixed 5 kb window |

**Site snapping.** `peak_index × bpp` pins every site to a multiple of `bpp`. A single upstream
indel shifts the phase of the downstream lattice, and then only the sites near a cell boundary jump
a whole cell. Over the 703,252 sites that pair up between HG002 v1.1 and v1.2 the true movement is
6 bp at the median, while a 200 bp map records "0 bp" for 95.0% of them and "≥200 bp" for 5.03%
with nothing in between, so 4.5% of sites carry an error of about 194 bp that does not exist.
maligner consumes those fragment lengths directly, which makes its answer sensitive to the
reference version: across v1.1 and v1.2 the call is unchanged for 67.0% of molecules on a binned
map against 96.2% on a snapped one. Snapping moves site coordinates only. The site set itself is
unchanged in count, order and neighbours, and the TIFFs are not affected.

**BBFIX_THRESH.** At 40 the backbone fix fires on 33 of the 3,164 molecules (1.0%), while roughly
6-8% of them have a trace that zigzags through a faint stretch.

**HDOM_TRIM_RULE.** A fixed 5 kb window can cut away real interior sites along with the terminal
ones. `outermost` cuts only up to just outside the outermost label cluster at each end.
