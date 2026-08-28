# Method specification

Everything needed to describe this pipeline in a paper: what each stage computes, every constant
it uses, and what has and has not been measured. Values are read from the code in this package.

Terminology: **red** = the AT-dye backbone intensity profile along the molecule; **green** = the
nick-label channel. A **reference** is one chromosome rendered as the same two channels at
`bpp` bases per pixel.

---

## 0. Input and output

Input per molecule: a two-frame 16-bit TIFF, frame 0 red, frame 1 green, one molecule per file.
Reference: one `{name}.tif` per chromosome plus `genome.maps`. Default `bpp = 200`.

Output per molecule: a ranked table of candidate loci, each with `ref_name`, `aligned_bp`
(reference start), `dir` (±1, molecule orientation), `scale`, and the scores below.

The dataset the numbers in §7 come from: 3,148 human molecules (HG002), median length ~150 kb,
aligned against 47 reference chromosomes (22 autosome pairs plus chrX and chrY).

---

## 1. Label detection

A random forest scores every candidate green peak in the molecule, replacing a fixed intensity
threshold. It is retrained rather than re-tuned when the imaging changes.

Candidate pool — deliberately over-inclusive, so faint true labels are not lost before scoring:

```
g = max over rows of the green frame              (1-D profile, length W)
find_peaks(g, height    = mean(g) + 0.02·(max(g) − mean(g)),
              prominence = 0.02·max(g),
              width      = 1,
              distance   = 3)
```

Twelve features per candidate, all scale-free or normalised so the model does not learn one
exposure level: `h_norm`, `h_over_med`, `prom_frac`, `fwhm`, `contrast_abs`, `contrast_ratio`,
`dleft`, `dright`, `density`, `pos_norm`, `voff`, `voff_dev`. Here `contrast_*` compare the peak
against the local background in a ±20 px window excluding ±3 px around the peak; `dleft`/`dright`
are distances to neighbouring candidates; `density` counts candidates within ±30 px; `voff` is
the peak's vertical offset from the traced red backbone, which separates labels on the DNA from
background speckle.

Classifier: `RandomForestClassifier`, 400 trees, unlimited depth. A candidate is called a label
when `P(label) ≥ 0.5`.

---

## 2. Straightening

Molecules are curved and drift across the field, so an image column is not a position along the
DNA. The molecule is resampled along its own backbone.

**Backbone trace.** Per column, the intensity-weighted centre of the red frame within ±4 px of
the previous column's centre, seeded at the brightest column and grown outwards in both
directions. Columns whose window carries no signal above `max(0.15·typical_peak, 0.02·max)` are
retried once with a ±20 px window, so a dim stretch or a sharp terminal bend does not truncate
the trace.

**Arc-length resampling.** The trace is smoothed over ±6 px, its arc length `S` is accumulated,
and the molecule is sampled at `ceil(S)+1` points evenly spaced *along the arc*, each sample
taken on the normal to the backbone. Both frames go through the same transform.

**Label carry-through.** Labels are detected on the raw image (§1) and mapped through the same
arc-length transform, then snapped to the nearest swept green peak within ±4 px. A label with no
peak within ±4 px keeps its mapped position, so the label count is preserved rather than warped.

**Green re-bake.** For scoring, the green frame is replaced by Gaussians of amplitude
`C_AMP = 846` at the label positions, each with the FWHM measured at that label in the swept
image (`σ = FWHM / 1.6651`). The displayed copy (`_orig.tif`) keeps the real green.

Measured on 2,180 molecules: arc length exceeds the backbone's horizontal span by a median of
0.2%, 90th percentile 0.5%, maximum 3.8%. **Straightening is not a length-recovery step** — its
effect is that the profile is sampled along the DNA rather than down image columns.

---

## 3. Candidate generation

Two independent searches over every reference, each returning its top 10. Their union is the
pool the re-rank scores. Both use the full label set (the terminal rule in §5 applies later).

**pcc.** The molecule's red profile is correlated against each reference at a grid of shifts,
directions and scales; the green profiles likewise. Candidates are ranked by

```
cc_rg2 = cc_r · cc_g²
```

so the label channel enters quadratically. This matters for interpretation: pcc is not a
red-only scorer, and the green term is what keeps many low-frequency red matches from winning.

**maligner.** `maligner_dp` aligns the molecule's inter-label fragment lengths against
`genome.maps`:

```
--max-alignments 10   --max-alignments-per-reference 3
--query-max-misses 3  --ref-max-misses 6
--min-query-rescaling 0.82   --max-query-rescaling 1.18
```

Candidates whose mapped reference span overlaps an assembly N-gap (≥100 bp run of N) are dropped
from every method's ranking, since neither DNA nor labels exist there.

---

## 4. Local re-ranking

Every candidate is re-scored. This is the step the method turns on.

**Windowed red correlation.** A whole-molecule Pearson correlation is dominated by the
low-frequency envelope: uneven illumination and bleaching along the molecule, and equally
low-frequency structure in the reference. Instead the molecule is cut into overlapping windows
and each window is correlated at its own best local shift:

| | value | in bases |
|---|---|---|
| `WIN` | 50 px | 10 kb |
| `STEP` | 25 px | 5 kb |
| `MAXS` | 15 px, searched in steps of 3 | ±3 kb |

```
cc_r^loc = mean over windows of  max over shifts s∈[−MAXS, MAXS]  corr(molecule[w], reference[w+s])
```

Letting each window find its own shift is what absorbs non-uniform stretch; taking the
correlation inside a window is what removes the envelope. Red uses the whole molecule — the
backbone at the ends is real signal.

**Label agreement.** Molecule labels are mapped into reference coordinates by
`shift + scale · x` (or the mirrored form when `dir = −1`) and greedily matched to reference
labels within `TOL_PX = 12` reference pixels (2.4 kb). `F1 = 2M / (2M + FP + FN)`.

**Offset and scale re-optimisation.** For each candidate, offset is searched over ±25 px in
steps of 5 px (±5 kb in 1 kb steps) and scale over an 8-point grid anchored at 0.905:

```
0.905  0.930  0.955  0.980  1.005  1.030  1.055  1.080      (range 0.905–1.100, step 0.025)
```

The pair maximising `cc_r^loc · F1` is kept. Selecting on `F1` alone lets a scale that matches
more labels win where the red profile disagrees.

**Ranking.** Candidates are ordered by `cc_r^loc · F1`.

`scale` is the stretch of the molecule relative to `bpp`: a molecule pixel covers `bpp · scale`
bases, so a molecule stretched beyond its contour length has `scale < 1`.

---

## 5. The terminal 5 kb

Arc-length resampling maps labels through the traced backbone `[x0, x1]`; a label outside that
span is clamped onto pixel 0 or `w−1`. The terminal 5 kb (`TRIM_PX = 25 px`) is where that
artifact lands, and a clamped label becomes a query fragment maligner must explain.

| step | terminal labels |
|---|---|
| pcc candidates | kept |
| maligner candidates (§3) | kept |
| `cc_r^loc × F1` re-rank | **dropped** |
| the maligner row reported per molecule | **dropped** — a second `maligner_dp` run |

Both searches therefore cover the same candidate pool, while the scores that decide the ranking
ignore the unreliable labels. The reported maligner pass uses `--max-alignments 3000` and keeps
the 10 best by `total_rescaled_score` (ascending is better). A molecule with fewer than two
surviving labels gets no trimmed pass and reports its full-label row.

---

## 6. Null model (h-score)

Optional (`--logp`). For each molecule, `K = 5000` decoy windows are drawn from every chromosome
*except* the one the molecule matched, each scored by the same `cc_r^loc × F1`. The observed
score is converted to a z-score against that distribution and then to `log₁₀ P`, corrected for
the number of independent windows a genome of this size offers a molecule of this length. More
negative is more significant. A reference holding only the matched chromosome has no population
to draw from and reports `nodecoy`.

---

## 7. Reference construction

At `bpp` bases per pixel and `n = ceil(len(seq)/bpp)` pixels, from `numpy.histogram` with edges
`arange(0, len+bpp, bpp)`:

```
red    = histogram of A/T positions            → gaussian_filter1d(σ = 1)          → rescale 5–255
green  = histogram of enzyme sites (motif and its reverse complement), binarised to
         presence per pixel                    → gaussian_filter1d(σ = 1) → sqrt() → rescale 5–255
```

The square root is not cosmetic. It makes overlapping sites combine sub-linearly — two sites in
adjacent pixels reach 1.27× an isolated site rather than the 1.78× plain summation gives — and
it widens the apparent kernel, since the square root of a σ=1 gaussian is a σ=√2 one. An
isolated site therefore prints a 9-pixel footprint `[7, 21, 63, 127, 162, 127, 63, 21, 7]`.

`genome.maps` is that green channel reduced to inter-site fragment lengths: the profile is
thresholded at `mean + 0.1·(max − mean)`, `find_peaks` gives the sites, and the line is
`name <TAB> length_bp <TAB> n_fragments <TAB> frag1 <TAB> …`.

Default motif `CTCGTG` (Nb.BssSI). Two references built with even slightly different parameters
disagree about which homologue a molecule came from, so a rebuild must be checked against the
reference the results were produced with.

---

## 8. What has been measured

On the 3,148-molecule set:

| observation | value |
|---|---|
| pcc rank 1 differs from the re-ranked rank 1 | **1,368 / 3,148 = 43.5%** |
| — by group | dense-label 7/7, sparse-label 61/71, B 46.4%, A 42.6%, XY 17.3% |
| the highest whole-molecule `cc_r` candidate is not the answer, and the answer beats it on both `cc_r^loc` and `F1` | **667 / 1,828 = 36.5%** (group A) |
| — decoy vs answer, medians | `cc_r` 0.47 vs 0.37; `F1` 0.49 vs 0.88 |
| — acrocentric chromosomes among those decoys | 47.4% against a 24.1% base rate (**1.97×**) |
| trap rate vs one-sided red intensity | 26.3% → 41.4% → **60.0%** as \|asymmetry\| rises through 0.2, 0.4, 0.6 |
| arc length over backbone span | median +0.2%, max +3.8% (n = 2,180) |
| straightening determinism | 7 of 8 molecules byte-identical on re-run; the eighth differed by 1 px |

Worked example, molecule `0312_5_6_merged_t028_t029_mol001_147kb` (155 kb, 15 labels):

| candidate | `cc_r` | `cc_r^loc` | `F1` |
|---|---|---|---|
| pcc rank 1, chr07_p 56.60 Mb | 0.55 | 0.21 | 0.60 |
| re-ranked rank 1, chr14_p 47.39 Mb | 0.09 | **0.53** | **0.97** |

For the pcc pick the reference is near-constant over the first 106 kb (std 2.6 on a 5–255
scale); reducing that alignment to a two-level step reproduces `cc_r = 0.54`, i.e. essentially
all of its whole-molecule correlation. For the re-ranked pick the same two-level envelope gives
−0.54: the envelope works against the locus whose features actually match.

---

## 9. Limits

- **The reference must match.** Rebuilding it with a nearly-correct recipe (red exact, green
  kernel exact, only the normalisation differing over ~4% of pixels) moved 1 of 7 chr21
  molecules to its homologue. With the exact reference, 8 of 8 reproduced.
- **"Correct" is the pipeline's own answer** in §8 unless stated otherwise. The label `F1` is
  independent of the red channel and is what makes those comparisons more than circular, but no
  external ground truth was used. Simulated molecules from known loci would settle it.
- **Detection threshold is fixed at 0.5.** On a molecule whose green background sits at 28% of
  its peak, ten candidates at 1.3–1.9× background pass at probabilities 0.50–0.75 while the
  twelve real labels score 0.99–1.00. The model separates them; the threshold does not.
- **A GPU is required.** The pcc scan runs through CuPy and exits without it.
- **One sample, one enzyme, one instrument.** Everything here is HG002 with Nb.BssSI at
  200 bp/pixel.
