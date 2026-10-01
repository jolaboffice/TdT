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

**Curated label lists.** Setting `HDOM_LABEL_EDIT` to a folder makes the pipeline read
`{molecule}_labeledit.json` from it, and use the positions in its `final` array in place of what
the classifier returned for that molecule. Positions are raw-frame pixels, the same frame the
detector reports, so straightening carries them onward exactly as it carries detected ones: a
hand-checked molecule differs from an automatic one only in where its labels came from. Molecules
with no such file are unaffected. The variable is unset by default; it exists so that a molecule
whose labels were checked by eye can be re-run without editing the image.

A sidecar is used in one of two ways, and the difference matters because it decides whether the
straightener's `snap` correction applies (section 2). A sidecar that only **pins** the detector's
own positions -- so that a reference change cannot also change the labels -- carries positions that
came out of the detector, and they are snapped like any detected label. A sidecar whose positions a
**person placed** is exempt from snap, because snap's peak list uses the detector's own 25%-of-max
floor: a label deliberately put on a fainter peak would be dragged onto whichever accepted peak lies
within 4 px, merging it with a neighbour and discarding the edit. Set `"hand_placed": true` or
`false` in the sidecar to state which it is. Without that field the file is read as a pin when
`pinned_to` holds a bare provenance tag (a version string, no spaces) and the `final` array still
holds as many labels as `n_detect` records; anything else is treated as hand-placed. Pinning is by
far the commoner case, so reading every sidecar as hand-placed would shift a label of nearly every
molecule in a run by a pixel or two.

---

## 2. Straightening

Molecules are curved and drift across the field, so an image column is not a position along the
DNA. The molecule is resampled along its own backbone.

**Backbone trace.** Per column, the intensity-weighted centre of the red frame within ±4 px of
the previous column's centre, seeded at the brightest column and grown outwards in both
directions. Columns whose window carries no signal above `max(0.15·typical_peak, 0.02·max)` are
retried once with a ±20 px window, so a dim stretch or a sharp terminal bend does not truncate
the trace.

**Arc-length resampling.** The trace is smoothed over ±6 px (`SWEEP_SM`), its arc length `S` is
accumulated, and the molecule is sampled at `ceil(S)+1` points evenly spaced *along the arc*, each
sample taken on the normal to the backbone. Both frames go through the same transform.

**Backbone repair (BBFIX).** Where the red signal is locally faint the tracer can step off the
strand onto a nearby speck, and arc-length resampling then bakes that excursion into the
straightened image as a bump the molecule never had. The repair is applied after straightening,
not before: the molecule is straightened once with the trace above, that result is tested, and only
a molecule that fails the test is straightened a second time with a repaired trace.

The repaired trace drops a column whose centroid sits more than `BBFIX_MAXSTEP = 2.0 px` from the
previous column's, rather than accepting it. Three median-filter passes (windows 101, 61 and 31
columns, widest first, so a narrow excursion cannot drag the median toward itself) then drop
surviving points more than 2.0 px from the local median, and a surviving island shorter than
`BBFIX_MINRUN = 10` columns with a rejected run at least as long on either side is dropped as
re-acquisition residue.

The test is a product of two quantities, both in pixels:

```
d99       99th percentile of |stock trace − repaired trace| over their common span
band_dev  99th percentile of |c − median₁₅₁(c)|, c = per-column intensity centroid of the
          straightened image the stock trace produced
repair when   d99 · band_dev ≥ BBFIX_THRESH = 40
```

Neither factor alone is sufficient. The two traces disagree on many molecules without the
disagreement reaching the image, and a band wanders for reasons other than a trace excursion; the
product is large only when a real excursion has been baked in. The test reads the straightening
this run produced rather than a stored result, so a re-run judges its own output.

Why the repair is gated rather than always on: a backbone's net vertical drift is small (median
0.004 px per column over 59 molecules, at most 0.034), but the largest *single*-column step in an
ordinary trace is itself about 2 px (median 2.09 px over the same molecules, 90th percentile 7.9).
`BBFIX_MAXSTEP = 2.0` therefore rejects legitimate columns on a large fraction of molecules, and
applying the repaired trace unconditionally would perturb many straightenings to fix a few.

Measured on 585 molecules drawn at random from five groups, the test fires on 17 (**2.9%**); the
product has median 0.82, 90th percentile 15.3, 99th percentile 62.1. A separate 30-molecule set was
run twice, once with the repair forced on for every molecule and once with it disabled. These 30
are not a random sample: the test selects 14 of them, far more than the 2.9% base rate. Among those 14 the
rank-1 `F1` rises in 5 and falls in **none** (mean +0.012, at most +0.06), no molecule changes
chromosome, 4 move to the other haplotype and 10 stay on the same reference with a median shift of
800 bp. Among the other 16, which the repair would also have altered but the test rejects, `F1`
rises in 6 and falls in 3 (−0.03 to +0.13). Separating those two groups is what the threshold buys.

`BBFIX=0` disables the pass, and the straightening is then exactly the trace described above.

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

`cc_rg2`, and the dual-channel scheme it scores — a sequence-specific label channel read together
with an A/T density profile — are from the DOM method this pipeline builds on:

> Jo K et al., *DOM: Dual Optical Mapping Integrating Sequence-Specific Barcodes and A/T Density Profiles*, Small Methods, 2026. doi:10.1002/smtd.70883

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

The window is `HDOM_TRIM_PX` if that is set, and `5000/bpp` otherwise. Everything reported here
uses the default; `HDOM_TRIM_PX=0` turns the trim off and scores every label.

---

## 6. Null model (h-score)

Optional (`--logp`). For each molecule, `K = 5000` decoy windows are drawn from every chromosome
*except* the one the molecule matched (both haplotypes of that chromosome are excluded), each
scored by `cc_r^loc × F1 × w_coh` with the freedom over the scale grid and both directions that a
real candidate also gets. That freedom is the one asymmetry between the two sides of `z`: a decoy is
re-optimised over the grid at its drawn position, while `S_obs` is taken at the scale and direction
the ranking settled on, which need not be a grid point. `μ̂` therefore sits a little above what the
row-score definition alone would give a decoy, and `log₁₀ P` is conservative by that much.
The draw is seeded from `md5(molecule name)`, so a re-run reproduces the same null. A reference
holding only the matched chromosome has no population to draw from and reports `nodecoy`.

Scored rows are the ranked `local` rows and the `H_local` homologue row, so the homologue is
measured against the same null as the loci it is being compared with. The pipeline recomputes this
null even when a `{base}_vsrandom.txt` sidecar already exists, because the re-ranking step
overwrites the `log(P)` column with the within-molecule `null_logp` of §4 on its way past; without
the recompute, the column read back would be that value rather than the `K = 5000` one.

The observed score becomes a p-value **parametrically, not by counting**. A normal distribution
is fitted to the null by moment matching, and `log₁₀ P` is read from its upper tail. `S_obs` is the
ranking score the candidate row carries, `cc_r^loc × F1 × w_coh`, read from the columns the table
shows, so `h_score` orders the candidates exactly as the score column does. With `S₁ … S_K` the
decoy scores:

```
μ̂       = (1/K) Σ Sₖ                     # np.mean
σ̂       = sqrt( (1/K) Σ (Sₖ − μ̂)² )      # np.std, i.e. divisor K, not K−1
z       = (S_obs − μ̂) / σ̂
log₁₀ P = log₁₀(1 − Φ(z))                # scipy.stats.norm.logsf(z) / ln(10)
log₁₀ E = log₁₀ P + log₁₀ N              # N = (total reference bp) / (molecule bp)
```

The hats matter: `μ̂` and `σ̂` are estimates from a finite Monte Carlo sample, not known
parameters. The 5000 draws supply those two numbers and nothing else. The empirical shape of the
null is discarded and no exceedance count is formed, so there is no zero-count case to handle.
`h_score` is the per-position `log₁₀ P`; `log₁₀ E` goes to the sidecar beside it. More negative is
more significant. `N` counts *non-overlapping* windows of the molecule's own length, while the
search tries every offset, so `log₁₀ E` understates the multiple-testing burden rather than
matching it.

**Note what this is not.** The tail is never fitted. It is implied by two moments taken over the
whole null, so the bulk of the distribution — where 99% of the mass sits — dictates the shape far
to the right. This is *not* tail extrapolation in the extreme-value sense (peaks-over-threshold
with a generalised Pareto fit, or block maxima with a GEV fit), and it should not be described as
such: those methods estimate the tail's own shape from the data and are easier to defend than
what is done here.

**The reported values lie far outside what 5000 draws can establish.** An empirical p-value from
`K` draws cannot go below `1/(K+1)`, i.e. `log₁₀ P = −3.70`, and the pipeline routinely reports
values well past it. Counting is not a usable alternative: over 6,392 molecules the real
placement outscores **every one** of its 5000 decoys in **97.8%** of cases, so the empirical
estimate collapses to the same `< 1/5001` for almost every molecule and cannot order them. The
fitted tail exists to recover that ordering. Observed `z` runs 3.39 to 9.34 over the 5–95% range
(median 5.97, maximum 29.91).

Sample size is not the limiting factor. Bootstrapping one null (A004, 4000 resamples) gives `μ̂`
to 0.40% and `σ̂` to 1.84%, and a 95% interval on `log₁₀ P` of [−24.91, −21.87], about ±1.5 on a
value of −23.28. That amplification is expected: for large `z` the normal tail behaves as
`log₁₀ P ≈ −z²/(2 ln 10)` while `z ∝ 1/σ̂`, so a relative error in `σ̂` enters `log₁₀ P` roughly
doubled. Raising `K` shrinks it slowly and does not touch the real problem — `K = 1000` gives
±1.9, `K = 5000` ±0.9, `K = 20000` ±0.4 by the same linearisation, while the error from assuming
normality at all is far larger (below).

The normal assumption does not hold, and the direction of the error is known:

| check | result |
|---|---|
| null's own maximum, in `z` | median **2.82** (5–95%: 2.17–4.11) |
| where a normal sample of 5000 would put it | ≈ **3.54** |
| molecules whose null maximum falls below that | **85.0%** (n = 6,392) |
| generalised Pareto fit to one null (A004, threshold `q0.90`) | shape **ξ = −0.17** |
| — fitted right endpoint vs the observed score | **0.323** vs **0.653** |

The upper tail is **thinner** than normal, with a finite right endpoint. The fitted normal
therefore assigns *more* probability far to the right than the data support, so the reported
`log₁₀ P` is conservative rather than optimistic. Read it as a separation statistic that orders
candidates on a log scale; read it as a calibrated probability only as far as the normal
assumption is granted. (The generalised Pareto fit in the table is a **check** on that assumption,
not the reported method: under it the observed score falls outside the fitted support entirely and
`log₁₀ P` goes to −∞, which is a bound rather than a value.)

The `sv_confidence` tool, by contrast, computes a genuinely empirical p-value with the add-one
correction, `(1 + #{null ≥ obs}) / (1 + n)` (`null_cross.py`), because its observed statistics
fall inside the null's range and counting works there.

---

## 6b. Image quality checks

Two per-molecule checks (`hdom/qc.py`), written to `{base}_qc.json` and turned into a banner on
the molecule's page. Neither looks at the alignment. They ask whether the image could support one,
which the alignment score cannot answer: a molecule whose labels sit barely above their own
background still produces a rank 1.

```
green_snr    median(height of detected green peaks) / median(green profile)     fail < 2.5
red_medroll  per-column backbone width -> std over a rolling ±30-column window
             -> median over columns                                            fail > 0.33
             (None, also a failure, when fewer than 70 columns are bright
              enough to trace: too short or too broken to measure)
```

`green_snr` is a ratio, not `(signal − background) / noise`. `red_medroll` is large when the
stained backbone thickens and thins along its length.

Nothing is filtered by default: the sidecar is written, a warning is printed, and the molecule is
aligned and displayed as usual. `--qc-filter` moves failures to `<output>/qc_failed/` instead.

The thresholds are empirical, set by looking at molecules already judged by eye rather than
derived, so they are a boundary drawn around known-good data. What they catch well is data that
does not resemble it:

| set | molecules | failed |
|---|---|---|
| this study (groups A, AR, B, C+, XY) | 3,079 | **8 (0.3%)** |
| a set imaged with an unoptimised green fluorophore | 39 | **37 (94.9%)** |

Here `green_snr` has median 6.6 and 5th percentile 3.4; there it clusters at 1.2–1.8, i.e. label
peaks only 30–80% brighter than their own background, at which point no height threshold separates
labels from background bumps and label detection over-reports several-fold. `red_medroll` failed
nothing in that set, so the checks also say *which* channel is at fault. On this study's data the
checks change nothing, which is the point of reporting both numbers: they exist to tell a new user
whether the pipeline is being applied to the kind of image it was built for.

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

That `0.1` is inherited rather than derived. It has no measurement behind it, was never tuned,
and is written out verbatim in five files (`make_reference.py`, and three times each in
`DOM_lib.py`, `DOM_lib_Y.py`, `DOM_lib_Q.py` and `hDOM_calc_Q.py`). Changing it in one place
only would put `genome.maps` out of step with the profiles the pcc stage correlates against.

Default motif `CTCGTG` (Nb.BssSI). Two references built with even slightly different parameters
disagree about which homologue a molecule came from, so a rebuild must be checked against the
reference the results were produced with.

---

## 8. What has been measured

The groups the numbers below are broken down by. A molecule is assigned by what its top two
candidates are, not by how well it scored:

| group | meaning |
|---|---|
| **A** | rank 1 and rank 2 are the two haplotypes of one chromosome, and maligner agrees with one of them |
| **B** | the same, but maligner disagrees or returned nothing |
| **XY** | rank 1 is a sex chromosome |
| **AR** | rank 1 and rank 2 are *different* acrocentric chromosomes, and rank 1 falls inside an acrocentric-specific region |

**AR** exists because the short arms of chromosomes 13, 14, 15, 21 and 22 share sequence between
*non-homologous* chromosomes, so a molecule from one of them has no unique locus to be assigned
to. The region is the union of the pseudo-homologous regions (PHRs) described by Guarracino et
al., the acrocentric α-satellite composite, and rDNA:

> Guarracino A, Buonaiuto S, de Lima LG, Potapova T, Rhie A, Koren S, Rubinstein B, Fischer C,
> Gerton JL, Phillippy AM, Colonna V, Garrison E and the Human Pangenome Reference Consortium,
> *Recombination between heterologous human acrocentric chromosomes*, Nature 617(7960):335–343,
> 2023. doi:10.1038/s41586-023-05976-y

PHR coordinates were obtained by aligning that paper's CHM13 sequences to HG002 with minimap2.
The union spans **43.1 Mb, 0.72%** of the 6.0 Gb diploid genome (PHR 0.663%, α-satellite
composite 0.055%, rDNA 0.041%, mostly inside the PHRs). Pan-genomic satellites (HSat1A/1B, bSat,
SST1) are excluded because they also occur on chromosomes 3, 4, 8 and Y and so are not
acrocentric-specific; defining AR by the whole short arm instead would cover 138.6 Mb (2.31%),
which does not match the published PHR extent. As a check, 24 of 3,166 molecules (0.95%) have
rank 1 inside the region against the 0.72% the genome fraction predicts, within Poisson 1.5σ.

The AR assignment itself is made by the classification tooling on the analysis server and is not
part of this package.

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
- **`BBFIX_THRESH = 40` is empirical, and its check set is not independent.** Like the §6b
  thresholds it is a boundary drawn around molecules already judged by eye, and the 30-molecule
  comparison in §2 that shows it separating help from harm is drawn from the same kind of case it
  was set on. What that comparison does establish is the direction: where the test fires, no
  molecule's `F1` fell. A held-out set was not used.
- **Detection threshold is fixed at 0.5.** On a molecule whose green background sits at 28% of
  its peak, ten candidates at 1.3–1.9× background pass at probabilities 0.50–0.75 while the
  twelve real labels score 0.99–1.00. The model separates them; the threshold does not.
- **`log₁₀ P` comes from a fitted normal, not from counting.** §6 gives the formula, the log
  base, and the measured direction of the error. The tail is not itself fitted, so the phrase
  "tail extrapolation" (which names a different, better-supported technique) does not apply. The
  same `norm.logsf` conversion is
  applied to three other score columns, and in all three the observation lies *inside* the
  population its own mean and sd are taken from: `pcc` (`hDOM_calc_Q.py`) and `local`/`mDOM` when
  `--logp` is not given (`null_logp` in `local_align_pipeline.py`, calibrated on that molecule's
  own three-to-six candidate scores). Those three order candidates within one molecule; they are
  not calibrated probabilities. Only `--logp` puts the observation outside its null.
- **The `pcc` column's null is neither decoy-only nor a sample of wrong positions.** Unlike
  `--logp`, `hDOM_calc_Q.py` excludes no chromosome: the matched chromosome and its homologue are
  in the population, and the score being converted is that population's maximum. Each element is
  one `(reference, scale)` pair already maximised over every shift and both directions, so the
  distribution is one of per-chromosome maxima rather than of wrong answers. Removing the matched
  chromosome and recomputing, over 20 runs (10 molecules, raw and straightened), inflates `sd` in
  18 of 20 (median 1.08×, max 1.94×), i.e. self-inclusion pulls `log₁₀ P` toward zero and is
  conservative — but by an amount that grows with match quality:

  | molecule | rank 1 | as computed | matched chr removed | `sd` ratio |
  |---|---|---|---|---|
  | `0316…mol001_135kb` | chr01_m | −23.73 | −77.52 | 1.82× |
  | `0327…mol006` | chr08_m | −6.60 | −23.53 | 1.94× |
  | `0315…mol003_176kb` | chr06_m | −7.71 | −12.09 | 1.26× |
  | `0128_Re5_37_38_2` | chr09_m | −3.34 | −3.62 | 1.04× |

  Confident molecules inflate their own `sd` and dock themselves several-fold in log units while
  weak ones near −3 barely move, so the column is **not comparable across molecules**. Those runs
  used 17–30 references, where the matched chromosome is 2–14% of the population; a full
  94-haplotype run dilutes it to about 2% and shrinks the effect without changing its direction.
  `cc_rg2 = cc_r · cc_g²` compounds this, `cc_g` being a cosine taken without mean subtraction.
- **A GPU is required.** The pcc scan runs through CuPy and exits without it.
- **One sample, one enzyme, one instrument.** Everything here is HG002 with Nb.BssSI at
  200 bp/pixel.

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
