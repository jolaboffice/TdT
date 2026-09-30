# Single-molecule SV confidence

Statistical-confidence scoring for structural variants (SVs) called on single
straightened hDOM molecules. Two independent tests:

- **Inversion** (pure red / A-T backbone channel). The breakpoint is localized by a
  fully-**joint hybrid** search: for each per-window **local stretch `gamma`** (correcting local
  contraction/stretch vs the reference beyond the single global scale), green re-localizes the
  center (label anchor) and red refines the symmetric width by a length-normalized score `z`, both
  scanned over the **±10 kb (ENDPOINT_R) AND label-preserving** band — the statistic is
  `T = max-z` over `(center, width, gamma)`. An inversion is supported when
  `red_d = corr(reversed) − corr(forward) > 0` and the empirical `p_T ≤ 0.05`, where `p_T` tests
  `T` against a **matched null-of-T** (null windows run through the identical full
  `(center, width, gamma)` ±10 kb hybrid, so the localization search and the extra stretch DoF are
  matched on both sides; Kulldorff scan-statistic). The verdict uses only the backbone channel; the
  sparse-label channel nominates and localizes the candidate but is excluded from the verdict to
  avoid circular analysis. Inversion self-FP ≈ 4.7–5.0% (calibrated at α).
- **Translocation / junction** (dual channel, intersection-union test). At a
  junction localized within ±10 kb, two statistics contrast each flank fitting
  its own locus against the opposite flank's alignment continued across the
  junction:
  - `cross_red  = Σ_half (own global_ccr − alt global_ccr)` — backbone Pearson correlation.
  - `cross_green = Σ_half (own galign − alt galign)` — effective matched-label count.
  A junction is supported when `max(p_red, p_green) ≤ 0.05`, both statistics are
  positive, and the two loci are distinct.

## Layout

| File | Role |
|------|------|
| `config.py` | Central paths; all overridable by environment variables. |
| `sv_la.py` | Core: molecule placement, region/whole scoring, inversion scan (green-argmax fallback), `analyze_inversion` (hybrid localization + matched null-of-T verdict). |
| `hybrid_localize.py` | Hybrid inversion localization (`hybrid_scan_A`), length-fair width score, and the matched null-of-T (`null_bank_T`). |
| `junction_search.py` | Junction localization over mark ±10 kb (whole-piece `global_ccr`). |
| `confidence.py` | Translocation cross statistics (`cross_red`, `cross_green`). |
| `null_cross.py` | Translocation null distribution and self false-positive rate. |
| `sv_la_transloc.py` | Sub-segment aligned-bp helper (`_cont_aligned_bp`). |
| `run_transloc_fixed.py` | Green label-overlap statistic (`_green_overlap`). |

## Requirements

Python 3.9 with the pinned packages in `requirements.txt`
(`pip install -r requirements.txt`). The shared analysis modules
(`local_align_pipeline`, `mdom_reanchor`, `DOM_lib_Y`, `Label_Detect_ML_based`)
are included with this code set and must be importable — put their directories on
`PYTHONPATH`, or point `SV_LA_CODE` / `SV_ML_CODE` / `SV_HOMOLOG_CODE` at them
(`config.py` prepends those to `sys.path` and errors clearly if the import fails).

## Paths / configuration

`config.py` centralises every path and lets each be overridden with an environment
variable. Three kinds:

- **Shared code** — `SV_LA_CODE`, `SV_ML_CODE`, `SV_HOMOLOG_CODE`: default to `../hdom`
  in this repository, so nothing needs setting. `SV_PROJECT_ROOT` (default: parent of
  this folder) roots the artifacts below.
- **Input data** — `SV_ROOT`, `SV_POOL_DIR`, `SV_FIG4`: the molecule images the analysis runs
  on, which are output of the mapping pipeline in `..`. **These have no default and must be
  set** — they are your own data, and a guessed path would either fail late or read the wrong
  images. Nothing here ships with the package:
  - `SV_ROOT` — molecule tifs + `matching_results.csv`; `analyze_inversion` default source (base_dir overrides).
  - `SV_POOL_DIR` — the reference-consistent normal-molecule pool that seeds the null bank.
  - `SV_FIG4` — translocation candidates: left/right fragments aligned separately.
- **Artifacts** (default under `SV_PROJECT_ROOT`, shipped with this package) —
  `SV_POOL_CACHE`, `SV_TNULL_CACHE`, `SV_FINAL_INV`, `SV_CALIB_LENGTHNULL`, `SVLA_TIFDIR`,
  `SV_NULL_DIR`. `SV_POOL_CACHE` is reused as-is so the nulls reproduce
  bit-for-bit (do not rebuild it from `SV_POOL_DIR`). Figure reproduction runs off these
  caches — the input-data paths are only needed to re-run the analysis from the images
  (and `SV_FIG4` for the translocation figures).

## Reproduce

Run from this folder (or set `SV_PROJECT_ROOT`).

```bash
# Translocation null distribution + self false-positive rate.
python null_cross.py 1000                      # -> transloc_iut_null_cross.json

# Translocation cross scores (5 reference cases).
python confidence.py
```

Programmatic entry points:

- Inversion: `sv_la.analyze_inversion(stem, a_bp, b_bp, group=..., dump_dir=...)`.
- Junction search: `python junction_search.py <stem> <Lchrom> <Lbp> <Ldir> <Lscale> <Rchrom> <Rbp> <Rdir> <Rscale> <mark_bp>`.

## Reproducibility

All randomness is seeded (the matched null-of-T is keyed on marked width and its
per-key seed is deterministic; the translocation null uses `RandomState(3)`), so
re-running yields identical numbers. The prebuilt pool cache (`SV_POOL_CACHE`), the
length calibration curve (`SV_CALIB_LENGTHNULL`), the null-of-T cache
(`SV_TNULL_CACHE`), and the stored per-molecule null-of-T (`diag.null_T` inside each
`*_svla.json`) are the fixed inputs that keep the reported p-values stable.
