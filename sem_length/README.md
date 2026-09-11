# SEM single-molecule DNA length measurement & alignment

Tools for turning negative-stain SEM image tiles of single DNA molecules into
nick-fragment size profiles and aligning them to a reference genome. Each
molecule yields a `.maps` fragment profile that `maligner_dp` aligns against
reference maps; reference maps can be smoothed to a lower resolution to compare
SEM against fluorescence-microscopy (FM) maps.

## Pipeline

```
 SEM tiles ─▶ pick_points.py ─▶ config_<folder>.json
                                     │
                                     ▼
                              stitching.py ─▶ <QD###>.tif        (stitched image)
                                     │
                                     ▼
              run_interactive.py (dna_length_measurement.py)
                 trace backbone + measure segments
                                     │
                    ┌────────────────┼─────────────────┐
                    ▼                ▼                 ▼
              <QD###>.csv     <QD###>.maps      <QD###>_points.json
              (nm/bp/px)      (fragment bp)     (click points)
                                     │
                                     ▼
                        maligner_dp  (DP alignment)
                                     │
                                     ▼
                              <QD###>_mal.txt   (ranked alignments)
                                     │
                        maligner_multirank_strip.R
                                     │
                                     ▼
                          <QD###>_rank<n>_strip.png
                                     │
                                     ▼
                   sem_view.py ─▶ <QD###>_view_ytb.html   (one page per molecule,
                                   <QD###>_matching_results.csv   served by index.php?file=)
                                     │
                                     ▼
                 sem_trace_panel.py ─▶ <QD###>_sem_trace.png
                   redraw the SEM panel: raw photo + lifted trace + QD arrows,
                   and resample the query / rank strips onto that photo's axis

 smooth_maps_batch.py:  genome_*.maps ─▶ genome_*_sm<N>.maps
   merge sub-threshold fragments to match FM resolution (for FM↔SEM comparison)
```

## Files

| File | Role |
|------|------|
| `pick_points.py` | GUI: click matching overlap points on consecutive SEM tiles; writes `config_<folder>.json`. |
| `stitching.py` | Gradient-NCC stitching of the tiles into one image, using the config. |
| `dna_length_measurement.py` | Core: DNA backbone tracing (least-cost path in a corridor) + per-segment length (nm/bp); `DNALengthMeasurer` (headless) and `InteractiveDNAMeasurer` (matplotlib GUI). Exports `.csv` / `.maps` / `_points.json`. |
| `run_interactive.py` | Launcher for the measurement GUI on one image. |
| `smooth_maps_batch.py` | Smooths `.maps` (merges fragments below a bp threshold via maligner's `smooth_left_right`) into `*_sm<N>.maps`; originals untouched. |
| `maligner_multirank_strip.R` | Draws the per-molecule rank figures: query and reference label strips on a shared bp axis, one panel per rank, with the fragment-size table. Needs R with `ggplot2`. |
| `sem_view.py` | Re-frames what the steps above already produced into the hDOM align-page layout: one `<QD###>_view_ytb.html` per molecule plus a `<QD###>_matching_results.csv` in the hDOM schema, so `index.php?file=<QD###>` renders the same table, navigation and marking UI as a fluorescence molecule. Nothing is re-aligned. Borrows the page style from `../hdom/hdom` (override with `HDOM_DIR`). |
| `sem_trace_panel.py` | Redraws one molecule's SEM panel as the raw photo with the trace and QD marks lifted clear of the DNA, and resamples the query / rank strips onto that photo's axis so the labels line up with the arrows. Run it per molecule after `sem_view.py`. |

## Viewing the result

`sem_view.py` writes plain files; serving them needs a web server with PHP and an `index.php`
that renders `?file=<molecule>` (pass one with `--template-index`). The hDOM distribution's
viewer is the one these pages are written for.

```
python3 sem_view.py --data <data dir> --out <view dir> --template-index <index.php>
for m in <view dir>/*_matching_results.csv; do
    python3 sem_trace_panel.py "$(basename "${m%_matching_results.csv}")" \
        --data <data dir> --out <view dir>
done
```

`sem_trace_panel.py` keeps a copy of each view as `<QD###>_view_ytb.html.bak_sem` on its first
run and rebuilds from that copy afterwards, so re-running it with different options is safe.

## Requirements

Python 3.9 with the packages in `requirements.txt` (`pip install -r requirements.txt`):
numpy, scipy, scikit-image, matplotlib, tifffile, pillow.

- The GUI tools (`pick_points.py`, `run_interactive.py`) need a display and an
  interactive matplotlib backend (Tk / `tkinter`). The measurement core runs
  headless with `matplotlib.use("Agg")`.
- **R with `ggplot2`** — only for `maligner_multirank_strip.R`, the rank figures (step 6).
  Nothing else in this folder needs R.
- **maligner** — the optical-map aligner, bundled in this repository under `../maligner`
  (not on PyPI). A Linux build is included, so there is normally nothing to do:
  - put the binary on `$PATH`: `export PATH="$PWD/../maligner/build/bin:$PATH"`;
  - `smooth_maps_batch.py` finds the bundled `malignpy` by itself.

  On other platforms, rebuild from `../maligner`'s own README.

## Configuration (paths)

No absolute paths are hard-coded. Reference maps (`genome_*.maps`) and the
per-molecule `data/<QD###>/` folders are yours and are passed to every tool as explicit path
arguments. `maligner_dp` is invoked as a separate program, so it has to be on `$PATH`.

| What | How it is located | Override |
|------|-------------------|----------|
| `maligner_dp` binary | `$PATH` | — |
| `malignpy` package (`smooth_maps_batch.py`) | `../maligner/build/bin`, resolved from this file | `MALIGNER_PY=dir1:dir2:…` (searched first) |

`smooth_maps_batch.py` falls back to the bundled copy after anything in `MALIGNER_PY`, and exits
with a clear message if neither imports.

## Usage

```bash
# 1. Pick tile overlap points (GUI) -> config_<folder>.json
python3 pick_points.py /path/to/tiles/

# 2. Stitch the tiles
python3 stitching.py /path/to/tiles/config_<folder>.json

# 3. Trace + measure one molecule (GUI); --pixel-nm / --stretch optional
python3 run_interactive.py <QD###>/<QD###>.tif --pixel-nm 2.7 --stretch 100

# 4. Align each molecule's .maps to a reference (maligner_dp from ../maligner, on $PATH)
#    The reference maps are deposited ready-built (see "Reference maps" below); build your
#    own only for a different assembly or enzyme:
#      make_insilico_map /path/to/chr01_p.fasta CTCGTG -o genome_unsmoothed
#      python3 smooth_maps_batch.py /path/to/maps -m 100
maligner_dp <QD###>/<QD###>.maps ../genome_smoothed_rm100.maps > <QD###>/<QD###>_mal.txt

# 5. Smooth ref (or query) maps to a lower resolution for FM<->SEM comparison
#    (uses the bundled malignpy; --dry-run to preview)
python3 smooth_maps_batch.py /path/to/folder -m 500 --dry-run

# 6. Rank figures: query vs smoothed/unsmoothed reference strips per molecule
Rscript maligner_multirank_strip.R /path/to/data genome_smoothed_rm100.maps genome_unsmoothed.maps
```

## Reference maps

The in silico Nb.BssSI maps for HG002 v1.1 used in the paper are `SGM_reference_maps_HG002v11`,
one file per smoothing threshold. They are deposited as data rather than kept in this repository,
so step 4 needs no `make_insilico_map` run of your own:

| file | fragments below … merged | for |
|---|---|---|
| `genome_unsmoothed.maps` | — | every predicted site |
| `genome_smoothed_rm50.maps` | 50 bp | |
| `genome_smoothed_rm100.maps` | 100 bp | **the map the paper aligns against** |
| `genome_smoothed_rm150.maps` | 150 bp | |
| `genome_FM_sm200.maps` | 200 bp | fluorescence-microscopy resolution, for FM↔SEM |
| `genome_FM_sm600.maps` | 600 bp | fluorescence-microscopy resolution, for FM↔SEM |

Each holds 46 sequences (chr1–22 in both haplotypes, plus `ChrX_m` and `ChrY_p`); `genome_FM_sm200.maps`
additionally carries `ChrM`, for 47.

**rm100** is the one used for alignment: two nick sites closer than 100 bp (about 34 nm) are a
single spot in the SEM image, so merging them in the reference keeps them from being counted as
false negatives. Genome-wide this merges 2.2% of fragments (808,027 → 789,874).

## Data layout

- Reference: `genome_*.maps` — one per smoothing threshold (`unsmoothed`,
  `smoothed_rm{50,100,150}`, `FM_sm{200,600}`).
- Per molecule: `data/<QD###>/` with the SEM tiles, `config_<QD###>.json`, the
  stitched `<QD###>.tif`, `<QD###>.csv`, `<QD###>.maps`, `<QD###>_points.json`,
  and `<QD###>_mal.txt`.
- `.maps` format (maligner): `id <TAB> total_bp <TAB> n_frags <TAB> frag1 frag2 …`.

## Notes

- The measurement tool is tracing/measurement only. Its optional maligner FP/FN
  marker overlay was removed from this package (it depended on a separate
  `maligner_fpfn` module outside this scope).
- Line endings are normalised to LF.
