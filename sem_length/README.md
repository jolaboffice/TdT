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

## Requirements

Python 3.9 with the packages in `requirements.txt` (`pip install -r requirements.txt`):
numpy, scipy, scikit-image, matplotlib, tifffile, pillow.

- The GUI tools (`pick_points.py`, `run_interactive.py`) need a display and an
  interactive matplotlib backend (Tk / `tkinter`). The measurement core runs
  headless with `matplotlib.use("Agg")`.
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
maligner_dp <QD###>/<QD###>.maps ../genome_smoothed_rm100.maps > <QD###>/<QD###>_mal.txt

# 5. Smooth ref (or query) maps to a lower resolution for FM<->SEM comparison
#    (uses the bundled malignpy; --dry-run to preview)
python3 smooth_maps_batch.py /path/to/folder -m 500 --dry-run
```

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
