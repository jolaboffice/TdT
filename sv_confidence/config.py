"""Central path configuration. Every path defaults to its on-disk location and can
be overridden with an environment variable.
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
# project root (parent of deploy/)
PROJECT_ROOT = os.environ.get('SV_PROJECT_ROOT', os.path.dirname(_HERE))


# Repository root, two levels up from this file, so the defaults below work wherever the
# checkout happens to sit.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _env(key, default):
    return os.environ.get(key, default)


# --- Shared analysis code (must be importable) -----------------------------
# local_align_pipeline / mdom_reanchor / DOM_lib_Y
# The scoring engine lives in this repository, one level up: sv_confidence imports
# local_align_pipeline, mdom_reanchor and DOM_lib_Y from the mapping pipeline's package
# rather than carrying its own copy, so both lines score with the same code.
LA_CODE = _env('SV_LA_CODE', os.path.join(REPO_ROOT, 'hdom', 'hdom'))
# homologue-comparison helpers
# Haplotype PAF liftover; the mapping pipeline ships it as paf_liftover.
HOMOLOG_CODE = _env('SV_HOMOLOG_CODE', os.path.join(REPO_ROOT, 'hdom', 'hdom'))
# ML nick-label detector (Label_Detect_ML_based)
ML_CODE = _env('SV_ML_CODE', os.path.join(REPO_ROOT, 'hdom', 'hdom'))

# --- Input data (molecule images from the hDOM pipeline; not shipped here) --
# molecule tifs + matching_results.csv; analyze_inversion default source ({ROOT}/<group>)
# Molecule images and the pool they are drawn from: these are data, not code, and have no
# in-repository default. Point them at your own hDOM output before running anything.
ROOT = _env('SV_ROOT', '')
# reference-consistent molecule pool (null bank source)
POOL_DIR = _env('SV_POOL_DIR', '')
# translocation candidates: left/right fragments aligned separately
FIG4 = _env('SV_FIG4', '')

# --- Data artifacts inside the project -------------------------------------
# prebuilt molecule pool cache; reused as-is (bit-exact nulls, do not rebuild)
POOL_CACHE = _env('SV_POOL_CACHE', os.path.join(PROJECT_ROOT, '.pool_cache.pkl'))
# straightened working dir: crop-source tifs + per-molecule svla.json dumps
SVLA_TIFDIR = _env('SVLA_TIFDIR', os.path.join(PROJECT_ROOT, 'straighten_work', 'out'))
SVLA_BASE = _env('SVLA_BASE', SVLA_TIFDIR)
# per-molecule inversion results
FINAL_INV = _env('SV_FINAL_INV', os.path.join(PROJECT_ROOT, 'straighten_work', 'final_inv'))
# length -> null(red_d) calibration curve for the width score
CALIB_LENGTHNULL = _env('SV_CALIB_LENGTHNULL', os.path.join(PROJECT_ROOT, 'calib_lengthnull.json'))
# matched null-of-T cache (per marked-width key); reused as-is (bit-exact)
TNULL_CACHE_DIR = _env('SV_TNULL_CACHE', os.path.join(PROJECT_ROOT, '.tnull_cache'))
# translocation null json directory
NULL_DIR = _env('SV_NULL_DIR', PROJECT_ROOT)


def ensure_sys_path():
    """Prepend the shared-code dirs to sys.path (idempotent)."""
    import sys
    for p in (ML_CODE, HOMOLOG_CODE, LA_CODE):
        if p and p not in sys.path:
            sys.path.insert(0, p)
