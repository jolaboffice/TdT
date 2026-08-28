#!/usr/bin/env python
"""
hDOM_mal.py
Generate maps files from TIFF images and run maligner_dp alignment.
"""

import os
import sys
import argparse
import subprocess
import shutil
import numpy as np
import tifffile as tiff
from scipy import signal
import DOM_lib_Y as lib

# --- maligner-only green peak refinement (Y) ---
# PCC uses the continuous green profile (peaks discarded), so it is unaffected.
# maligner turns peaks into discrete labels, so noise peaks become false labels.
# Detect on the RAW green profile with prominence (reject shoulders) + width
# (reject sharp spikes) + distance (resolution floor). height = absolute intensity floor.
PEAK_PROMINENCE_FRAC = 0.10   # min prominence as fraction of g.max
PEAK_MIN_WIDTH_PX = 2         # reject sub-PSF-width spikes
PEAK_MIN_DISTANCE_PX = 4      # resolution floor


def refined_green_peaks(data_file):
    """Refined nick-label peaks for maligner (prominence + width + distance)."""
    green_raw = np.max(tiff.imread(data_file)[1], axis=0).astype(float)
    thres = green_raw.mean() + (green_raw.max() - green_raw.mean()) * 0.1
    peaks, _ = signal.find_peaks(
        green_raw,
        height=thres,
        prominence=PEAK_PROMINENCE_FRAC * green_raw.max(),
        width=PEAK_MIN_WIDTH_PX,
        distance=PEAK_MIN_DISTANCE_PX,
    )
    return peaks

# Default values used across all DOM scripts
DEFAULT_MAPS_FOLDER = 'ref_map'
DEFAULT_BPP = 200
DEFAULT_SHIFT_WINDOW = 50
DEFAULT_START_BPP = 181
DEFAULT_END_BPP = 240
DEFAULT_MALIGNER_PATH = 'maligner_dp'
DEFAULT_SEQ_FOR_R = 'A'
DEFAULT_SEQ_FOR_G = 'CTTAAG'


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate maps from TIFF images and run maligner_dp alignment'
    )
    parser.add_argument('data_file', help='Input TIFF data file')
    parser.add_argument('--ref-map', dest='ref_map', default=None,
                        help='Reference maps file (default: ref_map/genome.maps, fallback: ref_map/huma.maps)')
    parser.add_argument('--bpp', type=int, default=DEFAULT_BPP,
                        help=f'Base pairs per pixel (default: {DEFAULT_BPP})')
    parser.add_argument('--maligner-path', dest='maligner_path', 
                        default=DEFAULT_MALIGNER_PATH,
                        help=f'Path to maligner_dp executable (default: {DEFAULT_MALIGNER_PATH})')
    return parser.parse_args()


def resolve_maligner_path(path_hint):
    """Resolve maligner_dp executable path."""
    # First try the hint as-is (absolute or relative path)
    if os.path.isfile(path_hint):
        return path_hint
    
    # Try in the script directory (maligner/build/bin/maligner_dp)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_path = os.path.join(script_dir, 'maligner', 'build', 'bin', 'maligner_dp')
    if os.path.isfile(local_path):
        return local_path
    
    # Try which/shutil.which (for PATH lookup)
    mal_path = shutil.which(path_hint)
    if mal_path and os.path.isfile(mal_path):
        return mal_path
    
    print(f"Error: maligner_dp executable '{path_hint}' not found", file=sys.stderr)
    print(f"  Tried: {path_hint}")
    print(f"  Tried: {local_path}")
    sys.exit(1)


def create_maps_from_tif(data_file, bpp):
    """Create .maps file from TIFF image."""
    data_dir = os.path.dirname(data_file)
    data_name = os.path.splitext(os.path.basename(data_file))[0]
    maps_file = os.path.join(data_dir, f"{data_name}_{bpp}.maps")
    
    # Get image dimensions
    _, w2, _ = lib.tif2RGimage(data_file)
    # Refined peaks (maligner-only fix); PCC path still uses tif2RGprofile.
    peaks = refined_green_peaks(data_file)

    # Convert peaks to bp units and create maps file
    lib.peaks2maps(data_file, w2 * bpp, peaks * bpp, maps_file)
    
    return maps_file


def run_maligner_dp(data_maps, ref_maps, mal_path):
    """Run maligner_dp and return output."""
    if not os.path.exists(data_maps):
        print(f"Error: Data maps file not found: {data_maps}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(ref_maps):
        print(f"Error: Reference maps file not found: {ref_maps}", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Moderately expanded parameters — enough candidate diversity for the
        # downstream cc_r-aware re-ranking to recover mal-suboptimal hits,
        # but small enough that the combine step's per-candidate cc_r/cc_g
        # computation completes in reasonable time. Maligner defaults:
        # --max-alignments 10, --max-alignments-per-reference 100,
        # --query-max-misses 2, --ref-max-misses 5, scaling 0.85..1.15.
        mal_cmd = [
            mal_path,
            "--max-alignments", "10",
            "--max-alignments-per-reference", "3",
            "--query-max-misses", "3",
            "--ref-max-misses", "6",
            "--min-query-rescaling", "0.82",
            "--max-query-rescaling", "1.18",
            data_maps, ref_maps,
        ]
        result = subprocess.run(
            mal_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"Warning: maligner_dp returned non-zero exit code ({result.returncode})", file=sys.stderr)
            if result.stderr:
                print(f"Error message: {result.stderr}", file=sys.stderr)
        
        return result.stdout
        
    except FileNotFoundError:
        print(f"Error: maligner_dp executable not found: {mal_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error running maligner_dp: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    args = parse_args()
    
    # Resolve paths
    data_file = os.path.abspath(args.data_file)
    if args.ref_map:
        ref_maps = os.path.abspath(args.ref_map) if not os.path.isabs(args.ref_map) else args.ref_map
    else:
        # Prefer ref_map folder inside data folder, then next to data folder, then fallback to script ref_map
        data_dir = os.path.dirname(data_file)
        parent_dir = os.path.dirname(data_dir)
        ref_map_dir = None
        for candidate in [
            os.path.join(data_dir, 'ref_map'),
            os.path.join(parent_dir, 'ref_map'),
        ]:
            if os.path.isdir(candidate):
                ref_map_dir = candidate
                break
        if ref_map_dir is None:
            ref_map_dir = os.path.join(parent_dir, 'ref_map')
        primary_ref = os.path.join(ref_map_dir, 'genome.maps')
        fallback_ref = os.path.join(ref_map_dir, 'huma.maps')
        if not os.path.exists(primary_ref) and not os.path.exists(fallback_ref):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            primary_ref = os.path.join(script_dir, 'ref_map', 'genome.maps')
            fallback_ref = os.path.join(script_dir, 'ref_map', 'huma.maps')
        ref_maps = primary_ref if os.path.exists(primary_ref) else fallback_ref
    
    # Check if data file exists
    if not os.path.exists(data_file):
        print(f"Error: Data file not found: {data_file}", file=sys.stderr)
        sys.exit(1)
    
    # Check if reference maps exists
    if not os.path.exists(ref_maps):
        print(f"Error: Reference maps file not found: {ref_maps}", file=sys.stderr)
        sys.exit(1)
    
    # Resolve maligner path
    mal_path = resolve_maligner_path(args.maligner_path)
    
    print(f"Data file: {data_file}")
    print(f"Reference maps: {ref_maps}")
    print(f"BPP: {args.bpp}")
    print(f"Maligner path: {mal_path}\n")
    
    # Create maps file from TIFF
    print("Creating maps file from TIFF image...")
    data_maps = create_maps_from_tif(data_file, args.bpp)
    print(f"Created maps file: {data_maps}\n")
    
    # Run maligner_dp
    print("Running maligner_dp...")
    output = run_maligner_dp(data_maps, ref_maps, mal_path)

    # Save output to .txt file next to the data file
    txt_path = data_file.replace(".tif", "_maligner_dp.txt")
    try:
        with open(txt_path, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Saved maligner_dp output to: {txt_path}\n")
    except Exception as e:
        print(f"Warning: Failed to save maligner_dp output to {txt_path}: {e}", file=sys.stderr)
    


if __name__ == '__main__':
    main()
