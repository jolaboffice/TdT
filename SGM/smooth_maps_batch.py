#!/usr/bin/env python
"""
Batch smooth all .maps files in a directory (including subdirectories).
Uses smooth_left_right from maligner (merge smallest fragment with shorter neighbor).
Creates new files with _sm{threshold}.maps suffix. Original files are not modified.

Usage:
    smooth_maps_batch.py /path/to/folder [-m 500] [--dry-run]
    smooth_maps_batch.py /path/to/folder -m 1000
"""

import argparse
import os
import sys
import glob

# maligner python lib: import if already on PYTHONPATH, else prepend MALIGNER_PY entries
# (os.pathsep-separated, like PYTHONPATH; malignpy needs its lib dir + core/common subdirs).
# The maligner Python library ships in this repository under maligner/, so the default needs
# no environment variable; MALIGNER_PY still wins if someone points it at another build.
_BUNDLED = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'maligner', 'build', 'bin')
for _d in (os.environ.get('MALIGNER_PY', '').split(os.pathsep) + [_BUNDLED]):
    if _d:
        sys.path.insert(0, os.path.expanduser(_d))
try:
    from malignpy.maps.MalignerMap import smooth_left_right
    from malignpy.maps.utils import gen_maps
except ImportError:
    sys.exit("Error: cannot import malignpy. Put the maligner lib (and its core/common "
             "subdirs) on PYTHONPATH, or set MALIGNER_PY (colon-separated).")


def smooth_file(maps_path, min_frag, dry_run=False):
    """Smooth a single maps file and save as new file with _sm{threshold}.maps suffix."""
    base, ext = os.path.splitext(maps_path)
    out_path = f"{base}_sm{min_frag}{ext}"

    results = []
    smoothed_maps = []
    with open(maps_path, 'r') as f:
        for m in gen_maps(f):
            original_n = len(m.frags)
            smoothed = smooth_left_right(m, min_frag, boundary_frags=True)
            new_n = len(smoothed.frags)
            merged = original_n - new_n
            results.append((m.mapId, original_n, new_n, merged))
            smoothed_maps.append(smoothed)

    if not dry_run and smoothed_maps:
        with open(out_path, 'w') as out:
            for sm in smoothed_maps:
                sm.write(out)

    return results, out_path


def main():
    parser = argparse.ArgumentParser(
        description='Batch smooth all .maps files in a directory tree (creates new files)')
    parser.add_argument('folder', help='Target folder (searches recursively)')
    parser.add_argument('-m', '--min-frag', type=int, default=500,
                        help='Minimum fragment size in bp (default: 500)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be changed without creating files')
    args = parser.parse_args()

    # Find all .maps files but exclude already-smoothed ones (_sm*.maps)
    all_maps = sorted(glob.glob(os.path.join(args.folder, '**', '*.maps'), recursive=True))
    maps_files = [f for f in all_maps if '_sm' not in os.path.basename(f)]
    if not maps_files:
        print(f'No .maps files found in {args.folder}')
        return

    print(f'Found {len(maps_files)} .maps files (min_frag={args.min_frag}bp)')
    if args.dry_run:
        print('[DRY RUN - no files will be created]')
    print()

    total_merged = 0
    files_changed = 0

    for maps_path in maps_files:
        rel_path = os.path.relpath(maps_path, args.folder)
        results, out_path = smooth_file(maps_path, args.min_frag, args.dry_run)
        for map_id, orig_n, new_n, merged in results:
            if merged > 0:
                total_merged += merged
                files_changed += 1
                rel_out = os.path.relpath(out_path, args.folder)
                print(f'  {rel_path}: {orig_n} → {new_n} frags ({merged} merged) → {rel_out}')

    print()
    print(f'Total: {files_changed}/{len(maps_files)} files changed, {total_merged} fragments merged')


if __name__ == '__main__':
    main()
