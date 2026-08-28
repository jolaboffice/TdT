#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SEM DNA length measurement — interactive GUI launcher.

Usage:
    python3 run_interactive.py <image.tif> [--pixel-nm 2.7] [--stretch 100]
    python3 run_interactive.py QD066/QD066.tif
"""

import argparse
from dna_length_measurement import InteractiveDNAMeasurer

DEFAULT_PIXEL_TO_NM = 2.7   # px -> nm at the acquisition magnification
DEFAULT_DNA_STRETCH = 100   # assumed DNA stretch %, 0-100


def main():
    p = argparse.ArgumentParser(description='SEM DNA path tracing / length measurement')
    p.add_argument('image_path', help='Path to .tif image (e.g. QD066/QD066.tif)')
    p.add_argument('--pixel-nm', type=float, default=DEFAULT_PIXEL_TO_NM,
                   help='pixel size in nm (default %(default)s)')
    p.add_argument('--stretch', type=float, default=DEFAULT_DNA_STRETCH,
                   help='DNA stretch percent 0-100 (default %(default)s)')
    args = p.parse_args()

    measurer = InteractiveDNAMeasurer(args.image_path)
    measurer.set_scale(pixel_to_nm=args.pixel_nm)
    measurer.set_dna_stretch(args.stretch)
    measurer.run_interactive()


if __name__ == "__main__":
    main()
