#!/usr/bin/env python3

# example usage: scale AM1.5 from 400 to 500nm by 0.5 and from 600 to 700nm by 2.0:
# cat ASTMG173.csv | ./spectralmod.py -b 600 -e 700 -f 2 | ./spectralmod.py -b 400 -e 500 -f 0.5 > ASTMG173_scaled.csv

import argparse
import fileinput
import sys

parser = argparse.ArgumentParser(description='Scale parts of a spectrum. Reads from sdtin, writes to stdout', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i", "--input", type=str, default='-', help="Input file (defaults to stdin)")
parser.add_argument("-b", "--begin", type=float, required=True, help="Scale range start (inclusive)")
parser.add_argument("-e", "--end", type=float, required=True, help="Scale range end (inclusive)")
parser.add_argument("-f", "--scale-factor", type=float, default=1, help="Scale factor")
parser.add_argument("--col", default=2, type=int, help="Column of interest in input file")

args = parser.parse_args()

for line in fileinput.input(files=(args.input,), encoding="utf-8"):
    try:
        splitted = line.strip().split(',')
        wl = float(splitted[0])
        val = float(splitted[args.col])
        if (wl >= args.begin) and (wl <= args.end):
            val = val * args.scale_factor
            splitted[args.col] = f"{val:.4E}"
        unsplit = ','.join(splitted)
        print(unsplit)
    except Exception as e:
        print(e, file=sys.stderr)
        print(line, end="")
