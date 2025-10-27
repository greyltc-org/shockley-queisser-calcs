#!/usr/bin/env python3

# example usage: scale AM1.5 from 400 to 500nm by 0.5 and from 600 to 700nm by 2.0:
# cat ASTMG173.csv | ./spectralmod.py -b 600 -e 700 -f 2 -t -j | ./spectralmod.py -b 400 -e 500 -f 0.5 -t -j > ASTMG173_scaled.tsv

import argparse
import fileinput
import sys

parser = argparse.ArgumentParser(
    description="A tool to manipulate parts of a spectrum. Reads from sdtin, writes to stdout",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument(
    "-i", "--input", type=str, default="-", help="Input file (defaults to stdin)"
)
parser.add_argument(
    "-b", "--begin", type=float, required=True, help="Scale range start (inclusive)"
)
parser.add_argument(
    "-e", "--end", type=float, required=True, help="Scale range end (inclusive)"
)
parser.add_argument("-f", "--scale-factor", type=float, default=1, help="Scale factor")
parser.add_argument("-o", "--offset", type=float, default=0, help="Offset")
parser.add_argument(
    "-t",
    "--tsv",
    default=False,
    action="store_true",
    help="Output in tsv instead of csv",
)
parser.add_argument(
    "-j",
    "--just",
    default=False,
    action="store_true",
    help="Just output the rescaled column",
)
parser.add_argument(
    "--col", default=2, type=int, help="Column of interest in input file"
)

args = parser.parse_args()

for line in fileinput.input(files=(args.input,), encoding="utf-8"):
    try:
        if "\t" in line:
            dlm = "\t"
        else:
            dlm = ","
        splitted = line.strip().split(dlm)
        wl = float(splitted[0])
        if len(splitted) == 2:
            val = float(splitted[1])
        else:
            val = float(splitted[args.col])
        if (wl >= args.begin) and (wl <= args.end):
            val = val * args.scale_factor
            val = val + args.offset
            splitted[args.col] = f"{val:.4E}"
        if args.tsv:
            dlm = "\t"
        if args.just:
            splitted = (splitted[0], f"{val:.4E}")
        unsplit = dlm.join(splitted)
        print(unsplit)
    except Exception as e:
        print(e, file=sys.stderr)
        print(line, end="")
