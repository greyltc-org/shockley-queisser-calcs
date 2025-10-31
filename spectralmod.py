#!/usr/bin/env python3

# example usage: add photons worth 2mA/cm^2 to AM1.5 in the region from 400 to 500nm and subtract 2mA/cm^2 of them from the 600 to 700nm one:
# cat ASTMG173.csv | ./spectralmod.py -b 600 -e 700 -j=-2 | ./spectralmod.py -b 400 -e 500 -j=2 > ASTMG173_modded.tsv

import argparse
import pandas as pd
from sys import stdin, stdout

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
parser.add_argument(
    "-o",
    "--offset",
    type=float,
    default=0,
    help="Raw offset volume (if -pi, units here might be photons/m^2/s). takes precidence over -j option",
)
parser.add_argument(
    "-j",
    "--j-offset",
    type=float,
    default=0,
    help="Apply offset in mA/cm^2. auto-sets --photon-input",
)
parser.add_argument(
    "-t",
    "--tsv",
    default=False,
    action="store_true",
    help="Output in tsv instead of csv",
)
parser.add_argument(
    "-1",
    "--one-col",
    default=False,
    action="store_true",
    help="Just output the rescaled column",
)
parser.add_argument(
    "-pi",
    "--photon-input",
    default=False,
    action="store_true",
    help="Convert power input to photons before doing any mods",
)
parser.add_argument(
    "-po",
    "--photon-output",
    default=False,
    action="store_true",
    help="Output photons",
)
parser.add_argument(
    "--col", default=2, type=int, help="Column of interest in input file"
)

args = parser.parse_args()

h = 6.62607004081e-34  # [m^2*kg/s] planck constant
cee = 299792458  # [m/s] speed of light
hc = h * cee  # [J*m]
q = 1.60217657e-19  # [C] or [A*s] elementary charge
ma_per_a = 1000  # mA per A
cm2_per_m2 = 10000  # cm^2 per m^2

width = args.end - args.begin

if args.input == "-":
    infile = stdin
else:
    infile = args.input

df = pd.read_csv(infile, sep=None, engine="python", header=None)

mask = (df[0] >= args.begin) & (df[0] <= args.end)
if df.shape[1] == 2:
    c = 1
else:
    c = args.col
E = (1e9 * hc) / df[0]

use_offset = 0
if args.j_offset:
    args.photon_input = True
    use_offset = cm2_per_m2 * args.j_offset / ma_per_a / q / width

if args.photon_input:
    df.loc[:, c] = df[c] / E

if args.offset:
    use_offset = args.offset / width

df.loc[mask, c] = df[c][mask] * args.scale_factor
df.loc[mask, c] = df[c][mask] + use_offset

if args.photon_input:
    if args.photon_output:
        pass
    else:
        df.loc[:, c] = df[c] * E
else:
    if args.photon_output:
        df.loc[:, c] = df[c] / E

if args.one_col:
    df = pd.concat((df[0], df[c]), axis="columns")
if args.tsv:
    dlm = "\t"
else:
    dlm = ","
df.to_csv(stdout, sep=dlm, float_format="%.4E", header=False, index=False)
