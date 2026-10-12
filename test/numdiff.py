#!/usr/bin/env python3
"""Compare two text files like diff, but allow a relative tolerance for
floating point numbers. Everything else, including integers, must match
exactly. Usage: numdiff.py --rtol 0.01 REFERENCE OUTPUT"""

import argparse
import re
import sys

FLOAT = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")

parser = argparse.ArgumentParser()
parser.add_argument("--rtol", type=float, required=True, help="relative tolerance for floating point numbers")
parser.add_argument("reference")
parser.add_argument("output")
args = parser.parse_args()


def close(a, b):
    return abs(a - b) <= args.rtol * max(abs(a), abs(b))


def same_line(ref, out):
    if FLOAT.sub("#", ref) != FLOAT.sub("#", out):
        return False
    return all(close(float(a), float(b)) for a, b in zip(FLOAT.findall(ref), FLOAT.findall(out)))


with open(args.reference) as f:
    ref_lines = f.read().splitlines()
with open(args.output) as f:
    out_lines = f.read().splitlines()

ok = len(ref_lines) == len(out_lines)
for ref, out in zip(ref_lines, out_lines):
    if not same_line(ref, out):
        print(f"< {ref}\n> {out}")
        ok = False
if len(ref_lines) != len(out_lines):
    print(f"different number of lines: {len(ref_lines)} vs {len(out_lines)}")
sys.exit(0 if ok else 1)
