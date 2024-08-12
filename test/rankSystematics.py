#!/usr/bin/env python3
from __future__ import absolute_import

import json
import os
from itertools import combinations
from optparse import OptionParser
from pprint import pprint
from sys import argv, exit, stderr, stdout

import ROOT

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append("-b-")

ROOT.gROOT.SetBatch(True)
argv.remove("-b-")

masses = list()
parser = OptionParser(usage="usage: %prog [options] dir \nrun with --help to get list of options")
parser.add_option(
    "--dir",
    dest="dir",
    default="systs",
    type="string",
    help="Output directory for results and intermediates.",
)
parser.add_option(
    "--masses",
    dest="masses",
    default=[120],
    type="string",
    action="callback",
    help="Which mass points to scan.",
    callback=lambda option, opt_str, value, parser: masses.extend([int(x) for x in value.split(",")]),
)
parser.add_option(
    "--X-keep-global-nuisances",
    dest="keepGlobalNuisances",
    action="store_true",
    default=False,
    help="When excluding nuisances, do not exclude standard nuisances that are correlated CMS-wide even if their effect is small.",
)
(options, args) = parser.parse_args()

if masses:
    options.masses = masses
    del masses

OWD = os.getcwd()
os.chdir(options.dir)

limits = json.load(open("limits.json", "r"))
limits.sort(key=lambda x: len(x[0]))

allIn = limits[0][1]
allOut = limits[-1][1]
nuisances = set(limits[-1][0])
variations = limits[1:-1]


def rel(test, ref, what="500"):
    return (test[what] - ref[what]) / ref[what]


results = [
    (
        tuple(limit[0]),
        {
            "obs": {
                "AllIn": rel(limit[1], allIn, "-1000"),
                "AllOut": rel(limit[1], allOut, "-1000"),
            },
            "exp": {
                "AllIn": rel(limit[1], allIn, "500"),
                "AllOut": rel(limit[1], allOut, "500"),
            },
        },
    )
    for limit in variations
]


def metric(x, ref="AllIn"):
    return max((abs(x["obs"][ref]), abs(x["exp"][ref])))


results.sort(key=lambda x: metric(x[1], "AllIn"), reverse=True)

Singles = [x for x in results if len(x[0]) == len(nuisances) - 1]
nMinusOnes = [x for x in results if len(x[0]) == 1]

nuisancesKept = dict([(tuple(v[0]), str(list(set(nuisances) - set(v[0]))[0])) for v in Singles])

SinglesLabelled = [(nuisancesKept[v[0]], v[1]) for v in Singles]

# pprint(nMinusOnesLabelled)
g = open("Kept1.json", "w")
g.write(json.dumps(SinglesLabelled, indent=2))
g.close()

g = open("Removed1.json", "w")
g.write(json.dumps(nMinusOnes, indent=2))
g.close()

os.chdir(OWD)
