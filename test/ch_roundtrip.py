#!/usr/bin/env python3
"""Parse datacards with CombineHarvester and write them again.

Usage: ch_roundtrip.py [--mass=MASS] CARD OUTPUT [[--mass=MASS] CARD OUTPUT ...]

The "-ch" tests read OUTPUT instead of CARD. MASS is put in the shapes
patterns of the next CARD in place of $MASS, as combine does; without it no
mass is put in. The shapes of OUTPUT are written next to it, to OUTPUT with
"_shapes.root" in place of ".txt"."""

import os
import re
import sys

import ROOT

import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch


def roundtrip(card, mass, out_card):
    out_dir = os.path.dirname(out_card)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    cb = ch.CombineHarvester()
    cb.SetVerbosity(0)
    cb.ParseDatacard(card, "", "", "", 0, mass)
    cb.WriteDatacard(out_card, re.sub(r"\.txt$", "", out_card) + "_shapes.root")


ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

args = sys.argv[1:]
mass = ""
while args:
    if args[0].startswith("--mass="):
        # The mass as text2workspace.py puts it in the shapes patterns
        mass = float(args.pop(0)[len("--mass=") :])
        mass = "%d" % mass if mass.is_integer() else str(mass)
        continue
    if len(args) < 2:
        sys.exit(__doc__)
    card, out_card = args.pop(0), args.pop(0)
    roundtrip(card, mass, out_card)
    mass = ""
