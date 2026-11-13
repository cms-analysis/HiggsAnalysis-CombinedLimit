#!/usr/bin/env python3
"""Parse datacards with CombineHarvester and write them again.

Usage: ch_roundtrip.py [--mass=MASS] [--split] CARD OUTPUT [[--mass=MASS] [--split] CARD OUTPUT ...]

The "-ch" tests read OUTPUT instead of CARD. MASS is put in the shapes
patterns of the next CARD in place of $MASS, as combine does; without it no
mass is put in. The shapes of OUTPUT are written next to it, to OUTPUT with
"_shapes.root" in place of ".txt".

With --split, which needs a MASS, the next CARD is written with ch::CardWriter
instead: OUTPUT is a directory, with a datacard for each bin in
OUTPUT/MASS/BIN.txt and the shapes of all of them in OUTPUT/common/shapes.root.
As in a script that makes datacards with CombineHarvester, only the signals
keep the mass: the data, the backgrounds and their systematics get the mass
"*", which CardWriter writes into the datacards of every mass."""

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


def split(card, mass, out_dir):
    cb = ch.CombineHarvester()
    cb.SetVerbosity(0)
    cb.ParseDatacard(card, "", "", "", 0, mass)

    def any_mass(obj):
        if not obj.signal():
            obj.set_mass("*")

    cb.ForEachObs(lambda obs: obs.set_mass("*"))
    cb.ForEachProc(any_mass)
    cb.ForEachSyst(any_mass)
    shapes = f"{out_dir}/common/shapes.root"
    cards = ch.CardWriter("$TAG/$MASS/$BIN.txt", "$TAG/common/shapes.root").WriteCards(out_dir, cb)
    # A datacard for each bin, with the objects of that bin only
    written = {str(name): sorted(str(b) for b in c.bin_set()) for name, c in cards}
    expected = {f"{out_dir}/{mass}/{b}.txt": [str(b)] for b in cb.bin_set()}
    if written != expected:
        sys.exit(f"CardWriter wrote {written}, expected {expected}")
    missing = [f for f in list(expected) + [shapes] if not os.path.isfile(f)]
    if missing:
        sys.exit(f"CardWriter did not write {missing}")


ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

args = sys.argv[1:]
mass, use_writer = "", False
while args:
    if args[0].startswith("--mass="):
        # The mass as text2workspace.py puts it in the shapes patterns
        mass = float(args.pop(0)[len("--mass=") :])
        mass = "%d" % mass if mass.is_integer() else str(mass)
        continue
    if args[0] == "--split":
        use_writer = True
        args.pop(0)
        continue
    if len(args) < 2 or (use_writer and not mass):
        sys.exit(__doc__)
    card, output = args.pop(0), args.pop(0)
    if use_writer:
        split(card, mass, output)
    else:
        roundtrip(card, mass, output)
    mass, use_writer = "", False
