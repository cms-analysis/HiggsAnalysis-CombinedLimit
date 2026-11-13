#!/usr/bin/env python3
"""Check that the CombineHarvester python interface can parse, modify and write a datacard."""

import os
import sys
import tempfile

import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch

card = sys.argv[1]

cb = ch.CombineHarvester()
cb.SetVerbosity(0)
cb.ParseDatacard(card, "")

assert set(cb.bin_set()) == {"htt_tt_9_13TeV"}
assert set(cb.cp().signals().process_set()) == {"ggH", "bbH"}
assert set(cb.cp().backgrounds().process_set()) == {"ZL", "TTT", "VVT", "ZTT", "jetFakes"}
assert abs(cb.GetObservedRate() - 3416.0) < 1e-6
bkg_rate = cb.cp().backgrounds().GetRate()
assert abs(bkg_rate - 3608.6656) < 1e-3, bkg_rate

# Add a systematic with the python-only AddSyst implementation
cb.cp().signals().AddSyst(cb, "test_syst", "lnN", ch.SystMap()(1.05))
assert "test_syst" in cb.syst_name_set()

# Round trip through the datacard writer
with tempfile.TemporaryDirectory() as tmp:
    out_card = os.path.join(tmp, "card.txt")
    cb.WriteDatacard(out_card, os.path.join(tmp, "shapes.root"))
    cb2 = ch.CombineHarvester()
    cb2.SetVerbosity(0)
    cb2.ParseDatacard(out_card, "")
    assert abs(cb2.cp().backgrounds().GetRate() - bkg_rate) < 1e-6 * bkg_rate
    assert set(cb2.syst_name_set()) == set(cb.syst_name_set())

print("All checks passed")
