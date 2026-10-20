#!/usr/bin/env python3
"""Check that CombineHarvester parses, evaluates and writes gmN systematics."""

import math
import os
import re
import sys
import tempfile
from optparse import OptionParser

import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch
from HiggsAnalysis.CombinedLimit.DatacardParser import addDatacardParserOptions, parseCard

card = sys.argv[1]
n_events = 4
alpha = 0.16


def parse(card_file):
    cb = ch.CombineHarvester()
    cb.SetVerbosity(0)
    cb.ParseDatacard(card_file, "")
    return cb


def syst_info(cb):
    """Return {(name, process): (type, value_u, n_events)} for all systematics."""
    res = {}
    cb.ForEachSyst(lambda s: res.__setitem__((s.name(), s.process()), (s.type(), s.value_u(), s.n_events())))
    return res


def combine_parse(card_file):
    """Read the card with combine's own datacard parser."""
    parser = OptionParser()
    addDatacardParserOptions(parser)
    options, _ = parser.parse_args([])
    with open(card_file) as f:
        return parseCard(f, options)


cb = parse(card)
systs = syst_info(cb)
assert set(cb.syst_name_set()) == {"lumi", "xs_ggH", "WW_norm", "xs_ggWW", "bg_others"}, cb.syst_name_set()
assert set(cb.syst_type_set()) == {"lnN", "gmN"}, cb.syst_type_set()
assert [key for key in systs if key[0] == "WW_norm"] == [("WW_norm", "qqWW")]
assert systs[("WW_norm", "qqWW")] == ("gmN", alpha, n_events), systs[("WW_norm", "qqWW")]
assert cb.cp().syst_type(["gmN"]).SetFromSysts(lambda s: s.process()) == {"qqWW"}

# The gmN parameter is the expected number of control region events, with its
# Poisson uncertainty
par = cb.GetParameter("WW_norm")
assert par.val() == n_events, par.val()
assert par.err_d() == -math.sqrt(n_events) and par.err_u() == math.sqrt(n_events), (par.err_d(), par.err_u())

# Rates follow the datacard and scale with the parameter; the uncertainty from
# the sideband is alpha * sqrt(N)
qqWW = cb.cp().process(["qqWW"])
assert math.isclose(qqWW.GetRate(), alpha * n_events), qqWW.GetRate()
assert math.isclose(qqWW.GetUncertainty(), alpha * math.sqrt(n_events)), qqWW.GetUncertainty()
bkg_rate = cb.cp().backgrounds().GetRate()
bkg_unc = cb.cp().backgrounds().GetUncertainty()
assert math.isclose(bkg_rate, 0.64 + 0.06 + 0.22), bkg_rate
assert bkg_unc > qqWW.GetUncertainty()
par.set_val(6)
assert math.isclose(qqWW.GetRate(), alpha * 6), qqWW.GetRate()
par.set_val(n_events)

# Round trip through the datacard writer, checked with CombineHarvester and
# with combine's own parser
with tempfile.TemporaryDirectory() as tmp:
    out_card = os.path.join(tmp, "card.txt")
    cb.WriteDatacard(out_card, os.path.join(tmp, "shapes.root"))

    with open(out_card) as f:
        card_text = f.read()
    assert re.search(r"^WW_norm\s+gmN 4\s+-\s+0\.16\s+-\s+-\s*$", card_text, re.M), card_text

    cb2 = parse(out_card)
    assert syst_info(cb2) == systs
    assert math.isclose(cb2.cp().backgrounds().GetRate(), bkg_rate)
    assert math.isclose(cb2.cp().backgrounds().GetUncertainty(), bkg_unc)

    dc = combine_parse(out_card)
    gmn = [s for s in dc.systs if s[0] == "WW_norm"]
    assert len(gmn) == 1, gmn
    name, nofloat, pdf, args, errline = gmn[0]
    assert pdf == "gmN" and args == [n_events], (pdf, args)
    assert errline["bin1"]["qqWW"] == alpha, errline
    assert all(v == 0 for p, v in errline["bin1"].items() if p != "qqWW"), errline
    assert {s[0]: s[2] for s in dc.systs} == {s[0]: s[2] for s in combine_parse(card).systs}

print("All checks passed")
