#!/usr/bin/env python3
"""Check that CombineHarvester parses shape systematics and can rename them after parsing."""

import math
import os
import sys
import tempfile

import ROOT

import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch

card = sys.argv[1]
old_name = "CMS_scale_t_1prong0pi0_13TeV"
new_name = "CMS_scale_t_1prong0pi0_2016"
syst_procs = {"TTT", "VVT", "ZTT", "ggH", "bbH"}


def bin_contents(hist):
    return [hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)]


def get_systs(cb, name):
    """Return {process: (type, value_d, value_u, shape_d, shape_u)} for the systematic `name`."""
    res = {}

    def fill(syst):
        res[syst.process()] = (
            syst.type(),
            syst.value_d(),
            syst.value_u(),
            bin_contents(syst.ShapeDAsTH1F()),
            bin_contents(syst.ShapeUAsTH1F()),
        )

    cb.cp().syst_name([name]).ForEachSyst(fill)
    return res


def assert_close(a, b, rel_tol=1e-6):
    if isinstance(a, (list, tuple)):
        assert len(a) == len(b), (a, b)
        for x, y in zip(a, b):
            assert_close(x, y, rel_tol)
    elif isinstance(a, str):
        assert a == b, (a, b)
    else:
        assert math.isclose(a, b, rel_tol=rel_tol, abs_tol=1e-9), (a, b)


def parse(card_file):
    cb = ch.CombineHarvester()
    cb.SetVerbosity(0)
    cb.ParseDatacard(card_file, mass="200")
    return cb


cb = parse(card)

# Shape systematics are parsed together with their up/down templates
assert set(cb.syst_type_set()) == {"lnN", "shape", "rateParam"}
assert len(cb.cp().syst_type(["shape"]).syst_name_set()) == 53
systs = get_systs(cb, old_name)
assert set(systs) == syst_procs, systs.keys()
for proc, (stype, value_d, value_u, shape_d, shape_u) in systs.items():
    assert stype == "shape", (proc, stype)
    assert value_d != 1.0 and value_u != 1.0, (proc, value_d, value_u)
    assert len(shape_d) == len(shape_u) == 18, proc
    assert shape_d != shape_u, proc

bkg_rate = cb.cp().backgrounds().GetRate()
bkg_unc = cb.cp().backgrounds().GetUncertainty()
sig_unc = cb.cp().signals().GetUncertainty()

# Rename the systematic and its parameter
cb.RenameSystematic(cb, old_name, new_name)
assert old_name not in cb.syst_name_set()
assert new_name in cb.syst_name_set()
assert not cb.GetParameter(old_name)
assert cb.GetParameter(new_name).name() == new_name
assert get_systs(cb, new_name) == systs
assert_close(cb.cp().backgrounds().GetRate(), bkg_rate)
assert_close(cb.cp().backgrounds().GetUncertainty(), bkg_unc)
assert_close(cb.cp().signals().GetUncertainty(), sig_unc)

# The renamed systematic is written out with its templates and parsed back
with tempfile.TemporaryDirectory() as tmp:
    out_card = os.path.join(tmp, "card.txt")
    out_shapes = os.path.join(tmp, "shapes.root")
    cb.WriteDatacard(out_card, out_shapes)

    with open(out_card) as f:
        card_text = f.read()
    assert new_name in card_text
    assert old_name not in card_text

    f = ROOT.TFile.Open(out_shapes)
    hist_names = {key.GetName() for key in f.Get("htt_tt_9_13TeV").GetListOfKeys()}
    f.Close()
    for proc in syst_procs:
        hist_proc = proc + "200" if proc in ("ggH", "bbH") else proc
        for shift in ("Up", "Down"):
            assert f"{hist_proc}_{new_name}{shift}" in hist_names, (hist_proc, shift)
    assert not any(old_name in name for name in hist_names)

    cb2 = parse(out_card)
    assert old_name not in cb2.syst_name_set()
    systs2 = get_systs(cb2, new_name)
    assert set(systs2) == syst_procs, systs2.keys()
    for proc in syst_procs:
        assert_close(systs2[proc], systs[proc], rel_tol=1e-5)
    assert_close(cb2.cp().backgrounds().GetRate(), bkg_rate, rel_tol=1e-5)
    assert_close(cb2.cp().backgrounds().GetUncertainty(), bkg_unc, rel_tol=1e-5)

print("All checks passed")
