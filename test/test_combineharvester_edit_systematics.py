#!/usr/bin/env python3
"""Check the ways to rename and change the systematic uncertainties of a
datacard with CombineHarvester that docs/part2/editingsystematics.md gives.
Each one is applied to a parsed datacard, and the datacard that is written is
read with combine's own datacard parser."""

import math
import os
import re
import subprocess
import tempfile
from optparse import OptionParser

import ROOT

import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch
from HiggsAnalysis.CombinedLimit.DatacardParser import addDatacardParserOptions, parseCard

ROOT.gROOT.SetBatch(True)

CARD = """\
imax 2
jmax 1
kmax *
---------------
bin         ch_a ch_b
observation 10   20
---------------
bin      ch_a ch_a ch_b ch_b
process  sig  bkg  sig  bkg
process  0    1    0    1
rate     1    10   2    20
---------------
lumi     lnN  1.02 1.02 1.02 1.02
eff      lnN  1.05 -    1.07 -
norm     lnN  -    1.10 -    1.20
bkg_stat lnN  -    1.04 -    1.03
exp group = lumi eff
"""

SYSTS = {
    "lumi": ("lnN", {("ch_a", "sig"): 1.02, ("ch_a", "bkg"): 1.02, ("ch_b", "sig"): 1.02, ("ch_b", "bkg"): 1.02}),
    "eff": ("lnN", {("ch_a", "sig"): 1.05, ("ch_b", "sig"): 1.07}),
    "norm": ("lnN", {("ch_a", "bkg"): 1.10, ("ch_b", "bkg"): 1.20}),
    "bkg_stat": ("lnN", {("ch_a", "bkg"): 1.04, ("ch_b", "bkg"): 1.03}),
}

SHAPE_CARD = """\
imax 2
jmax 1
kmax *
shapes * * hists.root $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC
bin         ch_a ch_b
observation -1   -1
bin     ch_a ch_a ch_b ch_b
process sig  bkg  sig  bkg
process 0    1    0    1
rate    -1   -1   -1   -1
scale   shape 1  -    1    -
"""


def ch_parse(card_file):
    cb = ch.CombineHarvester()
    cb.SetVerbosity(0)
    cb.ParseDatacard(card_file, "")
    return cb


def combine_parse(card_file):
    """Read the card with combine's own datacard parser."""
    parser = OptionParser()
    addDatacardParserOptions(parser)
    options, _ = parser.parse_args([])
    with open(card_file) as f:
        return parseCard(f, options)


def check(cb, out_card, changes, groups={"exp": {"lumi", "eff"}}):
    """Write the datacard and compare what combine reads from it with the
    systematics of CARD after these changes: None for a name removes it."""
    cb.WriteDatacard(out_card)
    dc = combine_parse(out_card)
    found = {name: (pdf, {(b, p): v for b in errline for p, v in errline[b].items() if v}) for name, _, pdf, _, errline in dc.systs}
    expected = {name: syst for name, syst in {**SYSTS, **changes}.items() if syst is not None}
    assert found.keys() == expected.keys(), (out_card, sorted(found), sorted(expected))
    for name, (pdf, effects) in expected.items():
        assert found[name][0] == pdf, (out_card, name, found[name][0])
        assert found[name][1].keys() == effects.keys(), (out_card, name, found[name][1])
        for key, value in effects.items():
            for got, want in zip(*([v] if isinstance(v, float) else v for v in (found[name][1][key], value)), strict=True):
                assert math.isclose(got, want, rel_tol=1e-5), (out_card, name, key, found[name][1][key])
    assert {g: set(names) for g, names in dc.groups.items()} == groups, (out_card, dc.groups)
    # Every systematic has its parameter in CombineHarvester
    assert set(expected) <= {p.name() for p in cb.GetParameters()}, (out_card, [p.name() for p in cb.GetParameters()])


with tempfile.TemporaryDirectory() as tmp:
    card = os.path.join(tmp, "card.txt")
    with open(card, "w") as f:
        f.write(CARD)

    def out(name):
        return os.path.join(tmp, name + ".txt")

    # --- the card as it is
    check(ch_parse(card), out("unchanged"), {})

    # --- rename a systematic everywhere
    cb = ch_parse(card)
    cb.RenameSystematic(cb, "lumi", "lumi_13TeV")
    check(cb, out("rename"), {"lumi": None, "lumi_13TeV": SYSTS["lumi"]}, groups={"exp": {"lumi_13TeV", "eff"}})

    # --- rename all the systematics that match a pattern
    cb = ch_parse(card)
    for name in cb.syst_name_set():
        if re.fullmatch("eff|norm", name):
            cb.RenameSystematic(cb, name, "CMS_" + name)
    check(cb, out("rename_pattern"), {"eff": None, "norm": None, "CMS_eff": SYSTS["eff"], "CMS_norm": SYSTS["norm"]}, groups={"exp": {"lumi", "CMS_eff"}})

    # --- rename for a part of the processes and bins: the old name keeps its
    # parameter, and the new one is in the same groups
    cb = ch_parse(card)
    cb.cp().bin(["ch_b"]).process(["sig"]).RenameSystematic(cb, "eff", "eff_b")
    check(
        cb,
        out("rename_part"),
        {"eff": ("lnN", {("ch_a", "sig"): 1.05}), "eff_b": ("lnN", {("ch_b", "sig"): 1.07})},
        groups={"exp": {"lumi", "eff", "eff_b"}},
    )
    assert cb.GetParameter("eff") and cb.GetParameter("eff_b")
    subprocess.run(["text2workspace.py", out("rename_part"), "-o", os.path.join(tmp, "rename_part.root")], check=True, stdout=subprocess.DEVNULL)

    # --- add a systematic, with one value, with values per bin, and asymmetric
    cb = ch_parse(card)
    cb.cp().process(["sig"]).AddSyst(cb, "theory", "lnN", ch.SystMap()(1.1))
    cb.cp().process(["bkg"]).AddSyst(cb, "fakes_$BIN", "lnN", ch.SystMap("bin")(["ch_a"], 1.2)(["ch_b"], (0.8, 1.3)))
    check(
        cb,
        out("add"),
        {
            "theory": ("lnN", {("ch_a", "sig"): 1.1, ("ch_b", "sig"): 1.1}),
            "fakes_ch_a": ("lnN", {("ch_a", "bkg"): 1.2}),
            "fakes_ch_b": ("lnN", {("ch_b", "bkg"): [0.8, 1.3]}),
        },
    )

    # --- remove a systematic from a bin, and everywhere
    cb = ch_parse(card)
    cb.FilterSysts(lambda s: s.name() == "norm" and s.bin() == "ch_a")
    cb.FilterSysts(lambda s: s.name() == "bkg_stat")
    check(cb, out("remove"), {"norm": ("lnN", {("ch_b", "bkg"): 1.20}), "bkg_stat": None})

    # --- change the values: one value, an asymmetric one, and the direction
    cb = ch_parse(card)
    cb.cp().bin(["ch_a"]).syst_name(["eff"]).ForEachSyst(lambda s: s.set_value_u(1.08))

    def set_asymm(s):
        s.set_asymm(True)
        s.set_value_d(0.9)
        s.set_value_u(1.15)

    cb.cp().bin(["ch_b"]).syst_name(["norm"]).ForEachSyst(set_asymm)
    cb.cp().syst_name(["bkg_stat"]).ForEachSyst(lambda s: s.set_value_u(1 / s.value_u()))
    check(
        cb,
        out("values"),
        {
            "eff": ("lnN", {("ch_a", "sig"): 1.08, ("ch_b", "sig"): 1.07}),
            "norm": ("lnN", {("ch_a", "bkg"): 1.10, ("ch_b", "bkg"): [0.9, 1.15]}),
            "bkg_stat": ("lnN", {("ch_a", "bkg"): 1 / 1.04, ("ch_b", "bkg"): 1 / 1.03}),
        },
    )

    # --- change the type
    cb = ch_parse(card)
    cb.cp().syst_name(["norm"]).ForEachSyst(lambda s: s.set_type("lnU"))
    check(cb, out("type"), {"norm": ("lnU", SYSTS["norm"][1])})

    # --- split a systematic into two
    cb = ch_parse(card)
    cb.cp().syst_name(["norm"]).ForEachSyst(lambda s: s.set_value_u(1.06))
    cb.cp().process(["bkg"]).AddSyst(cb, "norm_theory", "lnN", ch.SystMap()(1.08))
    check(
        cb,
        out("split"),
        {"norm": ("lnN", {("ch_a", "bkg"): 1.06, ("ch_b", "bkg"): 1.06}), "norm_theory": ("lnN", {("ch_a", "bkg"): 1.08, ("ch_b", "bkg"): 1.08})},
    )

    # --- merge a systematic into another, in quadrature
    cb = ch_parse(card)
    stat = {}
    cb.cp().syst_name(["bkg_stat"]).ForEachSyst(lambda s: stat.update({(s.bin(), s.process()): s.value_u()}))

    def add_stat(s):
        if (s.bin(), s.process()) in stat:
            s.set_value_u(math.exp(math.hypot(math.log(s.value_u()), math.log(stat[s.bin(), s.process()]))))

    cb.cp().syst_name(["norm"]).ForEachSyst(add_stat)
    cb.FilterSysts(lambda s: s.name() == "bkg_stat")
    merged = {key: math.exp(math.hypot(math.log(value), math.log(SYSTS["bkg_stat"][1][key]))) for key, value in SYSTS["norm"][1].items()}
    check(cb, out("merge"), {"norm": ("lnN", merged), "bkg_stat": None})

    # --- groups, to freeze the parameters with --freezeNuisanceGroups
    cb = ch_parse(card)
    cb.SetGroup("bkg", ["norm", "bkg_.*"])
    cb.RemoveGroup("exp", ["eff"])
    check(cb, out("groups"), {}, groups={"exp": {"lumi"}, "bkg": {"norm", "bkg_stat"}})

    # --- freeze: a frozen rateParam is written as a "nuisance edit freeze"
    # line, and so is any line added with AddDatacardLineAtEnd
    cb = ch_parse(card)
    cb.cp().process(["bkg"]).AddSyst(cb, "r_bkg", "rateParam", ch.SystMap()(1.0))
    cb.GetParameter("r_bkg").set_frozen(True)
    cb.AddDatacardLineAtEnd("nuisance edit freeze lumi")
    cb.WriteDatacard(out("frozen"))
    with open(out("frozen")) as f:
        text = f.read()
    assert re.search(r"^r_bkg\s+rateParam", text, re.M), text
    assert combine_parse(out("frozen")).frozenNuisances == {"r_bkg", "lumi"}, text

    # --- rename a shape systematic for one bin: its histograms are written
    # with the new name
    f = ROOT.TFile.Open(os.path.join(tmp, "hists.root"), "RECREATE")
    for chan in ("ch_a", "ch_b"):
        f.mkdir(chan).cd()
        for name, rate, shift in (("data_obs", 60, 0), ("bkg", 50, 0), ("sig", 5, 0), ("sig_scaleUp", 5, 1), ("sig_scaleDown", 5, -1)):
            h = ROOT.TH1F(name, name, 2, 0, 2)
            h.SetBinContent(1, 0.5 * rate + shift)
            h.SetBinContent(2, 0.5 * rate - shift)
            h.Write()
    f.Close()
    shape_card = os.path.join(tmp, "shape_card.txt")
    with open(shape_card, "w") as f:
        f.write(SHAPE_CARD)
    cb = ch_parse(shape_card)
    cb.cp().bin(["ch_b"]).RenameSystematic(cb, "scale", "scale_b")
    cb.WriteDatacard(out("shape_rename"), os.path.join(tmp, "shape_rename.root"))
    dc = combine_parse(out("shape_rename"))
    found = {name: {(b, p) for b in errline for p, v in errline[b].items() if v} for name, _, _, _, errline in dc.systs}
    assert found == {"scale": {("ch_a", "sig")}, "scale_b": {("ch_b", "sig")}}, found
    subprocess.run(["text2workspace.py", out("shape_rename"), "-o", os.path.join(tmp, "shape_rename_ws.root")], check=True, stdout=subprocess.DEVNULL)
    f = ROOT.TFile.Open(os.path.join(tmp, "shape_rename_ws.root"))
    w = f.Get("w")
    assert w.var("scale") and w.var("scale_b")
    f.Close()

print("All checks passed")
