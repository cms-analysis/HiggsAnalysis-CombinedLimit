#!/usr/bin/env python3
"""Check that CombineHarvester keeps the discrete, flatParam and "nuisance edit"
lines of a datacard, matches autoMCStats patterns, and rejects the lines it
can't represent instead of dropping them."""

import os
import re
import subprocess
import sys
import tempfile
from optparse import OptionParser

import ROOT

import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch
from HiggsAnalysis.CombinedLimit.DatacardParser import addDatacardParserOptions, parseCard

ROOT.gROOT.SetBatch(True)

PARAMETRIC_CARD = """\
imax 1
jmax 1
kmax *
---------------
shapes sig      cat0 ws.root w:sig
shapes bkg      cat0 ws.root w:bkg
shapes data_obs cat0 ws.root w:data
---------------
bin         cat0
observation -1
---------------
bin      cat0 cat0
process  sig  bkg
process  0    1
rate     10   1
---------------
lumi      lnN  1.02  -
xs_sig    lnN  1.05  -
pdf_index discrete
p1        flatParam
bkg_norm  flatParam
nuisance edit freeze xs_sig
nuisance edit rename sig cat0 lumi lumi_13TeV
"""

COUNTING_CARD = """\
imax 3
jmax 1
kmax *
---------------
bin         ch_a ch_b other
observation 10   20   30
---------------
bin      ch_a ch_a ch_b ch_b other other
process  sig  bkg  sig  bkg  sig   bkg
process  0    1    0    1    0     1
rate     1    10   2    20   3     30
---------------
lumi     lnN  1.02 1.02 1.02 1.02 1.02 1.02
{extra}
"""


def make_workspace(path):
    """Signal on top of a background whose shape is a RooMultiPdf inside a
    RooAddPdf, so the index is not the parameter of a top-level RooMultiPdf."""
    w = ROOT.RooWorkspace("w")
    x = ROOT.RooRealVar("CMS_mass", "CMS_mass", 125, 100, 180)
    mean = ROOT.RooRealVar("mean", "mean", 125)
    sigma = ROOT.RooRealVar("sigma", "sigma", 2)
    sig = ROOT.RooGaussian("sig", "sig", x, mean, sigma)
    p1 = ROOT.RooRealVar("p1", "p1", -0.03, -0.1, 0)
    p2 = ROOT.RooRealVar("p2", "p2", -0.01, -0.1, 0)
    exp1 = ROOT.RooExponential("exp1", "exp1", x, p1)
    exp2 = ROOT.RooExponential("exp2", "exp2", x, p2)
    pdf_index = ROOT.RooCategory("pdf_index", "pdf_index")
    multipdf = ROOT.RooMultiPdf("multipdf", "multipdf", pdf_index, ROOT.RooArgList(exp1, exp2))
    flat = ROOT.RooUniform("flat", "flat", ROOT.RooArgSet(x))
    frac = ROOT.RooRealVar("frac", "frac", 0.9)
    bkg = ROOT.RooAddPdf("bkg", "bkg", ROOT.RooArgList(multipdf, flat), ROOT.RooArgList(frac))
    bkg_norm = ROOT.RooRealVar("bkg_norm", "bkg_norm", 1000, 0, 10000)
    ROOT.RooRandom.randomGenerator().SetSeed(1234)
    data = exp1.generate(ROOT.RooArgSet(x), 1000)
    data.SetName("data")
    for obj in (sig, bkg, bkg_norm, data):
        w.Import(obj, ROOT.RooFit.Silence())
    w.writeToFile(path)


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


def write_card(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def ch_rejects(card_file, match):
    """Check that CombineHarvester refuses to parse the card, with an error
    containing match. This runs in a separate process: depending on the
    platform, an exception from the C++ library reaches python or aborts."""
    script = (
        "import sys\n"
        "import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch\n"
        "cb = ch.CombineHarvester()\n"
        "cb.SetVerbosity(0)\n"
        "cb.ParseDatacard(sys.argv[1], '')\n"
    )
    res = subprocess.run([sys.executable, "-c", script, card_file], capture_output=True, text=True)
    out = res.stdout + res.stderr
    assert res.returncode != 0, f"{card_file} was accepted:\n{out}"
    assert match in out, f"expected {match!r} in:\n{out}"


def raises(func, match):
    """Check that func() raises an exception whose message contains match."""
    try:
        func()
    except Exception as e:
        assert match in str(e), f"expected {match!r} in: {e}"
        return
    raise AssertionError(f"no exception raised, expected {match!r}")


with tempfile.TemporaryDirectory() as tmp:
    # --- discrete, flatParam and nuisance edit lines survive a round trip
    in_dir = os.path.join(tmp, "in")
    os.mkdir(in_dir)
    make_workspace(os.path.join(in_dir, "ws.root"))
    card = write_card(os.path.join(in_dir, "card.txt"), PARAMETRIC_CARD)

    cb = ch_parse(card)
    # The edits are not applied by CombineHarvester, only kept
    assert set(cb.syst_name_set()) == {"lumi", "xs_sig"}, cb.syst_name_set()

    texts = []
    for name in ("out1", "out2"):
        out_dir = os.path.join(tmp, name)
        os.mkdir(out_dir)
        out_card = os.path.join(out_dir, "card.txt")
        cb.WriteDatacard(out_card, os.path.join(out_dir, "ws.root"))
        with open(out_card) as f:
            texts.append(f.read())
        # Parse what was just written, so that out2 is a second round trip
        cb = ch_parse(out_card)
    text = texts[0]
    assert re.search(r"^pdf_index\s+discrete$", text, re.M), text
    assert re.search(r"^p1\s+flatParam$", text, re.M), text
    assert re.search(r"^bkg_norm\s+flatParam$", text, re.M), text
    assert text.count("nuisance edit freeze xs_sig\n") == 1, text
    assert text.count("nuisance edit rename sig cat0 lumi lumi_13TeV\n") == 1, text
    assert texts[1] == texts[0], "second round trip changed the datacard"

    out_card = os.path.join(tmp, "out1", "card.txt")
    dc = combine_parse(out_card)
    assert dc.discretes == ["pdf_index"], dc.discretes
    assert set(dc.flatParamNuisances) == {"p1", "bkg_norm"}, dc.flatParamNuisances
    assert {s[0] for s in dc.systs} == {"lumi_13TeV", "xs_sig"}, [s[0] for s in dc.systs]
    assert dc.frozenNuisances == {"xs_sig"}, dc.frozenNuisances

    # combine builds a model from the written card, with the discrete index
    # and the flatParam attribute in place
    out_ws = os.path.join(tmp, "out1", "card.root")
    subprocess.run(["text2workspace.py", out_card, "-o", out_ws], check=True, stdout=subprocess.DEVNULL)
    f = ROOT.TFile.Open(out_ws)
    w = f.Get("w")
    assert [p.GetName() for p in w.genobj("discreteParams")] == ["pdf_index"]
    assert w.var("p1").getAttribute("flatParam")
    f.Close()

    # A card for another category of the same workspace doesn't get these lines
    cb.cp().process(["sig"]).WriteDatacard(os.path.join(tmp, "sig_only.txt"), os.path.join(tmp, "sig_only.root"))
    with open(os.path.join(tmp, "sig_only.txt")) as f:
        sig_text = f.read()
    assert "discrete" not in sig_text and "flatParam" not in sig_text, sig_text

    # --- autoMCStats bin patterns
    card = write_card(os.path.join(tmp, "counting_stats.txt"), COUNTING_CARD.format(extra="ch_* autoMCStats 5 1 1\nother autoMCStats 0"))
    cb = ch_parse(card)
    assert set(cb.GetAutoMCStatsBins()) == {"ch_a", "ch_b", "other"}, cb.GetAutoMCStatsBins()
    assert set(combine_parse(card).binParFlags) == {"ch_a", "ch_b", "other"}
    cb = ch_parse(write_card(card, COUNTING_CARD.format(extra="ch_* autoMCStats 5 1 1")))
    assert set(cb.GetAutoMCStatsBins()) == {"ch_a", "ch_b"}, cb.GetAutoMCStatsBins()

    # --- inline comments and values beyond the number of processes on the
    # nuisance lines are ignored
    card = write_card(
        os.path.join(tmp, "counting_comments.txt"),
        COUNTING_CARD.format(extra="xs lnN 1.1 - 1.2 - 1.3 -  # cross sections\nWW gmN 4 - 0.2 - 0.2 - 0.2 0.1"),
    )
    cb = ch_parse(card)
    ch_systs = {}
    cb.ForEachSyst(lambda s: ch_systs.__setitem__((s.name(), s.bin(), s.process()), s.value_u()))
    combine_systs = {(name, b, p): v for name, _, _, _, errline in combine_parse(card).systs for b in errline for p, v in errline[b].items() if v != 0}
    assert ch_systs == combine_systs, (ch_systs, combine_systs)

    # --- a constrained rateParam gets one param line, whatever the order of
    # the param and rateParam lines
    for order, extra in (
        ("rateparam_first", "r_bkg rateParam * bkg 1 [0,5]\nr_bkg param 1 0.1"),
        ("param_first", "r_bkg param 1 0.1\nr_bkg rateParam * bkg 1 [0,5]"),
    ):
        card = write_card(os.path.join(tmp, f"counting_{order}.txt"), COUNTING_CARD.format(extra=extra))
        out_card = os.path.join(tmp, f"counting_{order}_out.txt")
        ch_parse(card).WriteDatacard(out_card)
        with open(out_card) as f:
            text = f.read()
        assert len(re.findall(r"^r_bkg\s+param\s", text, re.M)) == 1, text
        assert [s[0] for s in combine_parse(out_card).systs].count("r_bkg") == 1, text

    # --- lines that CombineHarvester can't represent are rejected
    bad = {
        "missing_bin": ("nobin autoMCStats 0", "no bin named nobin", True),
        "short_line": ("xs lnN 1.1 - 1.1 -", "expected 6 values, one per process, found 4", True),
        "short_gmN": ("WW gmN 4 - 0.2 - 0.2 -", "expected 6 values, one per process, found 5", True),
        "unknown_type": ("xs lnX 1.1 - 1.1 - 1.1 -", "unknown systematic type lnX", True),
        "short_param": ("mu param 1", "param lines need a mean and a sigma", True),
        # Valid for combine, but CombineHarvester has no representation for it
        "unif": ("u unif 0 1 1 - 1 - 1 -", "systematics of type unif are not supported", False),
    }
    for name, (line, match, combine_rejects) in bad.items():
        card = write_card(os.path.join(tmp, name + ".txt"), COUNTING_CARD.format(extra=line))
        ch_rejects(card, match)
        if combine_rejects:
            raises(lambda: combine_parse(card), "")
        else:
            combine_parse(card)

print("All checks passed")
