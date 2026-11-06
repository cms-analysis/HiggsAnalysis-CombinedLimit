"""Pre-fit Asimov for RooHistPdf shapes on an observable with variable-width bins.

make:  write shapes_RooHistPdf_variable_binning.root, with nominal bin contents 10, 20, 30, 40, 50
       on the bin edges 0, 1, 2, 4, 6, 10, shape-systematic templates, and data_obs equal to the nominal.
check: generate the Asimov dataset with combine at shape = 0, +1 and -1 and verify that it
       reproduces the corresponding template (normalised to the datacard rate) in every bin.
"""

import array
import math
import subprocess
import sys

import ROOT

EDGES = [0, 1, 2, 4, 6, 10]
CONTENTS = [10.0, 20.0, 30.0, 40.0, 50.0]
RATE = 150.0
SHIFTS = {"data_obs": 0.0, "sig": 0.0, "sig_shapeUp": 0.1, "sig_shapeDown": -0.1}


def template(shift):
    return [c * (1 + shift * (2 - i) / 2) for i, c in enumerate(CONTENTS)]


def make():
    x = ROOT.RooRealVar("x", "x", EDGES[0], EDGES[-1])
    x.setBinning(ROOT.RooBinning(len(EDGES) - 1, array.array("d", EDGES)))
    w = ROOT.RooWorkspace("w")
    keep = []
    for name, shift in SHIFTS.items():
        h = ROOT.TH1D("h_" + name, "", len(EDGES) - 1, array.array("d", EDGES))
        for i, content in enumerate(template(shift)):
            h.SetBinContent(i + 1, content)
        dh = ROOT.RooDataHist(name if name == "data_obs" else "dh_" + name, "", ROOT.RooArgList(x), h)
        keep += [h, dh]
        obj = dh if name == "data_obs" else ROOT.RooHistPdf(name, "", ROOT.RooArgSet(x), dh)
        getattr(w, "import")(obj)
    w.writeToFile("shapes_RooHistPdf_variable_binning.root")


def check(toy_file, shape):
    # at shape = +-1 the vertical morphing returns the up/down template, normalised to the rate
    ref = template({0: 0.0, 1: SHIFTS["sig_shapeUp"], -1: SHIFTS["sig_shapeDown"]}[shape])
    expected = [c * RATE / sum(ref) for c in ref]
    centers = [(lo + hi) / 2 for lo, hi in zip(EDGES, EDGES[1:])]
    f = ROOT.TFile.Open(toy_file)
    asimov = f.Get("toys/toy_asimov")
    weights = {}
    for i in range(asimov.numEntries()):
        x = asimov.get(i).find("x").getVal()
        weights[x] = weights.get(x, 0.0) + asimov.weight()
    coordinates = sorted(weights)
    got = [weights[x] for x in coordinates]
    print("shape = %g: Asimov bin contents" % shape, got, "expected:", expected)
    if len(coordinates) != len(centers) or any(not abs(x - c) <= 1e-6 for x, c in zip(coordinates, centers)):
        sys.exit("Asimov dataset is not at the bin centers %s: %s" % (centers, coordinates))
    if any(not (math.isfinite(g) and abs(g - c) <= 1e-6 * c) for g, c in zip(got, expected)):
        sys.exit("Asimov dataset does not reproduce the RooHistPdf bin contents")


def run(workspace):
    for shape in (0, 1, -1):
        cmd = "combine -M GenerateOnly %s -t -1 --saveToys --expectSignal 1 --setParameters shape=%d -n VariableBinning%d" % (workspace, shape, shape)
        subprocess.check_call(cmd.split())
        check("higgsCombineVariableBinning%d.GenerateOnly.mH120.123456.root" % shape, shape)


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make()
    else:
        run(sys.argv[2])
