"""Pre-fit Asimov for RooHistPdf shapes on an observable with variable-width bins.

make:  write shapes_RooHistPdf_variable_binning.root, with bin contents 10, 20, 30, 40, 50
       on the bin edges 0, 1, 2, 4, 6, 10, and data_obs equal to the nominal shape.
check: verify that the Asimov dataset saved by combine reproduces those bin contents.
"""

import array
import sys

import ROOT

EDGES = [0, 1, 2, 4, 6, 10]
CONTENTS = [10.0, 20.0, 30.0, 40.0, 50.0]


def make():
    x = ROOT.RooRealVar("x", "x", EDGES[0], EDGES[-1])
    x.setBinning(ROOT.RooBinning(len(EDGES) - 1, array.array("d", EDGES)))
    w = ROOT.RooWorkspace("w")
    keep = []
    for name, shift in (("data_obs", 0.0), ("sig", 0.0), ("sig_shapeUp", 0.1), ("sig_shapeDown", -0.1)):
        h = ROOT.TH1D("h_" + name, "", len(EDGES) - 1, array.array("d", EDGES))
        for i, content in enumerate(CONTENTS):
            h.SetBinContent(i + 1, content * (1 + shift * (2 - i) / 2))
        dh = ROOT.RooDataHist(name if name == "data_obs" else "dh_" + name, "", ROOT.RooArgList(x), h)
        keep += [h, dh]
        obj = dh if name == "data_obs" else ROOT.RooHistPdf(name, "", ROOT.RooArgSet(x), dh)
        getattr(w, "import")(obj)
    w.writeToFile("shapes_RooHistPdf_variable_binning.root")


def check(toy_file):
    f = ROOT.TFile.Open(toy_file)
    asimov = f.Get("toys/toy_asimov")
    weights = {}
    for i in range(asimov.numEntries()):
        x = asimov.get(i).find("x").getVal()
        weights[x] = weights.get(x, 0.0) + asimov.weight()
    got = [weights[x] for x in sorted(weights)]
    print("Asimov bin contents:", got, "expected:", CONTENTS)
    if len(got) != len(CONTENTS) or any(abs(g - c) > 1e-6 * c for g, c in zip(got, CONTENTS)):
        sys.exit("Asimov dataset does not reproduce the RooHistPdf bin contents")


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make()
    else:
        check(sys.argv[2])
