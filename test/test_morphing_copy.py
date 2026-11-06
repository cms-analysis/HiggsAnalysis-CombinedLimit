#!/usr/bin/env python3
"""Check that copies of FastVerticalInterpHistPdf(2) follow their own nuisance parameter.

The copy (a RooWorkspace copy or a deep snapshot) is scanned in theta while the
original stays at theta = 0, and must reproduce the original's values. The
morphing coefficient is either theta itself or a scaled nuisance 0.5*theta, as
text2workspace builds for a shape systematic with a size other than 1.
"""

import math
import os
import sys
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

PDFS = ("morph_v1", "morph_v2")
THETAS = (-1.0, 0.0, 1.0, 2.0)


def make_workspace(path, coef):
    """Write a workspace with the same morphing in both classes.

    coef is "theta", "scaled_first" (0.5*theta added to the workspace before the
    pdfs) or "scaled_with_pdf" (0.5*theta enters only through the pdf import).
    """
    w = ROOT.RooWorkspace("w")
    w.factory("x[0,4]")
    w.var("x").setBins(4)
    w.factory("theta[0,-5,5]")
    if coef == "scaled_first":
        w.factory("prod::c(0.5, theta)")
        c = w.function("c")
    elif coef == "scaled_with_pdf":
        k = ROOT.RooConstVar("k", "", 0.5)
        c = ROOT.RooProduct("c", "", ROOT.RooArgList(k, w.var("theta")))
    else:
        c = w.var("theta")

    # nominal is flat, "up" tilts towards high x, "down" towards low x
    hists = []
    for tag, ramp in (("nom", 0.0), ("up", 0.6), ("dn", -0.2)):
        h = ROOT.TH1F("h_" + tag, "", 4, 0.0, 4.0)
        for i in range(1, 5):
            h.SetBinContent(i, 1.0 + ramp * i)
        hists.append(h)

    # RooArgList keeps bare pointers, so the inputs must stay alive until the import
    keep = []
    pdfs = ROOT.RooArgList()
    for h in hists:
        dh = ROOT.RooDataHist("d" + h.GetName(), "", ROOT.RooArgList(w.var("x")), h)
        p = ROOT.RooHistPdf("p" + h.GetName(), "", ROOT.RooArgSet(w.var("x")), dh)
        keep += [dh, p]
        pdfs.add(p)
    tlist = ROOT.TList()
    for h in hists:
        tlist.Add(h)
    v1 = ROOT.FastVerticalInterpHistPdf("morph_v1", "", w.var("x"), pdfs, ROOT.RooArgList(c), 1.0, 1)
    v2 = ROOT.FastVerticalInterpHistPdf2("morph_v2", "", w.var("x"), tlist, ROOT.RooArgList(c), 1.0, 1)
    for pdf in (v1, v2):
        getattr(w, "import")(pdf, ROOT.RooFit.RecycleConflictNodes())
    w.writeToFile(path)


def scan(pdf, x, theta):
    """pdf value in the first bin for each theta in THETAS"""
    out = []
    for t in THETAS:
        theta.setVal(t)
        x.setVal(0.5)
        out.append(pdf.getVal())
    return out


def copies(w, name):
    """(x, theta, pdf) of a workspace copy and of a deep snapshot of the pdf"""
    cw = ROOT.RooWorkspace(w)
    snap = ROOT.RooArgSet(w.pdf(name)).snapshot(True)
    return [
        ("workspace copy", cw, (cw.var("x"), cw.var("theta"), cw.pdf(name))),
        ("deep snapshot", snap, (snap.find("x"), snap.find("theta"), snap.find(name))),
    ]


failures = []
with tempfile.TemporaryDirectory() as tmp:
    for coef in ("theta", "scaled_first", "scaled_with_pdf"):
        path = os.path.join(tmp, coef + ".root")
        make_workspace(path, coef)
        for name in PDFS:
            f = ROOT.TFile.Open(path)
            w = f.Get("w")
            ref = scan(w.pdf(name), w.var("x"), w.var("theta"))
            f.Close()
            # "cold": the original was never evaluated before the copy, "warm": it was
            for warm in (False, True):
                f = ROOT.TFile.Open(path)
                w = f.Get("w")
                w.var("theta").setVal(0.0)
                if warm:
                    w.pdf(name).getVal()
                # owner holds the copied objects and must stay alive while they are used
                for how, owner, (x, theta, pdf) in copies(w, name):
                    pdf.getVal()  # evaluate the copy at theta = 0 first
                    got = scan(pdf, x, theta)
                    if not all(math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=1e-9) for a, b in zip(got, ref)):
                        failures.append(f"{name}, {coef}, {how}, {'warm' if warm else 'cold'}: got {got}, expected {ref}")
                f.Close()

for failure in failures:
    print("FAILED:", failure)
if failures:
    sys.exit(1)
print("All checks passed")
