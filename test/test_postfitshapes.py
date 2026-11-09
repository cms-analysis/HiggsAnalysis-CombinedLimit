#!/usr/bin/env python3
"""Check the shapes of PostFitShapesFromWorkspace against FitDiagnostics.

Usage: test_postfitshapes.py CARD MASS [TEXT2WORKSPACE_ARGS ...]

The workspace of CARD is fitted with FitDiagnostics --saveShapes, and
PostFitShapesFromWorkspace makes the pre- and post-fit shapes from the same
workspace and fit result. Both are computed from the same model, so the
contents of all the histograms must agree. Both take the post-fit uncertainty
of a bin as the RMS of its yields around the post-fit one, in samples of the
fit covariance matrix, but with different random numbers: they only agree
within the spread of NSAMPLES samples. The pre-fit uncertainties are not
compared: PostFitShapesFromWorkspace adds the shifts of the parameters in
quadrature, while FitDiagnostics samples them."""

import os
import subprocess
import sys
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)

# Number of samples of the fit covariance matrix, in both programs
NSAMPLES = 1000
# Relative tolerance of the contents, and of the post-fit uncertainties. An RMS
# of N samples is uncertain by about 1/sqrt(2N), so the ratio of two of them by
# 1/sqrt(N), 3% for 1000 samples: the largest difference over the 1100 bins of
# the five tests is 10%
RTOL_CONTENT = 2e-3
RTOL_POSTFIT_ERROR = 0.2

# Names of the totals in the output of PostFitShapesFromWorkspace and FitDiagnostics
TOTALS = {"TotalProcs": "total", "TotalBkg": "total_background", "TotalSig": "total_signal"}

card, mass = sys.argv[1], sys.argv[2]
t2w_args = sys.argv[3:]
failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)


def run(cmd, cwd):
    res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    if res.returncode != 0:
        sys.exit(f"{' '.join(cmd)} failed:\n{res.stdout[-3000:]}")


def contents(h, per_width):
    """Bin contents and errors of h, multiplied by the bin widths if per_width"""
    w = [h.GetBinWidth(i) if per_width else 1.0 for i in range(1, h.GetNbinsX() + 1)]
    return [h.GetBinContent(i) * w[i - 1] for i in range(1, h.GetNbinsX() + 1)], [h.GetBinError(i) * w[i - 1] for i in range(1, h.GetNbinsX() + 1)]


def graph_contents(g, binning):
    """Numbers of events in the data graph of FitDiagnostics, which has a
    point per bin of the histogram binning"""
    ys = g.GetY()
    return [ys[i] * binning.GetBinWidth(i + 1) for i in range(g.GetN())]


def close(a, b, rtol, atol):
    return abs(a - b) <= rtol * max(abs(a), abs(b)) + atol


def compare(name, ref, val, rtol, atol):
    """Compare the bins of val to the first bins of ref. FitDiagnostics pads
    the bins of a channel up to those of the longest channel, so the bins of
    ref after those of val must be empty"""
    check(len(ref) >= len(val), f"{name}: {len(val)} bins, FitDiagnostics has {len(ref)}")
    for i, (a, b) in enumerate(zip(ref, val)):
        check(close(a, b, rtol, atol), f"{name}, bin {i + 1}: {b:.6g}, FitDiagnostics {a:.6g}")
    for i in range(len(val), len(ref)):
        check(abs(ref[i]) <= atol, f"{name}, bin {i + 1}: missing, FitDiagnostics {ref[i]:.6g}")


def check_channel(fd_chan, pf_dir, name, stage, sums):
    """Compare the histograms of a channel at one stage, prefit or postfit.
    The shapes of FitDiagnostics are divided by the bin widths; for TH1
    shapes the bins are numbered, so their widths are 1"""
    total = fd_chan.Get("total")
    # Contents, and uncertainties ten times this, are zero below atol
    atol = 1e-6 * max(contents(total, True)[0])
    procs = [
        k.GetName()
        for k in fd_chan.GetListOfKeys()
        if k.GetName() not in TOTALS.values()
        and k.GetName() != "data"
        and fd_chan.Get(k.GetName()).InheritsFrom("TH1")
        and not fd_chan.Get(k.GetName()).InheritsFrom("TH2")
    ]
    for proc in procs + list(TOTALS):
        h = pf_dir.Get(proc)
        if not h:
            check(False, f"no histogram {name}/{proc}")
            continue
        ref_vals, ref_errs = contents(fd_chan.Get(TOTALS.get(proc, proc)), True)
        vals, errs = contents(h, False)
        compare(f"{name}/{proc}", ref_vals, vals, RTOL_CONTENT, atol)
        if stage == "postfit":
            compare(f"{name}/{proc} uncertainty", ref_errs, errs, RTOL_POSTFIT_ERROR, 1e-3 * max(ref_errs) + 10 * atol)
        if proc in TOTALS:
            sums[proc] = sums.get(proc, 0.0) + sum(vals)
    data = pf_dir.Get("data_obs")
    check(bool(data), f"no histogram {name}/data_obs")
    if data:
        compare(f"{name}/data_obs", graph_contents(fd_chan.Get("data"), total), contents(data, False)[0], RTOL_CONTENT, 1e-9)


card = os.path.abspath(card)
# Everything runs in a temporary directory, as combine writes its outputs to
# the working directory
with tempfile.TemporaryDirectory() as tmp:
    run(["text2workspace.py", card, "-o", "ws.root", "-m", mass] + t2w_args, tmp)
    run(["combine", "-M", "FitDiagnostics", "ws.root", "-m", mass, "--saveShapes", "--saveWithUncertainties", "--numToysForShapes", str(NSAMPLES)], tmp)
    pfs = ["PostFitShapesFromWorkspace", "-w", "ws.root", "-d", card, "-m", mass, "-o", "shapes.root", "-f", "fitDiagnosticsTest.root:fit_s"]
    run(pfs + ["--postfit", "--total-shapes", "--samples", str(NSAMPLES)], tmp)
    fd = ROOT.TFile.Open(os.path.join(tmp, "fitDiagnosticsTest.root"))
    pf = ROOT.TFile.Open(os.path.join(tmp, "shapes.root"))

    channels = sorted(k.GetName() for k in fd.Get("shapes_prefit").GetListOfKeys() if k.GetClassName() == "TDirectoryFile")
    check(channels, "FitDiagnostics saved no shapes")
    for stage, fd_dir in (("prefit", "shapes_prefit"), ("postfit", "shapes_fit_s")):
        sums = {}
        for chan in channels:
            pf_dir = pf.Get(f"{chan}_{stage}")
            check(bool(pf_dir), f"no directory {chan}_{stage}")
            if pf_dir:
                check_channel(fd.Get(f"{fd_dir}/{chan}"), pf_dir, f"{chan}_{stage}", stage, sums)
        # --total-shapes: the totals of all the channels together
        for proc in list(TOTALS) + ["data_obs"]:
            h = pf.Get(f"{stage}/{proc}")
            check(bool(h), f"no histogram {stage}/{proc}")
            if h and proc in TOTALS:
                expected = sums.get(proc, 0.0)
                check(
                    close(h.Integral(), expected, RTOL_CONTENT, 1e-6 * max(expected, 1.0)),
                    f"{stage}/{proc}: {h.Integral():.6g}, sum of the channels {expected:.6g}",
                )
    fd.Close()
    pf.Close()

if failures:
    print("\n".join(failures))
    sys.exit(f"{len(failures)} checks failed for {card}")
print("All checks passed")
