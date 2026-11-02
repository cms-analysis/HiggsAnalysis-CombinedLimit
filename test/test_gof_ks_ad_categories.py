#!/usr/bin/env python3
"""Test that the KS and AD goodness-of-fit tests handle several categories correctly.

The datacard lists its channels in an order (zz, aa, mm) that differs from the alphabetical
order in which RooAbsData::split() returns them. The channels are independent (no shared
parameters, r frozen), so the value stored in each per-category branch of a combined run must
equal the statistic of that channel run on its own, and the total must equal their sum.

This catches comparing the data of one channel with the model of another, and storing a
channel's value in another channel's branch.
"""

import subprocess

import ROOT

ROOT.gROOT.SetBatch(True)

PREFIX = "gof_categories_"
# Datacard order, deliberately not alphabetical
CHANNELS = ["ch_zz", "ch_aa", "ch_mm"]
NBINS = 20
# Each channel gets its own background shape and its own data distortion, so that using the
# wrong channel's data or model changes the statistic a lot
BKG_SHAPE = {
    "ch_zz": lambda x: 200.0 * (1.2 - x),
    "ch_aa": lambda x: 200.0 * (0.2 + x),
    "ch_mm": lambda x: 120.0,
}
DISTORTION = {"ch_zz": 0.0, "ch_aa": 0.8, "ch_mm": 0.35}
FIT_OPTIONS = ["--setParameters", "r=1", "--freezeParameters", "r"]
REL_TOLERANCE = 1e-3


def run(cmd):
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr)
    assert res.returncode == 0, "Command failed with exit code {ret}: {cmd}".format(ret=res.returncode, cmd=" ".join(cmd))
    return res


def make_inputs():
    rng = ROOT.TRandom3(42)
    shapes = PREFIX + "shapes.root"
    f = ROOT.TFile(shapes, "RECREATE")
    for ch in CHANNELS:
        bkg = ROOT.TH1F(ch + "_bkg", "", NBINS, 0, 1)
        sig = ROOT.TH1F(ch + "_sig", "", NBINS, 0, 1)
        data = ROOT.TH1F(ch + "_data_obs", "", NBINS, 0, 1)
        for i in range(1, NBINS + 1):
            x = bkg.GetBinCenter(i)
            b = BKG_SHAPE[ch](x)
            bkg.SetBinContent(i, b)
            sig.SetBinContent(i, 5.0)
            data.SetBinContent(i, rng.Poisson(b * (1 + DISTORTION[ch] * (x - 0.5)) + 5.0))
        for h in (bkg, sig, data):
            h.Write()
        with open(PREFIX + ch + ".txt", "w") as card:
            card.write(
                "\n".join(
                    [
                        "imax 1",
                        "jmax 1",
                        "kmax 1",
                        "---",
                        "shapes * {ch} {shapes} {ch}_$PROCESS".format(ch=ch, shapes=shapes),
                        "shapes data_obs {ch} {shapes} {ch}_data_obs".format(ch=ch, shapes=shapes),
                        "---",
                        "bin {ch}".format(ch=ch),
                        "observation -1",
                        "---",
                        "bin      {ch} {ch}".format(ch=ch),
                        "process  sig  bkg",
                        "process  0    1",
                        "rate     -1   -1",
                        "---",
                        "norm_{ch} lnN - 1.10".format(ch=ch),
                        "",
                    ]
                )
            )
    f.Close()

    comb = run(["combineCards.py"] + ["{ch}={prefix}{ch}.txt".format(ch=ch, prefix=PREFIX) for ch in CHANNELS])
    with open(PREFIX + "comb.txt", "w") as card:
        card.write(comb.stdout)

    for name in CHANNELS + ["comb"]:
        run(["text2workspace.py", PREFIX + name + ".txt", "-o", PREFIX + name + ".root"])


def goodness_of_fit(name, algo):
    """Run GoodnessOfFit and return the total and any per-category values."""
    tag = ".{prefix}{name}.{algo}".format(prefix=PREFIX, name=name, algo=algo)
    run(["combine", "-M", "GoodnessOfFit", PREFIX + name + ".root", "--algo", algo, "-n", tag] + FIT_OPTIONS)
    f = ROOT.TFile.Open("higgsCombine{tag}.GoodnessOfFit.mH120.root".format(tag=tag))
    tree = f.Get("limit")
    tree.GetEntry(0)
    values = {"total": tree.limit}
    for ch in CHANNELS:
        if tree.GetBranch(ch):
            values[ch] = getattr(tree, ch)
    f.Close()
    return values


def close(a, b):
    return abs(a - b) <= REL_TOLERANCE * max(abs(a), abs(b), 1e-12)


make_inputs()

for algo in ["KS", "AD"]:
    alone = {ch: goodness_of_fit(ch, algo)["total"] for ch in CHANNELS}
    combined = goodness_of_fit("comb", algo)
    for ch in CHANNELS:
        assert ch in combined, "{algo}: no branch {ch} in the combined output".format(algo=algo, ch=ch)
        print("{algo} {ch}: combined branch {c:.6g}, channel alone {a:.6g}".format(algo=algo, ch=ch, c=combined[ch], a=alone[ch]))
        assert close(combined[ch], alone[ch]), "{algo}: branch {ch} holds {c:.6g}, but {ch} alone gives {a:.6g}".format(
            algo=algo, ch=ch, c=combined[ch], a=alone[ch]
        )
    total = sum(alone.values())
    print("{algo} total: combined {c:.6g}, sum of channels {s:.6g}".format(algo=algo, c=combined["total"], s=total))
    assert close(combined["total"], total), "{algo}: total {c:.6g} differs from the sum of channels {s:.6g}".format(algo=algo, c=combined["total"], s=total)

print("All checks passed.")
