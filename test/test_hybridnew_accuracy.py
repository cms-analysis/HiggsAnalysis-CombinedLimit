#!/usr/bin/env python3
"""Test that HybridNew reaches the relative accuracy on the limit (issue #1232)."""

import subprocess

import ROOT

CARD = "hybridnew_small_limit_counting.txt"
R_REL_ACC = 0.05
# Factor two of slack; pre-fix was 1.0
TOLERANCE = 2.0

cmd = "combine {card} -M HybridNew --LHCmode LHC-limits --rMax 1 -T 500 --clsAcc 0.02 --rRelAcc {acc}".format(card=CARD, acc=R_REL_ACC)
ret = subprocess.call(cmd.split(" "))
assert ret == 0, "combine failed with exit code {ret}".format(ret=ret)

f = ROOT.TFile.Open("higgsCombineTest.HybridNew.mH120.root")
assert f and not f.IsZombie(), "Could not open the combine output file"

tree = f.Get("limit")
assert tree is not None and tree.GetEntries() > 0, "No limit found in the combine output"
tree.GetEntry(0)
limit, limit_err = tree.limit, tree.limitErr

assert limit > 0, "Expected a positive upper limit, got {limit}".format(limit=limit)
relative = limit_err / limit
assert (
    relative < TOLERANCE * R_REL_ACC
), "HybridNew stopped at r < {limit} +/- {limit_err} ({relative:.1%} relative), " "although --rRelAcc {acc} was requested".format(
    limit=limit, limit_err=limit_err, relative=relative, acc=R_REL_ACC
)

# Both switched off leaves nothing to converge to
both_off = subprocess.run(
    "combine {card} -M HybridNew --LHCmode LHC-limits --rAbsAcc 0 --rRelAcc 0".format(card=CARD).split(" "),
    capture_output=True,
    text=True,
)
assert both_off.returncode != 0, "combine accepted --rAbsAcc 0 --rRelAcc 0"
assert "cannot both be zero" in both_off.stdout + both_off.stderr, "Unexpected error message for --rAbsAcc 0 --rRelAcc 0:\n" + both_off.stdout + both_off.stderr

print("All checks passed.")
