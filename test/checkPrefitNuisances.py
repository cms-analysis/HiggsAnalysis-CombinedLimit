#!/usr/bin/env python3
"""Check the pre-fit nuisance parameter information written by FitDiagnostics.

Usage: checkPrefitNuisances.py fitDiagnostics.root [name=expected_error ...]

Prints the pre-fit fit result status and the pre-fit values and uncertainties
of all nuisance parameters, and exits with a non-zero code if the fit result is
missing or not valid, or if the uncertainty of a parameter differs from the
expected one by more than 2 %.
"""

import sys

import ROOT

ok = True

tfile = ROOT.TFile.Open(sys.argv[1])

res = tfile.Get("nuisances_prefit_res")
if not res:
    print("ERROR: nuisances_prefit_res not found")
    ok = False
else:
    print(f"nuisances_prefit_res: status={res.status()} covQual={res.covQual()}")
    if res.status() != 0 or res.covQual() != 3:
        print("ERROR: the pre-fit nuisance parameter fit result is not valid")
        ok = False

snapshot = tfile.Get("nuisances_prefit")
for var in snapshot:
    print(f"{var.GetName()}: {var.getVal():.3f} +/- {var.getError():.3f} constant={int(var.isConstant())}")

for spec in sys.argv[2:]:
    name, expected = spec.split("=")
    expected = float(expected)
    var = snapshot.find(name)
    if not var:
        print(f"ERROR: {name} not found in nuisances_prefit")
        ok = False
    elif abs(var.getError() - expected) > 0.02 * max(1.0, expected):
        print(f"ERROR: pre-fit uncertainty of {name} is {var.getError()}, expected {expected}")
        ok = False

sys.exit(0 if ok else 1)
