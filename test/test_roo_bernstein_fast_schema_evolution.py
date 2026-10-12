#!/usr/bin/env python3
"""Test that a RooBernsteinFast<N> object written with the old v1 schema
(inheriting from RooAbsPdf, c0=1 implicit) is correctly evolved on read-back:
  - the object inherits from RooBernstein (v2 class hierarchy)
  - the PDF evaluates to a finite, positive value (proxies correctly wired)
"""

import math
import ROOT

f = ROOT.TFile.Open("test_bernstein_fast_v1.root")
assert f and not f.IsZombie(), "Could not open test_bernstein_fast_v1.root"

ws = f.Get("w")
assert ws is not None, "Workspace 'w' not found"

pdf = ws.pdf("pdf")
assert pdf is not None, "PDF 'pdf' not found in workspace 'w'"

assert pdf.InheritsFrom("RooBernstein"), f"Expected 'pdf' to inherit from RooBernstein after schema evolution, got {pdf.ClassName()}"

val = pdf.getVal()
assert math.isfinite(val) and val > 0, f"pdf.getVal() returned invalid value: {val}"

print("All checks passed.")
