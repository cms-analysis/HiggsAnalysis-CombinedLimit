#!/usr/bin/env python3
"""Check that text2workspace.py compares the observation in the datacard with
the yield of the data histogram, also when the data histogram is empty (#312).

Run it in an empty directory: it writes the input histograms and datacards
there, runs text2workspace.py on each case and checks whether it fails."""

import subprocess
import sys

import ROOT

CARD = """imax 1
jmax 1
kmax 0
---
shapes data_obs * shapes.root {data}
shapes *        * shapes.root $PROCESS
---
bin         ch1
observation {obs}
---
bin         ch1  ch1
process     sig  bkg
process     0    1
rate        -1   -1
"""

# (name, data histogram, observation in the card, should text2workspace.py fail?)
CASES = [
    ("empty_data_nonzero_obs", "data_empty", "50", True),
    ("empty_data_zero_obs", "data_empty", "0", False),
    ("empty_data_obs_from_shape", "data_empty", "-1", False),
    ("wrong_data", "data_wrong", "50", True),
    ("matching_data", "data_ok", "50", False),
]


def make_histograms():
    f = ROOT.TFile("shapes.root", "RECREATE")
    for name, total in [("sig", 5), ("bkg", 50), ("data_empty", 0), ("data_wrong", 20), ("data_ok", 50)]:
        h = ROOT.TH1F(name, name, 5, 0, 5)
        for i in range(1, h.GetNbinsX() + 1):
            h.SetBinContent(i, total / 5.0)
        h.Write()
    f.Close()


def main():
    make_histograms()
    failures = []
    for name, data, obs, should_fail in CASES:
        card = f"observed_norm_{name}.txt"
        with open(card, "w") as f:
            f.write(CARD.format(data=data, obs=obs))
        res = subprocess.run(["text2workspace.py", card, "-o", f"observed_norm_{name}.root"], capture_output=True, text=True)
        mismatch = "Mismatch in normalizations for observed data" in res.stderr
        failed = res.returncode != 0
        ok = (failed and mismatch) if should_fail else not failed
        print(f"{name:28s} expected {'failure' if should_fail else 'success':7s}: {'OK' if ok else 'WRONG'}")
        if not ok:
            failures.append(name)
            sys.stdout.write(res.stdout + res.stderr)
    if failures:
        print("FAILED: " + ", ".join(failures))
        sys.exit(1)


if __name__ == "__main__":
    main()
