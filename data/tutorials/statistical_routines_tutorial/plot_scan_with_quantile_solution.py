"""Plot the 1D likelihood scan together with the 68% quantile of the test
statistic computed from toys at each scanned r value.

It should be run after producing the toy files for each r value and the likelihood scan file using the commands in the tutorial.
For the former, the command is something like:

for r_val in $(seq -2 0.08 6); do
    combine -M MultiDimFit datacard.txt --rMin -10 --rMax 10 --algo fixed \
        --fixedPointPOIs r=${r_val} --setParameters r=${r_val} \
        -t 1000 --toysFrequentist -n ".r_${r_val}"
done

The script:
  1. Loops over higgsCombine.r_<value>.MultiDimFit.mH120.123456.root files
     and computes the 0.68 quantile of 2*deltaNLL for each r value.
  2. Reads the likelihood scan from higgsCombineTest.MultiDimFit.mH120.root.
  3. Plots both curves and finds their crossings via interpolation.
"""

import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from scipy.interpolate import interp1d

from ROOT import TFile


def get_quantile_for_file(filepath, n_sigma=1):
    """Compute the test statistic cutoff at the given sigma level.

    Follows the same logic as get_quantile.py.
    """
    quantile_val = 2 * st.norm().cdf(-n_sigma)

    f = TFile(filepath, "READ")
    limit = f.Get("limit")
    if not limit:
        return None

    m2nll_vals = []
    for i in range(limit.GetEntries()):
        limit.GetEntry(i)
        if limit.quantileExpected < 0:
            continue
        m2nll_vals.append(2 * limit.deltaNLL)

    f.Close()

    if len(m2nll_vals) == 0:
        return None

    return np.quantile(m2nll_vals, 1 - quantile_val)


def read_scan(filepath, poi="r"):
    """Read the likelihood scan from a MultiDimFit output file.

    Returns sorted arrays of (poi_values, 2*deltaNLL).
    """
    f = TFile(filepath, "READ")
    limit = f.Get("limit")

    r_vals = []
    dnll_vals = []
    for i in range(limit.GetEntries()):
        limit.GetEntry(i)
        if limit.quantileExpected < -1.5:
            continue
        r_vals.append(getattr(limit, poi))
        dnll_vals.append(2 * limit.deltaNLL)

    f.Close()

    r_vals = np.array(r_vals)
    dnll_vals = np.array(dnll_vals)

    order = np.argsort(r_vals)
    return r_vals[order], dnll_vals[order]


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # --- Step 1: compute 0.68 quantile for each r value from toy files ---
    pattern = os.path.join(base_dir, "higgsCombine.r_*.MultiDimFit.mH120.123456.root")
    toy_files = sorted(glob.glob(pattern))

    r_points = []
    quantile_points = []
    for fpath in toy_files:
        basename = os.path.basename(fpath)
        match = re.search(r"higgsCombine\.r_([-\d.eE+]+)\.MultiDimFit", basename)
        if not match:
            continue
        r_val = float(match.group(1))
        q = get_quantile_for_file(fpath, n_sigma=1)
        if q is not None:
            r_points.append(r_val)
            quantile_points.append(q)

    r_points = np.array(r_points)
    quantile_points = np.array(quantile_points)

    order = np.argsort(r_points)
    r_points = r_points[order]
    quantile_points = quantile_points[order]

    print(f"Computed quantile for {len(r_points)} r values")

    # --- Step 2: read the likelihood scan ---
    scan_file = os.path.join(base_dir, "higgsCombineTest.MultiDimFit.mH120.root")
    scan_r, scan_dnll = read_scan(scan_file)
    print(f"Read scan with {len(scan_r)} points")

    # --- Step 3: interpolate the quantile curve ---
    quantile_interp = interp1d(r_points, quantile_points, kind="cubic")

    # --- Step 4: find crossings ---
    # Restrict to the overlap region
    r_min = max(r_points.min(), scan_r.min())
    r_max = min(r_points.max(), scan_r.max())

    # Interpolate both curves on a fine common grid
    r_fine = np.linspace(r_min, r_max, 5000)
    scan_interp = interp1d(scan_r, scan_dnll, kind="cubic")
    scan_fine = scan_interp(r_fine)
    quantile_fine = quantile_interp(r_fine)

    diff = scan_fine - quantile_fine
    crossings = []
    for i in range(len(diff) - 1):
        if diff[i] * diff[i + 1] < 0:
            # Linear interpolation for the crossing
            r_cross = r_fine[i] - diff[i] * (r_fine[i + 1] - r_fine[i]) / (diff[i + 1] - diff[i])
            crossings.append(r_cross)

    # Find the best-fit point (minimum of the scan)
    bestfit = scan_r[np.argmin(scan_dnll)]

    # Sort crossings into lo/hi relative to best fit
    crossings_lo = sorted([c for c in crossings if c < bestfit])
    crossings_hi = sorted([c for c in crossings if c >= bestfit])
    lo = crossings_lo[-1] if crossings_lo else None
    hi = crossings_hi[0] if crossings_hi else None

    err_lo = bestfit - lo if lo is not None else None
    err_hi = hi - bestfit if hi is not None else None

    label_parts = [f"r = {bestfit:.3f}"]
    if err_hi is not None:
        label_parts.append(f"+{err_hi:.3f}")
    if err_lo is not None:
        label_parts.append(f"-{err_lo:.3f}")
    bestfit_label = " ".join(label_parts)

    print(f"Best fit: {bestfit_label}")
    print(f"Found {len(crossings)} crossing(s):")
    for c in crossings:
        print(f"  r = {c:.4f}")

    # --- Step 5: plot ---
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(
        scan_r,
        scan_dnll,
        "k-",
        linewidth=2,
        label="Likelihood scan ($2\\Delta\\mathrm{NLL}$)",
    )
    ax.plot(
        r_fine,
        quantile_fine,
        "-",
        color="0.4",
        linewidth=1.5,
        label="68% quantile (toys)",
    )
    ax.axhline(
        1.0,
        color="lightblue",
        linewidth=1.5,
        linestyle="-",
        label="$2\\Delta\\mathrm{NLL} = 1$",
    )

    for c in crossings:
        y_cross = float(scan_interp(c))
        ax.axvline(c, color="gray", linestyle=":", alpha=0.7)
        ax.plot(c, y_cross, "o", color="red", markersize=8, zorder=5)

    # Annotate best-fit value with uncertainties (stacked like plot1DScan)
    hi_str = f"+{err_hi:.3f}" if err_hi is not None else ""
    lo_str = f"-{err_lo:.3f}" if err_lo is not None else ""
    bestfit_text = f"$r = {bestfit:.3f}^{{{hi_str}}}_{{{lo_str}}}$"
    ax.text(
        0.80,
        0.95,
        bestfit_text,
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    ax.set_xlabel("r", fontsize=13)
    ax.set_ylabel("$2\\Delta\\mathrm{NLL}$", fontsize=13)
    ax.set_title("Likelihood scan vs toy-based 68% quantile")
    ax.legend(fontsize=11, loc="upper left")
    ax.set_ylim(bottom=0)

    outpath = os.path.join(base_dir, "scan_vs_quantile.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {outpath}")

    # Also save as pdf
    outpath_pdf = os.path.join(base_dir, "scan_vs_quantile.pdf")
    fig.savefig(outpath_pdf, bbox_inches="tight")
    print(f"Saved plot to {outpath_pdf}")


if __name__ == "__main__":
    main()
