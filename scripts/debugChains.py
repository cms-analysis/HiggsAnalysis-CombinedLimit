#!/usr/bin/env python3
import matplotlib.pyplot as plt

# import HiggsAnalysis.CombinedLimit.util.plotting as plot
import argparse
import sys
import numpy as np
import ROOT

# recent versions of numpy complain but can ignore
import warnings

warnings.filterwarnings("ignore", message="The value of the smallest subnormal for")

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()
parser.add_argument("input", help="""Input root file""")
parser.add_argument(
    "--output",
    "-o",
    default="",
    help="""Name of the output
    plot without file extension""",
)
parser.add_argument(
    "--param",
    "-p",
    default="r",
    help="""Name of the parameter to plot""",
)
parser.add_argument(
    "--range",
    "-r",
    nargs=2,
    type=float,
    help="""Range of the parameter to plot""",
)
parser.add_argument(
    "--nbins",
    type=int,
    default=100,
    help="""Number of bins to use in the histogram""",
)
parser.add_argument(
    "--CL",
    type=float,
    default=0.95,
    help="""Confidence level to use for the interval""",
)
parser.add_argument(
    "--burnInFraction",
    "-b",
    type=float,
    default=0.5,
    help="""Fraction of the chain to use as burn-in""",
)
parser.add_argument(
    "--chainsF",
    type=float,
    default=0.1,
    help="""Fraction of the chains to plot in the trace plot (doesn't effect the histogram)""",
)
parser.add_argument(
    "--mode",
    default="interval",
    help="""Choose from interval (default)/upperlim/lowerlim""",
)

modes = ["upperlim", "lowerlim", "interval"]

args = parser.parse_args()

if args.mode not in modes:
    print(f"ERROR: for --mode, must pick from {modes}")
    sys.exit(0)


def weighted_percentile(data, weights, perc):
    # Source - https://stackoverflow.com/a/61343915
    # Posted by imbr, modified by community. See post 'Timeline' for change history
    # Retrieved 2026-03-03, License - CC BY-SA 4.0
    # based off https://en.wikipedia.org/wiki/Percentile#Definition_of_the_Weighted_Percentile_method

    ix = np.argsort(data)
    data = data[ix]  # sort data
    weights = weights[ix]  # sort weights
    cdf = (np.cumsum(weights) - 0.5 * weights) / np.sum(weights)  # 'like' a CDF function
    return np.interp(perc, cdf, data)


def findInterval(arr, weights, CL, mode="interval"):
    # have an array of values and a CL
    # find the interval that contains CL of them

    if mode == "interval":
        left_q = (1 - CL) / 2
        right_q = 1 - left_q

    elif mode == "upperlim":
        left_q = 0
        right_q = CL

    elif mode == "lowerlim":
        left_q = 1 - CL
        right_q = 1

    q1 = weighted_percentile(arr, weights, left_q)
    q2 = weighted_percentile(arr, weights, right_q)
    return [q1, q2]


fi_MCMC = ROOT.TFile.Open(args.input)

# array to hold the values of the parameter in the chain
param_value_chunks = []
param_weight_chunks = []

graphs = []
all_graphs = []
kept_chain = 0

average_chain_length = 0
total_chains = 0

for j, k in enumerate(fi_MCMC.Get("toys").GetListOfKeys()):
    mychain = k.ReadObj().GetAsDataSet()
    lenchain = mychain.numEntries()
    average_chain_length += lenchain
    burnin = int(lenchain * args.burnInFraction)
    start_idx = burnin + 1
    total_chains += 1

    keep_gr = np.random.uniform(0, 1) < args.chainsF or j == 0
    chain_vals = np.empty(lenchain, dtype=float)
    chain_weights = np.empty(lenchain, dtype=float)

    for i in range(0, lenchain):
        row = mychain.get(i)
        chain_vals[i] = row.getRealValue(args.param)
        chain_weights[i] = mychain.weight()

    if start_idx < lenchain:
        param_value_chunks.append(chain_vals[start_idx:])
        param_weight_chunks.append(chain_weights[start_idx:])

    graphs.append(keep_gr)
    if keep_gr:
        kept_chain += 1

    all_graphs.append(chain_vals)

if param_value_chunks:
    param_values = np.concatenate(param_value_chunks)
    param_weights = np.concatenate(param_weight_chunks)
else:
    param_values = np.array([])
    param_weights = np.array([])

average_chain_length = float(average_chain_length) / (total_chains)

plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(1, 2, figsize=(14, 5))

param_values = np.array(param_values, dtype=float)
ax[0].hist(param_values, density=True, color="black", bins=args.nbins, range=args.range, weights=param_weights, histtype="step")
ax[0].set_xlabel(args.param)
ax[0].set_ylabel("Posterior probability density")

interval = findInterval(param_values, param_weights, args.CL, args.mode)
print(f"Average chain length: {average_chain_length:.1f}")
print(f"Number of chains: {j + 1}")
print(f"Burn-in fraction: {args.burnInFraction:.2f} (average burn-in length: {average_chain_length * args.burnInFraction:.1f} entries)")

if args.mode == "interval":
    label = f"{args.CL * 100:.1f}% CL interval: {interval[0]:.3f} < {args.param} < {interval[1]:.3f}"
elif args.mode == "upperlim":
    label = f"{args.CL * 100:.1f}% CL interval: {args.param} < {interval[1]:.3f}"
elif args.mode == "lowerlim":
    label = f"{args.CL * 100:.1f}% CL interval: {args.param} > {interval[0]:.3f}"

print(label)
ax[0].axvline(interval[0], color="red", linestyle="--", label=label)
ax[0].axvline(interval[1], color="red", linestyle="--")
# put legend in the upper left corner, outside the plot
ax[0].legend(loc="upper left", bbox_to_anchor=(0.01, 1.14))

for k, gr in enumerate(all_graphs):
    if graphs[k]:
        ax[1].plot(np.arange(len(gr)), gr, color="black", marker=None, linestyle="-", linewidth=0.2, alpha=0.4)

ax[1].set_ylabel(args.param)
ax[1].axvline(args.burnInFraction * average_chain_length, color="blue", linestyle="--", label="Burn-in fraction")
ax[1].set_xlabel("Chain index")
ax[1].set_title(f"Trace plot of {kept_chain} chains / {j + 1} chains")

# make a sliding  average and 68% interval plot on top of the trace plot
# this should be across the graphs and take ~5% of the average chain length as the window size
window_size = int(average_chain_length * 0.05)
num_windows = int(average_chain_length / window_size)
running_avg = np.empty(num_windows)
running_avg_upper = np.empty(num_windows)
running_avg_lower = np.empty(num_windows)
window_centers = np.empty(num_windows)

for i in range(num_windows):
    window_vals = []
    window_weights = []
    for gr in all_graphs:
        if (i + 1) * window_size > len(gr):
            continue
        window_vals.append(gr[i * window_size : (i + 1) * window_size])
        window_weights.append(np.ones(window_size))  # equal weights for the running average
    window_vals = np.concatenate(window_vals)
    window_weights = np.concatenate(window_weights)
    running_avg[i] = np.average(window_vals, weights=window_weights)
    interval = findInterval(window_vals, window_weights, 0.68, mode="interval")
    running_avg_lower[i] = interval[0]
    running_avg_upper[i] = interval[1]
    window_center = (i * window_size) + window_size / 2
    window_centers[i] = window_center

ax[1].plot(window_centers, running_avg, color="red", marker=None, linestyle="-", linewidth=2, label="Sliding average")
ax[1].fill_between(window_centers, running_avg_lower, running_avg_upper, color="red", alpha=0.5, label="68% interval")
ax[1].legend(loc="upper right")

if args.range:
    ax[1].set_ylim(args.range[0], args.range[1])

if not args.output:
    args.output = args.param
plt.savefig(args.output + ".pdf")
plt.savefig(args.output + ".png")
print(f"Saved output as {args.output}.pdf/png")
