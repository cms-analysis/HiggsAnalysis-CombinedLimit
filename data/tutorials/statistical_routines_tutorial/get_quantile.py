from ROOT import TFile
import numpy as np
import scipy.stats as st

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--input", default="higgsCombineTest.MultiDimFit.mH120.123456.root", help="input root file with toy results from MultiDimFit")
parser.add_argument("--q0", action="store_true", help="use the q_0 test statistic rather than the profile likelihood ratio.")
args = parser.parse_args()

n_sigma = 1
quantile_val = 2 * st.norm().cdf(-n_sigma)  # Get the quantile corresponding to the N sigma interval

f = TFile(args.input, "READ")
limit = f.Get("limit")
n_entries = limit.GetEntries()

m2nll_vals = []
r_vals = []
last_toy_num = -1
for i in range(n_entries):
    limit.GetEntry(i)
    if limit.quantileExpected < 0:
        if args.q0:
            r_vals.append(limit.r)
        continue

    m2nll_vals.append(2 * limit.deltaNLL)


test_stat_vals = m2nll_vals
if args.q0:
    test_stat_vals = np.where(np.array(r_vals) > 0, test_stat_vals, 0)

test_stat_cutoff = np.quantile(test_stat_vals, 1 - quantile_val)
t_stat_name = "q0" if args.q0 else "-2*deltaNLL"
print(f"This point is rejected at the {n_sigma} sigma level if the test stat {t_stat_name}  > {test_stat_cutoff}")
