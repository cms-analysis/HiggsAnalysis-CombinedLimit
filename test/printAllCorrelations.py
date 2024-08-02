#!/usr/bin/env python3

from __future__ import absolute_import, print_function

import argparse

from six.moves import range

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input mlfit file")
parser.add_argument("-m", "--max", type=int, default=10, help="max number to print")

args = parser.parse_args()

fin = ROOT.TFile(args.input.split(":")[0])
rfr = fin.Get(args.input.split(":")[1])
# rfr.Print()

all_pars = rfr.floatParsFinal()
correlations = []
for i in range(all_pars.getSize()):
    par = all_pars.at(i)
    for j in range(all_pars.getSize()):
        if j < i:
            par_corr = all_pars.at(j)
            correlations.append(
                (
                    par.GetName(),
                    par_corr.GetName(),
                    rfr.correlation(par.GetName(), par_corr.GetName()),
                )
            )

correlations.sort(key=lambda x: abs(x[2]), reverse=True)
# print correlations

print("covQual: %i" % rfr.covQual())
for i in range(min(len(correlations), args.max)):
    print("%i %-50s %-50s  %+.3f" % (i, correlations[i][0], correlations[i][1], correlations[i][2]))
