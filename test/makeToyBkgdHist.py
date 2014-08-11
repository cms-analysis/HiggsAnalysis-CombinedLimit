########
# Author: Jessica Brinson (jbrinson@cern.ch), The Ohio State Unversity
# Date created: 7 Aug, 2014
##############
# Script to generate root file with a histogram containing
# the simulated number of background  events for a given number of toys
##############
# Usage: First run Higgs Combine tool, for example:
# combine datacard.txt -S 0 -M GenerateOnly -t 10000 --saveToys
#
# The -M GenerateOnly option generates only toys without running the limits.
# Greatly decreases the processing time: for a single-bin counting experiment,
# can generate 10K toys in ~1 min
#
# Can use the --toysFreqentist to generate frequentist toys.  See:
# https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/572/1/2/1/1/1.html 
#
##############
# Usage: python parseCombine.py name_of_output_root_file_from_combine_command ntoys
##############
# 
#
#!/usr/bin/env python

import time
import os
import sys
import math
import copy
import re
import subprocess
import glob

from ROOT import TF1, TFile, TH2F, gROOT, gStyle,TH1F, TCanvas, TString, TLegend, TPaveLabel, TH2D, TPave, Double

if len(sys.argv) < 2:
    print "Must specify the root file that you wish to parse."
    sys.exit()
    
if len(sys.argv) < 3:
    print "Must specify the number of toys that you ran over."
    sys.exit()
    
inputFile = str(sys.argv[1])
nToys = int(sys.argv[2])

file = TFile.Open(inputFile)
file.cd()

for key in  file.GetListOfKeys():
    if (key.GetClassName() != "TDirectoryFile"):
        continue
    rootDirectory = key.GetName()

#loop over the toys and extract nToyBkgd
nToyBkgdList = []
for i in range (1,nToys+1):
    toyTemp = file.Get(rootDirectory + "/toy_" + str(i)).Clone()
    nToyBkgd = toyTemp.get(0).getRealValue("n_obs_binMyChan")
    nToyBkgdList.append(nToyBkgd)

# fill histogram with the number of simulated background events for each toy
h = TH1F("nToyBkgd", ";background yield;number of toys", 100, float(min(nToyBkgdList)), float(max(nToyBkgdList)))
for item in nToyBkgdList:
    h.Fill(float(item))

    outputFileName = (inputFile).rstrip(".root")
    outputFile = TFile(outputFileName + "_histToyNBkgd.root", "RECREATE")
h.Write()
outputFile.Write()
outputFile.Close()

print "Finished writing output to " + outputFile.GetName() 
