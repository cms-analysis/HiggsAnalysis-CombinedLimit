#!/usr/bin/env python3
# Script for batch sumbission of Asymptotic Limit Calculator
from __future__ import absolute_import, print_function

import commands
import glob
import os
import sys
from optparse import OptionParser

import numpy

parser = OptionParser()
parser.add_option("-w", "--workspace", dest="workspace")
parser.add_option("-m", "--mass", dest="mass", default=120.0, type="float")
parser.add_option("-n", "--npoints", dest="npoints", default=20, type="int")
parser.add_option("-r", "--range", dest="range", nargs=2, type="float", default=(0.0, 2.0))
parser.add_option("", "--nCPU", dest="nCPU", type="int", default=1)
parser.add_option("-o", "", dest="outputname", default="")
parser.add_option("-O", "", dest="options", default="")

(options, args) = parser.parse_args()

ws = options.workspace
mass = options.mass
outputname = options.outputname
print(options.range)
step = (options.range[1] - options.range[0]) / options.npoints
points = numpy.arange(options.range[0], options.range[1] + step, step)

workingDir = os.getcwd()
wslocation = workingDir + "/" + ws

jobname = "limitgrid_%.1f%s" % (mass, outputname)
f = open("%s.sh" % (jobname), "w")

# Create job script
f.write("#!/bin/bash\n")
f.write("set -x\n")
f.write("cd %s\n" % workingDir)
f.write("eval `scramv1 runtime -sh`\n")

f.write("cd -\n")
f.write("cp -p $CMSSW_BASE/lib/${CMSSW_ARCH}/libHiggsAnalysisCombinedLimit.{so,dylib} . \n")
f.write("mkdir scratch\n")
f.write("cd scratch\n")
f.write("cp -p %s . \n" % (wslocation))

# Loop over number of points, and add wait after npoints/nCPU jobs
ext = " & " if options.nCPU > 1 else " "
for i, po in enumerate(points):
    f.write("combine -M Asymptotic %s -m %.1f --singlePoint %f -n %.1f %s %s \n" % (ws, mass, po, po, ext, options.options))
    if i > 0 and i % options.nCPU == 0:
        f.write("wait\n")
print(i)
if options.nCPU > i:
    f.write("wait\n")
f.write("hadd -f grid_%.1f%s.root higgsCombine* \n" % (mass, outputname))
f.write("cp -f grid_%.1f%s.root %s \n" % (mass, outputname, workingDir))
f.write("echo 'DONE' \n")

# make it executable
os.system("chmod 755 %s.sh" % (jobname))
