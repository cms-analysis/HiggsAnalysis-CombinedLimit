#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit, modules
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
argv.remove( '-b-' )

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ModelTools import *
from HiggsAnalysis.CombinedLimit.ShapeTools import *
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

parser = OptionParser(usage="usage: %prog [options] datacard.txt -o output \nrun with --help to get list of options")
addDatacardParserOptions(parser)
parser.add_option("-P", "--physics-model", dest="physModel", default="HiggsAnalysis.CombinedLimit.PhysicsModel:defaultModel",  type="string", help="Physics model to use. It should be in the form (module name):(object name)")
parser.add_option("--PO", "--physics-option", dest="physOpt", default=[],  type="string", action="append", help="Pass a given option to the physics model (can specify multiple times)")
parser.add_option("", "--dump-datacard", dest="dumpCard", default=False, action='store_true',  help="Print to screen the DataCard as a python config and exit")
parser.add_option("--just-check-physics-model", dest="justCheckPhysicsModel", default=False, action='store_true',  help="Just check if the physics model is ok, without building the workspace.")
parser.add_option("--remove-multipdf", dest="removeMultiPdf", default=False, action='store_true',  help="Swap multipdf pdfs with their current index pdf")
(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_usage()
    exit(1)

options.fileName = args[0]
if options.fileName.endswith(".gz"):
    import gzip
    file = gzip.open(options.fileName, "rb")
    options.fileName = options.fileName[:-3]
else:
    file = open(options.fileName, "r")

## Parse text file 
DC = parseCard(file, options)

if options.dumpCard:
    DC.print_structure()
    exit()

## Load tools to build workspace
MB = None
if DC.hasShapes:
    MB = ShapeBuilder(DC, options)
else:
    MB = CountingModelBuilder(DC, options)

## Load physics model
(physModMod, physModName) = options.physModel.split(":")
__import__(physModMod)
mod = modules[physModMod]
physics = getattr(mod, physModName)
if mod     == None: raise RuntimeError, "Physics model module %s not found" % physModMod
if physics == None or not isinstance(physics, PhysicsModelBase): 
    raise RuntimeError, "Physics model %s in module %s not found, or not inheriting from PhysicsModelBase" % (physModName, physModMod)
physics.setPhysicsOptions(options.physOpt)
## Attach to the tools, and run
MB.setPhysics(physics)
MB.doModel(justCheckPhysicsModel=options.justCheckPhysicsModel)
