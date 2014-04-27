#!/usr/bin/env python
import math
from sys import argv
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

from HiggsAnalysis.CombinedLimit.DatacardParser import Datacard
from HiggsAnalysis.CombinedLimit.ModelTools import ModelBuilder
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from HiggsAnalysis.CombinedLimit.PhysicsModel import SM_HIGG_DECAYS as SM_HIGGS_DECAYS

parser = OptionParser(usage="usage: %prog [options] rvrf.root -o output \nrun with --help to get list of options")
parser.add_option("-v", "--verbose",  dest="verbose",  default=0,  type="int",      help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")
parser.add_option("-m", "--mass",     dest="mass",     default=125.0,  type="float",  help="Higgs mass to use.")
parser.add_option("-d" ,"--decays",   dest="decays",   default=None,  type="string", help="decays to include (comma separated, if unspecified use all)")
parser.add_option("-f" ,"--format",   dest="fmt",   default="5.3f",  type="string", help="Format string (default: 5.3f)")
parser.add_option("-k" ,"--kind",     dest="kind",   default="lnN",  type="string", help="Kind of output: lnN for log-normal kappas, % for percent, eps for bare values ")
parser.add_option("--mp", "--merge-param", dest="mergeParam", default=False, action="store_true", help="Merge parametric uncertainties into a single entry")
parser.add_option("-t", "--threshold",     dest="negligible",     default=-1.0,  type="float",  help="Threshold below which to ignore a number")
(options, args) = parser.parse_args()

decays = [ x.strip() for x in options.decays.split(",") ] if options.decays != None else SM_HIGGS_DECAYS

## set up some dummy options for the sake of the ModelBuilder
options.bin = True
options.fileName = "dummy.txt"
options.out      = "dummy.root"
options.cexpr = False
## and create a model builder
DC = Datacard()
MB = ModelBuilder(DC, options)
MB.doVar("MH[%g,110,140]" % options.mass)
SMH = SMHiggsBuilder(MB)

widthUncertainties = {}; widthUncertaintiesKeys = []
for line in open(SMH.brpath+"/WidthUncertainties_126GeV.txt"):
    if widthUncertaintiesKeys == []:
        widthUncertaintiesKeys = line.split()[1:]
    else:
        fields = line.split()
        widthUncertainties[fields[0]] = dict([(k,0.01*float(v)) for (k,v) in zip(widthUncertaintiesKeys, fields[1:])]) 

BRs = {}
for d in SM_HIGGS_DECAYS:
    SMH.makeBR(d)
    BRs[d] = MB.out.function("SM_BR_"+d).getVal()

THU_GROUPS = [
   ('hvv' , [ 'hww', 'hzz' ] ),
   ('hqq' , [ 'hbb', 'hcc', 'hss' ] ),
   ('hll' , [ 'htt', 'hmm' ] ),
   ('hgg' , [ 'hgg' ] ),
   ('hzg' , [ 'hzg' ] ),
   ('hgluglu' , [ 'hgluglu' ] ),
]

uncertainties = []
for s in widthUncertaintiesKeys[:-1]:
    key = "HiggsDecayWidth_"+s
    val = {}
    for d in decays:
        effect = (1.0 - BRs[d]) * widthUncertainties[d][s]
        for d2 in SM_HIGGS_DECAYS:
            if d2 == d: continue
            effect -= BRs[d2] * widthUncertainties[d2][s]
        val[d] = effect
    uncertainties.append( (key, val) )
if options.mergeParam:
    key = "HiggsDecayWidth_param"
    val = {}
    for d in decays:
        effect = math.sqrt(sum([v[d]**2 for (k,v) in uncertainties]))
        if d in ["hbb","hcc"]: effect *= -1
        val[d] = effect
    uncertainties = [ (key,val) ]
for s,decs in THU_GROUPS:
    key = "HiggsDecayWidth_"+s+"THU"
    val = {}
    for d in decays:
        effect = 0;
        if d in decs:
            #print "decay %s takes effect of %sTHU at numerator" % (d,s)
            effect += (1.0 - BRs[d]) * widthUncertainties[d]['thu']
        for d2 in SM_HIGGS_DECAYS:
            if d2 != d and d2 in decs: 
                #print "decay %s takes effect of %sTHU at denominator" % (d2,s)
                effect -= BRs[d2] * widthUncertainties[d2]['thu']
        val[d] = effect
    uncertainties.append( (key, val) )

print "%-30s   " % "Nuisance", "  ".join(["%-7s" % d for d in decays])
print "%-30s   " % "-----------------", "  ".join(["%-7s" % "-----" for d in decays])
for k,v in uncertainties:
    print "%-30s  " % k, 
    for d in decays:
        effect = v[d]
        if abs(effect) < options.negligible:
            print " %-7s" % (" -"),
        else:
            if options.kind == "lnN": effect = math.exp(effect)
            elif options.kind == "%": effect *= 100.0
            print " %-7s" % (("%"+options.fmt) % effect),
    print ""

