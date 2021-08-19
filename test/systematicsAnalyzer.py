#!/usr/bin/env python
import re
import os.path
from math import *
from optparse import OptionParser
import sys 

parser = OptionParser()
parser.add_option("-f", "--format",  type="string",   dest="format", default="html", help="Format for output number (choose html or brief)")
parser.add_option("-m", "--mass",    dest="mass",     default=0,  type="float",  help="Higgs mass to use. Will also be written in the Workspace as RooRealVar 'MH'.")
parser.add_option("-p", "--process",    dest="process",     default=None,  type="string",  help="Higgs process to use. Will also be written in the Workspace as RooRealVar 'MH'.")
parser.add_option("-D", "--dataset", dest="dataname", default="data_obs",  type="string",  help="Name of the observed dataset")
parser.add_option("-s", "--search", "--grep", dest="grep", default=[], action="append",  type="string",  help="Selection of nuisance parameters (regexp, can be used multiple times)")
parser.add_option("-a", "--all", dest="all", default=False,action='store_true',  help="Report all nuisances (default is only lnN)")
parser.add_option("", "--t2w", dest="t2w", default=False,action='store_true',  help="Run text2workspace -- only effects datacards with param lines")
parser.add_option("", "--noshape", dest="noshape", default=False,action='store_true',  help="Counting experiment only (alternatively, build a shape analysis from combineCards.py -S card.txt > newcard.txt )")
(options, args) = parser.parse_args()
options.stat = False
options.bin = True # fake that is a binary output, so that we parse shape lines
options.out = "tmp.root"
options.fileName = args[0]
options.cexpr = False
options.fixpars = False
options.libs = []
options.verbose = 0
options.poisson = 0
options.nuisanceFunctions = []
options.nuisanceGroupFunctions = []
options.noOptimizePdf=False
options.optimizeBoundNuisances=False
options.useHistPdf=False
options.optimizeExistingTemplates=False
options.packAsymPows=False
options.noBOnly=True
options.moreOptimizeSimPdf="none"
options.doMasks=False
options.defMorph = "shape"
options.nuisancesToRescale = ""
options.nuisancesToExclude = []
options.noJMax = True
options.allowNoSignal = True
options.modelparams = [] 
options.optimizeMHDependency = False 
options.optimizeTemplateBins=False
options.forceNonSimPdf = False
options.physModel = "HiggsAnalysis.CombinedLimit.PhysicsModel:defaultModel"
# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
import sys
sys.argv = [ '-b-']
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ShapeTools     import *
if options.fileName.endswith(".gz"):
    import gzip
    file = gzip.open(options.fileName, "rb")
    options.fileName = options.fileName[:-3]
else:
    file = open(options.fileName, "r")
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

DC = parseCard(file, options)

if not DC.hasShapes: DC.hasShapes = True
MB = ShapeBuilder(DC, options)
if not options.noshape: 
  MB.prepareAllShapes()

MODELBUILT=False
def buildModel():
  ## Load physics model
  (physModMod, physModName) = options.physModel.split(":")
  __import__(physModMod)
  mod = sys.modules[physModMod]
  physics = getattr(mod, physModName)
  ## Attach to the tools, and run
  MB.setPhysics(physics)
  MB.doModel(False)
  return True

def commonStems(list, sep="_"):
    hits = {}
    for item in list:
        base = ""
        for token in item.split(sep):
            base = base + sep + token if base else token
            if base not in hits: hits[base] = 0
            hits[base] += 1
    veto = {}
    for k,v in hits.iteritems():
        pieces = k.split(sep)
        for i in xrange(1, len(pieces)):
            k2 = "_".join(pieces[:-i])
            if hits[k2] == v: 
                veto[k2] = True
            else:
                veto[k] = True
    ret = []
    for k,v in hits.iteritems():
       if k not in veto: ret.append((k,v))
    ret.sort()
    return ret 

def addTo(dic,key,l):
  if key in dic.keys(): 
    dic[key].append(l)
  else:
    dic[key]=[l]

if options.t2w: 
	buildModel()
	MODELBUILT=True
	options.all=True

report = {}; errlines = {}; outParams = {}; check_list = {}; 
seen_systematics = []
for (lsyst,nofloat,pdf,pdfargs,errline) in DC.systs:
    if ("rateParam" in pdf) or ("discrete" in pdf): 
         if options.all: outParams[lsyst]=[pdf,pdfargs]
    if not options.t2w and ("param" in pdf): outParams[lsyst]=[pdf,pdfargs]

    if not options.all and pdf != "lnN": continue
    if not options.t2w and "param" in pdf : continue 
    if "param" in pdf:
      if not lsyst in seen_systematics: 
        if not len(errline): errline = {b:{p:0 for p in DC.exp[b].iterkeys() } for b in DC.bins}
      else: errline=errlines[lsyst]
    types = []
    minEffect, maxEffect = 999.0, 1.0
    processes = {}
    channels  = []
    errlines[lsyst] = errline
    vals = []
    if "param" in pdf and not  MODELBUILT:continue 
    for b in DC.bins:
        numKeysFound = 0
        channels.append(b)
        for p in DC.exp[b].iterkeys():
	    if lsyst in check_list.keys():
	      if [p,b] in check_list[lsyst]: continue
            if (not pdf=="param") and errline[b][p] == 0: continue
            if pdf == "gmN":
               numKeysFound+=1
               minEffect = pdfargs[0] 
               maxEffect = pdfargs[0]
               processes[p] = True
	       addTo(check_list,lsyst,[p,b])
            elif pdf == "param":
	       if not MODELBUILT: continue 
               formula = "n_exp_final_bin%s_proc_%s"%(b,p)
               if not MB.out.function(formula) : formula = "n_exp_bin%s_proc_%s"%(b,p)
	       if not MB.out.function(formula) : sys.exit("No formula %s"%formula)
               if not (MB.out.function(formula).getParameters(ROOT.RooArgSet())).contains(MB.out.var(lsyst)): continue
               centralVal = float(pdfargs[0])
               if "/" in pdfargs[1]: 
                 minError, maxError = float(pdfargs[1].split("/")[0]),float(pdfargs[1].split("/")[1])
               else: 
                 minError, maxError =  -1*float(pdfargs[1]),float(pdfargs[1])
               MB.out.var(lsyst).setVal(centralVal)
               centralNorm =  MB.out.function(formula).getVal()
	       if centralNorm<=0: continue
               MB.out.var(lsyst).setVal(centralVal+minError)
               lowNorm =  MB.out.function(formula).getVal()
               MB.out.var(lsyst).setVal(centralVal+maxError)
               highNorm =  MB.out.function(formula).getVal()
	       errlines[lsyst][b][p] = "%.3f/%.3f (param) "%(lowNorm/centralNorm, highNorm/centralNorm)
	       vals.append(lowNorm/centralNorm) 
	       vals.append(highNorm/centralNorm) 
               MB.out.var(lsyst).setVal(centralVal)
               numKeysFound+=1
	       types.append(pdf)
               processes[p] = True
	       addTo(check_list,lsyst,[p,b])

            elif "shape" in pdf and MB.isShapeSystematic(b,p,lsyst):

	      if errline[b][p]==0: continue
              systShapeName = lsyst
              #vals = []
              if (lsyst,b,p) in DC.systematicsShapeMap.keys(): systShapeName = DC.systematicsShapeMap[(lsyst,b,p)]

              objU,objD,objC = MB.getShape(b,p,systShapeName+"Up"), MB.getShape(b,p,systShapeName+"Down"), MB.getShape(b,p)
             
              if objC.InheritsFrom("TH1"): valU,valD,valC =  objU.Integral(), objD.Integral(), objC.Integral()
              elif objC.InheritsFrom("RooDataHist"): valU,valD,valC =  objU.sumEntries(), objD.sumEntries(), objC.sumEntries()
              if valC!=0: 
                  errlines[lsyst][b][p] = "%.3f/%.3f (shape)"%(valD/valC,valU/valC)
                  vals.append(valD/valC)
                  vals.append(valU/valC)
              else: 
                  errlines[lsyst][b][p] = "NAN/NAN"
                  vals.append(1.)
                  vals.append(1.)
              numKeysFound+=1
              types.append(pdf)
              processes[p] = True
	      addTo(check_list,lsyst,[p,b])

            else:
                vals.extend(errline[b][p] if type(errline[b][p]) == list else [ errline[b][p] ])
                numKeysFound+=1
                types.append(pdf)
                processes[p] = True
	        addTo(check_list,lsyst,[p,b])
            for val in vals:
                if val < 1: val = 1.0/val
                minEffect = min(minEffect, val)
                maxEffect = max(maxEffect, val)
        if numKeysFound == 0 : channels.remove(b)
    #if no effect just skip 
    if not len(vals): continue 
    channelsShort = commonStems(channels)
    types = set(types)
    types = ",".join(types)
    if lsyst in report.keys():
       report[lsyst]['channels'].extend(channelsShort)
       report[lsyst]['bins'].extend(channels)
       report[lsyst]['processes'].extend(processes)
       
       report[lsyst]['channels'] = list(set(report[lsyst]['channels']))
       report[lsyst]['bins'] = list(set(report[lsyst]['bins']))
       report[lsyst]['processes'] = list(sorted(set(report[lsyst]['processes'])))
       
       report[lsyst]['effect'] = ["%5.3f"%(min(float(report[lsyst]['effect'][0]),minEffect)),"%5.3f"%(max(float(report[lsyst]['effect'][1]),maxEffect))]
       if types not in report[lsyst]['types']: report[lsyst]['types']+=","+types
    else: report[lsyst] = { 'channels':channelsShort, 'bins' : channels, 'processes': sorted(processes.keys()), 'effect':["%5.3f"%minEffect,"%5.3f"%maxEffect], 'types':types }
    seen_systematics.append(lsyst)

# Get list
names = report.keys() 
if "brief" in options.format:
    names = [ k for (k,v) in report.iteritems()  ]
if options.process:
    names = [ k for k in names if any(p for p in report[k]['processes'] if re.match(options.process, p)) ]
if options.grep:
    names = [ n for n in names if any(p for p in options.grep if re.match(p,n)) ]

# alphabetic sort
names.sort()
# now re-sort by category (preserving alphabetic sort inside)
namesCommon = [ n for n in names if re.match(r"(pdf_|QCD|lumi|UE|BR).*", n) ]
namesCMS1   = [ n for n in names if re.match(r"CMS_(eff|scale|fake|res|trig).*", n) ]
namesCMS2   = [ n for n in names if re.match(r"CMS_.*", n) and n not in namesCMS1 ]
namesRest   = [ n for n in names if n not in namesCommon and n not in namesCMS1 and n not in namesCMS2 ]
names = namesCommon + namesCMS1 + namesCMS2 + namesRest

if "html" in options.format:
    print """
<html>
<head>
<style type="text/css">
body { font-family: 'Consolas', 'Courier New', courier, monospace; font-size: small; }
td, th { border-bottom: 1px solid black; padding: 1px 1em; vertical-align: top; }
td.channDetails { font-size: x-small; }
</style>
<script type="text/javascript">
function toggleChann(id) {
    if (document.getElementById(id+"_chann_toggle").innerHTML == "[+]") {
        document.getElementById(id+"_chann").style = "";
        document.getElementById(id+"_chann_toggle").innerHTML = "[-]";
    } else {
        document.getElementById(id+"_chann").style = "display: none";
        document.getElementById(id+"_chann_toggle").innerHTML = "[+]";
    }
}
</script>
<title>Nuisance Report</title>
</head><body>
<h1>Nuisance Report</h1>
All numbers shown report the +/- 1-sigma variation in the yield for each affected channel/process. The Range shows the minimum and maximum effects across all channels/processes.
%s 
<table>
<tr><th>Nuisance (types)</th><th colspan="2">Range</th><th>Processes</th><th>Channels</th></tr>
"""%("You didn't run with the option --t2w so param types will only show the line from the datacard" if not options.t2w else "")
    for nuis in names:
        val = report[nuis]
        print "<tr><td><a name=\"%s\"><b>%s</b></a></td>" % (nuis,nuis+"  ("+val['types']+")")
        #print "<td>%5.3f</td><td>%5.3f</td>" % ( val['effect'][0],val['effect'][1] )
        print "<td>%s</td><td>%s</td>" % ( val['effect'][0],val['effect'][1] )
        print "<td>",", ".join(val['processes']), "</td>"
        print "<td>%s" % ", ".join(["%s(%d)" % (k,v) for (k,v) in sorted(val['channels'])]),
        print "<a id=\"%s_chann_toggle\" href=\"#%s\" onclick=\"toggleChann(&quot;%s&quot;)\">[+]</a></td>" % (nuis,nuis,nuis)
        print "</tr>"
        print "<tr id=\"%s_chann\" style=\"display: none\">" % nuis
        print "\t<td colspan=\"5\"><table class=\"channDetails\">" 
        for x in sorted(val["bins"]): print "\t\t<tr><td>%s</td><td>%s</td></li>" % (x, ", ".join(["%s(%s)"%(k,v) for (k,v) in errlines[nuis][x].iteritems() if v != 0]))
        print "\t</table></td>"
        print "</tr>\n"
    for x in outParams.keys():
        print "\t\t<tr><td><b>%s(%s)</b></td><td>%s</td></li>" % (x,  outParams[x][0] , ", ".join([a for a in outParams[x][1]]))
        print "</tr>\n"
    print """
</table>
</body>
</html>"""
else:
    if "brief" in options.format:
        print "%-50s  [%5s, %5s]   %-40s  %-30s" % ("   NUISANCE (TYPE)", " MIN", " MAX", "PROCESSES", "CHANNELS(#SUBCHANNELS)" )
        print "%-50s  %14s   %-40s  %-30s" % ("-"*50, "-"*14, "-"*30, "-"*30)
        for nuis in names:
            val = report[nuis]
            print "%-50s (%s)  [%s, %s]   %-40s  %-30s" % ( nuis,val['types'], val['effect'][0],val['effect'][1], 
                                                                ",".join(val['processes']),
                                                                ",".join(["%s(%d)" % (k,v) for (k,v) in sorted(val['channels'])]))
            
