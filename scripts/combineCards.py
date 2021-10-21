#!/usr/bin/env python
import re
from sys import argv,exit
import os.path
from pprint import pprint
from optparse import OptionParser
parser = OptionParser(
    usage="%prog [options] [label=datacard.txt | datacard.txt]",
    epilog="The label=datacard.txt syntax allows to specify the label that channels from datacard.txt will have in the combined datacard. To combine cards with different energies one can use dc_7TeV=datacard7.txt dc_8TeV=datacard8.txt (avoid using labels starting with numbers)."
    )
parser.add_option("-s", "--stat",   dest="stat",          default=False, action="store_true", help="Drop all systematics")
parser.add_option("-S", "--force-shape", dest="shape",    default=False, action="store_true", help="Treat all channels as shape analysis. Useful for mixed combinations")
parser.add_option("-a", "--asimov", dest="asimov",  default=False, action="store_true", help="Replace observation with asimov dataset. Works only for counting experiments")
parser.add_option("-P", "--prefix", type="string", dest="fprefix", default="",  help="Prefix this to all file names")
parser.add_option("--xc", "--exclude-channel", type="string", dest="channelVetos", default=[], action="append", help="Exclude channels that match this regexp; can specify multiple ones")
parser.add_option("--ic", "--include-channel", type="string", dest="channelIncludes", default=[], action="append", help="Only include channels that match this regexp; can specify multiple ones")
parser.add_option("--X-no-jmax",  dest="noJMax", default=False, action="store_true", help="FOR DEBUG ONLY: Turn off the consistency check between jmax and number of processes.")
parser.add_option("--xn-file", "--exclude-nuisances-from-file", type="string", dest="nuisVetoFile", help="Exclude all the nuisances in this file")
parser.add_option("--en-file", "--edit-nuisances-from-file", type="string", dest="editNuisFile", help="edit the nuisances in this file")

(options, args) = parser.parse_args()
options.bin = True # fake that is a binary output, so that we parse shape lines
options.nuisancesToExclude = []
options.verbose = 0
options.allowNoSignal = True
options.allowNoBackground = True
options.evaluateEdits = False

if options.nuisVetoFile:
    for line in open(options.nuisVetoFile,"r"):
        options.nuisancesToExclude.append(re.compile(line.strip()))

from HiggsAnalysis.CombinedLimit.DatacardParser import *

obsline = []; obskeyline = [] ;
keyline = []; expline = []; systlines = {}
signals = []; backgrounds = []; shapeLines = []
paramSysts = {}; flatParamNuisances = {}; discreteNuisances = {}; groups = {}; rateParams = {}; rateParamsOrder = set();
extArgs = {}; binParFlags = {}
nuisanceEdits = [];

def compareParamSystLines(a,b):
  if float(a[0])!=float(b[0]): return False
  if "/" in a[1]:
    if "/" not in b[1] : return False
    a1,a2 = a[1].split("/")
    b1,b2 = b[1].split("/")
    if float(a1)!=float(b1): return False
    if float(a2)!=float(b2): return False
  else:
    if float(a[1])!=float(a[1]): return False
  return True


cmax = 5 # column width
if not args:
    raise RuntimeError, "No input datacards specified."
for ich,fname in enumerate(args):
    label = "ch%d" % (ich+1)
    if "=" in fname: (label,fname) = fname.split("=")
    fname = options.fprefix+fname
    dirname = os.path.dirname(fname)
    if fname.endswith(".gz"):
        import gzip
        file = gzip.open(fname, "rb")
        fname = fname[:-3]
    else:
        file = open(fname, "r")
    DC = parseCard(file, options)
    singlebin = (len(DC.bins) == 1)
    if label == ".":
        label=DC.bins[0] if singlebin else "";
    elif not singlebin:
        label += "_";
    for b in DC.bins:
        bout = label if singlebin else label+b
        b_in  = label if singlebin else b
        if isVetoed(b_in,options.channelVetos): continue
        if not isIncluded(b_in,options.channelIncludes): continue
        obskeyline.append(bout)
        for (p,e) in DC.exp[b].items(): # so that we get only self.DC.processes contributing to this bin
            if DC.isSignal[p] == False: continue
            #print "in DC.exp.items:b,p", b,p
            expline.append("%s" % FloatToString(e)) if (e == 0 or e > 1e-3) else expline.append("%s" % FloatToStringScientific(e))
            keyline.append((bout, p, DC.isSignal[p]))
        for (p,e) in DC.exp[b].items(): # so that we get only self.DC.processes contributing to this bin
            if DC.isSignal[p]: continue
            #print "in DC.exp.items:b,p", b,p
            expline.append("%s" % FloatToString(e)) if (e == 0 or e > 1e-3) else expline.append("%s" % FloatToStringScientific(e))
            keyline.append((bout, p, DC.isSignal[p]))
    # systematics
    for (lsyst,nofloat,pdf,pdfargs,errline) in DC.systs:
        systeffect = {}
        if pdf == "param":
            if paramSysts.has_key(lsyst):
               #if paramSysts[lsyst] != pdfargs:
	       if not compareParamSystLines(paramSysts[lsyst],pdfargs) : raise RuntimeError, "Parameter uncerainty %s mismatch between cards, %g != %g" % lsyst
            else:
                paramSysts[lsyst] = pdfargs
            continue
        for b in DC.bins:
            bout = label if singlebin else label+b
            b_in  = label if singlebin else b
            if isVetoed(b_in,options.channelVetos): continue
            if not isIncluded(b_in,options.channelIncludes): continue
            if not systeffect.has_key(bout): systeffect[bout] = {}
            for p in DC.exp[b].keys(): # so that we get only self.DC.processes contributing to this bin
                r = str(errline[b][p]);
                if type(errline[b][p]) == list: r = "%s/%s" % (FloatToString(errline[b][p][0]), FloatToString(errline[b][p][1]))
                elif type in ("lnN",'gmM'): r = "%s" % FloatToString(errline[b][p])
                if errline[b][p] == 0: r = "-"
                if len(r) > cmax: cmax = len(r) # get max col length, as it's more tricky to do it later with a map
                systeffect[bout][p] = r
        if systlines.has_key(lsyst):
            (otherpdf, otherargs, othereffect, othernofloat) = systlines[lsyst]
            if otherpdf != pdf:
                if (pdf == "lnN" and otherpdf.startswith("shape")):
                    if systlines[lsyst][0][-1] != '?': systlines[lsyst][0] += '?'
                    for b,v in systeffect.items(): othereffect[b] = v;
                elif (pdf.startswith("shape") and otherpdf == "lnN"):
                    if pdf[-1] != '?': pdf += '?'
                    systlines[lsyst][0] = pdf
                    for b,v in systeffect.items(): othereffect[b] = v;
                elif (pdf == otherpdf+"?") or (pdf+"?" == otherpdf):
                    systlines[lsyst][0] = pdf.replace("?","")+"?"
                    for b,v in systeffect.items(): othereffect[b] = v;
                else:
                    raise RuntimeError, "File %s defines systematic %s as using pdf %s, while a previous file defines it as using %s" % (fname,lsyst,pdf,otherpdf)
            else:
                if pdf == "gmN" and int(pdfargs[0]) != int(otherargs[0]):
                    raise RuntimeError, "File %s defines systematic %s as using gamma with %s events in sideband, while a previous file has %s" % (fname,lsyst,pdfargs[0],otherargs[0])
                for b,v in systeffect.items(): othereffect[b] = v;
        else:
            pdfargs = [ str(x) for x in pdfargs ]
            systlines[lsyst] = [pdf,pdfargs,systeffect,nofloat]
    # flat params
    for K in DC.flatParamNuisances.iterkeys():
        flatParamNuisances[K] = True
    for K in DC.extArgs.keys():
        extArgs[K] = DC.extArgs[K]
    for K in DC.binParFlags.iterkeys():
        tbin = label if singlebin else label+K
        binParFlags[tbin] = DC.binParFlags[K]
    # rate params
    for K in DC.rateParams.iterkeys():
        tbin,tproc = K.split("AND")[0],K.split("AND")[1]
        b_in = tbin
        tbin = label if singlebin else label+tbin
        if isVetoed(b_in,options.channelVetos): continue
        if not isIncluded(b_in,options.channelIncludes): continue
        nK = tbin+"AND"+tproc
        rateParams[nK] = DC.rateParams[K]
        rateParamsOrder.update(DC.rateParamsOrder)
    # discrete nuisance
    for K in DC.discretes:
        if discreteNuisances.has_key(K): raise RuntimeError, "Cannot currently correlate discrete nuisances across categories. Rename %s in one."%K
        else: discreteNuisances[K] = True
    # put shapes, if available
    if len(DC.shapeMap):
        for b in DC.bins:
            bout = label if singlebin else label+b
            b_in  = label if singlebin else b
            if isVetoed(b_in,options.channelVetos): continue
            if not isIncluded(b_in,options.channelIncludes): continue
            p2sMap  = DC.shapeMap[b]   if DC.shapeMap.has_key(b)   else {}
            p2sMapD = DC.shapeMap['*'] if DC.shapeMap.has_key('*') else {}
            for p, x in p2sMap.items():
                xrep = [xi.replace("$CHANNEL",b) for xi in x]
                if xrep[0] != 'FAKE' and dirname != '' and not xrep[0].startswith("/"):
                    xrep[0] = dirname+"/"+xrep[0]
                shapeLines.append((p,bout,xrep))
            for p, x in p2sMapD.items():
                if p2sMap.has_key(p): continue
                xrep = [xi.replace("$CHANNEL",b) for xi in x]
                if xrep[0] != 'FAKE' and dirname != '' and not xrep[0].startswith("/"):
                    xrep[0] = dirname+"/"+xrep[0]
                shapeLines.append((p,bout,xrep))
    elif options.shape:
        for b in DC.bins:
            bout = label if singlebin else label+b
            shapeLines.append(('*',bout,['FAKE']))
    # combine observations, but remove line if any of the datacards doesn't have it
    if len(DC.obs) == 0:
        obsline = None
    elif obsline != None:
        for b in DC.bins:
            bout = label if singlebin else label+b
            b_in  = label if singlebin else b
            if isVetoed(b_in,options.channelVetos): continue
            if not isIncluded(b_in,options.channelIncludes): continue
            obsline += [FloatToString(DC.obs[b])];
    #get the groups - keep nuisances in a set so that they are never repetitions
    for groupName,nuisanceNames in DC.groups.iteritems():
        if groupName in groups:
            groups[groupName].update(set(nuisanceNames))
        else:
            groups[groupName] = set(nuisanceNames)

    # Finally report nuisance edits propagated to end of card
    for editline in DC.nuisanceEditLines:
      if len(editline)==2: nuisanceEdits.append("%s %s"%(editline[0]," ".join(editline[1])))
      elif len(editline)==4 and not editline[3]: nuisanceEdits.append(" ".join(editline[0:3]))
      else:
        tmp_chan = editline[2]
        tmp_proc = editline[1]
        if tmp_chan == "*": # all channels
          tmp_chan = "%s(%s)"%(label,"|".join(c for c in DC.bins)) if len (DC.bins)>1 else label
	  if "ifexists" not in editline[3]: editline[3].append("ifexists")
	else: tmp_chan = label+tmp_chan
        if tmp_proc == "*":
          tmp_proc = "(%s)"%("|".join(p for p in DC.processes))
	  if "ifexists" not in editline[3]: editline[3].append("ifexists")
        nuisanceEdits.append("%s %s %s %s"%(editline[0],tmp_proc,tmp_chan," ".join(editline[3])))


bins = []
check_processes = {}
process_errors = []
for (b,p,s) in keyline:
    if b not in bins: bins.append(b)
    if p not in check_processes:
        check_processes[p] = (s,b)
    else:
        if check_processes[p][0] != s:
            process_errors.append(" - process %s %s a signal in %s and %s a signal in %s" % (p,
                            "is" if s else "is not", b,
                            "is" if check_processes[p][0] else "is not", check_processes[p][1]))
    if s:
        if p not in signals: signals.append(p)
    else:
        if p not in backgrounds: backgrounds.append(p)

if process_errors:
    raise RuntimeError("ERROR: mismatch between process signal labels:\n%s" % ("\n".join(process_errors)))

print "Combination of", "  ".join(args)
print "imax %d number of bins" % len(bins)
print "jmax %d number of processes minus 1" % (len(signals) + len(backgrounds) - 1)
print "kmax %d number of nuisance parameters" % (len(systlines) + len(paramSysts))
print "-" * 130

if shapeLines:
    chmax = max([max(len(p),len(c)) for p,c,x in shapeLines]);
    cfmt = "%-"+str(chmax)+"s ";
    shapeLines.sort( lambda x,y : cmp(x[0],y[0]) if x[1] == y[1] else cmp(x[1],y[1]) )
    for (process,channel,stuff) in shapeLines:
        print "shapes", cfmt % process, cfmt % channel, ' '.join(stuff);
    print "-" * 130

if obsline:
    cmax = max([cmax]+[len(l) for l in obskeyline]+[len(x) for x in obsline])
    cfmt = "%-"+str(cmax)+"s";
    print "bin         ", "  ".join([cfmt % x for x in obskeyline])
    print "observation ", "  ".join([cfmt % x for x in obsline])

print "-" * 130

pidline = []; signals = []; backgrounds = []
tmpsignals = [];
for (b,p,s) in keyline:
    if s:
        if p not in tmpsignals: tmpsignals.append(p)
for (b,p,s) in keyline:
    if s:
        if p not in signals: signals.append(p)
        pidline.append(signals.index(p)-len(tmpsignals)+1)
    else:
        if p not in backgrounds: backgrounds.append(p)
        pidline.append(1+backgrounds.index(p))
cmax = max([cmax]+[max(len(p),len(b)) for p,b,s in keyline]+[len(e) for e in expline])
hmax = max([10] + [len("%-12s[nofloat]  %s %s" % (l,p,a)) for l,(p,a,e,nf) in systlines.items()])
cfmt  = "%-"+str(cmax)+"s"; hfmt = "%-"+str(hmax)+"s  ";
print hfmt % "bin",     "  ".join([cfmt % p for p,b,s in keyline])
print hfmt % "process", "  ".join([cfmt % b for p,b,s in keyline])
print hfmt % "process", "  ".join([cfmt % x for x in pidline])
print hfmt % "rate",    "  ".join([cfmt % x for x in expline])

print "-" * 130

sysnamesSorted = systlines.keys(); sysnamesSorted.sort()
for name in sysnamesSorted:
    (pdf,pdfargs,effect,nofloat) = systlines[name]
    if nofloat: name += "[nofloat]"
    systline = []
    for b,p,s in keyline:
        try:
            systline.append(effect[b][p])
        except KeyError:
            systline.append("-");
    print hfmt % ("%-21s   %s  %s" % (name, pdf, " ".join(pdfargs))), "  ".join([cfmt % x for x in systline])
for (pname, pargs) in paramSysts.items():
    print "%-12s  param  %s" %  (pname, " ".join(pargs))

for pname in flatParamNuisances.iterkeys():
    print "%-12s  flatParam" % pname
for pname in rateParams.iterkeys():
    for pk in range(len(rateParams[pname])):
     print "%-12s  rateParam %s"% (rateParams[pname][pk][0][0],pname.replace("AND"," ")),
     for p in rateParams[pname][pk][0][1:-1]: print p,
     print rateParams[pname][pk][1],
     print "\n",
for dname in discreteNuisances.iterkeys():
    print "%-12s  discrete" % dname
for ext in extArgs.iterkeys():
    print "%s" % ' '.join(extArgs[ext])
for groupName,nuisanceNames in groups.iteritems():
    nuisances = ' '.join(nuisanceNames)
    print '%(groupName)s group = %(nuisances)s' % locals()
for bpf in binParFlags.iterkeys():
    if len(binParFlags[bpf]) == 1:
      print "%s autoMCStats %g" % (bpf,binParFlags[bpf][0])
    if len(binParFlags[bpf]) == 2:
      print "%s autoMCStats %g %i" % (bpf,binParFlags[bpf][0], binParFlags[bpf][1])
    if len(binParFlags[bpf]) == 3:
      print "%s autoMCStats %g %i %i" % (bpf,binParFlags[bpf][0], binParFlags[bpf][1], binParFlags[bpf][2])

nuisanceEdits = set(nuisanceEdits)
nuisanceEdits_lengths = [ [len(e.split()),e] for e in nuisanceEdits ]
nuisanceEdits_lengths = sorted(nuisanceEdits_lengths,reverse=True)
nuisanceEdits = [e[1] for e in nuisanceEdits_lengths]
nuisanceEdits_lengths = 0

for edit in nuisanceEdits:
    print "nuisance edit ", edit

if options.editNuisFile:
    file = open(options.editNuisFile, "r")
    str = file.read();
    print str
