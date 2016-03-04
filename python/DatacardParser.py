import re
from sys import stderr

globalNuisances = re.compile('(lumi|pdf_(qqbar|gg|qg)|QCDscale_(ggH|qqH|VH|ggH1in|ggH2in|VV)|UEPS|FakeRate|CMS_(eff|fake|trigger|scale|res)_([gemtjb]|met))')

def addDatacardParserOptions(parser):
    parser.add_option("-s", "--stat",   dest="stat",    default=False, action="store_true", help="keep only statistical uncertainties, no systematics") 
    parser.add_option("-f", "--fix-pars", dest="fixpars",default=False, action="store_true", help="fix all floating parameters of the pdfs except for the POI") 
    parser.add_option("-c", "--compiled", dest="cexpr", default=False, action="store_true", help="use compiled expressions (not suggested)")
    parser.add_option("-a", "--ascii",    dest="bin",   default=True, action="store_false", help="produce a Workspace in a rootfile in an HLF file (legacy, unsupported)")
    parser.add_option("-b", "--binary",   dest="bin",   default=True, action="store_true",  help="produce a Workspace in a rootfile (default)")
    parser.add_option("-o", "--out",      dest="out",   default=None,  type="string", help="output file (if none, it will print to stdout). Required for binary mode.")
    parser.add_option("-v", "--verbose",  dest="verbose",  default=0,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")
    parser.add_option("-m", "--mass",     dest="mass",     default=0,  type="float",  help="Higgs mass to use. Will also be written in the Workspace as RooRealVar 'MH'.")
    parser.add_option("-D", "--dataset",  dest="dataname", default="data_obs",  type="string",  help="Name of the observed dataset")
    parser.add_option("-L", "--LoadLibrary", dest="libs",  type="string" , action="append", help="Load these libraries")
    parser.add_option("--poisson",  dest="poisson",  default=0,  type="int",    help="If set to a positive number, binned datasets wih more than this number of entries will be generated using poissonians")
    parser.add_option("--default-morphing",  dest="defMorph", type="string", default="shape2N", help="Default template morphing algorithm (to be used when the datacard has just 'shape')")
    parser.add_option("--no-b-only","--for-fits",    dest="noBOnly", default=False, action="store_true", help="Do not save the background-only pdf (saves time)")
    parser.add_option("--no-optimize-pdfs",    dest="noOptimizePdf", default=False, action="store_true", help="Do not save the RooSimultaneous as RooSimultaneousOpt and Gaussian constraints as SimpleGaussianConstraint")
    parser.add_option("--optimize-simpdf-constraints",    dest="moreOptimizeSimPdf", default=False, action="store_true", help="Deeper optimization of RooSimultaneous: add the constraints only at the end (RooFit-incompatible!)")
    #parser.add_option("--use-HistPdf",  dest="useHistPdf", type="string", default="always", help="Use RooHistPdf for TH1s: 'always' (default), 'never', 'when-constant' (i.e. not when doing template morphing)")
    parser.add_option("--use-HistPdf",  dest="useHistPdf", type="string", default="never", help="Use RooHistPdf for TH1s: 'always', 'never' (default), 'when-constant' (i.e. not when doing template morphing)")
    parser.add_option("--X-exclude-nuisance", dest="nuisancesToExclude", type="string", action="append", default=[], help="Exclude nuisances that match these regular expressions.")
    parser.add_option("--X-rescale-nuisance", dest="nuisancesToRescale", type="string", action="append", nargs=2, default=[], help="Rescale by this factor the nuisances that match these regular expressions (the rescaling is applied to the sigma of the gaussian constraint term).")
    parser.add_option("--X-force-no-simpdf",  dest="forceNonSimPdf", default=False, action="store_true", help="FOR DEBUG ONLY: Do not produce a RooSimultaneous if there is just one channel (note: can affect performance)")
    parser.add_option("--X-no-check-norm",  dest="noCheckNorm", default=False, action="store_true", help="FOR DEBUG ONLY: Turn off the consistency check between datacard norms and shape norms. Will give you nonsensical results if you have shape uncertainties.")
    parser.add_option("--X-no-jmax",  dest="noJMax", default=False, action="store_true", help="FOR DEBUG ONLY: Turn off the consistency check between jmax and number of processes.")
    parser.add_option("--X-allow-no-signal",  dest="allowNoSignal", default=False, action="store_true", help="Allow parsing datacards that contain only backgrounds")
    parser.add_option("--X-allow-no-background",  dest="allowNoBackground", default=False, action="store_true", help="Allow parsing datacards that contain only signal")
    #parser.add_option("--X-optimize-templates",  dest="optimizeExistingTemplates", default=False, action="store_true", help="Optimize templates on the fly (relevant for HZZ)")
    #parser.add_option("--X-optimize-bound-nusances",  dest="optimizeBoundNuisances", default=False, action="store_true", help="Flag nuisances to have a different implementation of bounds")
    parser.add_option("--X-no-optimize-templates",  dest="optimizeExistingTemplates", default=True, action="store_false", help="Don't optimize templates on the fly (relevant for HZZ)")
    parser.add_option("--X-no-optimize-bound-nusances",  dest="optimizeBoundNuisances", default=True, action="store_false", help="Don't flag nuisances to have a different implementation of bounds")
    parser.add_option("--X-no-optimize-bins",  dest="optimizeTemplateBins", default=True, action="store_false", help="Don't optimize template bins")


from HiggsAnalysis.CombinedLimit.Datacard import Datacard
from HiggsAnalysis.CombinedLimit.NuisanceModifier import doEditNuisance

def isVetoed(name,vetoList):
    for pattern in vetoList:
        if not pattern: continue 
        if re.match(pattern,name): return True
    return False

def isIncluded(name,includeList):
    if not len(includeList): return True
    for pattern in includeList:
        if not pattern: continue 
        if re.match(pattern,name): return True
    return False

def addRateParam(lsyst,f,ret):

    if len(f) > 6: raise RuntimeError, "Directives of type 'rateParam' can only have a single channel. name rateParam channel process [init / expression vars]. repeat if rate is to affect multiple channels"
    if len(f)==5  : tmp_exp = [lsyst,f[4],0]
    elif len(f)==6: tmp_exp = [lsyst,f[4],f[5],1]
    if ("%sAND%s"%(f[2],f[3])) in ret.rateParams.keys(): ret.rateParams["%sAND%s"%(f[2],f[3])].append(tmp_exp)
    else: ret.rateParams["%sAND%s"%(f[2],f[3])] = [tmp_exp]

def parseCard(file, options):
    if type(file) == type("str"):
        raise RuntimeError, "You should pass as argument to parseCards a file object, stream or a list of lines, not a string"
    ret = Datacard()
    ret.discretes=[]
    ret.groups={}
    #
    nbins      = -1; 
    nprocesses = -1; 
    nuisances  = -1;
    binline = []; processline = []; sigline = []
    shapesUseBin = False
    lineNumber = None
    try:
        for lineNumber,l in enumerate(file):
            f = l.split();
            if len(f) < 1: continue
            if f[0] == "imax": 
                nbins = int(f[1]) if f[1] != "*" else -1
            if f[0] == "jmax": 
                nprocesses = int(f[1])+1 if f[1] != "*" else -1
            if f[0] == "kmax": 
                nuisances = int(f[1]) if f[1] != "*" else -1
            if f[0] == "shapes":
                if not options.bin: raise RuntimeError, "Can use shapes only with binary output mode"
                if len(f) < 4: raise RuntimeError, "Malformed shapes line"
                if not ret.shapeMap.has_key(f[2]): ret.shapeMap[f[2]] = {}
                if ret.shapeMap[f[2]].has_key(f[1]): raise RuntimeError, "Duplicate definition for process '%s', channel '%s'" % (f[1], f[2])
                ret.shapeMap[f[2]][f[1]] = f[3:]
                if "$CHANNEL" in l: shapesUseBin = True
                if f[2] != "*":     shapesUseBin = True
            if f[0] == "Observation" or f[0] == "observation": 
                ret.obs = [ float(x) for x in f[1:] ]
                if nbins == -1: nbins = len(ret.obs)
                if len(ret.obs) != nbins: raise RuntimeError, "Found %d observations but %d bins have been declared" % (len(ret.obs), nbins)
                if binline != []:
                    if len(binline) != len(ret.obs): raise RuntimeError, "Found %d bins (%s) but %d bins have been declared" % (len(ret.bins), ret.bins, nbins)
                    ret.bins = binline
                    ret.obs = dict([(b,ret.obs[i]) for i,b in enumerate(ret.bins)])
                    binline = []
            if f[0] == "bin": 
                binline = []
                for b in f[1:]:
                    if re.match("[0-9]+", b):
                        if shapesUseBin: stderr.write("Warning: Bin %(b)s starts with a digit. Will call it 'bin%(b)s' but this may break shapes.\n" % locals())
                        b = "bin"+b
                        # TODO Here should be some patching of the shapes names in order to not get errors later.
                    binline.append(b)
            if f[0] == "process": 
                if processline == []: # first line contains names
                    processline = f[1:]
                    if len(binline) != len(processline): raise RuntimeError, "'bin' line has a different length than 'process' line."
                    continue
                sigline = f[1:] # second line contains ids
                if re.match("-?[0-9]+", processline[0]) and not re.match("-?[0-9]+", sigline[0]):
                    (processline,sigline) = (sigline,processline)
                if len(sigline) != len(processline): raise RuntimeError, "'bin' line has a different length than 'process' line."
                hadBins = (len(ret.bins) > 0)
                for i,b in enumerate(binline):
                    p = processline[i];
                    s = (int(sigline[i]) <= 0) # <=0 for signals, >0 for backgrounds
                    ret.keyline.append((b, processline[i], s))
                    if hadBins:
                        if b not in ret.bins: raise RuntimeError, "Bin %s not among the declared bins %s" % (b, ret.bins)
                    else:
                        if b not in ret.bins: ret.bins.append(b)
                    if p not in ret.processes: ret.processes.append(p)
                if nprocesses == -1: nprocesses = len(ret.processes)
                if nbins      == -1: nbins      = len(ret.bins)
                if not options.noJMax:
                    if nprocesses != len(ret.processes): raise RuntimeError, "Found %d processes (%s), declared jmax = %d" % (len(ret.processes),ret.processes,nprocesses)
                if nbins      != len(ret.bins):      raise RuntimeError, "Found %d bins (%s), declared imax = %d" % (len(ret.bins),ret.bins,nbins)
                ret.exp = dict([(b,{}) for b in ret.bins])
                ret.isSignal = dict([(p,None) for p in ret.processes])
                if ret.obs != [] and type(ret.obs) == list: # still as list, must change into map with bin names
                    ret.obs = dict([(b,ret.obs[i]) for i,b in enumerate(ret.bins)])
                for (b,p,s) in ret.keyline:
                    if ret.isSignal[p] == None: 
                        ret.isSignal[p] = s
                    elif ret.isSignal[p] != s:
                        raise RuntimeError, "Process %s is declared as signal in some bin and as background in some other bin" % p
                ret.signals = [p for p,s in ret.isSignal.items() if s == True]
                if len(ret.signals) == 0 and not options.allowNoSignal: raise RuntimeError, "You must have at least one signal process (id <= 0)"
            if f[0] == "rate":
                if processline == []: raise RuntimeError, "Missing line with process names before rate line" 
                if sigline == []:     raise RuntimeError, "Missing line with process id before rate line" 
                if len(f[1:]) != len(ret.keyline): raise RuntimeError, "Malformed rate line: length %d, while bins and process lines have length %d" % (len(f[1:]), len(ret.keyline))
                for (b,p,s),r in zip(ret.keyline,f[1:]):
                    ret.exp[b][p] = float(r)
                break # rate is the last line before nuisances
        # parse nuisances   
        for lineNumber,l in enumerate(file):
            if l.startswith("--"): continue
            l  = re.sub("\\s*#.*","",l)
            l = re.sub("(?<=\\s)-+(\\s|$)"," 0\\1",l);
            f = l.split();
            if len(f) <= 1: continue
            nofloat = False
            lsyst = f[0]; pdf = f[1]; args = []; numbers = f[2:];
            if lsyst.endswith("[nofloat]"):
              lsyst = lsyst.replace("[nofloat]","")
              nofloat = True
            if options.nuisancesToExclude and isVetoed(lsyst, options.nuisancesToExclude):
                if options.verbose > 0: stderr.write("Excluding nuisance %s selected by a veto pattern among %s\n" % (lsyst, options.nuisancesToExclude))
                if nuisances != -1: nuisances -= 1
                continue
            if re.match("[0-9]+",lsyst): lsyst = "theta"+lsyst
            if pdf == "lnN" or pdf == "lnU" or pdf == "gmM" or pdf == "trG" or pdf.startswith("shape"):
                pass # nothing special to do
            elif pdf == "gmN":
                args = [int(f[2])]; numbers = f[3:];
            elif pdf == "unif":
                args = [float(f[2]), float(f[3])]; numbers = f[4:];
            elif pdf == "dFD" or pdf == "dFD2":
                args = [float(f[2])];  numbers = f[3:];
            elif pdf == "param":
                # for parametric uncertainties, there's no line to account per bin/process effects
                # just assume everything else is an argument and move on
                args = f[2:]
                if len(args) <= 1: raise RuntimeError, "Uncertainties of type 'param' must have at least two arguments (mean and sigma)"
                ret.systs.append([lsyst,nofloat,pdf,args,[]])
                continue
            elif pdf == "flatParam":
                ret.flatParamNuisances[lsyst] = True
                #for flat parametric uncertainties, code already does the right thing as long as they are non-constant RooRealVars linked to the model
                continue
            elif pdf == "rateParam":
	        if f[3]=="*" and f[2]=="*": # all channels 
		  for c in ret.processes: 
		   for b in ret.bins:
		    f_tmp = f[:]
		    f_tmp[2]=b
		    f_tmp[3]=c
	            addRateParam(lsyst,f_tmp,ret)
	        elif f[3]=="*": # all channels 
		  for c in ret.processes: 
		    f_tmp = f[:]; f_tmp[3]=c
	            addRateParam(lsyst,f_tmp,ret)
	        elif f[2]=="*": # all channels 
		  for b in ret.bins:
		    f_tmp = f[:]; f_tmp[2]=b
	            addRateParam(lsyst,f_tmp,ret)
		else : addRateParam(lsyst,f,ret)
                continue
            elif pdf=="discrete":
                args = f[2:]
                ret.discretes.append(lsyst)
                continue
            elif pdf=="edit":
                if nuisances != -1: nuisances = -1
                if options.verbose > 1: print "Before edit: \n\t%s\n" % ("\n\t".join( [str(x) for x in ret.systs] ))
                if options.verbose > 1: print "Edit command: %s\n" % numbers
                doEditNuisance(ret, numbers[0], numbers[1:])
                if options.verbose > 1: print "After edit: \n\t%s\n" % ("\n\t".join( [str(x) for x in ret.systs] ))
                continue
            elif pdf=="group":
                # This is not really a pdf type, but a way to be able to name groups of nuisances together
                groupName = lsyst
                groupNuisances = numbers

                if not groupNuisances:
                    raise RuntimeError, "Syntax error for group '%s': empty line after 'group'." % groupName

                defToks = ('=','+=')
                defTok = groupNuisances.pop(0)
                if defTok not in defToks:
                    raise RuntimeError, "Syntax error for group '%s': first thing after 'group' is not '[+]=' but '%s'." % (groupName,defTok)
                
                if groupName not in ret.groups:
                    if defTok=='=':
                        ret.groups[groupName] = set(groupNuisances)
                    else:
                        raise RuntimeError, "Cannot append to group '%s' as it was not yet defined." % groupName                                                                                                    
                else:
                    if defTok=='+=' :
                        ret.groups[groupName].update( set(groupNuisances) )
                    else:
                        raise RuntimeError, "Will not redefine group '%s'. It previously contained '%s' and you now wanted it to contain '%s'." % (groupName,ret.groups[groupName],groupNuisances)                        

                continue
            else:
                raise RuntimeError, "Unsupported pdf %s" % pdf
            if len(numbers) < len(ret.keyline): raise RuntimeError, "Malformed systematics line %s of length %d: while bins and process lines have length %d" % (lsyst, len(numbers), len(ret.keyline))
            errline = dict([(b,{}) for b in ret.bins])
            nonNullEntries = 0 
            for (b,p,s),r in zip(ret.keyline,numbers):
                if "/" in r: # "number/number"
                    if (pdf not in ["lnN","lnU"]) and ("?" not in pdf): raise RuntimeError, "Asymmetric errors are allowed only for Log-normals"
                    errline[b][p] = [ float(x) for x in r.split("/") ]
                    for v in errline[b][p]:
                        if v <= 0.00: raise ValueError('Found "%s" in the nuisances affecting %s for %s. This would lead to NANs later on, so please fix it.'%(r,p,b))
                else:
                    errline[b][p] = float(r)
                    #values of 0.0 are treated as 1.0; scrap negative values.
                    if pdf not in ["trG", "dFD", "dFD2"] and errline[b][p] < 0: raise ValueError('Found "%s" in the nuisances affecting %s in %s. This would lead to NANs later on, so please fix it.'%(r,p,b))
                # set the rate to epsilon for backgrounds with zero observed sideband events.
                if pdf == "gmN" and ret.exp[b][p] == 0 and float(r) != 0: ret.exp[b][p] = 1e-6
            ret.systs.append([lsyst,nofloat,pdf,args,errline])
    except Exception, ex:
        if lineNumber != None:
            msg = "Error reading line %d" % (lineNumber + 1)
            if hasattr(file,'name'):
                msg += " of file " + file.name

            msg += ": " + ex.args[0]
            ex.args = (msg, ) + ex.args[1:]

        raise


    # check if there are bins with no rate
    for b in ret.bins:
        np_bin = sum([(ret.exp[b][p] != 0) for (b1,p,s) in ret.keyline if b1 == b])
        ns_bin = sum([(ret.exp[b][p] != 0) for (b1,p,s) in ret.keyline if b1 == b and s == True])
        nb_bin = sum([(ret.exp[b][p] != 0) for (b1,p,s) in ret.keyline if b1 == b and s != True])
        if np_bin == 0: raise RuntimeError, "Bin %s has no processes contributing to it" % b
        if ns_bin == 0 and not options.allowNoSignal: stderr.write("Warning: Bin %s has no signal processes contributing to it\n" % b)
        if nb_bin == 0 and not options.allowNoBackground: raise RuntimeError, "Bin %s has no background processes contributing to it" % b
    # cleanup systematics that have no effect to avoid zero derivatives
    syst2 = []
    for lsyst,nofloat,pdf,args,errline in ret.systs:
        nonNullEntries = 0 
        if pdf == "param" or pdf=="discrete" or pdf=="rateParam": # this doesn't have an errline
            syst2.append((lsyst,nofloat,pdf,args,errline))
            continue
        for (b,p,s) in ret.keyline:
            r = errline[b][p]
            nullEffect = (r == 0.0 or (pdf == "lnN" and r == 1.0))
            if not nullEffect and ret.exp[b][p] != 0: nonNullEntries += 1 # is this a zero background?
        if nonNullEntries != 0: syst2.append((lsyst,nofloat,pdf,args,errline))
        elif nuisances != -1: nuisances -= 1 # remove from count of nuisances, since qe skipped it
    ret.systs = syst2
    # remove them if options.stat asks so
    if options.stat: 
        nuisances = 0
        ret.systs = []
    # check number of nuisances
    if nuisances == -1: 
        nuisances = len(ret.systs)
    elif len(ret.systs) != nuisances: 
        raise RuntimeError, "Found %d systematics, expected %d" % (len(ret.systs), nuisances)
    # set boolean to know about shape
    ret.hasShapes = (len(ret.shapeMap) > 0)
    # return result
    return ret
