import re
import sys
from math import log,exp,hypot


def fullmatch(regex, line):
    """Performs regex match requiring the full string to match

    Args:
        regex (str or precompiled regex): The pattern to match
        line (str): The line of text to match

    Returns:
        bool: true if a complete match was found, false otherwise
    """
    if isinstance(regex, str):
        compiled = re.compile(regex)
    else:
        compiled = regex
    m = compiled.match(line)
    if m:
        start, stop = m.span()
        if stop - start == len(line):
            return True
    return False


def quadratureAdd(pdf, val1, val2, context=None):
    if type(val1) == list and len(val1) != 2: raise RuntimeError("{} is a list of length != 2".format(val1))
    if type(val2) == list and len(val2) != 2: raise RuntimeError("{} is a list of length != 2".format(val2))

    if type(val1) == list and type(val2) == list:
        return [ quadratureAdd(pdf, val1[0], val2[0], context), quadratureAdd(pdf, val1[1], val2[1], context) ]
    elif type(val1) == list and type(val2) == float:
        return [ quadratureAdd(pdf, val1[0], 1.0/val2, context), quadratureAdd(pdf, val1[1], val2, context) ]
    elif type(val2) == list and type(val1) == float:
        return [ quadratureAdd(pdf, 1.0/val1, val2[0], context), quadratureAdd(pdf, val1, val2[1], context) ]
    if pdf in [ "lnN", "lnU" ]:
        if log(val1) * log(val2) < 0:
            raise RuntimeError, "Can't add in quadrature nuisances of pdf %s with values %s, %s that go in different directions (at %s)" % (pdf, val1, val2, context)
        ret = exp(hypot(abs(log(val1)),abs(log(val2))))
        return ret if val1 > 1 else 1.0/ret
    else:
        raise RuntimeError, "Quadrature add not implemented for pdf %s (at %s)" % (pdf, context)

def doAddNuisance(datacard, args):
    if len(args) < 5:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit add process channel name pdf value [ options ]"
    (process, channel, name, pdf, value) = args[:5]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel.replace("+","\+"))
    opts = args[5:]
    found = False
    errline = dict([(b,dict([(p,0) for p in datacard.exp[b]])) for b in datacard.bins])
    for lsyst,nofloat,pdf0,args0,errline0 in datacard.systs:
        if lsyst == name:
            if pdf != pdf0: raise RuntimeError, "Can't add nuisance %s with pdf %s ad it already exists as %s" % (name,pdf,pdf0)
            found = True
            errline = errline0
    if not found:
        datacard.systs.append([name,False,pdf,[],errline])
    if isinstance(value, (int, float, list)):
        pass
    elif "/" in value:
        value = [ float(x) for x in value.split("/") ]
    else:
        value = float(value)
    foundChann, foundProc = False, False
    for b in errline.keys():
        if channel == "*" or fullmatch(cchannel, b):
            foundChann = True
            for p in datacard.exp[b]:
                if process == "*" or fullmatch(cprocess, p):
                    foundProc = True
                    if value in [ 0.0, 1.0 ]:
                        pass   #do nothing, there's nothing to add
                    elif p in errline[b] and errline[b][p] not in [ 0.0, 1.0 ]:
                        if "addq" in opts:
                            errline[b][p] = quadratureAdd(pdf, errline[b][p], value, context="nuisance edit add, args = %s" % args)
                        elif "overwrite" in opts:
                            errline[b][p] = value
                        else:
                            raise RuntimeError, "Can't add nuisance with args = %s for bin %s, process %s: found existing non-null value %s, and no option 'addq' or 'overwrite' given" % (args,b,p,errline[b][p])
                    else:
                        errline[b][p] = value
    if not foundChann:
        sys.stderr.write("Warning: no channel found matching nuisance edit add with args = %s\n" % args)
    if not foundProc:
        sys.stderr.write("Warning: no process found matching nuisance edit add with args = %s\n" % args)

def doDropNuisance(datacard, args):
    if len(args) < 3:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit drop process channel name [ options ]"
    (process, channel, name) = args[:3]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel.replace("+","\+"))
    opts = args[3:]
    foundProc = False
    for lsyst,nofloat,pdf,args0,errline in datacard.systs:
        if fullmatch(name,lsyst):
            for b in errline.keys():
                if channel == "*" or fullmatch(cchannel, b):
                    #if channel != "*": foundProc = False
                    for p in datacard.exp[b]:
                        if process == "*" or fullmatch(cprocess, p):
                            foundProc = True
                            errline[b][p] = 0
            #if channel != "*" and foundProc == False:
            #    if "ifexists" not in opts:
            #        raise RuntimeError, "Error: nuisance edit drop %s found nothing in channel %s" % (args, channel)
            #    else:
            #        sys.stderr.write("Warning1: nuisance edit drop %s found nothing in channel %s\n" % (args, channel))
    if not foundProc and channel != "*":
        if "ifexists" not in opts:
            raise RuntimeError, "Error: nuisance edit drop %s found nothing" % (args)
        else:
            sys.stderr.write("Warning2: nuisance edit drop %s found nothing\n" % (args))


def doRenameNuisance(datacard, args):
    if len(args) == 2: # newname oldname 
      nuisanceID = i = -1
      (oldname, newname) = args[:2]
      for lsyst,nofloat,pdf0,args0,errline0 in (datacard.systs[:]):
        i+=1
        if lsyst == oldname : # found the nuisance
	  nuisanceID = i
	  if pdf0 == "flatParam" :
            raise RuntimeError, "Error: Cannot use nuisance edit rename with flatParam type nuisances currently - you should rename the parameter in your input workspace."	  
	  if pdf0 != "param": 
            raise RuntimeError, "Missing arguments: the syntax is: nuisance edit rename process channel oldname newname"	  
          for lsyst2,nofloat2,pdf02,args02,errline02 in (datacard.systs[:]):
	    if lsyst2 == newname:
	     if pdf02 != "param":
	      if (args0[0]) not in ["0.0","0.","0"] or (args0[1]) not in ["1.0","1.","1"] : raise RuntimeError, "Can't rename nuisance %s with Gaussian pdf G(%s,%s) to name %s which already exists with G(0,1) constraint!" % (lsyst,args0[0],args0[1],lsyst2)
             else: 
	      if (args0[0])!=(args02[0])  or float(args0[1])!=float(args02[1]) : raise RuntimeError, "Can't rename nuisance %s with Gaussian pdf G(%s,%s) to name %s which already exists with Gaussian pdf G(%s,%s) constraint!" % (lsyst,args0[0],args0[1],lsyst2,args02[0],args02[1])
	  break
      if nuisanceID >-1 : 
      	datacard.systs[nuisanceID][0]=newname 
	datacard.systematicsParamMap[oldname]=newname
      else: raise RuntimeError, "No nuisance parameter found with name %s in the datacard"%oldname 
      return

    if len(args) < 4:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit rename process channel oldname newname"
    (process, channel, oldname, newname) = args[:4]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel.replace("+","\+"))
    opts = args[4:]
    foundChann, foundProc = False, False
    for lsyst,nofloat,pdf0,args0,errline0 in datacard.systs[:]:
        lsystnew = lsyst
        if fullmatch(oldname, lsyst):
            lsystnew = newname
        if lsystnew != lsyst:
            if pdf0=="param" : raise RuntimeError, "Incorrect syntax. Cannot specify process and channel for %s with pdf %s. Use 'nuisance edit rename oldname newname'"% (lsyst,pdf0)
            found = False
            errline2 = dict([(b,dict([(p,0) for p in datacard.exp[b]])) for b in datacard.bins])
            for lsyst2,nofloat2,pdf2,args2,errline2b in datacard.systs:
                if lsyst2 == lsystnew:
                    found = True
                    errline2 = errline2b
                    if pdf2 != pdf0 and pdf2 not in ['lnN']: raise RuntimeError, "Can't rename nuisance %s with pdf %s to name %s which already exists as %s" % (lsyst,pdf0,lsystnew,pdf2)
            if not found:
                datacard.systs.append([lsystnew,nofloat,pdf0,args0,errline2])
            for b in errline0.keys():
                if channel == "*" or fullmatch(cchannel, b):
                    foundChann = True
                    if channel != "*": foundProc = False
                    for p in datacard.exp[b].keys():
                        if process == "*" or fullmatch(cprocess, p):
                            foundProc = True
                            if errline0[b][p] in [0.0]:
                                continue
    			    if "shape" in pdf0 : datacard.systematicsShapeMap[newname,b,p]=oldname
                            if p in errline0[b] and errline2[b][p] not in [ 0.0, 1.0 ]:
                                if "addq" in opts:
                                    errline2[b][p] = quadratureAdd(pdf0, errline0[b][p], errline2[b][p], context="nuisance edit rename, args = %s" % args)
                                elif "overwrite" in opts:
                                    errline2[b][p] = errline0[b][p]
                                else:
                                    raise RuntimeError, "Can't rename nuisance with args = %s for bin %s, process %s: found existing non-null value %s, and no option 'addq' or 'overwrite' given" % (args,b,p,errline2[b][p])
                            else:
                                errline2[b][p] = errline0[b][p]
                            errline0[b][p] = 0
                    if channel != "*" and not foundProc:
                        if "ifexists" not in opts:
                            raise RuntimeError, "Error: nuisance edit rename %s found nothing in channel %s" % (args, channel)
                        else:
                            sys.stderr.write("Warning: nuisance edit rename %s found nothing in channel %s\n" % (args, channel))
    if not foundProc and channel != "*":
        if "ifexists" not in opts:
            raise RuntimeError, "Error: no process found matching nuisance edit rename with args = %s (and option 'ifexist' not specified)\n" % args
        else:
            sys.stderr.write("Warning: no process found matching nuisance edit rename with args = %s\n" % args)
    if not foundChann:
        sys.stderr.write("Warning: no channel found matching nuisance edit rename with args = %s\n" % args)

def doChangeNuisancePdf(datacard, args):
    if len(args) < 2:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit changepdf name newpdf [ options ]"
    (name, newpdf) = args[:2]
    found = False
    for i,(lsyst,nofloat,pdf,args0,errline) in enumerate(datacard.systs):
        if fullmatch(name,lsyst):
            found = True; ok = False
            if newpdf == pdf: continue
            if pdf in [ "lnN", "lnU"]:
                if newpdf in ["lnN","lnU"]:
                    ok = True
                elif newpdf in [ "trG", "unif", "dFD", "dFD2"]:
                    ok = True
                    for b in errline.keys():
                        for p in errline[b].keys():
                            errline[b][p] = log(errline[b][p]) if errline[b][p] not in [0., 1.] else 0.
            elif pdf in ["trG", "unif", "dFD", "dFD2"]:
                if newpdf in [ "trG", "unif", "dFD", "dFD2"]:
                    ok = True
                elif newpdf in ["lnN","lnU"]:
                    ok = True
                    for b in errline.keys():
                        for p in errline[b].keys():
                            errline[b][p] = exp(errline[b][p]) if errline[b][p] not in [0., 1.] else 0.
            elif "shape" in pdf:
                if "shape" in newpdf:
                    ok = True
            if ok:
                datacard.systs[i][2] = newpdf
            else:
                raise RuntimeError, "I can't convert pdf %s from pdf %s to %s" % (lsyst,pdf,newpdf)
    if not found:
        sys.stderr.write("Warning: no pdf found for changepdf with args %s\n" % args)

def doMergeNuisance(datacard, args):
    if len(args) < 4:
        raise RuntimeError("Missing arguments: the syntax is: nuisance edit merge process channel name1 name2 [ options ]")
    (process, channel, name1, name2) = args[:4]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel.replace("+","\+"))
    opts = args[4:]
    foundProc = False

    for lsyst2,nofloat2,pdf2,args02,errline2 in datacard.systs:
        if fullmatch(name2, lsyst2):
            for b in errline2.keys():
                if channel == "*" or fullmatch(cchannel, b):
                    for p in datacard.exp[b]:
                        if process == "*" or fullmatch(cprocess, p):
                            foundProc = True
                            doAddNuisance(datacard, [p+"$", b+"$", name1, pdf2, errline2[b][p], "addq"])
                            errline2[b][p] = 0

    if not foundProc and channel != "*":
        if "ifexists" not in opts:
            raise RuntimeError("Error: nuisance edit merge %s found nothing" % (args))
        else:
            sys.stderr.write("Warning2: nuisance edit merge %s found nothing\n" % (args))


def doSplitNuisance(datacard, args):
    if len(args) < 7:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit split process channel oldname newname1 newname2 value1 value2"
    (process, channel, oldname, newname1, newname2, value1, value2) = args[:7]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel.replace("+","\+"))
    opts = args[7:]
    foundProc = False
    for lsyst,nofloat,pdf,args0,errline in datacard.systs:
        if fullmatch(oldname,lsyst):
            for b in errline.keys():
                if channel == "*" or fullmatch(cchannel, b):
                    for p in datacard.exp[b]:
                        if process == "*" or fullmatch(cprocess, p):
                            foundProc = True
                            if errline[b][p] not in [0., 1.]:
                                doAddNuisance(datacard, [p, b, newname1, pdf, value1, "overwrite"])
                                doAddNuisance(datacard, [p, b, newname2, pdf, value2, "overwrite"])
                            errline[b][p] = 0

    if not foundProc and channel != "*":
        if "ifexists" not in opts:
            raise RuntimeError, "Error: nuisance edit split %s found nothing" % (args)
        else:
            sys.stderr.write("Warning2: nuisance edit split %s found nothing\n" % (args))

def doFreezeNuisance(datacard, args):
    if len(args) < 1:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit freeze name [ifexists] (name can be a pattern)"
    pat = re.compile("^"+args[0]+"$")
    opts = args[1:]
    found = []

    # first check in the list of paramters as flatParam, rateParam or discretes not included in datacard.systs (smaller usually)
    for lsyst in datacard.flatParamNuisances.keys()+list(datacard.rateParamsOrder)+datacard.discretes +datacard.extArgs.keys():
         if fullmatch(pat,lsyst):
            datacard.frozenNuisances.add(lsyst)
            found.append(lsyst)

    if not found: 
      for lsyst,nofloat,pdf,args0,errline in datacard.systs:
        if fullmatch(pat,lsyst):
            datacard.frozenNuisances.add(lsyst)
            found.append(lsyst)

        
    # Warn user/exit  
    if not found:
        if "ifexists" not in opts:
            raise RuntimeError, "Error: nuisance edit freeze %s found nothing" % args[0]
        else:
            sys.stderr.write("Warning2: nuisance edit freeze %s found nothing\n" % args[0])

def doFlipNuisance(datacard,args):
    if len(args) < 3:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit flip process channel name [options: ifexists, p2n (only positive to negative), n2p (only negative to positive)]"
    (process, channel, name) = args[:3]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel.replace("+","\+"))
    opts = args[3:]
    if "n2p" not in opts and "p2n" not in opts:
        raise RuntimeError, "Error: nuisance edit flip %s missed option n2p and/or p2n" % (args)
    foundProc = False
    for lsyst,nofloat,pdf,args0,errline in datacard.systs:
        if fullmatch(name,lsyst):
            if pdf not in ["lnN"]:
                raise RuntimeError, "Error: nuisance edit flip %s currently not support pdftype %s" % (args,pdf)
            for b in errline.keys():
                if channel == "*" or fullmatch(cchannel, b):
                    for p in datacard.exp[b]:
                        if process == "*" or fullmatch(cprocess, p):
                            foundProc = True
                            if errline[b][p] not in [0., 1.]:
                                if type(errline[b][p]) is list:
                                    if errline[b][p][0] < 1 :
                                        if "p2n" in opts:
                                            errline[b][p][0] = 1./errline[b][p][0]
                                            errline[b][p][1] = 1./errline[b][p][1]
                                    elif errline[b][p][0] > 1:
                                        if "n2p" in opts:
                                            errline[b][p][0] = 1./errline[b][p][0]
                                            errline[b][p][1] = 1./errline[b][p][1]
                                else:
                                    if errline[b][p] > 1:
                                        if "p2n" in opts:  errline[b][p] = 1./errline[b][p]
                                    elif errline[b][p] < 1:
                                        if "n2p" in opts:  errline[b][p] = 1./errline[b][p]

    if not foundProc and channel != "*":
        if "ifexists" not in opts:
            raise RuntimeError, "Error: nuisance edit flip %s found nothing" % (args)
        else:
            sys.stderr.write("Warning2: nuisance edit flip %s found nothing\n" % (args))

def doEditNuisance(datacard, command, args):
    if command == "add":
        doAddNuisance(datacard, args)
    elif command == "drop":
        doDropNuisance(datacard, args)
    elif command == "rename":
        doRenameNuisance(datacard, args)
    elif command == "changepdf":
        doChangeNuisancePdf(datacard, args)
    elif command == "merge":
        doMergeNuisance(datacard, args)
    elif command == "split":
        doSplitNuisance(datacard, args)
    elif command == "freeze":
        doFreezeNuisance(datacard, args)
    elif command == "flip":
        doFlipNuisance(datacard, args)
    else:
        raise RuntimeError, "Error, unknown nuisance edit command %s (args %s)" % (command, args)
        
