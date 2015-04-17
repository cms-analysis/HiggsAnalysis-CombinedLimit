import re
import sys
from math import log,exp,hypot

def quadratureAdd(pdf, val1, val2, context=None):
    if type(val1) == list and type(val2) == list:
        return [ quadratureAdd(pdf, val1[0], val2[0], context), quadratureAdd(pdf, val1[0], val2[0], context) ]
    elif type(val1) == list and type(val2) == float:
        return [ quadratureAdd(pdf, val1[0], 1.0/val2, context), quadratureAdd(pdf, val1[0], val2, context) ]
    elif type(val2) == list and type(val1) == float:
        return [ quadratureAdd(pdf, 1.0/val1, val2[0], context), quadratureAdd(pdf, val1, val2[0], context) ]
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
    if channel != "*": cchannel = re.compile(channel)
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
    if "/" in value:
        value = [ float(x) for x in value.split("/") ]
    else:
        value = float(value)
    foundChann, foundProc = False, False
    for b in errline.keys():
        if channel == "*" or cchannel.search(b):
            foundChann = True
            for p in datacard.exp[b]:
                if process == "*" or cprocess.search(p):
                    foundProc = True
                    if p in errline[b] and errline[b][p] not in [ 0.0, 1.0 ]:
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
    if channel != "*": cchannel = re.compile(channel)
    opts = args[3:]
    foundProc = False
    for lsyst,nofloat,pdf,args0,errline in datacard.systs:
        if re.match(name,lsyst):
            for b in errline.keys():
                if channel == "*" or cchannel.search(b):
                    #if channel != "*": foundProc = False
                    for p in datacard.exp[b]:
                        if process == "*" or cprocess.search(p):
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
    if len(args) < 4:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit rename process channel oldname newname"
    (process, channel, oldname, newname) = args[:4]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel)
    opts = args[5:]
    for lsyst,nofloat,pdf0,args0,errline0 in datacard.systs[:]:
        lsystnew = re.sub(oldname,newname,lsyst)
        if lsystnew != lsyst:
            found = False
            errline2 = dict([(b,dict([(p,0) for p in datacard.exp[b]])) for b in datacard.bins])
            for lsyst2,nofloat2,pdf2,args2,errline2b in datacard.systs:
                if lsyst2 == lsystnew:
                    found = True
                    errline2 = errline2b
                if pdf2 != pdf0: raise RuntimeError, "Can't rename nuisance %s with pdf %s to name %s which already exists as %s" % (lsyst,pdf0,lsystnew,pdf2)
            if not found:
                datacard.systs.append([lsystnew,nofloat,pdf0,args0,errline2])
            foundChann, foundProc = False, False
            for b in errline0.keys():
                if channel == "*" or cchannel.search(b):
                    foundChann = True
                    if channel != "*": foundProc = False
                    for p in datacard.exp[b].keys():
                        if process == "*" or cprocess.search(p):
                            foundProc = True
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
        if re.match(name,lsyst):
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

def doSplitNuisance(datacard, args):
    if len(args) < 7:
        raise RuntimeError, "Missing arguments: the syntax is: nuisance edit split process channel oldname newname1 newname2 value1 value2"
    (process, channel, oldname, newname1, newname2, value1, value2) = args[:7]
    if process != "*": cprocess = re.compile(process)
    if channel != "*": cchannel = re.compile(channel)
    opts = args[8:]
    foundProc = False
    for lsyst,nofloat,pdf,args0,errline in datacard.systs:
        if re.match(oldname,lsyst):
            for b in errline.keys():
                if channel == "*" or cchannel.search(b):
                    for p in datacard.exp[b]:
                        if process == "*" or cprocess.search(p):
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

def doEditNuisance(datacard, command, args):
    if command == "add":
        doAddNuisance(datacard, args)
    elif command == "drop":
        doDropNuisance(datacard, args)
    elif command == "rename":
        doRenameNuisance(datacard, args)
    elif command == "changepdf":
        doChangeNuisancePdf(datacard, args)
    elif command == "split":
        doSplitNuisance(datacard, args)
    else:
        raise RuntimeError, "Error, unknown nuisance edit command %s (args %s)" % (command, args)
