import re
from sys import argv, stdout, stderr, exit

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

if len(argv) == 0: raise RuntimeError, "Usage: mlfitNormsToText.py [ -u ] mlfit.root";

errors = False
if len(argv) > 2 and argv[1] == "-u": 
    errors = True
    argv[1] = argv[2];
file = ROOT.TFile.Open(argv[1]);
prefit = file.Get("norm_prefit")
fit_s = file.Get("norm_fit_s")
fit_b = file.Get("norm_fit_b")
if prefit == None: stderr.write("Missing fit_s in %s. Did you run MaxLikelihoodFit in a recent-enough version of combine and with --saveNorm?\n" % file);
if fit_s  == None: raise RuntimeError, "Missing fit_s in %s. Did you run MaxLikelihoodFit with --saveNorm?" % file;
if fit_b  == None: raise RuntimeError, "Missing fit_b in %s. Did you run MaxLikelihoodFit with --saveNorm?" % file;

iter = fit_s.createIterator()
while True:
    norm_s = iter.Next()
    if norm_s == None: break;
    norm_b = fit_b.find(norm_s.GetName())
    norm_p = prefit.find(norm_s.GetName()) if prefit else None
    m = re.match(r"(\w+)/(\w+)", norm_s.GetName());
    if m == None: m = re.match(r"n_exp_(?:final_)?(?:bin)+(\w+)_proc_(\w+)", norm_s.GetName());
    if m == None: raise RuntimeError, "Non-conforming object name %s" % norm_s.GetName()
    if norm_b == None: raise RuntimeError, "Missing normalization %s for background fit" % norm_s.GetName()
    if prefit and norm_p and errors:
        print "%-30s %-30s %7.3f +/- %7.3f   %7.3f +/- %7.3f  %7.3f +/- %7.3f" % (m.group(1), m.group(2), norm_p.getVal(), norm_p.getError(), norm_s.getVal(), norm_s.getError(), norm_b.getVal(), norm_b.getError())
    else:
        if errors:
            print "%-30s %-30s %7.3f +/- %7.3f  %7.3f +/- %7.3f" % (m.group(1), m.group(2), norm_s.getVal(), norm_s.getError(), norm_b.getVal(), norm_b.getError())
        else:
            print "%-30s %-30s %7.3f %7.3f" % (m.group(1), m.group(2), norm_s.getVal(), norm_b.getVal())
