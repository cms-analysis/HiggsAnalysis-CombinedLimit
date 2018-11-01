import ROOT 
from sys import argv, stdout, stderr, exit
import datetime
from optparse import OptionParser

hasHelp = False
for X in ("-h", "-?", "--help"):
    if X in argv:
        hasHelp = True
        argv.remove(X)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
#ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
argv.remove( '-b-' )
if hasHelp: argv.append("-h")

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("", "--printValueOnly",  dest="printValueOnly", default=False, action='store_true', help="Just print the default value of the normalisation.")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

file_in = ROOT.TFile(args[0])
ws = file_in.Get("w")

def find_chan_proc(name): 
 chan = norm_name[norm_name.find("_bin")+len("_bin"):norm_name.find("_proc")]
 if not "proc" in name: return chan,""
 proc = norm_name[norm_name.find("_proc_")+len("_proc_"):]
 return chan,proc

chan_procs = {}

all_norms   = ws.allFunctions().selectByName("n_exp_final*")
norm_it = all_norms.createIterator()
for i in range(all_norms.getSize()): 
 norm = norm_it.Next()
 norm_name = norm.GetName()
 chan,proc = find_chan_proc(norm_name)
 if len(proc.strip())==0: continue #ignore summations

 if chan in chan_procs: chan_procs[chan].append([proc,norm,1])
 else: chan_procs[chan]= [[proc,norm,1]]

# now look for cases where there wasn't a "final"
all_norms   = ws.allFunctions().selectByName("n_exp*")
norm_it = all_norms.createIterator()
for i in range(all_norms.getSize()): 
 norm = norm_it.Next()
 norm_name = norm.GetName()
 chan,proc = find_chan_proc(norm_name)
 if len(proc.strip())==0: continue #ignore summations
 
 if chan in chan_procs: 
  if proc in [chan_procs[chan][i][0] for i in range(len(chan_procs[chan]))] : continue 
  chan_procs[chan].append([proc,norm,0,1])
 else: chan_procs[chan]= [[proc,norm,0,1]]

# Finally look for the simplest case where the normalisation is a constant fixed number 
all_norms   = ws.allVars().selectByName("n_exp*")
norm_it = all_norms.createIterator()
for i in range(all_norms.getSize()): 
 norm = norm_it.Next()
 norm_name = norm.GetName()
 chan,proc = find_chan_proc(norm_name)
 if len(proc.strip())==0: continue #ignore summations
 
 if chan in chan_procs: 
  if proc in [chan_procs[chan][i][0] for i in range(len(chan_procs[chan]))] : continue 
  chan_procs[chan].append([proc,norm,0,0])
 else: chan_procs[chan]= [[proc,norm,0,0]]

# Now print out information 
for chan in chan_procs.keys(): 
 print "---------------------------------------------------------------------------"
 print "---------------------------------------------------------------------------"
 print "Channel - %s "%chan 
 chanInfo  = chan_procs[chan]
 for proc in chanInfo :
  print "---------------------------------------------------------------------------"
  print "  Top-level normalisation for process %s -> %s"%(proc[0],proc[1].GetName())
  print "  -------------------------------------------------------------------------"
  if options.printValueOnly: print "  default value = ",proc[1].getVal()
  else: 
    if proc[2]: 
      proc_norm_var = ws.function("n_exp_bin%s_proc_%s"%(chan,proc[0]))
      proc[1].Print()
      if (proc_norm_var.Class().GetName()==ROOT.ProcessNormalization().Class().GetName()):
        print " ... is a product, which contains ",proc_norm_var.GetName()
        proc_norm_var.dump()
      else: 
        print " ... is a product, which contains ",proc_norm_var.GetName()
        #proc_norm_var = ws.var("n_exp_bin%s_proc_%s"%(chan,proc[0]))
	proc_norm_var.Print()
    else:
      if proc[3]: 
        if (proc[1].Class().GetName()==ROOT.ProcessNormalization().Class().GetName()): proc[1].dump()
	else: 
	  proc[1].Print()
      	  print " ... is a constant (formula)"
      else: 
        proc[1].Print()
      	print " ... is a constant (RooRealVar)"
    print "  -------------------------------------------------------------------------"
    print "  default value = ",proc[1].getVal()
