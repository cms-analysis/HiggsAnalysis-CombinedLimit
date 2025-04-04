import numpy as np
import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="input ws file")
parser.add_argument("--backend", help='set evaluation backend, default is "combine"')
parser.add_argument("--write-debug-macro", action="store_true", help='write debug macro in case of the "codegen" backend')


args = parser.parse_args()

with ROOT.TFile.Open(args.input) as f:
    ws = f.Get("w")

global_observables = ws.set("globalObservables")
constrain = ws.set("nuisances")

pdf = ws["model_s"]
data = ws["data_obs"]

# To use the evaluation backends of RooFit
ROOT.gInterpreter.Declare("#include <HiggsAnalysis/CombinedLimit/interface/Combine.h>")

ROOT.Combine.nllBackend_ = args.backend

# Change verbosity
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)

ROOT.gInterpreter.Declare("#include <HiggsAnalysis/CombinedLimit/interface/CombineMathFuncs.h>")

nll = pdf.createNLL(data, Constrain=constrain, GlobalObservables=global_observables, Offset="initial")

ROOT.gInterpreter.Declare(
    """
void writeDebugMacro(RooAbsReal & real, std::string const &name)
{
   static_cast<RooFit::Experimental::RooFuncWrapper&>(real).writeDebugMacro(name);
}
"""
)

if args.backend == "codegen" and args.write_debug_macro:
    ROOT.writeDebugMacro(nll, "codegen")

cfg = ROOT.RooMinimizer.Config()
minim = ROOT.RooMinimizer(nll, cfg)
minim.setEps(1.0)
minim.setStrategy(0)
minim.minimize("Minuit2", "")
minim.save().Print()
