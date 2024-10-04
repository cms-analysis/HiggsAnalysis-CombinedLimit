import numpy as np
import ROOT
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="input ws file")
parser.add_argument("--cpu-backend", action="store_true", help="set EvalBackend=cpu, default is EvalBackend=codegen")

args = parser.parse_args()

with ROOT.TFile.Open(args.input) as f:
    ws = f.Get("w")

global_observables = ws.set("globalObservables")
constrain = ws.set("nuisances")

pdf = ws["model_s"]
data = ws["data_obs"]

# To use the evaluation backends of RooFit
pdf.useCombineNLL(False)

# Change verbosity
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)

ROOT.gInterpreter.Declare("#include <HiggsAnalysis/CombinedLimit/interface/CombineMathFuncs.h>")

if not args.cpu_backend:
    nll = pdf.createNLL(data, Constrain=constrain, GlobalObservables=global_observables, EvalBackend="codegen")
else:
    nll = pdf.createNLL(data, Constrain=constrain, GlobalObservables=global_observables, EvalBackend="cpu", Offset="initial")

cfg = ROOT.RooMinimizer.Config()
minim = ROOT.RooMinimizer(nll, cfg)
minim.setEps(1.0)
minim.setStrategy(0)
minim.minimize("Minuit2", "")
minim.save().Print()
