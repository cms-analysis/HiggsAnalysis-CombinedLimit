import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--model", help="One of our three models, RPV, SYY, or SHH", default="RPV")
parser.add_argument("--masses", nargs='+', help="Masses to make plots for")
parser.add_argument("--year", help="2016, 2017, or Combo", default="2016")
parser.add_argument("--data", type=int, default=0, help="Data (0), pseudo-data (1), or pseudo-data with signal (2)")
parser.add_argument("--poi", type=str, default="r", help="Default is r: can pass in a parameter that the profile was performed for")
parser.add_argument("--inFile", type=str, default="")
parser.add_argument("--basedir", help="Directory that contains the output Fit_data_year directories", default="OutputForFreezing")
args=parser.parse_args()

import sys 
sys.argv.append( '-b' ) #setting batch mode for ROOT
import ROOT


year = args.year
model = args.model
data = "Data"
if args.data == 1:
    data = "pseudoData"
elif args.data == 2:
    data = "pseudoDataS"
poi = args.poi
inFile = args.inFile
basedir = args.basedir

for mass in args.masses:
    model_mass = "%s_%s" % (model, mass)

    inputfilename = "%s/Fit_%s_%s/output-files/%s_%s/higgsCombineSCAN_r_wSig.MultiDimFit.mH%s.MODEL%s.root" % (basedir, data, year, model_mass, year, mass, model)
    if inFile != "": inputfilename = inFile         
    inputfile = ROOT.TFile.Open(inputfilename)

    limit = inputfile.Get("limit")
    limit.Draw("2*deltaNLL:"+poi)

    graph = ROOT.gROOT.FindObject("Graph").Clone()
    graph.Sort()
    graph.SetLineWidth(3)
    graph.SetLineColor(1)
    graph.SetMaximum(8)
    graph.SetMinimum(ROOT.TMath.MinElement(graph.GetN(),graph.GetY()))
    
    canvas = ROOT.TCanvas("c","c", 600, 600)
    canvas.cd()
    canvas.SetBottomMargin(0.12)
    graph.Draw("AC")
    graph.GetXaxis().SetTitleSize(0.05)
    graph.GetXaxis().SetLabelSize(0.04)
    graph.GetXaxis().SetTitle("parameter "+poi)
    graph.GetYaxis().SetTitle("-2#Delta lnL")
    graph.GetYaxis().SetTitleSize(0.05)
    graph.GetYaxis().SetLabelSize(0.04)
    graph.SetTitle("Profile scan for %s in %s %s" % (model_mass, year, data))

    outputdir = "%s/Fit_%s_%s/output-files/%s_%s" % (basedir, data, year, model_mass, year)
    if inFile != "": outputdir = "."         
    canvas.SaveAs(outputdir + "/profilescan_%s_%s_%s_%s.png" % (model_mass, year, data, poi))
