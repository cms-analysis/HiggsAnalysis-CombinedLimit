#!/usr/bin/env python3
from __future__ import absolute_import
from __future__ import print_function
import HiggsAnalysis.CombinedLimit.util.plotting as plot
import CombineHarvester.CombineTools.maketable as maketable
import ROOT
import math
import argparse
import json
from six.moves import range

ROOT.gROOT.SetBatch(ROOT.kTRUE)
parser = argparse.ArgumentParser()
parser.add_argument("--file", "-f", help="named input file")
parser.add_argument("--process", help="The process on which a limit has been calculated. [gg#phi, bb#phi]", default="gg#phi")
parser.add_argument("--custom_y_range", help="Fix y axis range", action="store_true", default=False)
parser.add_argument("--y_axis_min", help="Fix y axis minimum", default=0.001)
parser.add_argument("--y_axis_max", help="Fix y axis maximum", default=100.0)
parser.add_argument("--custom_x_range", help="Fix x axis range", action="store_true", default=False)
parser.add_argument("--x_axis_min", help="Fix x axis minimum", default=90.0)
parser.add_argument("--x_axis_max", help="Fix x axis maximum", default=1000.0)
parser.add_argument("--verbosity", "-v", help="verbosity", default=0)
parser.add_argument("--log", help="Set log range for x and y axis", action="store_true", default=False)
parser.add_argument("--expected_only", help="Plot expected only", action="store_true", default=False)
parser.add_argument("--outname", "-o", help="Output plot name", default="mssm_limit")
parser.add_argument("--title", help="Plot title", default="19.8 fb^{-1} (8 TeV)")
# parser.add_argument('--table_vals', help='Amount of values to be written in a table for different masses', default=10)
args = parser.parse_args()


# Store the mass list convert from json file or directly via tgraphs
graph_obs = ROOT.TGraph()
graph_minus2sigma = ROOT.TGraph()
graph_minus1sigma = ROOT.TGraph()
graph_exp = ROOT.TGraph()
graph_plus1sigma = ROOT.TGraph()
graph_plus2sigma = ROOT.TGraph()

if ".root" in args.file:
    file = ROOT.TFile(args.file, "r")
    if not args.expected_only:
        graph_obs = plot.SortGraph(file.Get("observed"))
    graph_minus2sigma = plot.SortGraph(file.Get("minus2sigma"))
    graph_minus1sigma = plot.SortGraph(file.Get("minus1sigma"))
    graph_exp = plot.SortGraph(file.Get("expected"))
    graph_plus1sigma = plot.SortGraph(file.Get("plus1sigma"))
    graph_plus2sigma = plot.SortGraph(file.Get("plus2sigma"))
    maketable.Tablefrom1DGraph(args.file, "mssm_limit_table.txt")
else:
    data = {}
    with open(args.file) as jsonfile:
        data = json.load(jsonfile)
    if not args.expected_only:
        graph_obs = plot.LimitTGraphFromJSON(data, "observed")
    graph_minus2sigma = plot.LimitTGraphFromJSON(data, "-2")
    graph_minus1sigma = plot.LimitTGraphFromJSON(data, "-1")
    graph_exp = plot.LimitTGraphFromJSON(data, "expected")
    graph_plus1sigma = plot.LimitTGraphFromJSON(data, "+1")
    graph_plus2sigma = plot.LimitTGraphFromJSON(data, "+2")
    maketable.TablefromJson(args.file, "mssm_limit_table.txt")

process_label = args.process

mass_list = []
for i in range(graph_exp.GetN()):
    mass_list.append(float(graph_exp.GetX()[i]))
mass_list = sorted(set(mass_list))
mass_bins = len(mass_list)
if int(args.verbosity) > 0:
    print("mass_list: ", mass_list, "Total number: ", mass_bins)

# Create canvas and TH1D
plot.ModTDRStyle(width=600, l=0.12)
ROOT.gStyle.SetFrameLineWidth(2)
c1 = ROOT.TCanvas()
axis = plot.makeHist1D("hist1d", mass_bins, graph_exp)
if process_label == "gg#phi":
    axis.GetYaxis().SetTitle("95% CL limit on #sigma#font[42]{(gg#phi)}#upoint#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} [pb]")
elif process_label == "bb#phi":
    axis.GetYaxis().SetTitle("95% CL limit on #sigma#font[42]{(bb#phi)}#upoint#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} [pb]")
else:
    exit("Currently process is not supported")
if args.custom_y_range:
    axis.GetYaxis().SetRangeUser(float(args.y_axis_min), float(args.y_axis_max))
axis.GetXaxis().SetTitle("m_{#phi} [GeV]")
if args.custom_x_range:
    axis.GetXaxis().SetRangeUser(float(args.x_axis_min), float(args.x_axis_max))
# Create two pads, one is just for the Legend
pad_leg = ROOT.TPad("pad_leg", "pad_leg", 0, 0.82, 1, 1)
pad_leg.SetFillStyle(4000)
pad_leg.Draw()
pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 1, 0.82)
pad_plot.SetFillStyle(4000)
pad_plot.Draw()
pads = [pad_leg, pad_plot]
pads[1].cd()
if args.log:
    pad_plot.SetLogx(1)
    pad_plot.SetLogy(1)
    axis.SetNdivisions(50005, "X")
    axis.GetXaxis().SetMoreLogLabels()
    axis.GetXaxis().SetNoExponent()
    axis.GetXaxis().SetLabelSize(0.040)
axis.Draw()

innerBand = plot.MakeErrorBand(graph_minus1sigma, graph_plus1sigma)
outerBand = plot.MakeErrorBand(graph_minus2sigma, graph_plus2sigma)

outerBand.SetLineWidth(1)
outerBand.SetLineColor(ROOT.kBlack)
#    if(injected) outerBand->SetFillColor(kAzure-9);
#    else if(BG_Higgs) outerBand->SetFillColor(kSpring+5);
outerBand.SetFillColor(ROOT.TColor.GetColor(252, 241, 15))
outerBand.Draw("3")

innerBand.SetLineWidth(1)
innerBand.SetLineColor(ROOT.kBlack)
#      if(injected) innerBand->SetFillColor(kAzure-4);
#      else if(BG_Higgs) innerBand->SetFillColor(kGreen+2);
innerBand.SetFillColor(ROOT.kGreen)
innerBand.Draw("3same")

graph_exp.SetLineColor(ROOT.kRed)
graph_exp.SetLineWidth(3)
graph_exp.SetLineStyle(1)
#  if(mssm_log){
#      expected->SetLineColor(kBlack);
#      expected->SetLineStyle(2);
#      }
graph_exp.Draw("L")

if not args.expected_only:
    graph_obs.SetMarkerColor(ROOT.kBlack)
    graph_obs.SetMarkerSize(1.0)
    graph_obs.SetMarkerStyle(20)
    graph_obs.SetLineWidth(3)
    graph_obs.Draw("PLsame")

pads[0].cd()
legend = plot.PositionedLegend(0.5, 0.9, 2, 0.03)
legend.SetNColumns(2)
legend.SetFillStyle(1001)
legend.SetTextSize(0.15)
legend.SetTextFont(62)
legend.SetHeader("95% CL Excluded:")
if not args.expected_only:
    legend.AddEntry(graph_obs, "Observed", "L")
legend.AddEntry(innerBand, "#pm 1#sigma Expected", "F")
legend.AddEntry(graph_exp, "Expected", "L")
legend.AddEntry(outerBand, "#pm 2#sigma Expected", "F")
legend.Draw("same")

plot.DrawCMSLogo(pads[1], "", "", 11, 0.045, 0.035, 1.2)
plot.DrawTitle(pads[1], "%s" % args.title, 3)
plot.FixOverlay()
c1.SaveAs("%s.pdf" % args.outname)
c1.SaveAs("%s.png" % args.outname)
