#!/usr/bin/env python3
from __future__ import absolute_import
from __future__ import print_function
import ROOT
import HiggsAnalysis.CombinedLimit.util.plotting as plot
import argparse

# import CombineHarvester.CombineTools.maketable as maketable

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+", help="""Input json files""")
parser.add_argument(
    "--output",
    "-o",
    default="limit",
    help="""Name of the output
    plot without file extension""",
)
parser.add_argument("--show", default="exp,obs")
# parser.add_argument(
#     '--debug-output', '-d', help="""If specified, write the
#     TGraphs into this output ROOT file""")
parser.add_argument("--x-title", default="m_{H} (GeV)", help="""Title for the x-axis""")
parser.add_argument("--y-title", default=None, help="""Title for the y-axis""")
parser.add_argument("--limit-on", default="#sigma/#sigma_{SM}", help="""Shortcut for setting the y-axis label""")
parser.add_argument("--cms-sub", default="Internal", help="""Text below the CMS logo""")
parser.add_argument(
    "--scenario-label",
    default="",
    help="""Scenario name to be drawn in top
    left of plot""",
)
parser.add_argument("--title-right", default="", help="""Right header text above the frame""")
parser.add_argument("--title-left", default="", help="""Left header text above the frame""")
parser.add_argument("--logy", action="store_true", help="""Draw y-axis in log scale""")
parser.add_argument("--logx", action="store_true", help="""Draw x-axis in log scale""")
parser.add_argument("--ratio-to", default=None)
parser.add_argument("--pad-style", default=None, help="""Extra style options for the pad, e.g. Grid=(1,1)""")
parser.add_argument("--auto-style", nargs="?", const="", default=None, help="""Take line colors and styles from a pre-defined list""")
parser.add_argument("--table_vals", help="Amount of values to be written in a table for different masses", default=10)
args = parser.parse_args()


def DrawAxisHists(pads, axis_hists, def_pad=None):
    for i, pad in enumerate(pads):
        pad.cd()
        axis_hists[i].Draw("AXIS")
        axis_hists[i].Draw("AXIGSAME")
    if def_pad is not None:
        def_pad.cd()


## Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()
ROOT.gStyle.SetNdivisions(510, "XYZ")  # probably looks better

canv = ROOT.TCanvas(args.output, args.output)

if args.ratio_to is not None:
    pads = plot.TwoPadSplit(0.30, 0.01, 0.01)
else:
    pads = plot.OnePad()

# Set the style options of the pads
for padx in pads:
    # Use tick marks on oppsite axis edges
    plot.Set(padx, Tickx=1, Ticky=1, Logx=args.logx)
    if args.pad_style is not None:
        settings = {x.split("=")[0]: eval(x.split("=")[1]) for x in args.pad_style.split(",")}
        print("Applying style options to the TPad(s):")
        print(settings)
        plot.Set(padx, **settings)

graphs = []
graph_sets = []

legend = plot.PositionedLegend(0.45, 0.10, 3, 0.015)
plot.Set(legend, NColumns=2)

axis = None

defcols = [
    ROOT.kGreen + 3,
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kBlack,
    ROOT.kYellow + 2,
    ROOT.kOrange + 10,
    ROOT.kCyan + 3,
    ROOT.kMagenta + 2,
    ROOT.kViolet - 5,
    ROOT.kGray,
]

deflines = [1, 2, 3]

if args.auto_style is not None:
    icol = {x: 0 for x in args.auto_style.split(",")}
    icol["default"] = 0
    iline = {}
    iline["default"] = 1
    for i, x in enumerate(args.auto_style.split(",")):
        iline[x] = i + 1

# Process each input argument
for src in args.input:
    splitsrc = src.split(":")
    file = splitsrc[0]
    # limit.json => Draw as full obs + exp limit band
    if len(splitsrc) == 1:
        graph_sets.append(plot.StandardLimitsFromJSONFile(file, args.show.split(",")))
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), list(graph_sets[-1].values())[0], True)
            DrawAxisHists(pads, axis, pads[0])
        plot.StyleLimitBand(graph_sets[-1])
        plot.DrawLimitBand(pads[0], graph_sets[-1], legend=legend)
        pads[0].RedrawAxis()
        pads[0].RedrawAxis("g")
        pads[0].GetFrame().Draw()

    # limit.json:X => Draw a single graph for entry X in the json file
    # 'limit.json:X:Title="Blah",LineColor=4,...' =>
    # as before but also apply style options to TGraph
    elif len(splitsrc) >= 2:
        settings = {}
        settings["Title"] = src
        if args.auto_style is not None:
            nm = "default"
            for x in icol.keys():
                if x in splitsrc[1]:
                    nm = x
            i = icol[nm]  # take the next default color...
            j = iline[nm]  # take the next default line style...
            settings["LineColor"] = defcols[i]
            settings["MarkerColor"] = defcols[i]
            settings["LineStyle"] = j
            icol[nm] = (i + 1) if (i + 1) < len(defcols) else 0
        graphs.append(plot.LimitTGraphFromJSONFile(file, splitsrc[1]))
        if len(splitsrc) >= 3:
            settings.update({x.split("=")[0]: eval(x.split("=")[1]) for x in splitsrc[2].split(",")})
        plot.Set(graphs[-1], **settings)
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graphs[-1], True)
            DrawAxisHists(pads, axis, pads[0])
        graphs[-1].Draw("PLSAME")
        legend.AddEntry(graphs[-1], "", "PL")


axis[0].GetYaxis().SetTitle("95%% CL limit on %s" % args.limit_on)
if args.y_title is not None:
    axis[0].GetYaxis().SetTitle(args.y_title)
axis[0].GetXaxis().SetTitle(args.x_title)
axis[0].GetXaxis().SetLabelOffset(axis[0].GetXaxis().GetLabelOffset() * 2)

if args.logy:
    axis[0].SetMinimum(0.1)  # we'll fix this later
    pads[0].SetLogy(True)
    # axis[0].GetYaxis().SetMoreLogLabels()
    # axis[0].SetNdivisions(50005, "X")

y_min, y_max = (plot.GetPadYMin(pads[0]), plot.GetPadYMax(pads[0]))
plot.FixBothRanges(pads[0], y_min if args.logy else 0, 0.05 if args.logy else 0, y_max, 0.25)

ratio_graph_sets = []
ratio_graphs = []

if args.ratio_to is not None:
    pads[1].cd()
    plot.SetupTwoPadSplitAsRatio(pads, axis[0], axis[1], "Ratio_{}", True, 0.1, 2.4)
    axis[1].SetNdivisions(506, "Y")
    splitsrc = args.ratio_to.split(":")
    ref = plot.LimitTGraphFromJSONFile(splitsrc[0], splitsrc[1])
    for gr_set in graph_sets:
        ratio_set = {}
        for key in gr_set:
            ratio_set[key] = plot.GraphDivide(gr_set[key], ref)
        ratio_graph_sets.append(ratio_set)
        plot.DrawLimitBand(pads[1], ratio_graph_sets[-1])
        pads[1].RedrawAxis()
        pads[1].RedrawAxis("g")
        pads[1].GetFrame().Draw()
    for gr in graphs:
        ratio_graphs.append(plot.GraphDivide(gr, ref))
        ratio_graphs[-1].Draw("LP")
    ry_min, ry_max = (plot.GetPadYMin(pads[1]), plot.GetPadYMax(pads[1]))
    plot.FixBothRanges(pads[1], ry_min, 0.1, ry_max, 0.1)


pads[0].cd()
if legend.GetNRows() == 1:
    legend.SetY1(legend.GetY2() - 0.5 * (legend.GetY2() - legend.GetY1()))
legend.Draw()

# line = ROOT.TLine()
# line.SetLineColor(ROOT.kBlue)
# line.SetLineWidth(2)
# plot.DrawHorizontalLine(pads[0], line, 1)

box = ROOT.TPave(pads[0].GetLeftMargin(), 0.81, 1 - pads[0].GetRightMargin(), 1 - pads[0].GetTopMargin(), 1, "NDC")
box.Draw()

legend.Draw()

plot.DrawCMSLogo(pads[0], "CMS", args.cms_sub, 11, 0.045, 0.035, 1.2, "", 0.8)
plot.DrawTitle(pads[0], args.title_right, 3)
plot.DrawTitle(pads[0], args.title_left, 1)

canv.Print(".pdf")
canv.Print(".png")
# maketable.TablefromJson(args.table_vals, args.file, "TablefromJson.txt")
