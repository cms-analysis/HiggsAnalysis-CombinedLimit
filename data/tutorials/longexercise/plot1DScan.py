#!/usr/bin/env python
import ROOT
import math
from functools import partial
import CombineHarvester.CombineTools.plotting as plot
import json
import argparse
import os.path

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot.ModTDRStyle(width=700, l=0.13)
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(0.7)

NAMECOUNTER = 0

def read(scan, param, files, ycut):
    goodfiles = [f for f in files if plot.TFileIsGood(f)]
    limit = plot.MakeTChain(goodfiles, 'limit')
    graph = plot.TGraphFromTree(limit, param, '2*deltaNLL', 'quantileExpected > -1.5')
    graph.SetName(scan)
    graph.Sort()
    plot.RemoveGraphXDuplicates(graph)
    plot.RemoveGraphYAbove(graph, ycut)
    # graph.Print()
    return graph


def Eval(obj, x, params):
    return obj.Eval(x[0])


def BuildScan(scan, param, files, color, yvals, ycut):
    graph = read(scan, param, files, ycut)
    bestfit = None
    for i in xrange(graph.GetN()):
        if graph.GetY()[i] == 0.:
            bestfit = graph.GetX()[i]
    graph.SetMarkerColor(color)
    spline = ROOT.TSpline3("spline3", graph)
    global NAMECOUNTER
    func = ROOT.TF1('splinefn'+str(NAMECOUNTER), partial(Eval, spline), graph.GetX()[0], graph.GetX()[graph.GetN() - 1], 1)
    NAMECOUNTER += 1
    func.SetLineColor(color)
    func.SetLineWidth(3)
    assert(bestfit is not None)
    crossings = {}
    cross_1sig = None
    cross_2sig = None
    other_1sig = []
    other_2sig = []
    val = None
    val_2sig = None
    for yval in yvals:
        crossings[yval] = plot.FindCrossingsWithSpline(graph, func, yval)
        for cr in crossings[yval]:
            cr["contains_bf"] = cr["lo"] <= bestfit and cr["hi"] >= bestfit
    for cr in crossings[yvals[0]]:
        if cr['contains_bf']:
            val = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
            cross_1sig = cr
        else:
            other_1sig.append(cr)
    if len(yvals) > 1:
        for cr in crossings[yvals[1]]:
            if cr['contains_bf']:
                val_2sig = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
                cross_2sig = cr
            else:
                other_2sig.append(cr)
    else:
        val_2sig = (0., 0., 0.)
        cross_2sig = cross_1sig
    return {
        "graph"     : graph,
        "spline"    : spline,
        "func"      : func,
        "crossings" : crossings,
        "val"       : val,
        "val_2sig": val_2sig,
        "cross_1sig" : cross_1sig,
        "cross_2sig" : cross_2sig,
        "other_1sig" : other_1sig,
        "other_2sig" : other_2sig
    }

parser = argparse.ArgumentParser()

parser.add_argument('main', help='Main input file for the scan')
parser.add_argument('--y-cut', type=float, default=7., help='Remove points with y > y-cut')
parser.add_argument('--y-max', type=float, default=8., help='y-axis maximum')
parser.add_argument('--output', '-o', help='output name without file extension', default='scan')
parser.add_argument('--POI', help='use this parameter of interest', default='r')
parser.add_argument('--translate', default=None, help='json file with POI name translation')
parser.add_argument('--main-label', default='Observed', type=str, help='legend label for the main scan')
parser.add_argument('--main-color', default=1, type=int, help='line and marker color for main scan')
parser.add_argument('--others', nargs='*', help='add secondary scans processed as main: FILE:LABEL:COLOR')
parser.add_argument('--breakdown', help='do quadratic error subtraction using --others')
parser.add_argument('--logo', default='CMS')
parser.add_argument('--logo-sub', default='Internal')
args = parser.parse_args()

print '--------------------------------------'
print  args.output
print '--------------------------------------'

fixed_name = args.POI
if args.translate is not None:
    with open(args.translate) as jsonfile:
        name_translate = json.load(jsonfile)
    if args.POI in name_translate:
        fixed_name = name_translate[args.POI]

yvals = [1., 4.]


main_scan = BuildScan(args.output, args.POI, [args.main], args.main_color, yvals, args.y_cut)

other_scans = [ ]
other_scans_opts = [ ]
if args.others is not None:
    for oargs in args.others:
        splitargs = oargs.split(':')
        other_scans_opts.append(splitargs)
        other_scans.append(BuildScan(args.output, args.POI, [splitargs[0]], int(splitargs[2]), yvals, args.y_cut))


canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()
main_scan['graph'].SetMarkerColor(1)
main_scan['graph'].Draw('AP')

axishist = plot.GetAxisHist(pads[0])

axishist.SetMaximum(args.y_max)
axishist.GetYaxis().SetTitle("- 2 #Delta ln L")
axishist.GetXaxis().SetTitle("%s" % fixed_name)

new_min = axishist.GetXaxis().GetXmin()
new_max = axishist.GetXaxis().GetXmax()
mins = []
maxs = []
for other in other_scans:
    mins.append(other['graph'].GetX()[0])
    maxs.append(other['graph'].GetX()[other['graph'].GetN()-1])

if len(other_scans) > 0:
    if min(mins) < main_scan['graph'].GetX()[0]:
        new_min = min(mins) - (main_scan['graph'].GetX()[0] - new_min)
    if max(maxs) > main_scan['graph'].GetX()[main_scan['graph'].GetN()-1]:
        new_max = max(maxs) + (new_max - main_scan['graph'].GetX()[main_scan['graph'].GetN()-1])
    axishist.GetXaxis().SetLimits(new_min, new_max)

for other in other_scans:
    if args.breakdown is not None:
        other['graph'].SetMarkerSize(0.4)
    other['graph'].Draw('PSAME')

line = ROOT.TLine()
line.SetLineColor(16)
# line.SetLineStyle(7)
for yval in yvals:
    plot.DrawHorizontalLine(pads[0], line, yval)
    if (len(other_scans) == 0):
        for cr in main_scan['crossings'][yval]:
            if cr['valid_lo']: line.DrawLine(cr['lo'], 0, cr['lo'], yval)
            if cr['valid_hi']: line.DrawLine(cr['hi'], 0, cr['hi'], yval)

main_scan['func'].Draw('same')
for other in other_scans:
    if args.breakdown is not None:
        other['func'].SetLineStyle(2)
        other['func'].SetLineWidth(2)
    other['func'].Draw('SAME')



box = ROOT.TBox(axishist.GetXaxis().GetXmin(), 0.625*args.y_max, axishist.GetXaxis().GetXmax(), args.y_max)
box.Draw()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

crossings = main_scan['crossings']
val_nom = main_scan['val']
val_2sig = main_scan['val_2sig']

textfit = '%s = %.3f{}^{#plus %.3f}_{#minus %.3f}' % (fixed_name, val_nom[0], val_nom[1], abs(val_nom[2]))


pt = ROOT.TPaveText(0.59, 0.82 - len(other_scans)*0.08, 0.95, 0.91, 'NDCNB')
pt.AddText(textfit)

if args.breakdown is None:
    for i, other in enumerate(other_scans):
        textfit = '#color[%s]{%s = %.3f{}^{#plus %.3f}_{#minus %.3f}}' % (other_scans_opts[i][2], fixed_name, other['val'][0], other['val'][1], abs(other['val'][2]))
        pt.AddText(textfit)


if args.breakdown is not None:
    pt.SetX1(0.50)
    if len(other_scans) >= 3:
        pt.SetX1(0.19)
        pt.SetX2(0.88)
        pt.SetY1(0.66)
        pt.SetY2(0.82)
    breakdown = args.breakdown.split(',')
    v_hi = [val_nom[1]]
    v_lo = [val_nom[2]]
    for other in other_scans:
        v_hi.append(other['val'][1])
        v_lo.append(other['val'][2])
    assert(len(v_hi) == len(breakdown))
    textfit = '%s = %.3f' % (fixed_name, val_nom[0])
    for i, br in enumerate(breakdown):
        if i < (len(breakdown) - 1):
            if (abs(v_hi[i+1]) > abs(v_hi[i])):
                print 'ERROR SUBTRACTION IS NEGATIVE FOR %s HI' % br
                hi = 0.
            else:
                hi = math.sqrt(v_hi[i]*v_hi[i] - v_hi[i+1]*v_hi[i+1])
            if (abs(v_lo[i+1]) > abs(v_lo[i])):
                print 'ERROR SUBTRACTION IS NEGATIVE FOR %s LO' % br
                lo = 0.
            else:
                lo = math.sqrt(v_lo[i]*v_lo[i] - v_lo[i+1]*v_lo[i+1])
        else:
            hi = v_hi[i]
            lo = v_lo[i]
        textfit += '{}^{#plus %.3f}_{#minus %.3f}(%s)' % (hi, abs(lo), br)
    pt.AddText(textfit)


pt.SetTextAlign(11)
pt.SetTextFont(42)
pt.Draw()

plot.DrawCMSLogo(pads[0], args.logo, args.logo_sub, 11, 0.045, 0.035, 1.2,  cmsTextSize = 1.)

legend_l = 0.69
if len(other_scans) > 0:
    legend_l = legend_l - len(other_scans) * 0.04
legend = ROOT.TLegend(0.15, legend_l, 0.45, 0.78, '', 'NBNDC')
if len(other_scans) >= 3:
        legend = ROOT.TLegend(0.46, 0.83, 0.95, 0.93, '', 'NBNDC')
        legend.SetNColumns(2)

legend.AddEntry(main_scan['func'], args.main_label, 'L')
for i, other in enumerate(other_scans):
    legend.AddEntry(other['func'], other_scans_opts[i][1], 'L')
legend.Draw()

save_graph = main_scan['graph'].Clone()
save_graph.GetXaxis().SetTitle('%s = %.3f %+.3f/%+.3f' % (fixed_name, val_nom[0], val_nom[2], val_nom[1]))
outfile = ROOT.TFile(args.output+'.root', 'RECREATE')
outfile.WriteTObject(save_graph)
outfile.Close()
canv.Print('.pdf')
canv.Print('.png')

