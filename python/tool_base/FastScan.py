#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import sys
import json
import ROOT
import re
import HiggsAnalysis.CombinedLimit.tool_base.utils as utils

from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase
import HiggsAnalysis.CombinedLimit.util.plotting as plot
from six.moves import range


class FastScan(CombineToolBase):
    description = "Calculate nuisance parameter impacts"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_intercept_args(self, group):
        CombineToolBase.attach_intercept_args(self, group)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument(
            "-w",
            "--workspace",
            required=True,
            help="Input ROOT file and workspace object name, in the format [file.root]:[name]. For workspaces produced by combine, the name is usually w.",
        )
        group.add_argument(
            "-d",
            "--data",
            help="By default reads data_obs from the input workspace. Alternative can be specified as [file.root]:[dataset name] or [file.root]:[wsp name]:[dataset name], where in both cases [dataset name] identifies an object inheriting from RooAbsData",
        )
        group.add_argument("-f", "--fitres", help="Optionally supply a RooFitResult to update the initial parameter values, format [file.root]:[RooFitResult]")
        group.add_argument("--match", help="Regular expression to only run for matching parameter names")
        group.add_argument("--no-match", help="Regular expression to skip certain parameter names")
        group.add_argument("-o", "--output", default="nll", help="Name of the output file, without the .pdf extension")
        group.add_argument("-p", "--points", default=200, type=int, help="Number of NLL points to sample in each scan")

    def run_method(self):
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        outfile = ROOT.TFile("%s.root" % self.args.output, "RECREATE")
        points = self.args.points
        file = ROOT.TFile(self.args.workspace.split(":")[0])
        wsp = file.Get(self.args.workspace.split(":")[1])
        mc = wsp.genobj("ModelConfig")
        pdf = mc.GetPdf()
        if self.args.data is None:
            data = wsp.data("data_obs")
        else:
            ws_d = self.args.data.split(":")
            print(">> Data: " + str(ws_d))
            f_d = ROOT.TFile(ws_d[0])
            if len(ws_d) == 2:
                data = f_d.Get(ws_d[1])
            else:
                data = f_d.Get(ws_d[1]).data(ws_d[2])
        ll = ROOT.RooLinkedList()
        nll = pdf.createNLL(data, ll)
        pars = pdf.getParameters(data)
        pars.Print()
        snap = pars.snapshot()
        # nll.setZeroPoint()
        nll.Print()
        if self.args.fitres is not None:
            fitfile = ROOT.TFile(self.args.fitres.split(":")[0])
            rfr = fitfile.Get(self.args.fitres.split(":")[1])
            snap = rfr.floatParsFinal()
        pars.assignValueOnly(snap)

        page = 0
        doPars = []

        for par in pars:
            if par.isConstant():
                continue
            if self.args.match is not None:
                if not re.match(self.args.match, par.GetName()):
                    continue
            if self.args.no_match is not None:
                if re.match(self.args.no_match, par.GetName()):
                    continue
            par.Print()
            if not (par.hasMax() and par.hasMin()):
                print("Parameter does not have an associated range, skipping")
                continue
            doPars.append(par)
        plot.ModTDRStyle(width=700, height=1000)
        for idx, par in enumerate(doPars):
            print("%s : (%i/%i)" % (par.GetName(), idx + 1, len(doPars)))
            nlld1 = nll.derivative(par, 1)
            nlld2 = nll.derivative(par, 2)
            xmin = par.getMin()
            xmax = par.getMax()
            gr = ROOT.TGraph(points)
            grd1 = ROOT.TGraph(points)
            grd2 = ROOT.TGraph(points)
            gr.SetName(par.GetName())
            grd1.SetName(par.GetName() + "_d1")
            grd2.SetName(par.GetName() + "_d2")
            w = (xmax - xmin) / float(points)
            for i in range(points):
                x = xmin + (float(i) + 0.5) * w
                par.setVal(x)
                gr.SetPoint(i, x, nll.getVal())
                grd1.SetPoint(i, x, nlld1.getVal())
                grd2.SetPoint(i, x, nlld2.getVal())
            plot.ReZeroTGraph(gr, True)
            # plot.RemoveGraphYAbove(gr, 2.)
            # gr.Print()
            outfile.cd()
            gr.Write()
            grd1.Write()
            grd2.Write()
            pars.assignValueOnly(snap)
            canv = ROOT.TCanvas(self.args.output, self.args.output)
            pads = plot.MultiRatioSplit([0.4, 0.3], [0.005, 0.005], [0.005, 0.005])
            pads[0].cd()
            plot.Set(gr, MarkerSize=0.5)
            gr.Draw("APL")
            axis1 = plot.GetAxisHist(pads[0])
            axis1.GetYaxis().SetTitle("NLL")
            pads[1].cd()
            plot.Set(grd1, MarkerSize=0.5)
            grd1.Draw("APL")
            axis2 = plot.GetAxisHist(pads[1])
            axis2.GetYaxis().SetTitle("NLL'")
            pads[2].cd()
            plot.Set(grd2, MarkerSize=0.5)
            grd2.Draw("APL")
            axis3 = plot.GetAxisHist(pads[2])
            axis3.GetYaxis().SetTitle("NLL''")
            plot.Set(
                axis3.GetXaxis(),
                Title=par.GetName(),
                TitleSize=axis3.GetXaxis().GetTitleSize() * 0.5,
                TitleOffset=axis3.GetXaxis().GetTitleOffset() * 2,
            )
            extra = ""
            if page == 0:
                extra = "("
            if page == len(doPars) - 1:
                extra = ")"
            print(extra)
            canv.Print(".pdf%s" % extra)
            page += 1

        outfile.Write()
