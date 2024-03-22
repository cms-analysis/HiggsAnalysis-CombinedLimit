#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import ROOT

import HiggsAnalysis.CombinedLimit.tool_base.utils as utils
from HiggsAnalysis.CombinedLimit.tool_base.opts import OPTS

from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase
from six.moves import range
import ctypes


class CovMatrix(CombineToolBase):
    description = "Build a fit covariance matrix from scan results"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("-i", "--input", nargs="+", default=[], help="The input file containing the MultiDimFit singles mode output")
        group.add_argument("-o", "--output", help="The output name in the format file:prefix")
        group.add_argument("-P", "--POIs", help="The params that were scanned (in scan order)")
        group.add_argument("--POIs-from-set", help="Extract from file:workspace:set instead")
        group.add_argument("--compare", help="Compare to RooFitResult")

    def run_method(self):
        POIs = []
        if self.args.POIs is not None:
            POIs = self.args.POIs.split(",")
        if self.args.POIs_from_set is not None:
            ws_in = self.args.POIs_from_set.split(":")
            print(ws_in)
            POIs = utils.list_from_workspace(ws_in[0], ws_in[1], ws_in[2])

        compare = self.args.compare is not None
        if compare:
            f_in = self.args.compare.split(":")
            f = ROOT.TFile(f_in[0])
            fitres = f.Get(f_in[1])
            fitres_cov = ROOT.TMatrixDSym(len(POIs))
            fitres_cov_src = fitres.covarianceMatrix()
            fitres_cor = ROOT.TMatrixDSym(len(POIs))
            fitres_cor_src = fitres.correlationMatrix()
            ipos = []
            for p in POIs:
                ipos.append(fitres.floatParsFinal().index(p))
            for i, ip in enumerate(POIs):
                for j, jp in enumerate(POIs):
                    fitres_cor[i][j] = ctypes.c_double(fitres_cor_src[ipos[i]][ipos[j]])
                    fitres_cov[i][j] = ctypes.c_double(fitres_cov_src[ipos[i]][ipos[j]])

        if compare:
            print("RooFitResult correlation matrix:")
            fitres_cor.Print()

        if compare:
            print("RooFitResult covariance matrix:")
            fitres_cov.Print()
        if self.args.output is not None:
            out = self.args.output.split(":")
            fout = ROOT.TFile(out[0], "RECREATE")
            prefix = out[1]
            if compare:
                fout.WriteTObject(fitres_cor, prefix + "_comp_cor")
                h_cor_compare = self.fix_TH2(ROOT.TH2D(fitres_cor), POIs)
                fout.WriteTObject(h_cor_compare, prefix + "_comp_h_cor")
                fout.WriteTObject(fitres_cov, prefix + "_comp_cov")
                h_cov_compare = self.fix_TH2(ROOT.TH2D(fitres_cov), POIs)
                fout.WriteTObject(h_cov_compare, prefix + "_comp_h_cov")

    def fix_TH2(self, h, labels):
        h_fix = h.Clone()
        for y in range(1, h.GetNbinsY() + 1):
            for x in range(1, h.GetNbinsX() + 1):
                h_fix.SetBinContent(x, y, h.GetBinContent(x, h.GetNbinsY() + 1 - y))
        for x in range(1, h_fix.GetNbinsX() + 1):
            h_fix.GetXaxis().SetBinLabel(x, labels[x - 1])
        for y in range(1, h_fix.GetNbinsY() + 1):
            h_fix.GetYaxis().SetBinLabel(y, labels[-y])
        return h_fix
