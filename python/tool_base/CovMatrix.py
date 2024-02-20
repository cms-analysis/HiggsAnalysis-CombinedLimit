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
    description = 'Build a fit covariance matrix from scan results'
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument('-i', '--input', nargs='+', default=[],
                           help='The input file containing the MultiDimFit singles mode output')
        group.add_argument(
            '-o', '--output', help='The output name in the format file:prefix')
        group.add_argument(
            '-P', '--POIs', help='The params that were scanned (in scan order)')
        group.add_argument(
            '--POIs-from-set', help='Extract from file:workspace:set instead')
        group.add_argument('--compare', help='Compare to RooFitResult')

    def run_method(self):
        POIs = []
        if args.POIs is not None:
            POIs = args.POIs.split(',')
        if args.POIs_from_set is not None:
            ws_in = args.POIs_from_set.split(':')
            print(ws_in)
            POIs = list_from_workspace(ws_in[0], ws_in[1], ws_in[2])
        # res = { }
        # if len(args.input) == 1:
        #   res.update(get_singles_results(args.input, POIs, POIs))
        # elif len(args.input) > 1:
        #   assert len(args.input) == len(POIs)
        #   for i in range(len(POIs)):
        #     res.update(get_singles_results(args.input[i], [POIs[i]], POIs))
        # for p in POIs:
        #   val = res[p][p]
        #   print '%s = %.3f -%.3f/+%.3f' % (p, val[1], val[1] - val[0], val[2] - val[1])
        # print res
        # cor = ROOT.TMatrixDSym(len(POIs))
        # cov = ROOT.TMatrixDSym(len(POIs))
        # for i,p in enumerate(POIs):
        # cor[i][i] = ctypes.c_double(1.) # diagonal correlation is 1
        # cov[i][i] = ctypes.c_double(pow((res[p][p][2] - res[p][p][0])/2.,2.)) # symmetrized error
        # for i,ip in enumerate(POIs):
        #   for j,jp in enumerate(POIs):
        #     if i == j: continue
        #     val_i = ((res[ip][jp][2] - res[ip][jp][0])/2.)/math.sqrt(cov[j][j])
        #     val_j = ((res[jp][ip][2] - res[jp][ip][0])/2.)/math.sqrt(cov[i][i])
        # correlation = (val_i+val_j)/2. # take average correlation?
        # correlation = min(val_i,val_j, key=abs) # take the max?
        #     cor[i][j] = correlation
        #     cor[j][i] = correlation
        #     covariance = correlation * math.sqrt(cov[i][i]) * math.sqrt(cov[j][j])
        #     cov[i][j] = covariance
        #     cov[j][i] = covariance
        compare = args.compare is not None
        if compare:
            f_in = args.compare.split(':')
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
                    fitres_cor[i][j] = ctypes.c_double(
                        fitres_cor_src[ipos[i]][ipos[j]])
                    fitres_cov[i][j] = ctypes.c_double(
                        fitres_cov_src[ipos[i]][ipos[j]])
        # print 'My correlation matrix:'
        # cor.Print()
        if compare:
            print('RooFitResult correlation matrix:')
            fitres_cor.Print()
        # print 'My covariance matrix:'
        # cov.Print()
        if compare:
            print('RooFitResult covariance matrix:')
            fitres_cov.Print()
        if args.output is not None:
            out = args.output.split(':')
            fout = ROOT.TFile(out[0], 'RECREATE')
            prefix = out[1]
            # fout.WriteTObject(cor, prefix+'_cor')
            # h_cor = self.fix_TH2(ROOT.TH2D(cor), POIs)
            # fout.WriteTObject(h_cor, prefix+'_h_cor')
            # fout.WriteTObject(cov, prefix+'_cov')
            # h_cov = self.fix_TH2(ROOT.TH2D(cov), POIs)
            # fout.WriteTObject(h_cov, prefix+'_h_cov')
            if compare:
                fout.WriteTObject(fitres_cor, prefix + '_comp_cor')
                h_cor_compare = self.fix_TH2(ROOT.TH2D(fitres_cor), POIs)
                fout.WriteTObject(h_cor_compare, prefix + '_comp_h_cor')
                fout.WriteTObject(fitres_cov, prefix + '_comp_cov')
                h_cov_compare = self.fix_TH2(ROOT.TH2D(fitres_cov), POIs)
                fout.WriteTObject(h_cov_compare, prefix + '_comp_h_cov')

    def fix_TH2(self, h, labels):
        h_fix = h.Clone()
        for y in range(1, h.GetNbinsY() + 1):
            for x in range(1, h.GetNbinsX() + 1):
                h_fix.SetBinContent(
                    x, y, h.GetBinContent(x, h.GetNbinsY() + 1 - y))
        for x in range(1, h_fix.GetNbinsX() + 1):
            h_fix.GetXaxis().SetBinLabel(x, labels[x - 1])
        for y in range(1, h_fix.GetNbinsY() + 1):
            h_fix.GetYaxis().SetBinLabel(y, labels[-y])
        return h_fix
