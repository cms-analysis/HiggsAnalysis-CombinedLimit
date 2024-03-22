#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import argparse
import os
import re
import sys
import json
import math
import itertools
import stat
import glob
import ROOT
from array import array
from multiprocessing import Pool
from numpy import matrix
from numpy.linalg import solve
import HiggsAnalysis.CombinedLimit.tool_base.utils as utils
from HiggsAnalysis.CombinedLimit.tool_base.opts import OPTS

from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase
from six.moves import range
import ctypes


class ImpactsFromScans(CombineToolBase):
    description = "Calculate nuisance parameter impacts"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_intercept_args(self, group):
        CombineToolBase.attach_intercept_args(self, group)
        group.add_argument("--name", "-n", default="Test")
        group.add_argument("-m", "--mass", required=True)

    def get_fixed_results(self, file, POIs):
        """Extracts the output from the MultiDimFit singles mode
        Note: relies on the list of parameters that were run (scanned) being correct"""
        res = {}
        f = ROOT.TFile(file)
        if f is None or f.IsZombie():
            return None
        t = f.Get("limit")
        for i, evt in enumerate(t):
            if i != 1:
                continue
            for POI in POIs:
                res[POI] = getattr(evt, POI)
        return res

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        # group.add_argument('--offset', default=0, type=int,
        #     help='Start the loop over parameters with this offset (default: %(default)s)')
        # group.add_argument('--advance', default=1, type=int,
        #     help='Advance this many parameters each step in the loop (default: %(default)s')
        group.add_argument("--input-json", help=("json file and dictionary containing the fit values, of form file:key1:key2.."))
        group.add_argument("--do-fits", action="store_true", help=("Actually do the fits"))
        group.add_argument("--cov-method", choices=["full", "asymm"], default="full")
        group.add_argument("--cor-method", choices=["full", "asymm", "average"], default="full")
        group.add_argument("--asymm-vals", default="")
        group.add_argument("--do-generic", action="store_true")

    def run_method(self):
        mass = self.args.mass
        self.put_back_arg("mass", "-m")
        in_json = self.args.input_json.split(":")
        with open(in_json[0]) as jsonfile:
            js = json.load(jsonfile)
        for key in in_json[1:]:
            js = js[key]
        POIs = sorted([str(x) for x in js.keys()])
        print(POIs)
        for POI in POIs:
            if not self.args.do_fits:
                break
            arg_str = "-M MultiDimFit --algo fixed --saveInactivePOI 1 --floatOtherPOIs 1 -P %s" % POI
            cmd_hi = arg_str + " -n %s --fixedPointPOIs %s=%f" % (self.args.name + ".%s.Hi" % POI, POI, js[POI]["Val"] + js[POI]["ErrorHi"])
            cmd_lo = arg_str + " -n %s --fixedPointPOIs %s=%f" % (self.args.name + ".%s.Lo" % POI, POI, js[POI]["Val"] + js[POI]["ErrorLo"])
            self.job_queue.append("combine %s %s" % (cmd_hi, " ".join(self.passthru)))
            self.job_queue.append("combine %s %s" % (cmd_lo, " ".join(self.passthru)))
        self.flush_queue()
        if self.args.do_fits:
            print(">> Re-run without --do-fits to harvest the results")
            return
        res = {}
        for POI in POIs:
            res[POI] = {}
            name_hi = "higgsCombine%s.%s.Hi.MultiDimFit.mH%s.root" % (self.args.name, POI, mass)
            name_lo = "higgsCombine%s.%s.Lo.MultiDimFit.mH%s.root" % (self.args.name, POI, mass)
            res_hi = self.get_fixed_results(name_hi, POIs)
            res_lo = self.get_fixed_results(name_lo, POIs)
            for fPOI in POIs:
                res[POI][fPOI] = [res_lo[fPOI], js[fPOI]["Val"], res_hi[fPOI]]
        # print res
        cor = ROOT.TMatrixDSym(len(POIs))
        cov = ROOT.TMatrixDSym(len(POIs))
        bf_vals = {x.split("=")[0]: float(x.split("=")[1]) for x in self.args.asymm_vals.split(",") if x != ""}

        xvars = []
        muvars = []
        covvars = []
        xvec = ROOT.RooArgList()
        mu = ROOT.RooArgList()
        for POI in POIs:
            xvars.append(ROOT.RooRealVar(POI, "", js[POI]["Val"], -100, 100))
            muvars.append(ROOT.RooRealVar(POI + "_In", "", js[POI]["Val"], -100, 100))
            muvars[-1].setConstant(True)
            xvec.add(xvars[-1])
            mu.add(muvars[-1])

        print("-----------------------------------------------------------")
        print("Diagonal Covariance")
        print("-----------------------------------------------------------")
        print("%-30s %-7s %-7s %-7s %-7s %-7s" % ("POI", "Val", "Sym", "Hi", "Lo", "(Hi-Lo)/(Hi+Lo)"))
        for i, p in enumerate(POIs):
            cor[i][i] = ctypes.c_double(1.0)  # diagonal correlation is 1
            d1 = res[p][p][1]
            d21 = res[p][p][2] - res[p][p][1]
            d10 = res[p][p][1] - res[p][p][0]
            d20 = (res[p][p][2] - res[p][p][0]) / 2.0
            vlo = js[p]["ValidErrorLo"]
            print("%-30s %+.3f  %+.3f  %+.3f  %+.3f %+.3f" % (p, d1, d20, d21, d10, (d21 - d10) / (d21 + d10)))
            covv = 0.0
            if self.args.cov_method == "full":
                covv = d20
            elif self.args.cov_method == "asymm":
                bf_val = 1.0
                for x in bf_vals:
                    if x in p:
                        bf_val = bf_vals[x]
                        print("Using %s=%g" % (x, bf_vals[x]))
                covv = d21 if bf_val >= d1 else d10
            if p == "mu_XS_ZH_BR_WW":
                covv = covv * 0.89
            if p == "mu_XS_ttHtH_BR_tautau":
                covv = covv * 1.2
            # if p == 'mu_XS_ttHtH_BR_tautau': covv = 6.3
            if not vlo:
                print("No ValidErrorLo, using d21")
                covv = d21
            print("Chosen: %+.3f" % covv)
            cov[i][i] = ctypes.c_double(pow(covv, 2.0))

            x1 = -1.0 * d10
            x2 = 0.0
            x3 = d21
            x4 = js[p]["2sig_ErrorHi"]
            y1 = d10 * d10
            y2 = d20 * d20
            y3 = d21 * d21
            y4 = (x4 / 2.0) * (x4 / 2.0)
            if not vlo and abs(d10) < 1e-4:
                x1 = -1.0 * d21
                y1 = d21 * d21
            print((x1, y1))
            print((x2, y2))
            print((x3, y3))
            print((x4, y4))

            mtx = matrix([[x1 * x1, x1, 1], [x3 * x3, x3, 1], [x4 * x4, x4, 1]])
            yvec = matrix([[y1], [y3], [y4]])
            # print mtx
            # print yvec
            xres = solve(mtx, yvec)
            # print xres
            covvars.append(
                ROOT.RooFormulaVar("cov%i" % i, "", "%g*(@0-%g)*(@0-%g)+%g*(@0-%g)+%g" % (xres[0], d1, d1, xres[1], d1, xres[2]), ROOT.RooArgList(xvars[i]))
            )
            # covvars.append(ROOT.RooFormulaVar('cov%i'%i,'', '%g' % (y2), ROOT.RooArgList()))
            covvars[-1].Print()

        print("-----------------------------------------------------------")
        print("Correlation")
        print("-----------------------------------------------------------")
        print(
            "%-30s %-30s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s"
            % ("i", "j", "Val_i", "Val_j", "ij_Sym", "ij_Hi", "ij_Lo", "ji_Sym", "ji_Hi", "ji_Lo", "Sym_Asym")
        )
        cors = []
        mvals = ROOT.RooArgList()
        mvals_store = []
        for i, ip in enumerate(POIs):
            for j, jp in enumerate(POIs):
                if i == j:
                    mvals_store.append(ROOT.RooFormulaVar("ele_%i_%i" % (i, j), "@0", ROOT.RooArgList(covvars[i])))
                    mvals.add(mvals_store[-1])
                    continue
                # Check the scans
                di_1 = res[ip][ip][1]
                di_21 = res[ip][ip][2] - res[ip][ip][1]
                di_10 = res[ip][ip][1] - res[ip][ip][0]
                di_20 = (res[ip][ip][2] - res[ip][ip][0]) / 2.0
                cj_21 = res[ip][jp][2] - res[ip][jp][1]
                cj_10 = res[ip][jp][1] - res[ip][jp][0]
                cj_20 = (res[ip][jp][2] - res[ip][jp][0]) / 2.0
                vi_lo = js[ip]["ValidErrorLo"]
                dj_1 = res[jp][jp][1]
                dj_21 = res[jp][jp][2] - res[jp][jp][1]
                dj_10 = res[jp][jp][1] - res[jp][jp][0]
                dj_20 = (res[jp][jp][2] - res[jp][jp][0]) / 2.0
                ci_21 = res[jp][ip][2] - res[jp][ip][1]
                ci_10 = res[jp][ip][1] - res[jp][ip][0]
                ci_20 = (res[jp][ip][2] - res[jp][ip][0]) / 2.0
                vj_lo = js[jp]["ValidErrorLo"]

                cij_20 = ci_20 / di_20
                cij_21 = ci_21 / (di_21 if (ci_21 >= 0 or not vi_lo) else di_10)
                cij_10 = ci_10 / (di_21 if (ci_21 < 0 or not vi_lo) else di_10)
                # cij_21 = ci_21/di_21
                # cij_10 = ci_10/di_21

                cji_20 = cj_20 / dj_20
                cji_21 = cj_21 / (dj_21 if (cj_21 >= 0 or not vj_lo) else dj_10)
                cji_10 = cj_10 / (dj_21 if (cj_21 < 0 or not vj_lo) else dj_10)
                # cji_21 = cj_21/dj_21
                # cji_10 = cj_10/dj_21

                a_20 = (cij_20 - cji_20) / ((cij_20 + cji_20) if (cij_20 + cji_20) != 0.0 else 1.0)

                a_i = (cij_21 - cij_10) / ((cij_21 + cij_10) if (cij_21 + cij_10) != 0.0 else 1.0)
                a_j = (cji_21 - cji_10) / ((cji_21 + cji_10) if (cji_21 + cji_10) != 0.0 else 1.0)

                max_c = max([abs(x) for x in [cij_20, cij_21, cij_10, cji_20, cji_21, cji_10]])

                line = "%-30s %-30s %+.3f  %+.3f | %+.3f  %+.3f  %+.3f  %+.3f |  %+.3f  %+.3f  %+.3f %+.3f |  %+.3f" % (
                    ip,
                    jp,
                    di_1,
                    dj_1,
                    cij_20,
                    cij_21,
                    cij_10,
                    a_i,
                    cji_20,
                    cji_21,
                    cji_10,
                    a_j,
                    a_20,
                )
                print(line)

                cors.append((line, max_c))

                val_i = 0.0
                val_j = 0.0
                if self.args.cor_method == "full":
                    val_i = cij_20
                    val_j = cji_20
                elif self.args.cor_method == "average":
                    val_i = (cij_21 + cij_10) / 2.0
                    val_j = (cji_21 + cji_10) / 2.0
                elif self.args.cor_method == "asymm":
                    bf_val_i = 1.0
                    bf_val_j = 1.0
                    for x in bf_vals:
                        if x in ip:
                            bf_val_i = bf_vals[x]
                            print("Using %s=%g for POI i" % (x, bf_vals[x]))
                        if x in jp:
                            bf_val_j = bf_vals[x]
                            print("Using %s=%g for POI j" % (x, bf_vals[x]))

                    val_i = cji_21 if bf_val_i >= di_1 else cji_10
                    val_j = cij_21 if bf_val_j >= dj_1 else cij_10
                if not vi_lo:
                    print("No ValidErrorLo for POI i, using d21")
                    val_i = cji_21
                if not vj_lo:
                    print("No ValidErrorLo for POI j, using d21")
                    val_j = cij_21
                print("Chosen: %+.3f for val_i" % val_i)
                print("Chosen: %+.3f for val_j" % val_j)

                correlation = (val_i + val_j) / 2.0  # take average correlation?
                # if ip == 'mu_XS_ttHtH_BR_WW' and jp == 'mu_XS_ttHtH_BR_tautau': correlation = correlation * 1.15
                # if jp == 'mu_XS_ttHtH_BR_WW' and ip == 'mu_XS_ttHtH_BR_tautau': correlation = correlation * 1.15
                # correlation = min(sorted([val_i, val_j],key=lambda x: abs(x), reverse=True))
                # correlation = min(val_i,val_j, key=abs) # take the max?
                cor[i][j] = correlation
                cor[j][i] = correlation
                covariance = correlation * math.sqrt(cov[i][i]) * math.sqrt(cov[j][j])
                cov[i][j] = covariance
                cov[j][i] = covariance
                mvals_store.append(ROOT.RooFormulaVar("ele_%i_%i" % (i, j), "%g*sqrt(@0)*sqrt(@1)" % (correlation), ROOT.RooArgList(covvars[i], covvars[j])))
                # mvals_store.append(ROOT.RooFormulaVar('ele_%i_%i'%(i,j),'%g'%(covariance),ROOT.RooArgList()))
                mvals.add(mvals_store[-1])
        cors.sort(key=lambda tup: tup[1], reverse=True)
        for tup in cors:
            print(tup[0])
        # cor.Print()
        fout = ROOT.TFile("covariance_%s.root" % self.args.name, "RECREATE")
        fout.WriteTObject(cor, "cor")
        h_cor = self.fix_TH2(ROOT.TH2D(cor), POIs)
        fout.WriteTObject(h_cor, "h_cor")
        fout.WriteTObject(cov, "cov")
        h_cov = self.fix_TH2(ROOT.TH2D(cov), POIs)
        fout.WriteTObject(h_cov, "h_cov")

        xvec.Print("v")
        mu.Print("v")
        if self.args.do_generic:
            pdf = ROOT.RooGenericMultiVarGaussian("pdf", "", xvec, mu, mvals)
        else:
            pdf = ROOT.RooMultiVarGaussian("pdf", "", xvec, mu, cov)
        dat = ROOT.RooDataSet("global_obs", "", ROOT.RooArgSet(mu))
        dat.add(ROOT.RooArgSet(mu))
        pdf.Print()
        dat.Print()
        # fitRes = pdf.fitTo(dat, ROOT.RooFit.Minimizer('Minuit2', 'Migrad'), ROOT.RooFit.Hesse(True), ROOT.RooFit.Save(True))
        # fitRes.Print('v')
        wsp = ROOT.RooWorkspace("w", "")
        getattr(wsp, "import")(pdf)
        getattr(wsp, "import")(dat)
        wsp.Write()
        fout.Close()

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


#          self.job_queue.append('combine -M MultiDimFit -n _initialFit_%(name)s_POI_%(poi)s --algo singles --redefineSignalPOIs %(poistr)s --floatOtherPOIs 1 --saveInactivePOI 1 -P %(poi)s %(pass_str)s --altCommit' % vars())
#      else:
#        self.job_queue.append('combine -M MultiDimFit -n _initialFit_%(name)s --algo singles --redefineSignalPOIs %(poistr)s %(pass_str)s --altCommit' % vars())
#      self.flush_queue()
#      sys.exit(0)
#    initialRes = utils.get_singles_results('higgsCombine_initialFit_%(name)s.MultiDimFit.mH%(mh)s.root' % vars(), poiList, poiList)
#    if len(named) > 0:
#      paramList = named
#    else:
#      paramList = utils.list_from_workspace(ws, 'w', 'ModelConfig_NuisParams')
#    print 'Have nuisance parameters: ' + str(len(paramList))
#    prefit = utils.prefit_from_workspace(ws, 'w', paramList)
#    res = { }
#    res["POIs"] = []
#    res["params"] = []
#    # for poi in poiList:
#    #   res["POIs"].append({"name" : poi, "fit" : initialRes[poi][poi]})
#
#    missing = [ ]
#    for param in paramList:
#      pres = { }
#      # print 'Doing param ' + str(counter) + ': ' + param
#      if self.args.doFits:
#        self.job_queue.append('combine -M MultiDimFit -n _paramFit_%(name)s_%(param)s --algo singles --redefineSignalPOIs %(param)s,%(poistr)s -P %(param)s --floatOtherPOIs 1 --saveInactivePOI 1 %(pass_str)s --altCommit' % vars())
#      else:
#        paramScanRes = get_singles_results('higgsCombine_paramFit_%(name)s_%(param)s.MultiDimFit.mH%(mh)s.root' % vars(), [param], poiList + [param])
#        if paramScanRes is None:
#          missing.append(param)
#          continue
#        pres.update({"name" : param, "fit" : paramScanRes[param][param], "prefit" : prefit[param]})
#        for p in poiList:
#          pres.update({p : paramScanRes[param][p], 'impact_'+p : (paramScanRes[param][p][2] - paramScanRes[param][p][0])/2.})
#        res['params'].append(pres)
#    self.flush_queue()
#    jsondata = json.dumps(res, sort_keys=True, indent=2, separators=(',', ': '))
#    print jsondata
#    if self.args.output is not None:
#      with open(args.output, 'w') as out_file:
#        out_file.write(jsondata)
#    if len(missing) > 0:
#      print 'Missing inputs: ' + ','.join(missing)
