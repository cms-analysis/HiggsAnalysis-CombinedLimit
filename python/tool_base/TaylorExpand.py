#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import os
import six.moves.cPickle as pickle
import math
import json
import itertools
import hashlib
import ROOT
from array import array
import numpy as np
from pprint import pprint
from functools import partial
from HiggsAnalysis.CombinedLimit.tool_base.opts import OPTS
from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase
import six
from six.moves import range
from six.moves import zip


def Eval(obj, x, params):
    return obj.Eval(x[0])


def Permutations(indicies):
    """Calculate the number of unique permutations of a set of indicies

    Args:
        indicies (list): A list of integer indicies, e.g [0, 1, 2, 3]

    Returns:
        int: number of unique permutations
    """
    uniques = sorted(set(indicies))
    freqs = [indicies.count(x) for x in uniques]
    n_perms = math.factorial(len(indicies))
    for freq in freqs:
        n_perms = n_perms / math.factorial(freq)
    return n_perms


def GenerateStencils(d, h, s):
    N = len(s)
    smatrix = np.zeros((N, N))
    dvec = np.zeros((N))
    for i in range(N):
        for j in range(N):
            smatrix[i, j] = pow(s[j], i)
        dvec[i] = math.factorial(d) if i == d else 0.0
    # print smatrix
    # print dvec
    res = (1.0 / pow(h, d)) * np.dot(np.linalg.inv(smatrix), dvec)
    return res


def GenerateDebugGraph(wsp, var_name, ingraph):
    nll = wsp.function("nll")
    var = wsp.var(var_name)
    outgraph = ingraph.Clone()
    for i in range(ingraph.GetN()):
        var.setVal(ingraph.GetX()[i])
        outgraph.GetY()[i] = nll.getVal() * 2.0
    return outgraph


# Calculate partial derivatives using finite differences
class ExpansionTerm:
    def __init__(self, parameter_values, derivatives, stencils):
        self.fnval = 0.0
        self.parameter_values = np.array(parameter_values, dtype=np.float32)
        self.derivatives = np.array(derivatives, dtype=np.uint8)
        unique_derivatives = np.array(sorted(set(derivatives)), dtype=np.uint8)
        derivative_frequency = np.zeros(len(unique_derivatives), dtype=np.uint8)
        for idx, i in enumerate(unique_derivatives):
            derivative_frequency[idx] = np.count_nonzero(self.derivatives == i)
        terms = list()
        # self.terms = list()
        if len(derivative_frequency):
            parameter = unique_derivatives[0]
            order = derivative_frequency[0]
            self.fundamental = False
            stencil = stencils[parameter][order]
            self.coeffs = np.zeros(len(stencil), dtype=np.float64)
            for i in range(len(stencil)):
                remaining_derivatives = np.array(list([a for a in self.derivatives if a != parameter]), dtype=np.uint8)
                # Make a copy of the current parameters and adjust the
                # value for the current stencil point
                new_parameter_values = np.array(self.parameter_values, dtype=np.float32)
                new_parameter_values[parameter] += stencil[i][0]
                # Add this to the list of terms
                self.coeffs[i] = stencil[i][1]
                terms.append(ExpansionTerm(new_parameter_values, remaining_derivatives, stencils))
            self.terms = np.array(terms)
        else:
            self.terms = []
            self.fundamental = True

    def FormattedPars(self):
        return tuple([float("%f" % p) for p in self.parameter_values])

    def Eval(self, with_permutations=False, with_factorial=False):
        if self.fundamental:
            return self.fnval
        else:
            summed = 0.0
            for i in range(len(self.terms)):
                summed += self.coeffs[i] * self.terms[i].Eval()
            if with_permutations:
                n_perms = Permutations(list(self.derivatives))
                summed *= float(n_perms)
            if with_factorial:
                summed *= 1.0 / float(math.factorial(len(self.derivatives)))
            return summed

    def Print(self, indent=0, coeff=None):
        sp = " " * indent
        extra = ""
        if self.fundamental:
            extra = " %s = %f" % (list(self.FormattedPars()), self.fnval)
        if coeff is None:
            print("%s%s%s" % (sp, self.derivatives, extra))
        else:
            print("%s%+.1f*%s%s" % (sp, coeff, self.derivatives, extra))
        for i in range(len(self.terms)):
            self.terms[i].Print(indent + 2, self.coeffs[i])

    def GatherFundamentalTerms(self, termlist):
        if self.fundamental:
            termlist.append(self)
        else:
            for t in self.terms:
                t.GatherFundamentalTerms(termlist)


class TaylorExpand(CombineToolBase):
    description = "Calculate nuisance parameter impacts"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)
        self.nll = None

    def attach_intercept_args(self, group):
        CombineToolBase.attach_intercept_args(self, group)
        group.add_argument("--name", "-n", default="Test")
        group.add_argument("-m", "--mass", required=True)
        group.add_argument("-d", "--datacard", required=True)

    def get_results(self, file):
        res = []
        f = ROOT.TFile(file)
        if f is None or f.IsZombie():
            return None
        t = f.Get("limit")
        for i, evt in enumerate(t):
            if i == 0:
                continue
            res.append(getattr(evt, "deltaNLL"))
        print(res)
        if len(res) == 0:
            print(file)
        return res

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("--config", help=("json configuration file"))
        group.add_argument("--order", type=int, default=2, help=("Taylor expand up to and including this order"))
        group.add_argument("--cross-order", type=int, default=2, help=("Taylor expand up to and including this order for the cross-order terms"))
        group.add_argument("--choose-POIs", default=None, help=("Explict list POIs to expand in"))
        group.add_argument("--do-fits", action="store_true", help=("Actually do the fits"))
        group.add_argument("--test-mode", type=int, default=0, help=("Test on the workspace"))
        group.add_argument("--test-args", type=str, default="", help=("List of comma separated args to be interpreted by the test-mode"))
        group.add_argument("--save", default=None, help=("Save results to a json file"))
        group.add_argument("--load", default=None, help=("Load results from a json file"))
        group.add_argument("--stencil-add", type=int, default=0, help=("Add additional points to each stencil"))
        group.add_argument("--stencil-min", type=int, default=3, help=("Minimum number of points in stencil"))
        group.add_argument("--drop-threshold", type=float, default=0.0, help=("Drop contributions below this threshold"))
        group.add_argument("--multiple", type=int, default=1, help=("Run multiple fixed points in one combine job"))
        group.add_argument("--workspace-bestfit", action="store_true", help=("Update the best-fit using the workspace snapshot"))
        group.add_argument("--linear-workaround", default=None, help=("Comma separated list of POIs that require special treatment due to a linear NLL"))

    def load_workspace(self, file, POIs, data="data_obs", snapshot="MultiDimFit"):
        if self.nll is not None:
            return
        print("Loading NLL...")
        self.infile = ROOT.TFile(file)
        self.loaded_wsp = self.infile.Get("w")
        # infile.Close()
        mc = self.loaded_wsp.genobj("ModelConfig")
        pdf = mc.GetPdf()
        data = self.loaded_wsp.data(data)
        ll = ROOT.RooLinkedList()
        self.nll = pdf.createNLL(data, ll)
        self.loaded_wsp.loadSnapshot("MultiDimFit")
        print("...NLL loaded")
        # nll.Print()
        self.nll0 = self.nll.getVal()
        self.wsp_vars = {}
        for POI in POIs:
            self.wsp_vars[POI] = self.loaded_wsp.var(POI)

    def get_snpashot_pois(self, file, POIs, snapshot="MultiDimFit"):
        infile = ROOT.TFile.Open(file)
        loaded_wsp = infile.Get("w")
        loaded_wsp.loadSnapshot("MultiDimFit")
        fit_vals = {}
        for POI in POIs:
            fit_vals[POI] = loaded_wsp.var(POI).getVal()
        return fit_vals

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

    def run_method(self):
        mass = self.args.mass
        dc = self.args.datacard
        self.put_back_arg("mass", "-m")
        self.put_back_arg("datacard", "-d")

        ######################################################################
        # Step 1 - establish parameter ranges
        ######################################################################
        with open(self.args.config) as jsonfile:
            cfg = json.load(jsonfile)
        if self.args.choose_POIs is None:
            POIs = sorted([str(x) for x in cfg.keys()])
        else:
            POIs = self.args.choose_POIs.split(",")

        Nx = len(POIs)
        print(">> Taylor expansion in %i variables up to order %i:" % (Nx, self.args.order))
        pprint(cfg)

        if self.args.workspace_bestfit:
            fitvals = self.get_snpashot_pois(dc, POIs)
            for POI, val in six.iteritems(fitvals):
                print(">> Updating POI best fit from %f to %f" % (cfg[POI]["BestFit"], val))
                cfg[POI]["BestFit"] = val

        xvec = np.zeros(Nx, dtype=np.float32)
        # hvec = []
        valvec = []
        for i, P in enumerate(POIs):
            valvec.append(cfg[P]["BestFit"])
            xvec[i] = cfg[P]["BestFit"]

        ######################################################################
        # Step 2 - generate stencils
        ######################################################################

        linear_POIs = []
        if self.args.linear_workaround is not None:
            linear_POIs = self.args.linear_workaround.split(",")

        stencils = {}
        validity = []

        do_cheb = False

        for i, P in enumerate(POIs):
            stencils[i] = {}
            if "StencilRange" in cfg[P]:
                s_min = cfg[P]["StencilRange"][0] - valvec[i]
                s_max = cfg[P]["StencilRange"][1] - valvec[i]
            elif "StencilSize" in cfg[P]:
                s_min = -cfg[P]["StencilSize"]
                s_max = +cfg[P]["StencilSize"]
            s_min *= 1.0
            s_max *= 1.0
            validity.append(cfg[P]["Validity"])
            for n in range(self.args.order + 1):
                if n == 0:
                    continue
                stencil_size = max(self.args.stencil_min, 1 + (((n + 1) / 2) * 2) + self.args.stencil_add)
                stencil = list()
                stencil_spacing = (s_max - s_min) / (stencil_size - 1)
                for s in range(stencil_size):
                    stencil.append(s_min + float(s) * stencil_spacing)

                if do_cheb:
                    cheb_list = []
                    a = stencil[0]
                    b = stencil[-1]
                    chebN = len(stencil)
                    for inode in range(1, chebN + 1):
                        cheb_list.append(0.5 * (a + b) + 0.5 * (b - a) * math.cos((((2.0 * inode) - 1.0) * math.pi) / (2.0 * chebN)))
                    cheb_list.sort()
                    stencil = cheb_list
                coefficients = GenerateStencils(n, 1, stencil)

                ## special case here for linear
                if n == 2 and P in linear_POIs:
                    ## First requirement is that s_min or s_max  == the best-fit
                    if abs(s_min) < 1e-6:
                        xprime = s_max
                        print(xprime)
                        stencil = [s_min, s_max]
                        coefficients = [0.0, 2.0 / (xprime * xprime)]
                    elif abs(s_max) < 1e-6:
                        xprime = s_min
                        stencil = [s_min, s_max]
                        coefficients = [2.0 / (xprime * xprime), 0.0]
                    else:
                        print("Special treatment of %s not viable, one stencil range endpoint must correspond to the best fit" % P)

                stencils[i][n] = list(zip(stencil, coefficients))

        pprint(stencils)

        cached_terms = {}
        cached_evals = {}

        can_skip = []

        drop_thresh = self.args.drop_threshold

        stats = {}

        graph_func = None
        test_mode_args = []
        test_mode_ws = ROOT.RooWorkspace()
        ## In this mode we are loading a TGraph from a scan
        if self.args.test_mode == 2:
            test_mode_args = self.args.test_args.split(",")
            graph_filename = test_mode_args[0]
            graph_name = test_mode_args[1]
            graph_file = ROOT.TFile(graph_filename)
            graph = graph_file.Get(graph_name)
            spline = ROOT.TSpline3("spline3", graph)
            graph_func = ROOT.TF1("splinefn", partial(Eval, spline), graph.GetX()[0], graph.GetX()[graph.GetN() - 1], 1)

        if self.args.test_mode == 3:
            test_mode_args = [self.args.test_args]
            print(test_mode_args)
            test_mode_ws.factory("expr::func(%s)" % test_mode_args[0])

        ######################################################################
        # Step 3 - load pre-existing terms and evals
        ######################################################################
        if self.args.load is not None:
            term_cachefile = self.args.load + "_terms.pkl"
            eval_cachefile = self.args.load + "_evals.pkl"
            if os.path.isfile(term_cachefile):
                with open(term_cachefile) as pkl_file:
                    cached_terms = pickle.load(pkl_file)
            if os.path.isfile(eval_cachefile):
                with open(eval_cachefile) as pkl_file:
                    cached_evals = pickle.load(pkl_file)

        for i in range(self.args.order + 1):
            print(">> Order %i" % i)
            if i == 0 or i == 1:
                continue

            evallist = []
            termlist = []
            to_check_list = []

            stats[i] = {}
            stats[i]["nTerms"] = 0
            stats[i]["nCachedTerms"] = 0
            stats[i]["nSmallTerms"] = 0
            stats[i]["nAllNewTerms"] = 0
            stats[i]["nActualNewTerms"] = 0
            stats[i]["nEvals"] = 0
            stats[i]["nUniqueEvals"] = 0
            for item in itertools.combinations_with_replacement(list(range(Nx)), i):
                if len(set(item)) != 1 and i > self.args.cross_order:
                    if item in cached_terms:
                        del cached_terms[item]
                    continue

                stats[i]["nTerms"] += 1

                # If already in the cache we can skip evaluating this term, but first check
                # if it's small enough to be used added to the list of 'can_skip' terms
                if item in cached_terms:
                    stats[i]["nCachedTerms"] += 1
                    to_check_list.append((item, cached_terms[item]))
                    continue

                stats[i]["nAllNewTerms"] += 1

                # Check if this new term can be skipped
                skip_term = False
                for skip_item in can_skip:
                    has_all_terms = True
                    for x, freq in six.iteritems(skip_item[1]):
                        if item.count(x) < freq:
                            has_all_terms = False
                            break
                    if has_all_terms:
                        # print 'Testing if term %s is skipable based on %s' % (str(item), str(skip_item))
                        perm_ratio = float(Permutations(item)) / float(Permutations(skip_item[0]))
                        fact_ratio = float(math.factorial(len(skip_item[0]))) / float(math.factorial(len(item)))
                        expected = cached_terms[skip_item[0]] * perm_ratio * fact_ratio
                        for index in item:
                            expected *= validity[index]
                        # print 'Original = %g, permutations ratio = %g, factorial ratio = %g, final = %g' % (
                        #     skip_item[2],
                        #     perm_ratio,
                        #     fact_ratio,
                        #     expected)
                        # print 'Original estimate was %g, scaling with missing %s = %g' % (skip_item[2], str(remainder_terms), expected)
                        if abs(expected) < drop_thresh:
                            skip_term = True
                            break
                if skip_term:
                    # print 'Skipping negligible: %s' % str(item)
                    continue

                stats[i]["nActualNewTerms"] += 1

                termlist.append(ExpansionTerm(xvec, item, stencils))
                termlist[-1].GatherFundamentalTerms(evallist)

            stats[i]["nEvals"] = len(evallist)
            unique_evallist = [x for x in set(x.FormattedPars() for x in evallist)]
            stats[i]["nUniqueEvals"] = len(unique_evallist)
            actual_evallist = [x for x in unique_evallist if x not in cached_evals]
            stats[i]["nActualUniqueEvals"] = len(actual_evallist)

            ## Prepare the inputs for the different test modes
            if len(actual_evallist) > 0 and self.args.test_mode == 1:
                self.load_workspace(dc, POIs)

            multicount = 0
            multivars = []
            for idx, vals in enumerate(actual_evallist):
                if self.args.multiple == 1:
                    set_vals = []
                    for POI, val in zip(POIs, vals):
                        set_vals.append("%s=%f" % (POI, val))
                        if self.args.test_mode == 1:
                            self.wsp_vars[POI].setVal(val)
                    set_vals_str = ",".join(set_vals)
                else:
                    multivars.append(vals)
                    if multicount == 0:
                        multiargs = []
                        for POI in POIs:
                            multiargs.append("%s=" % POI)
                    for ival, val in enumerate(vals):
                        multiargs[ival] += "%f" % val
                        if multicount < (self.args.multiple - 1) and idx < (len(actual_evallist) - 1):
                            multiargs[ival] += ","
                    if multicount == (self.args.multiple - 1) or idx == (len(actual_evallist) - 1):
                        set_vals_str = ":".join(multiargs)
                        multicount = 0
                    else:
                        multicount += 1
                        continue

                hash_id = hashlib.sha1(set_vals_str).hexdigest()
                filename = "higgsCombine.TaylorExpand.%s.MultiDimFit.mH%s.root" % (hash_id, mass)
                arg_str = "-M MultiDimFit -n .TaylorExpand.%s --algo fixed --redefineSignalPOIs %s --fixedPointPOIs " % (hash_id, ",".join(POIs))
                arg_str += set_vals_str

                if self.args.do_fits:
                    if self.args.test_mode == 0 and not os.path.isfile(filename):
                        self.job_queue.append("combine %s %s" % (arg_str, " ".join(self.passthru)))
                    if self.args.test_mode == 1:
                        if idx % 10000 == 0:
                            print("Done %i/%i NLL evaluations..." % (idx, len(actual_evallist)))
                        cached_evals[vals] = self.nll.getVal() - self.nll0
                    if self.args.test_mode == 2:
                        # Divide by 2 here because the graph is already 2*deltaNLL
                        cached_evals[vals] = graph_func.Eval(vals[0]) / 2.0
                    if self.args.test_mode == 3:
                        # Divide by 2 here because the graph is already 2*deltaNLL
                        test_mode_ws.var("x").setVal(vals[0])
                        cached_evals[vals] = test_mode_ws.function("func").getVal()
                else:
                    if self.args.test_mode == 0:
                        if self.args.multiple == 1:
                            cached_evals[vals] = self.get_results(filename)[0]
                        else:
                            results = self.get_results(filename)
                            for xidx, xvals in enumerate(multivars):
                                cached_evals[xvals] = results[xidx]
                            multivars = list()

            if self.args.do_fits and len(self.job_queue):
                njobs = len(self.job_queue)
                self.flush_queue()
                pprint(stats[i])
                print("Number of jobs = %i" % njobs)
                print("Raw number of evaluations: %i" % len(evallist))
                return

            for x in evallist:
                x.fnval = cached_evals[x.FormattedPars()]

            for term in termlist:
                item = tuple(term.derivatives)
                term.Print()
                cached_terms[item] = term.Eval(with_permutations=True, with_factorial=True)
                print("Raw term: %f" % term.Eval(with_permutations=False, with_factorial=False))
                to_check_list.append((item, cached_terms[item]))

            for item, estimate in to_check_list:
                for index in item:
                    estimate *= validity[index]
                if abs(estimate) < drop_thresh:
                    can_skip.append((item, {x: item.count(x) for x in set(item)}, estimate))
                    stats[i]["nSmallTerms"] += 1

            pprint(stats[i])

        if self.args.save is not None:
            term_cachefile = self.args.save + "_terms.pkl"
            eval_cachefile = self.args.save + "_evals.pkl"
            jsondata = pickle.dumps(cached_terms)
            with open(term_cachefile, "w") as out_file:
                out_file.write(jsondata)
            jsondata = pickle.dumps(cached_evals)
            with open(eval_cachefile, "w") as out_file:
                out_file.write(jsondata)

        # Build the taylor expansion object
        xvars = []
        x0vars = []
        xvec = ROOT.RooArgList()
        x0vec = ROOT.RooArgList()
        for i, POI in enumerate(POIs):
            xvars.append(ROOT.RooRealVar(POI, "", cfg[POI]["BestFit"], cfg[POI]["OutputRange"][0], cfg[POI]["OutputRange"][1]))
            x0vars.append(ROOT.RooRealVar(POI + "_In", "", cfg[POI]["BestFit"], -100, 100))
            x0vars[-1].setConstant(True)
            xvec.add(xvars[-1])
            x0vec.add(x0vars[-1])

        te_terms = ROOT.vector("double")()
        pos = 0
        te_tracker = ROOT.vector("std::vector<int>")()

        save_cov_matrix = False
        if save_cov_matrix:
            hessian = ROOT.TMatrixDSym(len(POIs))
            cov_matrix = ROOT.TMatrixDSym(len(POIs))
            cor_matrix = ROOT.TMatrixDSym(len(POIs))
        sorted_terms = []
        for i in range(self.args.order + 1):
            sorted_tmp = []
            for tracker, val in six.iteritems(cached_terms):
                if len(tracker) == i:
                    sorted_tmp.append((tracker, val))
                    if i == 2 and save_cov_matrix:
                        multi = 1.0
                        if tracker[0] == tracker[1]:
                            multi = 2.0
                        hessian[int(tracker[0])][int(tracker[1])] = multi * val
                        hessian[int(tracker[1])][int(tracker[0])] = multi * val
            sorted_terms.extend(sorted(sorted_tmp, key=lambda x: x[0]))
        if save_cov_matrix:
            # hessian.Print()
            cov_matrix = hessian.Clone()
            cov_matrix.Invert()
            # cov_matrix.Print()
            cor_matrix = cov_matrix.Clone()
            for i in range(len(POIs)):
                for j in range(len(POIs)):
                    print(i, j, cor_matrix[i][j], (cov_matrix[i][i]), (cov_matrix[j][j]))
                    cor_matrix[i][j] = cor_matrix[i][j] / (math.sqrt(cov_matrix[i][i]) * math.sqrt(cov_matrix[j][j]))
            # cor_matrix.Print()
            fout = ROOT.TFile("covariance.root", "RECREATE")
            fout.WriteTObject(cor_matrix, "cor")
            h_cor = self.fix_TH2(ROOT.TH2D(cor_matrix), POIs)
            fout.WriteTObject(h_cor, "h_cor")
            fout.WriteTObject(cov_matrix, "cov")
            h_cov = self.fix_TH2(ROOT.TH2D(cov_matrix), POIs)
            fout.WriteTObject(h_cov, "h_cov")
        for tracker, val in sorted_terms:
            # Check if this is a really big value
            # if abs(val) > 1E9:
            #     print '%i -- %s --> %s: %f <ERROR LARGE VALUE>' % (pos, tracker, tracker, val)
            #     val = 0.

            # Push back the terms
            te_terms.push_back(val)
            i_tracker = ROOT.vector("int")(len(tracker))
            for idx, tr in enumerate(tracker):
                i_tracker[idx] = tr
            te_tracker.push_back(i_tracker)

            # Print it
            print("%i -- %s --> %s: %f" % (pos, tracker, tracker, val))

        # print 'NEGLIGIBLE TERMS AT ORDER %i: %i/%i' % (n, n_below, n_terms)

        nllfn = ROOT.RooTaylorExpansion("nll", "", xvec, x0vec, te_tracker, te_terms)

        fout = ROOT.TFile("taylor_expansion.root", "RECREATE")
        # xvec.Print('v')
        # x0vec.Print('v')
        wsp = ROOT.RooWorkspace("w", "")
        getattr(wsp, "import")(nllfn, ROOT.RooCmdArg())
        wsp.Write()
        fout.Close()

        if self.args.test_mode == 2:
            testgraph = GenerateDebugGraph(wsp, POIs[0], graph)
            fgout = ROOT.TFile("test_mode_%s.root" % (POIs[0]), "RECREATE")
            testgraph.SetName("main")
            testgraph.Write()
            fgout.Close()

        print("%-10s %-20s %-20s %-20s %-20s %-20s" % ("Order", "nTerms", "nCachedTerms", "nSmallTerms", "nAllNewTerms", "nActualNewTerms"))
        print("-" * 94)
        for i in stats:
            print(
                "%-10i %-20i %-20i %-20i %-20i %-20i"
                % (i, stats[i]["nTerms"], stats[i]["nCachedTerms"], stats[i]["nSmallTerms"], stats[i]["nAllNewTerms"], stats[i]["nActualNewTerms"])
            )

        print("\n%-10s %-20s %-20s %-20s" % ("Order", "nEvals", "nUniqueEvals", "nActualUniqueEvals"))
        print("-" * 74)
        for i in stats:
            print("%-10i %-20i %-20i %-20i" % (i, stats[i]["nEvals"], stats[i]["nUniqueEvals"], stats[i]["nActualUniqueEvals"]))
