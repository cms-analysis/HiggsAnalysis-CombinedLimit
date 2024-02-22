#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import ROOT
import json
import os
import pprint
from collections import defaultdict
from array import array

import HiggsAnalysis.CombinedLimit.tool_base.utils as utils
import HiggsAnalysis.CombinedLimit.util.plotting as plot

# from HiggsAnalysis.CombinedLimit.tool_base.opts import OPTS

from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase
import six


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


class PrintFit(CombineToolBase):
    description = "Print the output of MultimDitFit"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("input", help="The input file")
        group.add_argument("--algo", help="The algo used in MultiDimFit", default="none")
        group.add_argument("-P", "--POIs", help="The params that were scanned (in scan order)")
        group.add_argument("--json", help="Write json output (format file.json:key1:key2..")

    def run_method(self):
        if self.args.json is not None:
            json_structure = self.args.json.split(":")
            assert len(json_structure) >= 1
            if os.path.isfile(json_structure[0]):
                with open(json_structure[0]) as jsonfile:
                    js = json.load(jsonfile)
            else:
                js = {}
            js_target = js
            if len(json_structure) >= 2:
                for key in json_structure[1:]:
                    js_target[key] = {}
                    js_target = js_target[key]
        POIs = self.args.POIs.split(",")
        if self.args.algo == "none":
            res = utils.get_none_results(self.args.input, POIs)
            for p in POIs:
                val = res[p]
                print("%-30s = %+.3f" % (p, val))
            if self.args.json is not None:
                for key, val in six.iteritems(res):
                    js_target[key] = {"Val": val}
                with open(json_structure[0], "w") as outfile:
                    json.dump(js, outfile, sort_keys=True, indent=4, separators=(",", ": "))
        elif self.args.algo == "singles":
            res = utils.get_singles_results(self.args.input, POIs, POIs)
            for p in POIs:
                val = res[p][p]
                print("%s = %.3f -%.3f/+%.3f" % (p, val[1], val[1] - val[0], val[2] - val[1]))
        elif self.args.algo == "fixed":
            res = utils.get_fixed_results(self.args.input, POIs)
            print("%-30s   bestfit :   fixed" % (""))
            for p in POIs:
                print("%-30s = %+.3f  :   %+.3f" % (p, res["bestfit"][p], res["fixedpoint"][p]))
            print("-" * 60)
            print("2*deltaNLL = %f, nPOIs = %i, p-value = %0.4f" % (2.0 * res["deltaNLL"], len(POIs), res["pvalue"]))

            # pprint.pprint(res)


class CollectLimits(CombineToolBase):
    description = "Aggregate limit output from combine"
    requires_root = True
    default_name = "limits.json"

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("input", nargs="+", default=[], help="The input files")
        group.add_argument(
            "-o",
            "--output",
            nargs="?",
            const="limits.json",
            default="limits.json",
            help="""The name of the output json file.
            When the --use-dirs option is set the enclosing directory name
            will be appended to the filename given here.""",
        )
        group.add_argument(
            "--use-dirs",
            action="store_true",
            help="""Use the directory structure to create multiple limit
                 outputs and to set the output file names""",
        )
        group.add_argument("--toys", action="store_true", help="""Collect toy values""")
        group.add_argument("--limit-err", action="store_true", help="""Also store the uncertainties on the limit""")

    def run_method(self):
        limit_sets = defaultdict(list)
        for filename in self.args.input:
            if not plot.TFileIsGood(filename):
                print(">> File %s is corrupt or incomplete, skipping" % filename)
                continue
            if self.args.use_dirs is False:
                limit_sets["default"].append(filename)
            else:
                label = "default"
                dirs = filename.split("/")
                # The last dir could be the mass, if so we ignore it and check the next
                if len(dirs) > 1:
                    if not isfloat(dirs[-2]):
                        label = dirs[-2]
                    elif len(dirs) > 2:
                        label = dirs[-3]
                limit_sets[label].append(filename)
        # print limit_sets

        for label, filenames in six.iteritems(limit_sets):
            js_out = {}
            for filename in filenames:
                if plot.TFileIsGood(filename):
                    file = ROOT.TFile(filename)
                    tree = file.Get("limit")
                    for evt in tree:
                        mh = str(evt.mh)
                        if mh not in js_out:
                            js_out[mh] = {}
                            if self.args.toys:
                                js_out[mh]["toys"] = {}
                                for limit in ["obs", "exp0", "exp-2", "exp-1", "exp+1", "exp+2"]:
                                    js_out[mh]["toys"][limit] = []
                        if self.args.toys:
                            if evt.iToy > 0:
                                if evt.quantileExpected == -1:
                                    js_out[mh]["toys"]["obs"].append(evt.limit)
                                elif abs(evt.quantileExpected - 0.5) < 1e-4:
                                    js_out[mh]["toys"]["exp0"].append(evt.limit)
                                elif abs(evt.quantileExpected - 0.025) < 1e-4:
                                    js_out[mh]["toys"]["exp-2"].append(evt.limit)
                                elif abs(evt.quantileExpected - 0.160) < 1e-4:
                                    js_out[mh]["toys"]["exp-1"].append(evt.limit)
                                elif abs(evt.quantileExpected - 0.840) < 1e-4:
                                    js_out[mh]["toys"]["exp+1"].append(evt.limit)
                                elif abs(evt.quantileExpected - 0.975) < 1e-4:
                                    js_out[mh]["toys"]["exp+2"].append(evt.limit)
                            elif evt.iToy == 0:
                                if evt.quantileExpected == -1:
                                    js_out[mh]["obs"].append(evt.limit)

                        else:
                            if evt.quantileExpected == -1:
                                js_out[mh]["obs"] = evt.limit
                                if self.args.limit_err:
                                    js_out[mh]["obs_err"] = evt.limitErr
                            elif abs(evt.quantileExpected - 0.5) < 1e-4:
                                js_out[mh]["exp0"] = evt.limit
                                if self.args.limit_err:
                                    js_out[mh]["exp0_err"] = evt.limitErr
                            elif abs(evt.quantileExpected - 0.025) < 1e-4:
                                js_out[mh]["exp-2"] = evt.limit
                                if self.args.limit_err:
                                    js_out[mh]["exp-2_err"] = evt.limitErr
                            elif abs(evt.quantileExpected - 0.160) < 1e-4:
                                js_out[mh]["exp-1"] = evt.limit
                                if self.args.limit_err:
                                    js_out[mh]["exp-1_err"] = evt.limitErr
                            elif abs(evt.quantileExpected - 0.840) < 1e-4:
                                js_out[mh]["exp+1"] = evt.limit
                                if self.args.limit_err:
                                    js_out[mh]["exp+1_err"] = evt.limitErr
                            elif abs(evt.quantileExpected - 0.975) < 1e-4:
                                js_out[mh]["exp+2"] = evt.limit
                                if self.args.limit_err:
                                    js_out[mh]["exp+2_err"] = evt.limitErr

            if self.args.toys:
                for mh in js_out.keys():
                    print("Expected bands will be taken from toys")
                    print(mh)
                    limits = sorted(js_out[mh]["toys"]["obs"])
                    # if mh == '160.0' or mh == '90.0' :
                    #    limits = [x for x in limits if x > 0.1]
                    quantiles = array("d", [0.025, 0.160, 0.5, 0.840, 0.975])
                    res = array("d", [0.0, 0.0, 0.0, 0.0, 0.0])
                    empty = array("i", [0])
                    ROOT.TMath.Quantiles(len(limits), len(quantiles), array("d", limits), res, quantiles, True, empty, 1)
                    print(res)
                    js_out[mh]["exp-2"] = res[0]
                    js_out[mh]["exp-1"] = res[1]
                    js_out[mh]["exp0"] = res[2]
                    js_out[mh]["exp+1"] = res[3]
                    js_out[mh]["exp+2"] = res[4]
            # print js_out
            jsondata = json.dumps(js_out, sort_keys=True, indent=2, separators=(",", ": "))
            # print jsondata
            if self.args.output is not None:
                outname = self.args.output.replace(".json", "_%s.json" % label) if self.args.use_dirs else self.args.output
                with open(outname, "w") as out_file:
                    print(">> Writing output %s from files:" % outname)
                    pprint.pprint(filenames, indent=2)
                    out_file.write(jsondata)


class CollectGoodnessOfFit(CombineToolBase):
    description = "Aggregate Goodness of Fit output from fit and toys"
    requires_root = True
    default_name = "gof.json"

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("--input", nargs="+", default=[], help="The input files")
        group.add_argument(
            "-o",
            "--output",
            nargs="?",
            const="gof.json",
            default="gof.json",
            help="""The name of the output json file.
            When the --use-dirs option is set the enclosing directory name
            will be appended to the filename given here.""",
        )
        group.add_argument(
            "--use-dirs",
            action="store_true",
            help="""Use the directory structure to create multiple limit
                 outputs and to set the output file names""",
        )

    def run_method(self):
        limit_sets = defaultdict(list)
        for filename in self.args.input:
            if not plot.TFileIsGood(filename):
                print(">> File %s is corrupt or incomplete, skipping" % filename)
                continue
            if not self.args.use_dirs:
                if "default" not in limit_sets:
                    limit_sets["default"] = ([], [])
                limit_sets["default"][0].append(filename)
            else:
                label = "default"
                dirs = filename.split("/")
                # The last dir could be the mass, if so we ignore it and check the next
                if len(dirs) > 1:
                    if not isfloat(dirs[-2]):
                        label = dirs[-2]
                    elif len(dirs) > 2:
                        label = dirs[-3]
                if label not in limit_sets:
                    limit_sets[label] = ([], [])
                limit_sets[label][0].append(filename)

        for label, (filenames, toyfiles) in six.iteritems(limit_sets):
            js_out = {}
            for filename in filenames:
                file = ROOT.TFile(filename)
                tree = file.Get("limit")
                adding_cat_branch = False
                branches = []
                for branch in tree.GetListOfBranches():
                    # Current logic says any branch after quantileExpected is a special
                    # GOF branch labelled according to category
                    if adding_cat_branch:
                        branches.append(branch.GetName())
                    if branch.GetName() == "quantileExpected":
                        adding_cat_branch = True
                # print branches
                failedToys = 0
                nEvts = tree.GetEntries()
                for evt in tree:
                    mh = str(evt.mh)
                    if mh not in js_out:
                        js_out[mh] = {}
                    if evt.quantileExpected != -1:
                        continue
                    if evt.iToy > 0 and evt.limit < -0.5:  # Exclude toys with negative test statistic
                        failedToys += 1
                        continue
                    if branches:
                        for branch in branches:
                            if branch not in js_out[mh]:
                                js_out[mh][branch] = {}
                                js_out[mh][branch]["toy"] = []
                            if evt.iToy <= 0:
                                js_out[mh][branch]["obs"] = [getattr(evt, branch)]
                            else:
                                js_out[mh][branch]["toy"].append(getattr(evt, branch))
                    else:
                        if "toy" not in js_out[mh]:
                            js_out[mh]["toy"] = []
                        if evt.iToy <= 0:
                            js_out[mh]["obs"] = [evt.limit]
                        else:
                            js_out[mh]["toy"].append(evt.limit)
                if failedToys > 0:
                    print(
                        ">> %i/%i toys have negative test statistic values, and are excluded. This might indicate a failure in the calculation within combine, or for the KS and AD tests, an undefined value in toys with zero events. Note that the resulting p-value could be biased."
                        % (failedToys, nEvts)
                    )
            for mh in js_out:
                if all([entry in js_out[mh] for entry in ["toy", "obs"]]):
                    js_out[mh]["p"] = float(len([toy for toy in js_out[mh]["toy"] if toy >= js_out[mh]["obs"][0]])) / len(js_out[mh]["toy"])
                else:
                    for branch in js_out[mh]:
                        js_out[mh][branch]["p"] = float(len([toy for toy in js_out[mh][branch]["toy"] if toy >= js_out[mh][branch]["obs"][0]])) / len(
                            js_out[mh][branch]["toy"]
                        )

            # print js_out
            jsondata = json.dumps(js_out, sort_keys=True, indent=2, separators=(",", ": "))
            # print jsondata
            if self.args.output is not None:
                outname = self.args.output.replace(".json", "_%s.json" % label) if self.args.use_dirs else self.args.output
                with open(outname, "w") as out_file:
                    print(">> Writing output %s from files:" % outname)
                    pprint.pprint(filenames, indent=2)
                    out_file.write(jsondata)
