from __future__ import absolute_import
from __future__ import print_function
import ROOT
import json
import itertools
import glob
import sys
import re
import zipfile
import os
import bisect
from math import floor
from array import array

import HiggsAnalysis.CombinedLimit.tool_base.utils as utils
from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase
import HiggsAnalysis.CombinedLimit.util.plotting as plot
import six
from six.moves import range


class AsymptoticGrid(CombineToolBase):
    description = "Calculate asymptotic limits on parameter grids"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_intercept_args(self, group):
        CombineToolBase.attach_intercept_args(self, group)
        group.add_argument("--setParameters", default=None)
        group.add_argument("--freezeParameters", default=None)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("config", help="json config file")

    def run_method(self):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch(ROOT.kTRUE)

        # This is what the logic should be:
        #  - get the list of model points
        #  - figure out which files are:
        #    - completely missing
        #    - there but corrupt, missing tree
        #    - ok
        #  - If we have anything in the third category proceed to produce output files
        #  - Anything in the first two gets added to the queue only if --doFits is specified

        # Step 1 - open the json config file
        with open(self.args.config) as json_file:
            cfg = json.load(json_file)
        # to do - have to handle the case where it doesn't exist
        points = []
        blacklisted_points = []
        for igrid in cfg["grids"]:
            assert len(igrid) == 3
            if igrid[2] == "":
                points.extend(itertools.product(utils.split_vals(igrid[0]), utils.split_vals(igrid[1])))
            else:
                blacklisted_points.extend(itertools.product(utils.split_vals(igrid[0]), utils.split_vals(igrid[1]), utils.split_vals(igrid[2])))
        POIs = cfg["POIs"]
        opts = cfg["opts"]

        # remove problematic points (points with NaN values)
        points_to_remove = []
        grids_to_remove = cfg.get("grids_to_remove", None)
        if grids_to_remove is not None:
            for igrid in grids_to_remove:
                assert len(igrid) == 2
                points_to_remove.extend(itertools.product(utils.split_vals(igrid[0]), utils.split_vals(igrid[1])))
        for p in points_to_remove:
            points.remove(p)

        # Have to merge some arguments from both the command line and the "opts" in the json file
        to_freeze = []
        to_set = []
        set_opt, opts = self.extract_arg("--setParameters", opts)
        if set_opt is not None:
            to_set.append(set_opt)
        freeze_opt, opts = self.extract_arg("--freezeParameters", opts)
        if freeze_opt is not None:
            to_freeze.append(freeze_opt)
        if hasattr(self.args, "setParameters") and self.args.setParameters is not None:
            to_set.append(self.args.setParameters)
        if hasattr(self.args, "freezeParameters") and self.args.freezeParameters is not None:
            to_freeze.append(self.args.freezeParameters)

        file_dict = {}
        for p in points:
            file_dict[p] = []

        for f in glob.glob("higgsCombine.%s.*.%s.*.AsymptoticLimits.mH*.root" % (POIs[0], POIs[1])):
            # print f
            rgx = re.compile(r"higgsCombine\.%s\.(?P<p1>.*)\.%s\.(?P<p2>.*)\.AsymptoticLimits\.mH.*\.root" % (POIs[0], POIs[1]))
            matches = rgx.search(f)
            p = (matches.group("p1"), matches.group("p2"))
            if p in file_dict:
                file_dict[p].append(f)

        for key, val in six.iteritems(file_dict):
            name = "%s.%s.%s.%s" % (POIs[0], key[0], POIs[1], key[1])
            print(">> Point %s" % name)
            if len(val) == 0:
                print("Going to run limit for point %s" % (key,))
                set_arg = ",".join(["%s=%s,%s=%s" % (POIs[0], key[0], POIs[1], key[1])] + to_set)
                freeze_arg = ",".join(["%s,%s" % (POIs[0], POIs[1])] + to_freeze)
                point_args = "-n .%s --setParameters %s --freezeParameters %s" % (name, set_arg, freeze_arg)
                cmd = " ".join(["combine -M AsymptoticLimits", opts, point_args] + self.passthru)
                self.job_queue.append(cmd)

        bail_out = len(self.job_queue) > 0
        self.flush_queue()

        if bail_out:
            print(">> New jobs were created / run in this cycle, run the script again to collect the output")
            sys.exit(0)

        xvals = []
        yvals = []
        zvals_m2s = []
        zvals_m1s = []
        zvals_exp = []
        zvals_p1s = []
        zvals_p2s = []
        zvals_obs = []
        for key, val in six.iteritems(file_dict):
            for filename in val:
                fin = ROOT.TFile(filename)
                if fin.IsZombie():
                    continue
                tree = fin.Get("limit")
                for evt in tree:
                    if abs(evt.quantileExpected + 1) < 0.01:
                        xvals.append(float(key[0]))
                        yvals.append(float(key[1]))
                        # print 'At point %s have observed CLs = %f' % (key, evt.limit)
                        zvals_obs.append(float(evt.limit))
                    if abs(evt.quantileExpected - 0.025) < 0.01:
                        # print 'At point %s have -2sigma CLs = %f' % (key, evt.limit)
                        zvals_m2s.append(float(evt.limit))
                    if abs(evt.quantileExpected - 0.16) < 0.01:
                        # print 'At point %s have -1sigma CLs = %f' % (key, evt.limit)
                        zvals_m1s.append(float(evt.limit))
                    if abs(evt.quantileExpected - 0.5) < 0.01:
                        # print 'At point %s have expected CLs = %f' % (key, evt.limit)
                        zvals_exp.append(float(evt.limit))
                    if abs(evt.quantileExpected - 0.84) < 0.01:
                        # print 'At point %s have +1sigma CLs = %f' % (key, evt.limit)
                        zvals_p1s.append(float(evt.limit))
                    if abs(evt.quantileExpected - 0.975) < 0.01:
                        # print 'At point %s have +2sigma CLs = %f' % (key, evt.limit)
                        zvals_p2s.append(float(evt.limit))
        for POI1, POI2, CLs in blacklisted_points:
            xvals.append(float(POI1))
            yvals.append(float(POI2))
            zvals_m2s.append(float(CLs))
            zvals_m1s.append(float(CLs))
            zvals_exp.append(float(CLs))
            zvals_p1s.append(float(CLs))
            zvals_p2s.append(float(CLs))
            zvals_obs.append(float(CLs))
        graph_m2s = ROOT.TGraph2D(len(zvals_m2s), array("d", xvals), array("d", yvals), array("d", zvals_m2s))
        graph_m1s = ROOT.TGraph2D(len(zvals_m1s), array("d", xvals), array("d", yvals), array("d", zvals_m1s))
        graph_exp = ROOT.TGraph2D(len(zvals_exp), array("d", xvals), array("d", yvals), array("d", zvals_exp))
        graph_p1s = ROOT.TGraph2D(len(zvals_p1s), array("d", xvals), array("d", yvals), array("d", zvals_p1s))
        graph_p2s = ROOT.TGraph2D(len(zvals_p2s), array("d", xvals), array("d", yvals), array("d", zvals_p2s))
        graph_obs = ROOT.TGraph2D(len(zvals_obs), array("d", xvals), array("d", yvals), array("d", zvals_obs))
        # h_bins = cfg['hist_binning']
        # hist = ROOT.TH2F('h_observed', '', h_bins[0], h_bins[1], h_bins[2], h_bins[3], h_bins[4], h_bins[5])
        # for i in xrange(1, hist.GetNbinsX()+1):
        #  for j in xrange(1, hist.GetNbinsY()+1):
        #    hist.SetBinContent(i, j, graph.Interpolate(hist.GetXaxis().GetBinCenter(i), hist.GetYaxis().GetBinCenter(j)))
        fout = ROOT.TFile("asymptotic_grid.root", "RECREATE")
        fout.WriteTObject(graph_m2s, "exp-2")
        fout.WriteTObject(graph_m1s, "exp-1")
        fout.WriteTObject(graph_exp, "exp0")
        fout.WriteTObject(graph_p1s, "exp+1")
        fout.WriteTObject(graph_p2s, "exp+2")
        fout.WriteTObject(graph_obs, "obs")
        # fout.WriteTObject(hist)
        fout.Close()
        # Next step: open output files
        # Fill TGraph2D with CLs, CLs+b


class HybridNewGrid(CombineToolBase):
    description = "Calculate toy-based limits on parameter grids"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_intercept_args(self, group):
        CombineToolBase.attach_intercept_args(self, group)
        group.add_argument("--setParameters", default=None)
        group.add_argument("--freezeParameters", default=None)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("config", help="json config file")
        group.add_argument("--cycles", default=0, type=int, help="Number of job cycles to create per point")
        group.add_argument("--output", action="store_true", help="Write CLs grids into an output file")
        group.add_argument("--from-asymptotic", default=None, help="JSON file which will be used to create a limit grid automatically")

    def GetCombinedHypoTest(self, files):
        if len(files) == 0:
            return None
        results = []
        for file in files:
            found_res = False
            f = ROOT.TFile(file)
            ROOT.gDirectory.cd("toys")
            for key in ROOT.gDirectory.GetListOfKeys():
                if ROOT.gROOT.GetClass(key.GetClassName()).InheritsFrom(ROOT.RooStats.HypoTestResult.Class()):
                    results.append(ROOT.gDirectory.Get(key.GetName()))
                    found_res = True
            f.Close()
            if not found_res:
                print(">> Warning, did not find a HypoTestResult object in file %s" % file)
        if (len(results)) > 1:
            for r in results[1:]:
                results[0].Append(r)
        ntoys = min(results[0].GetNullDistribution().GetSize(), results[0].GetAltDistribution().GetSize())
        if ntoys == 0:
            print(">> Warning, HypoTestResult from file(s) %s does not contain any toy results, did something go wrong in your fits?" % "+".join(files))
        return results[0]

    def ValidateHypoTest(self, hyp_res, min_toys, max_toys, contours, signif, cl, output=False, verbose=False, precomputed=None, feldman_cousins=False):
        results = {}

        if hyp_res is None and precomputed is None:
            return (False, {"ntoys": 0})

        ntoys = None

        if hyp_res is not None:
            # We will take the number of toys thrown as the minimum of the number of b-only or s+b toys
            if feldman_cousins:
                # For F-C we expect only s+b toys
                ntoys = hyp_res.GetAltDistribution().GetSize()
            else:
                ntoys = min(hyp_res.GetNullDistribution().GetSize(), hyp_res.GetAltDistribution().GetSize())
            print(">>> Number of b toys %i" % (hyp_res.GetNullDistribution().GetSize()))
            print(">>> Number of s+b toys %i" % (hyp_res.GetAltDistribution().GetSize()))

        if precomputed is not None:
            ntoys = precomputed["ntoys"]

        results["ntoys"] = ntoys

        if verbose:
            print(">>> Toys completed: %i [min=%i, max=%i]" % (ntoys, min_toys, max_toys))

        # If we're not going to write the CLs grids out and we fail the ntoys criteria then we
        # don't need to waste time getting all the CLs values. Can just return the results dict as-is.
        # 1st test - have we thrown at least the minimum number of toys?
        if ntoys < min_toys and not output:
            return (False, results)
        # 2nd test - have we thrown the maximum (or more) toys?
        if ntoys >= max_toys and not output:
            return (True, results)

        if hyp_res is not None:
            # 3rd test - are we > X sigma away from the exclusion CLs? This must be true for all the
            # contours we're interested in
            btoys = sorted([x for x in hyp_res.GetNullDistribution().GetSamplingDistribution()])
            # save the real observed test stat, we'll change it in this
            # loop to get the expected but we'll want to restore it at the end
            q_obs = hyp_res.GetTestStatisticData()

        crossing = 1 - cl
        signif_results = {}

        if verbose:
            print(">>> CLs target is a significance of %.1f standard deviations from %.3f" % (signif, crossing))

        for contour in contours:
            # Start by assuming this contour passes, we'll set it to False if it fails
            signif_results[contour] = True

            # If this is an expected contour we will extract the quantile from the name
            if "exp" in contour:
                quantile = ROOT.Math.normal_cdf(float(contour.replace("exp", "")))
                if verbose:
                    print(">>> Checking the %s contour at quantile=%f" % (contour, quantile))
                if hyp_res is not None:
                    # Get the stat statistic value at this quantile by rounding to the nearest b-only toy
                    testStat = btoys[int(min(floor(quantile * len(btoys) + 0.5), len(btoys) - 1))]
                    hyp_res.SetTestStatisticData(testStat)
            elif contour == "obs":
                if verbose:
                    print(">>> Checking the %s contour" % contour)
            else:
                raise RuntimeError("Contour %s not recognised" % contour)

            if hyp_res is not None:
                # Currently assume we always want to use CLs, should provide option
                # for CLs+b at some point
                if not feldman_cousins:
                    CLs = hyp_res.CLs()
                    CLsErr = hyp_res.CLsError()
                else:
                    # For simplicity label CLs+b the same as CLs when using FC mode...
                    CLs = hyp_res.CLsplusb()
                    CLsErr = hyp_res.CLsplusbError()
                testStatObs = hyp_res.GetTestStatisticData()
            if precomputed is not None:
                CLs = precomputed[contour][0]
                CLsErr = precomputed[contour][1]
                testStatObs = precomputed[contour][3]
            if ntoys == 0:
                CLsErr = 0  # If there were no toys then ClsError() will return inf
            dist = 0.0
            if CLsErr == 0.0:
                if verbose:
                    print(">>>> CLs = %g +/- %g (infinite significance), will treat as passing" % (CLs, CLsErr))
                dist = 999.0
            else:
                dist = abs(CLs - crossing) / CLsErr
                if verbose:
                    print(">>>> CLs = %g +/- %g, reached %.1f sigma signifance" % (CLs, CLsErr, dist))
                if dist < signif:
                    signif_results[contour] = False
            results[contour] = (CLs, CLsErr, dist, testStatObs)
            if hyp_res is not None:
                # Set the observed test statistic back to the real data value
                hyp_res.SetTestStatisticData(q_obs)

        # Now do the full logic of the validation and return
        all_ok = ntoys >= min_toys  # OK if min toys passes
        for key, val in six.iteritems(signif_results):
            all_ok = all_ok and val  # still OK if all contour significances pass
        all_ok = all_ok or (ntoys >= max_toys)  # Always OK if we've reached the maximum
        results["ok"] = all_ok
        return (all_ok, results)

    def run_method(self):
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch(ROOT.kTRUE)

        # Open the json config file
        with open(self.args.config) as json_file:
            cfg = json.load(json_file)

        # Set all the parameter values locally using defaults if necessary
        grids = cfg["grids"]
        grids_to_remove = cfg.get("grids_to_remove", None)
        POIs = cfg["POIs"]
        opts = cfg["opts"]
        toys_per_cycle = cfg["toys_per_cycle"]
        zipname = cfg.get("zipfile", None)
        statfile = cfg.get("statusfile", None)
        contours = cfg.get("contours", ["obs", "exp-2", "exp-1", "exp0", "exp+1", "exp+2"])
        min_toys = cfg.get("min_toys", 500)
        max_toys = cfg.get("max_toys", 5000)
        signif = cfg.get("signif", 3.0)
        cl = cfg.get("CL", 0.95)
        verbose = cfg.get("verbose", False)
        make_plots = cfg.get("make_plots", False)
        # Write CLs values into the output even if current toys do not pass validation
        incomplete = cfg.get("output_incomplete", False)
        outfile = cfg.get("output", "hybrid_grid.root")
        from_asymptotic_settings = cfg.get("from_asymptotic_settings", dict())
        feldman_cousins = cfg.get("FC", False)
        # NB: blacklisting not yet implemented for this method

        # Have to merge some arguments from both the command line and the "opts" in the json file
        to_freeze = []
        to_set = []
        set_opt, opts = self.extract_arg("--setParameters", opts)
        if set_opt is not None:
            to_set.append(set_opt)
        freeze_opt, opts = self.extract_arg("--freezeParameters", opts)
        if freeze_opt is not None:
            to_freeze.append(freeze_opt)
        if hasattr(self.args, "setParameters") and self.args.setParameters is not None:
            to_set.append(self.args.setParameters)
        if hasattr(self.args, "freezeParameters") and self.args.freezeParameters is not None:
            to_freeze.append(self.args.freezeParameters)

        points = []
        blacklisted_points = []

        # For the automatic grid for the "from_asymptotic option" we should fix the format specifier for
        # the grid points, as the numerical precision of a given point may change once the grid spacing is
        # modified. By default we let split_vals do it's thing however
        fmt_spec = None

        # In this mode we're doing a classic limit search vs MH instead of a 2D grid.
        # Most of the same code can be used however. First we'll use the json file containing the
        # asymptotic limits to create a new grid from scratch.
        if self.args.from_asymptotic is not None:
            grids = []
            bound_vals = None
            bound_pars = []
            fmt_spec = "%.5g"
            with open(self.args.from_asymptotic) as limit_json:
                limits = json.load(limit_json)
            for m in limits.keys():
                limit_vals = [x for x in limits[m].values()]
                max_limit = max(limit_vals)
                min_limit = min(limit_vals)
                # print (min_limit, max_limit)
                width = max_limit - min_limit
                max_limit += width * 0.3
                min_limit = max(0.0, min_limit - width * 0.3)
                nsteps = from_asymptotic_settings.get("points", 100)
                step_width = (max_limit - min_limit) / nsteps
                grids.append([m, "%g:%g|%g" % (min_limit, max_limit, step_width), ""])
                boundlist_file = from_asymptotic_settings.get("boundlist", "")
                if boundlist_file:
                    with open(boundlist_file) as json_file:
                        bnd = json.load(json_file)
                    bound_pars = list(bnd.keys())
                    print("Found bounds for parameters %s" % ",".join(bound_pars))
                    bound_vals = {}
                    for par in bound_pars:
                        bound_vals[par] = list()
                        for mass, bounds in six.iteritems(bnd[par]):
                            bound_vals[par].append((float(mass), bounds[0], bounds[1]))
                        bound_vals[par].sort(key=lambda x: x[0])
                # print (min_limit, max_limit)
            # sys.exit(0)

        for igrid in grids:
            assert len(igrid) == 3
            if igrid[2] == "":
                points.extend(itertools.product(utils.split_vals(igrid[0], fmt_spec=fmt_spec), utils.split_vals(igrid[1], fmt_spec=fmt_spec)))
            else:
                blacklisted_points.extend(itertools.product(utils.split_vals(igrid[0]), utils.split_vals(igrid[1]), utils.split_vals(igrid[2])))

        # In between cycles of toys we may find there's something wrong with some of the points in the grid and therefore want to remove them:
        points_to_remove = []
        if grids_to_remove is not None:
            for igrid in grids_to_remove:
                assert len(igrid) == 2
                points_to_remove.extend(itertools.product(utils.split_vals(igrid[0]), utils.split_vals(igrid[1])))

        for p in points_to_remove:
            points.remove(p)

        # This dictionary will keep track of the combine output files for each model point
        file_dict = {}
        for p in points:
            file_dict[p] = {}

        # The regex we will use to identify output files and extract POI values
        rgx = re.compile(r"higgsCombine\.%s\.(?P<p1>.*)\.%s\.(?P<p2>.*)\.HybridNew\.mH.*\.(?P<toy>.*)\.root" % (POIs[0], POIs[1]))

        stats = {}
        if statfile and os.path.isfile(statfile):
            with open(statfile) as stat_json:
                stats = json.load(stat_json)

        # Can optionally copy output root files into a zip archive
        # If the user has specified a zipfile we will first
        # look for output files in this archive before scanning the
        # current directory
        if zipname:
            # Open the zip file in append mode, this should also
            # create it if it doesn't exist
            zipf = zipfile.ZipFile(zipname, "a")
            for f in zipf.namelist():
                matches = rgx.search(f)
                p = (matches.group("p1"), matches.group("p2"))
                seed = int(matches.group("toy"))
                if p in file_dict:
                    if seed not in file_dict[p]:
                        # For each model point have a dictionary keyed on the seed number
                        # with a value pointing to the file in the archive in the format
                        # ROOT expects: "zipfile.zip#higgsCombine.blah.root"
                        file_dict[p][seed] = zipname + "#" + f

        # Now look for files in the local directory
        for f in glob.glob("higgsCombine.%s.*.%s.*.HybridNew.mH*.root" % (POIs[0], POIs[1])):
            matches = rgx.search(f)
            p = (matches.group("p1"), matches.group("p2"))
            seed = int(matches.group("toy"))
            if p in file_dict:
                # Don't add this file to the list if its seed number is already
                # a value in the dict.
                if seed not in file_dict[p]:
                    # If we're using the zipfile we'll add this now and
                    # then delete it from the local directory
                    # But: only in the file is good, we don't want to pollute the zip
                    # file with incomplete or failed jobs
                    if zipname and plot.TFileIsGood(f):
                        zipf.write(f)  # assume this throws if it fails
                        print("Adding %s to %s" % (f, zipname))
                        file_dict[p][seed] = zipname + "#" + f
                        os.remove(f)
                    else:  # otherwise just add the file to the dict in the normal way
                        file_dict[p][seed] = f

        if zipname:
            zipf.close()

        # These lists will keep track of the CLs values which we will use
        # to create the output TGraph2Ds
        output_x = []
        output_y = []
        output_data = {}
        output_ntoys = []
        output_clserr = {}
        output_signif = {}
        # One list of Z-values per contour
        for contour in contours:
            output_data[contour] = []
            output_clserr[contour] = []
            output_signif[contour] = []

        # Also keep track of the number of model points which have met the
        # CLs criteria
        total_points = 0
        complete_points = 0

        for key, val in six.iteritems(file_dict):
            status_changed = True
            total_points += 1
            status_key = ":".join(key)
            name = "%s.%s.%s.%s" % (POIs[0], key[0], POIs[1], key[1])

            # First check if we use the status json
            all_files = list(val.values())
            status_files = []
            files = []

            if status_key in stats:
                status_files = stats[status_key]["files"]
                if set(all_files) == set(status_files):
                    print("For point %s, no files have been updated" % name)
                    status_changed = False
                    files = all_files
                else:
                    files = [x for x in val.values() if plot.TFileIsGood(x)]
                    if set(files) == set(status_files) and len(files) < len(all_files):
                        print("For point %s, new files exist but they are not declared good" % name)
                        status_changed = False
            else:
                files = [x for x in val.values() if plot.TFileIsGood(x)]

            # Merge the HypoTestResult objects from each file into one
            res = None
            precomputed = None
            if status_key in stats and not status_changed and stats[status_key]["ntoys"] > 0:
                precomputed = stats[status_key]
            else:
                res = self.GetCombinedHypoTest(files)

            # Do the validation of this model point
            #
            ok, point_res = self.ValidateHypoTest(
                res,
                min_toys=min_toys,
                max_toys=max_toys,
                contours=contours,
                signif=signif,
                cl=cl,
                output=self.args.output,
                verbose=verbose,
                precomputed=precomputed,
                feldman_cousins=feldman_cousins,
            )

            print(">> Point %s [%i toys, %s]" % (name, point_res["ntoys"], "DONE" if ok else "INCOMPLETE"))

            stats[status_key] = {"files": files, "ntoys": point_res["ntoys"]}
            for cont in contours:
                if cont in point_res:
                    stats[status_key][cont] = point_res[cont]

            if ok:
                complete_points += 1

            # Make plots of the test statistic distributions if requested
            if res is not None and make_plots:
                self.PlotTestStat(res, "plot_" + name, opts=cfg["plot_settings"], poi_vals=(float(key[0]), float(key[1])), point_info=point_res)

            # Add the resulting CLs values to the output arrays. Normally just
            # for the model points that passed the validation criteria, but if "output_incomplete"
            # has been set to true then we'll write all model points where at least one HypoTestResult
            # is present
            if (res is not None or precomputed is not None) and (ok or incomplete) and self.args.output:
                output_x.append(float(key[0]))
                output_y.append(float(key[1]))
                output_ntoys.append(point_res["ntoys"])
                for contour in contours:
                    output_data[contour].append(point_res[contour][0])
                    output_clserr[contour].append(point_res[contour][1])
                    output_signif[contour].append(point_res[contour][2])

            # Do the job cycle generation if requested
            if not ok and self.args.cycles > 0:
                print(">>> Going to generate %i job(s) for point %s" % (self.args.cycles, key))
                # Figure out the next seed numbers we need to run by finding the maximum seed number
                # so far
                done_cycles = list(val.keys())
                new_idx = max(done_cycles) + 1 if len(done_cycles) > 0 else 1
                new_cycles = list(range(new_idx, new_idx + self.args.cycles))

                print(">>> Done cycles: " + ",".join(str(x) for x in done_cycles))
                print(">>> New cycles: " + ",".join(str(x) for x in new_cycles))

                # Build to combine command. Here we'll take responsibility for setting the name and the
                # model parameters, making sure the latter are frozen
                if not feldman_cousins:
                    set_arg = ",".join(["%s=%s,%s=%s" % (POIs[0], key[0], POIs[1], key[1])] + to_set)
                    freeze_arg = ",".join(["%s,%s" % (POIs[0], POIs[1])] + to_freeze)
                    point_args = "-n .%s --setParameters %s --freezeParameters %s" % (name, set_arg, freeze_arg)
                else:
                    single_point_arg = ".".join(["%s=%s,%s=%s" % (POIs[0], key[0], POIs[1], key[1])])
                    if len(to_set) > 0 and len(to_freeze) > 0:
                        point_args = "-n .%s --singlePoint %s --setParameters %s --freezeParameters %s" % (name, single_point_arg, to_set, to_freeze)
                    elif len(to_set) > 0:
                        point_args = "-n .%s --singlePoint %s --setParameters %s" % (name, single_point_arg, to_set)
                    elif len(to_freeze) > 0:
                        point_args = "-n .%s --singlePoint %s --freezeParameters %s" % (name, single_point_arg, to_freeze)
                    else:
                        point_args = "-n .%s --singlePoint %s " % (name, single_point_arg)

                if self.args.from_asymptotic:
                    mval = key[0]
                    command = []
                    for par in bound_pars:
                        # The (mass, None, None) is just a trick to make bisect_left do the comparison
                        # with the list of tuples in bound_var[par]. The +1E-5 is to avoid float rounding
                        # issues
                        lower_bound = bisect.bisect_left(bound_vals[par], (float(mval) + 1e-5, None, None))
                        # If lower_bound == 0 this means we are at or below the lowest mass point,
                        # in which case we should increase by one to take the bounds from this lowest
                        # point
                        if lower_bound == 0:
                            lower_bound += 1
                        command.append("%s=%g,%g" % (par, bound_vals[par][lower_bound - 1][1], bound_vals[par][lower_bound - 1][2]))
                    if len(command) > 0:
                        point_args += " --setParameterRanges %s" % (":".join(command))
                    # print per_mass_point_args
                    point_args += " --singlePoint %s" % key[1]
                    point_args += " -m %s" % mval
                # Build a command for each job cycle setting the number of toys and random seed and passing through any other
                # user options from the config file or the command line
                for idx in new_cycles:
                    cmd = " ".join(["combine -M HybridNew", opts, point_args, "-T %i" % toys_per_cycle, "-s %i" % idx] + self.passthru)
                    self.job_queue.append(cmd)

        print(">> %i/%i points have completed and require no further toys" % (complete_points, total_points))
        self.flush_queue()

        # Create and write output CLs TGraph2Ds here
        # TODO: add graphs with the CLs errors, the numbers of toys and whether or not the point passes
        if self.args.output and not self.args.from_asymptotic:
            fout = ROOT.TFile(outfile, "RECREATE")
            for c in contours:
                graph = ROOT.TGraph2D(len(output_data[c]), array("d", output_x), array("d", output_y), array("d", output_data[c]))
                graph.SetName(c)
                fout.WriteTObject(graph, c)
                # Also write a Graph with the CLsErr
                graph = ROOT.TGraph2D(len(output_clserr[c]), array("d", output_x), array("d", output_y), array("d", output_clserr[c]))
                graph.SetName("clsErr_" + c)
                fout.WriteTObject(graph, "clsErr_" + c)
                # And a Graph with the significance
                graph = ROOT.TGraph2D(len(output_signif[c]), array("d", output_x), array("d", output_y), array("d", output_signif[c]))
                graph.SetName("signif_" + c)
                fout.WriteTObject(graph, "signif_" + c)
            graph = ROOT.TGraph2D(len(output_ntoys), array("d", output_x), array("d", output_y), array("d", output_ntoys))
            graph.SetName("ntoys" + c)
            fout.WriteTObject(graph, "ntoys")
            fout.Close()

        if self.args.output and self.args.from_asymptotic:
            # Need to collect all the files for each mass point and hadd them:
            files_by_mass = {}
            for key, val in six.iteritems(file_dict):
                if key[0] not in files_by_mass:
                    files_by_mass[key[0]] = list()
                files_by_mass[key[0]].extend(list(val.values()))
            for m, files in six.iteritems(files_by_mass):
                gridfile = "higgsCombine.gridfile.%s.%s.%s.root" % (POIs[0], m, POIs[1])
                self.job_queue.append("hadd -f %s %s" % (gridfile, " ".join(files)))
                for exp in ["", "0.025", "0.160", "0.500", "0.840", "0.975"]:
                    self.job_queue.append(
                        " ".join(
                            [
                                "combine -M HybridNew --rAbsAcc 0",
                                opts,
                                "--grid %s" % gridfile,
                                "-n .final.%s.%s.%s" % (POIs[0], m, POIs[1]),
                                "-m %s" % (m),
                                ("--expectedFromGrid %s" % exp) if exp else "--noUpdateGrid",
                            ]
                            + self.passthru
                        )
                    )
                self.flush_queue()

        if statfile:
            with open(statfile, "w") as stat_out:
                stat_json = json.dumps(stats, sort_keys=True, indent=2, separators=(",", ": "))
                stat_out.write(stat_json)

    def PlotTestStat(self, result, name, opts, poi_vals, point_info=None):
        sign = -1.0
        if opts["one_sided"]:
            sign = 1.0
        null_vals = [x * sign * 2.0 for x in result.GetNullDistribution().GetSamplingDistribution()]
        alt_vals = [x * sign * 2.0 for x in result.GetAltDistribution().GetSamplingDistribution()]
        if len(null_vals) == 0 or len(alt_vals) == 0:
            print(">> Errror in PlotTestStat for %s, null and/or alt distributions are empty")
            return
        plot.ModTDRStyle()
        canv = ROOT.TCanvas(name, name)
        pad = plot.OnePad()[0]
        min_val = min(min(alt_vals), min(null_vals))
        max_val = max(max(alt_vals), max(null_vals))
        min_plot_range = min_val - 0.05 * (max_val - min_val)
        if opts["one_sided"]:
            min_plot_range = 0.0
            pad.SetLogy(True)
        max_plot_range = max_val + 0.05 * (max_val - min_val)
        hist_null = ROOT.TH1F("null", "null", 40, min_plot_range, max_plot_range)
        hist_alt = ROOT.TH1F("alt", "alt", 40, min_plot_range, max_plot_range)
        for val in null_vals:
            hist_null.Fill(val)
        for val in alt_vals:
            hist_alt.Fill(val)
        hist_alt.SetLineColor(ROOT.TColor.GetColor(4, 4, 255))
        hist_alt.SetFillColor(plot.CreateTransparentColor(ROOT.TColor.GetColor(4, 4, 255), 0.4))
        hist_alt.GetXaxis().SetTitle("-2 #times ln(^{}L_{%s}/^{}L_{%s})" % (opts["alt_label"], opts["null_label"]))
        hist_alt.GetYaxis().SetTitle("Pseudo-experiments")
        hist_alt.Draw()
        hist_null.SetLineColor(ROOT.TColor.GetColor(252, 86, 11))
        hist_null.SetFillColor(plot.CreateTransparentColor(ROOT.TColor.GetColor(254, 195, 40), 0.4))
        hist_null.Draw("SAME")
        val_obs = result.GetTestStatisticData() * sign * 2.0
        obs = ROOT.TArrow(val_obs, 0, val_obs, hist_alt.GetMaximum() * 0.3, 0.05, "<-|")
        obs.SetLineColor(ROOT.kRed)
        obs.SetLineWidth(3)
        obs.Draw()
        # exp_line = ROOT.TLine()
        # plot.Set(exp_line, LineStyle=2, LineColor=ROOT.kRed, LineWidth=1)
        # if point_info is not None:
        #     for exp in ['exp-2', 'exp-1', 'exp0', 'exp+1', 'exp+2']:
        #         if exp in point_info:
        #             exp_line.DrawLine(2*sign*point_info[exp][3], 0, 2*sign*point_info[exp][3], hist_alt.GetMaximum() * 0.3)
        plot.FixTopRange(pad, plot.GetPadYMax(pad), 0.25)
        leg = plot.PositionedLegend(0.22, 0.2, 3, 0.02)
        leg.AddEntry(hist_alt, opts["alt_label"], "F")
        leg.AddEntry(hist_null, opts["null_label"], "F")
        leg.AddEntry(obs, "Observed", "L")
        leg.Draw()
        plot.DrawCMSLogo(pad, "CMS", opts["cms_subtitle"], 0, 0.15, 0.035, 1.2)
        pt_l = ROOT.TPaveText(0.23, 0.75, 0.33, 0.9, "NDCNB")
        pt_l.AddText("Model:")
        pt_l.AddText("Toys:")
        pt_l.AddText("CLs+b:")
        pt_l.AddText("CLb:")
        pt_l.AddText("CLs:")
        # if point_info is not None:
        #     for exp in ['exp-2', 'exp-1', 'exp0', 'exp+1', 'exp+2']:
        #         pt_l.AddText(exp)
        plot.Set(pt_l, TextAlign=11, TextFont=62, BorderSize=0)
        pt_l.Draw()
        pt_r = ROOT.TPaveText(0.33, 0.75, 0.63, 0.9, "NDCNB")
        pt_r.AddText("%s [%s = %.1f, %s = %.1f]" % (opts["model_label"], opts["poi_labels"][0], poi_vals[0], opts["poi_labels"][1], poi_vals[1]))
        pt_r.AddText(
            "%i (%s) + %i (%s)" % (result.GetNullDistribution().GetSize(), opts["null_label"], result.GetAltDistribution().GetSize(), opts["alt_label"])
        )
        pt_r.AddText("%.3f #pm %.3f" % (result.CLsplusb(), result.CLsplusbError()))
        pt_r.AddText("%.3f #pm %.3f" % (result.CLb(), result.CLbError()))
        pt_r.AddText("%.3f #pm %.3f" % (result.CLs(), result.CLsError()))
        # if point_info is not None:
        #     for exp in ['exp-2', 'exp-1', 'exp0', 'exp+1', 'exp+2']:
        #         pt_r.AddText('%.3f #pm %.3f' % (point_info[exp][0], point_info[exp][1]))
        plot.Set(pt_r, TextAlign=11, TextFont=42, BorderSize=0)
        pt_r.Draw()
        pad.GetFrame().Draw()
        pad.RedrawAxis()
        for fmt in opts["formats"]:
            canv.SaveAs(fmt)
