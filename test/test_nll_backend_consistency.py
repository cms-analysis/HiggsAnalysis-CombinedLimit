#!/usr/bin/env python3
"""Analytical NLL regressions and algorithmic consistency checks for various backends."""

import argparse
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

import ROOT


class Model:
    def __init__(self, layout="extra", observations=(0.0, 2.0), multipdf=False, signal=False):
        """
        A two-channel model (channels a and b) with a single global nuisance theta and a single POI r.
        The nuisance is constrained by a unit Gaussian, and the channels are simple Gaussians

        The layout argument controls how the constraint is stored in the model:
            - "extra": the constraint is added as an extra constraint to the RooSimultaneous
            - "embedded": the constraint is embedded in channel a's PDF

        The multipdf argument controls whether channel a uses a RooMultiPdf with two components:
            - The first component is a Gaussian with mean theta (or r + theta if signal=True)
            - The second component is a Gaussian with mean theta and a different width (alt_sigma)

        The signal argument controls whether the means of the Gaussians depend on the POI r:
            - If signal=True, channel a's mean is r + theta and channel b's mean is r - theta
            - If signal=False, both channels' means are just theta
        """
        self.x = ROOT.RooRealVar("x", "x", 0, -20, 20)
        self.theta = ROOT.RooRealVar("theta", "theta", 0, -7, 7)
        self.r = ROOT.RooRealVar("r", "r", 1, -5, 5)
        self.global_obs = ROOT.RooRealVar("theta_In", "theta_In", 0)
        self.global_obs.setConstant(True)
        self.sigma = ROOT.RooConstVar("sigma", "sigma", 1)
        self.constraint = ROOT.RooGaussian("theta_Pdf", "", self.theta, self.global_obs, self.sigma)
        self.mean_a = self.theta
        self.mean_b = self.theta
        if signal:
            # The means of the Gaussians depend on the POI r
            self.mean_a = ROOT.RooFormulaVar("mean_a", "@0+@1", ROOT.RooArgList(self.r, self.theta))
            self.mean_b = ROOT.RooFormulaVar("mean_b", "@0-@1", ROOT.RooArgList(self.r, self.theta))
        self.gauss_a = ROOT.RooGaussian("gauss_a", "", self.x, self.mean_a, self.sigma)
        self.gauss_b = ROOT.RooGaussian("gauss_b", "", self.x, self.mean_b, self.sigma)
        self.shape_a = self.gauss_a
        if multipdf:
            # Use a RooMultiPdf for channel a to toggle between two Gaussians
            self.index = ROOT.RooCategory("pdfindex", "")
            self.alt_sigma = ROOT.RooRealVar("alt_sigma", "", 1, 0.5, 2)
            self.alt = ROOT.RooGaussian("alternative", "", self.x, self.mean_a, self.alt_sigma)
            self.multi = ROOT.RooMultiPdf("envelope", "", self.index, ROOT.RooArgList(self.gauss_a, self.alt))
            self.index.setIndex(0)
            self.shape_a = self.multi
        self.yield_a = ROOT.RooConstVar("yield_a", "", 1)
        self.yield_b = ROOT.RooConstVar("yield_b", "", 1)
        self.a = ROOT.RooAddPdf("pdf_a", "", ROOT.RooArgList(self.shape_a), ROOT.RooArgList(self.yield_a))
        self.b = ROOT.RooAddPdf("pdf_b", "", ROOT.RooArgList(self.gauss_b), ROOT.RooArgList(self.yield_b))
        self.channel_a = self.a
        if layout == "embedded":
            # Embed the constraint in channel a's PDF
            self.prod = ROOT.RooProdPdf("pdf_a_constrained", "", ROOT.RooArgList(self.a, self.constraint))
            self.channel_a = self.prod
        self.channel = ROOT.RooCategory("channel", "")
        self.channel.defineType("a")
        self.channel.defineType("b")
        self.pdf = ROOT.RooSimultaneousOpt("model_s", "", self.channel)
        self.pdf.addPdf(self.channel_a, "a")
        self.pdf.addPdf(self.b, "b")
        if layout == "extra":
            self.pdf.addExtraConstraint(self.constraint)

        # Dummy 'anchor' variable to force Combine to keep the constraint-factorisation
        # path active when testing ignoreConstraint.
        self.anchor = ROOT.RooRealVar("anchor", "", 0)
        self.anchor.setConstant(True)
        self.anchor_constraint = ROOT.RooGaussian("anchor_Pdf", "", self.anchor, self.global_obs, self.sigma)
        self.pdf.addExtraConstraint(self.anchor_constraint)
        self.mask_a = ROOT.RooRealVar("mask_a", "", 0, 0, 1)
        self.mask_b = ROOT.RooRealVar("mask_b", "", 0, 0, 1)
        self.mask_a.setConstant(True)
        self.mask_b.setConstant(True)
        self.pdf.addChannelMasks(ROOT.RooArgList(self.mask_a, self.mask_b))
        self.observables = ROOT.RooArgSet(self.x, self.channel)
        self.data = ROOT.RooDataSet("data_obs", "", self.observables)
        for label, value in zip(("a", "b"), observations):
            self.channel.setLabel(label)
            self.x.setVal(value)
            self.data.add(self.observables)
        self.nuisances = ROOT.RooArgSet(self.theta)

    def _nll(self, backend, nuisances=None):
        """
        Create the NLL using a particular backend
        """
        previous_backend = str(ROOT.Combine.nllBackend())
        try:
            ROOT.Combine.setNllBackend(backend)
            return ROOT.CombineUtils.combineCreateNLL(self.pdf, self.data, self.nuisances if nuisances is None else nuisances, False)
        finally:
            ROOT.Combine.setNllBackend(previous_backend)

    def _delta(self, nll, value=1.0):
        """
        deltaNLL for a given value of theta, relative to a baseline at theta=0
        """
        self.theta.setVal(0)
        baseline = float(nll.getVal())
        self.theta.setVal(value)
        result = float(nll.getVal()) - baseline
        self.theta.setVal(0)
        return result

    def write_workspace(self, path):
        """
        Dump the model somewhere
        """
        workspace = ROOT.RooWorkspace("w", "w")
        getattr(workspace, "import")(self.pdf)
        getattr(workspace, "import")(self.data)
        mc = ROOT.RooStats.ModelConfig("ModelConfig", workspace)
        mc.SetPdf(workspace.pdf("model_s"))
        mc.SetParametersOfInterest(ROOT.RooArgSet(workspace.var("r")))
        mc.SetNuisanceParameters(ROOT.RooArgSet(workspace.var("theta")))
        mc.SetGlobalObservables(ROOT.RooArgSet(workspace.var("theta_In")))
        mc.SetObservables(ROOT.RooArgSet(workspace.var("x"), workspace.cat("channel")))
        getattr(workspace, "import")(mc)
        workspace.writeToFile(str(path))


class LikelihoodBoundaryTest(unittest.TestCase):
    def assert_close(self, actual, expected, message="", tolerance=1e-8):
        self.assertTrue(math.isfinite(actual), f"nonfinite value: {actual}")
        self.assertAlmostEqual(actual, expected, delta=tolerance, msg=message)

    def extra_constraints(self):
        """
        Add extra constraints to the model and verify that they enter the likelihood as expected.
        """
        model = Model(observations=(0, 0))
        nll = model._nll(ARGS.backend)
        for theta in (-3, -1, 1, 3):
            with self.subTest(theta=theta):
                self.assert_close(
                    model._delta(nll, theta), 3 * 0.5 * theta**2, "extra Gaussian penalty must enter the likelihood"  # 0.5 from the two channels + theta
                )

    def constraint_only_prefit(self):
        """
        Test the constraint-only likelihood used by FitDiagnostics for the prefit nuisance fit.
        """
        theta = ROOT.RooRealVar("theta", "", 0, -7, 7)
        mean = ROOT.RooConstVar("mean", "", 0)
        sigma = ROOT.RooConstVar("sigma", "", 1)
        constraint = ROOT.RooGaussian("constraint", "", theta, mean, sigma)
        category = ROOT.RooCategory("dummy", "")
        pdf = ROOT.RooSimultaneousOpt("prefit", "", category)
        pdf.addExtraConstraint(constraint)
        data = ROOT.RooDataSet("empty", "", ROOT.RooArgSet(category))
        nuisances = ROOT.RooArgSet(theta)
        previous_backend = str(ROOT.Combine.nllBackend())
        try:
            ROOT.Combine.setNllBackend(ARGS.backend)
            nll = ROOT.CombineUtils.combineCreateNLL(pdf, data, nuisances, True)
        finally:
            ROOT.Combine.setNllBackend(previous_backend)
        baseline = float(nll.getVal())
        theta.setVal(1)
        self.assert_close(
            float(nll.getVal()) - baseline, 0.5, "prefit constraint must survive generic NLL creation"  # No channels, so only the constraint contributes
        )

    def multipdf_penalty(self):
        """
        RooMultiPdf adds a penalty to the likelihood for more complex models.
        Test if that penalty is correctly included in the NLL
        """
        model = Model(layout="none", multipdf=True)
        nll = model._nll(ARGS.backend)
        model.index.setIndex(0)
        baseline = float(nll.getVal())
        first_penalty = float(model.multi.getCorrection())
        model.index.setIndex(1)

        # Check that the penalty matches the expected (one 'unit' for one extra parameter)
        penalty = float(model.multi.getCorrection()) - first_penalty
        self.assert_close(penalty, 0.5, "fixture must differ by one complexity parameter")
        with self.subTest(component="complexity penalty"):
            # Check that the penalty is included in the NLL
            self.assert_close(float(nll.getVal()) - baseline, penalty, "nested RooMultiPdf penalty must enter the likelihood")

        # Verify that the selected shape is evaluated too, not just the penalty
        model.alt_sigma.setVal(1.5)
        with self.subTest(component="shape and penalty"):
            self.assert_close(float(nll.getVal()) - baseline, math.log(1.5) + penalty)

    def ignore_constraint_poi_selection(self):
        """
        Check that the ignoreConstraint attribute is honored a nuisance is redefined as a POI
        """
        model = Model(layout="embedded")
        model.constraint.setAttribute("ignoreConstraint", True)
        selected_nuisances = ROOT.RooArgSet()
        nll = model._nll(ARGS.backend, selected_nuisances)
        # Ignoring the constraint, moving theta from 0 to 1 should change the NLL by -1
        self.assert_close(model._delta(nll), -1)

    def ignore_constraint_selected_nuisance(self):
        """
        Test that ignoreConstraint is *still* honoured when the nuisance is included in
        the Constrain set.
        """
        model = Model(layout="embedded")
        model.constraint.setAttribute("ignoreConstraint", True)
        nll = model._nll(ARGS.backend)
        self.assert_close(model._delta(nll), -1, "ignoreConstraint must exclude a selected nuisance constraint")

    def ignored_extra_constraint(self):
        """
        Test that a constraint stored in extraConstraints is excluded if ignoreConstraint is set
        """
        model = Model(observations=(0, 0))
        model.constraint.setAttribute("ignoreConstraint", True)
        nll = model._nll(ARGS.backend)
        self.assert_close(model._delta(nll), 1)

    def _run_algorithm(self, method, layout, directory, extra_args):
        """
        Helper to allow running Combine methods on a `Model()` instance
        """
        path = directory / "workspace.root"
        model = Model(layout=layout, signal=True)
        model.write_workspace(path)
        command = [
            ARGS.combine,
            "-M",
            method,
            str(path),
            "-n",
            ".boundary",
            "--noMCbonly",
            "1",
            "--nllbackend",
            ARGS.backend,
            "--cminDefaultMinimizerStrategy",
            "0",
            "--cminDefaultMinimizerTolerance",
            "0.001",
            "-v",
            "0",
        ] + extra_args
        result = subprocess.run(command, cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
        (directory / "combine.log").write_text(result.stdout)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def _saved_result(self, path, key, log):
        """
        Check that the Combine output file exists and has the expected items
        """
        self.assertTrue(path.exists(), f"missing fresh {path}\n{log}")
        file = ROOT.TFile.Open(str(path))
        self.assertTrue(bool(file) and not file.IsZombie(), log)
        result = file.Get(key)
        keys = [item.GetName() for item in file.GetListOfKeys()]
        self.assertTrue(bool(result), f"missing {key}; available objects: {keys}\n{log}")
        self.assertEqual(result.status(), 0, f"{key} status={result.status()}\n{log}")
        parameter = result.floatParsFinal().find("theta")
        self.assertTrue(bool(parameter), f"theta missing in {key}\n{log}")
        values = (float(parameter.getVal()), float(parameter.getError()))
        file.Close()
        return values

    def fitdiagnostics_prefit(self):
        """
        Build a Model with theta in extraConstraints and run FitDiagnostics
        Internally, FitDiagnostics switches to the ROOT default backend.

        Checks that the prefit nuisance fit is consistent with the expected unit Gaussian constraint.
        """
        with tempfile.TemporaryDirectory(prefix="fitdiagnostics-boundary-") as temp:
            directory = Path(temp)
            log = self._run_algorithm("FitDiagnostics", "extra", directory, ["--skipBOnlyFit", "--skipSBFit", "--out", str(directory)])
            value, error = self._saved_result(directory / "fitDiagnostics.boundary.root", "nuisances_prefit_res", log)
            self.assert_close(value, 0, tolerance=5e-3)
            self.assert_close(error, 1, "prefit must retain the unit Gaussian constraint", tolerance=5e-3)

    def channel_compatibility_constraints(self):
        """
        Test that ChannelCompatibilityCheck preserves nuisance constraints when it builds the alternative model.
        """
        for layout in ("embedded", "extra"):
            with self.subTest(layout=layout), tempfile.TemporaryDirectory(prefix="compatibility-boundary-") as temp:
                directory = Path(temp)
                log = self._run_algorithm("ChannelCompatibilityCheck", layout, directory, ["--saveFitResult", "--runMinos", "0"])
                path = directory / "higgsCombine.boundary.ChannelCompatibilityCheck.mH120.root"
                value, error = self._saved_result(path, "fit_alternate", log)
                self.assert_close(value, 0, tolerance=1e-2)
                self.assert_close(error, 1, "alternative model must retain nuisance constraints", tolerance=1e-2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["combine", "cpu", "codegen"], required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--combine", required=True)
    ARGS = parser.parse_args()
    if ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") < 0:
        raise RuntimeError("Cannot load the Combine library")
    ROOT.gInterpreter.AddIncludePath(str(Path(ARGS.combine).resolve().parents[1]))
    if not ROOT.gInterpreter.Declare("#include <HiggsAnalysis/CombinedLimit/interface/Combine.h>"):
        raise RuntimeError("Cannot declare the Combine interface")
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    if not hasattr(LikelihoodBoundaryTest, ARGS.case) or ARGS.case.startswith("_"):
        parser.error(f"unknown case: {ARGS.case}")
    suite = unittest.TestSuite([LikelihoodBoundaryTest(ARGS.case)])
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)
