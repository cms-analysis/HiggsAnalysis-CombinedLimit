// Tests of the "codegen" evaluation backend on the template histograms with
// bin-by-bin uncertainties (CMSHistSum and CMSHistErrorPropagator over
// CMSHistFunc processes): the likelihood values against the "cpu" backend and
// the generated gradients against finite differences, with and without the
// analytic Barlow-Beeston minimization. With the minimization, the gradient
// has to be the one of the profiled likelihood also for bins whose closed-form
// solution sits at the boundary of zero expected yield, see
// CMSHistBarlowBeestonGradientTerm.
//
// The workspaces come from the CTests cmshistsum (CMSHistSum), cmshistfunc
// (CMSHistErrorPropagator) and cmshist_poisson_bins (both, with an autoMCStats
// threshold that gives Poisson and per-process Gaussian bin parameters).

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "TFile.h"
#include "TRandom3.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"

#include "../interface/CMSHistErrorPropagator.h"
#include "../interface/CMSHistHelpers.h"
#include "../interface/CMSHistSum.h"
#include "../interface/ChannelSumNLL.h"
#include "../interface/Combine.h"
#include "../interface/CombineUtils.h"

#include <gtest/gtest.h>

namespace {

  struct Model {
    std::unique_ptr<TFile> file;
    RooWorkspace *w = nullptr;
    RooAbsPdf *pdf = nullptr;
    RooAbsData *data = nullptr;
    RooArgSet nuisances;
    std::vector<RooAbsReal const *> hists;  // the CMSHistSum and CMSHistErrorPropagator nodes
  };

}  // namespace

static Model openModel(std::string const &fileName) {
  Model m;
  m.file.reset(TFile::Open(fileName.c_str(), /*option=*/"READ"));
  if (!m.file || m.file->IsZombie())
    throw std::runtime_error("could not open " + fileName);
  m.w = dynamic_cast<RooWorkspace *>(m.file->Get("w"));
  auto *mc = m.w ? dynamic_cast<RooStats::ModelConfig *>(m.w->genobj("ModelConfig")) : nullptr;
  if (!mc)
    throw std::runtime_error("no ModelConfig in " + fileName);
  m.pdf = mc->GetPdf();
  m.data = m.w->data("data_obs");
  if (!m.pdf || !m.data)
    throw std::runtime_error("no pdf or data in " + fileName);
  m.nuisances.add(*mc->GetNuisanceParameters());
  for (RooAbsArg *func : m.w->allFunctions()) {
    if (dynamic_cast<CMSHistSum *>(func) || dynamic_cast<CMSHistErrorPropagator *>(func))
      m.hists.push_back(dynamic_cast<RooAbsReal *>(func));
  }
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  m.w->saveSnapshot(/*name=*/"start", m.w->allVars());
  return m;
}

static std::unique_ptr<RooAbsReal> makeNLL(Model const &m, std::string const &backend) {
  Combine::setNllBackend(backend);
  return combineCreateNLL(*m.pdf, *m.data, &m.nuisances, /*offset=*/true);
}

// Moves the floating parameters to a random point: the signal strength to
// [0.5, 1.5], the nuisance parameters with a symmetric range of at least
// [-3, 3] (the morphing and bin-by-bin parameters) to [-1.6, 1.6], which
// reaches beyond the smooth region of the morphing, and the others by up to
// 10 %.
static void randomPoint(RooAbsReal const &nll, TRandom3 &rng) {
  RooArgSet params;
  nll.getParameters(/*observables=*/nullptr, params, /*stripDisconnected=*/false);
  for (RooAbsArg *arg : params) {
    auto *var = dynamic_cast<RooRealVar *>(arg);
    if (!var || var->isConstant())
      continue;
    if (std::string(var->GetName()) == "r") {
      var->setVal(rng.Uniform(0.5, 1.5));
    } else if (var->getMin() <= -3. && var->getMax() >= 3.) {
      var->setVal(rng.Uniform(-1.6, 1.6));
    } else {
      var->setVal(var->getVal() * rng.Uniform(0.9, 1.1));
    }
  }
}

// Compares the gradient with central finite differences of the value, for
// every floating parameter. With the analytic Barlow-Beeston minimization
// active, the value is the profiled likelihood.
static void checkGradient(RooAbsReal &nll, double relTolerance, std::string const &label) {
  ASSERT_TRUE(nll.hasGradient()) << label;
  RooArgSet params;
  nll.getParameters(/*observables=*/nullptr, params, /*stripDisconnected=*/false);
  std::vector<double> grad(params.size(), 0.);
  nll.gradient(grad.data());
  std::size_t i = 0;
  for (RooAbsArg *arg : params) {
    auto *var = dynamic_cast<RooRealVar *>(arg);
    const std::size_t index = i++;
    if (!var || var->isConstant())
      continue;
    const double v0 = var->getVal();
    const double h = 1e-4 * (1. + std::abs(v0));
    var->setVal(v0 + h);
    const double up = nll.getVal();
    var->setVal(v0 - h);
    const double down = nll.getVal();
    var->setVal(v0);
    const double fd = (up - down) / (2 * h);
    EXPECT_NEAR(grad[index], fd, relTolerance * (1. + std::abs(fd))) << label << ": parameter " << var->GetName();
  }
}

// Whether any analytically profiled bin parameter of the histograms sits at
// the boundary of zero expected yield.
static bool anyBoundaryBin(Model const &m) {
  std::vector<double> slopes;
  std::vector<double> ratios;
  for (RooAbsReal const *hist : m.hists) {
    cmsHistVisit(*hist, [&](auto const &h) { h.barlowBeestonGradientTerms(slopes, ratios); });
    for (double slope : slopes) {
      if (slope != 0.)
        return true;
    }
  }
  return false;
}

// A workspace and whether its analytic Barlow-Beeston minimization can hit
// the boundary of zero expected yield: only bins with one parameter for the
// total are profiled analytically, and with the autoMCStats threshold of the
// cmshist_poisson_bins workspaces the bins without data get parameters per
// process instead.
struct Workspace {
  const char *name;  // the test name
  const char *file;
  bool boundaryBins;
};

class CMSHistCodegen : public ::testing::TestWithParam<Workspace> {};

// The codegen likelihood agrees with the cpu one and its gradient with finite
// differences, with all bin parameters floating.
TEST_P(CMSHistCodegen, ValuesAndGradient) {
  Model m = openModel(GetParam().file);
  ASSERT_FALSE(m.hists.empty());
  auto nllCpu = makeNLL(m, "cpu");
  auto nllCodegen = makeNLL(m, "codegen");
  ASSERT_TRUE(nllCpu && nllCodegen);
  ASSERT_TRUE(dynamic_cast<cacheutils::ChannelSumNLL *>(nllCodegen.get()));
  TRandom3 rng(/*seed=*/1);
  for (int point = 0; point < 3; ++point) {
    m.w->loadSnapshot(/*name=*/"start");
    randomPoint(*nllCodegen, rng);
    const double cpu = nllCpu->getVal();
    const double codegen = nllCodegen->getVal();
    EXPECT_NEAR(codegen, cpu, 1e-9 * std::abs(cpu)) << "point " << point;
    checkGradient(*nllCodegen, /*relTolerance=*/1e-5, "point " + std::to_string(point));
  }
}

// With the analytic Barlow-Beeston minimization, the gradient is the one of
// the profiled likelihood, also where the closed-form solution of a bin
// without data sits at the boundary of zero expected yield. The random points
// must contain such bins, otherwise the test would not exercise that case.
TEST_P(CMSHistCodegen, ProfiledGradient) {
  Model m = openModel(GetParam().file);
  auto nll = makeNLL(m, "codegen");
  auto *sum = dynamic_cast<cacheutils::ChannelSumNLL *>(nll.get());
  ASSERT_TRUE(sum);
  TRandom3 rng(/*seed=*/1);
  int pointsWithBoundaryBin = 0;
  for (int point = 0; point < 3; ++point) {
    m.w->loadSnapshot(/*name=*/"start");
    sum->setAnalyticBarlowBeeston(/*flag=*/true);
    randomPoint(*nll, rng);
    nll->getVal();  // runs the closed-form minimization at this point
    if (anyBoundaryBin(m))
      ++pointsWithBoundaryBin;
    checkGradient(*nll, /*relTolerance=*/1e-4, "profiled, point " + std::to_string(point));
    sum->setAnalyticBarlowBeeston(/*flag=*/false);
  }
  EXPECT_EQ(pointsWithBoundaryBin > 0, GetParam().boundaryBins);
}

INSTANTIATE_TEST_SUITE_P(Workspaces,
                         CMSHistCodegen,
                         ::testing::Values(Workspace{"histsum", "cmshistsum.root", true},
                                           Workspace{"errprop", "ws_template-analysis.root", true},
                                           Workspace{"poisson_bins_histsum", "cmshist_poisson_bins_histsum.root", false},
                                           Workspace{"poisson_bins_errprop", "cmshist_poisson_bins_errprop.root", false}),
                         [](::testing::TestParamInfo<Workspace> const &info) { return std::string(info.param.name); });
