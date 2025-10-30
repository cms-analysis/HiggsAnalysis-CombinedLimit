#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
#include <optional>
#include <string>

#include "TFile.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

#include "../../interface/Combine.h"
#include "../../interface/CascadeMinimizer.h"
#include "../../interface/RooMultiPdfCombine.h"
#include "../../interface/utils.h"
#include "../../interface/ProfilingTools.h"

#include <gtest/gtest.h>

namespace {

bool loadSnapshotIfExists(RooWorkspace &ws, const char *name) {
  if (!name)
    return false;
  if (!ws.getSnapshot(name))
    return false;
  return ws.loadSnapshot(name);
}

void configureCascadeMinimizerState(RooWorkspace &ws, RooStats::ModelConfig &mc) {
  // Adapted from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/src/Combine.cc#L1146-L1203
  auto &cfg = CascadeMinimizerGlobalConfigs::O();

  cfg.parametersOfInterest = RooArgList();
  if (const RooArgSet *poi = mc.GetParametersOfInterest()) {
    for (RooAbsArg *arg : *poi) {
      auto *var = dynamic_cast<RooRealVar *>(arg);
      if (!var)
        continue;
      cfg.parametersOfInterest.add(*var);
    }
  }

  cfg.nuisanceParameters = RooArgList();
  if (const RooArgSet *nuis = mc.GetNuisanceParameters()) {
    for (RooAbsArg *arg : *nuis) {
      auto *var = dynamic_cast<RooRealVar *>(arg);
      if (!var)
        continue;
      if (!var->isConstant())
        cfg.nuisanceParameters.add(*var);
    }
  }

  cfg.allFloatingParameters = RooArgList();
  RooArgSet allVars(ws.allVars());
  for (RooAbsArg *arg : allVars) {
    auto *var = dynamic_cast<RooRealVar *>(arg);
    if (!var)
      continue;
    if (!var->isConstant())
      cfg.allFloatingParameters.add(*var);
  }

  cfg.pdfCategories = RooArgList();
  cfg.allRooMultiPdfParams = RooArgList();
  cfg.allRooMultiPdfs = RooArgList();

  if (auto *discreteParameters = dynamic_cast<RooArgSet *>(ws.genobj("discreteParams"));
      discreteParameters != nullptr) {
    for (RooAbsArg *arg : *discreteParameters) {
      auto *cat = dynamic_cast<RooCategory *>(arg);
      if (!cat)
        continue;
      if (!cfg.pdfCategories.containsInstance(*cat))
        cfg.pdfCategories.add(*cat);
    }
  } else if (runtimedef::get("ADD_DISCRETE_FALLBACK")) {
    RooArgSet categories(ws.allCats());
    for (RooAbsArg *arg : categories) {
      auto *cat = dynamic_cast<RooCategory *>(arg);
      if (!cat)
        continue;
      if (std::string(cat->GetName()).find("pdfindex") == std::string::npos)
        continue;
      if (!cfg.pdfCategories.containsInstance(*cat))
        cfg.pdfCategories.add(*cat);
    }
  }

  if (cfg.pdfCategories.getSize() > 0) {
    // Mirrors https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/src/Combine.cc#L1205-L1230
    RooArgSet clients;
    utils::getClients(cfg.pdfCategories, ws.allPdfs(), clients);
    for (RooAbsArg *arg : clients) {
      auto *multi = dynamic_cast<RooMultiPdf *>(arg);
      if (!multi)
        continue;
      if (!cfg.allRooMultiPdfs.containsInstance(*multi))
        cfg.allRooMultiPdfs.add(*multi);

      std::unique_ptr<RooArgSet> pdfParams{multi->getParameters((const RooArgSet *)nullptr)};
      for (RooAbsArg *a : *pdfParams) {
        auto *var = dynamic_cast<RooRealVar *>(a);
        if (!var || var->isConstant())
          continue;
        if (!cfg.allRooMultiPdfParams.containsInstance(*var))
          cfg.allRooMultiPdfParams.add(*var);
      }
    }
  }
}
enum class BackendComparison { Cpu, Codegen };

struct FitSummary {
  std::map<std::string, double> poiValues;
  std::map<std::string, double> parameterErrors;
  std::unique_ptr<RooFitResult> fitResult;
  double minNll = std::numeric_limits<double>::quiet_NaN();
};

std::optional<FitSummary> runFit(RooWorkspace &workspace,
                                 RooStats::ModelConfig &modelConfig,
                                 RooAbsPdf &pdf,
                                 RooAbsData &data,
                                 const RooArgSet &poiSet,
                                 const std::string &backend,
                                 std::string *error = nullptr) {
  std:: cout << "Running fit with backend " << backend << std::endl;
  if (!loadSnapshotIfExists(workspace, "clean")) {
    for (RooAbsArg *arg : poiSet) {
      if (auto *var = dynamic_cast<RooRealVar *>(arg))
        var->setConstant(false);
    }
  }

  configureCascadeMinimizerState(workspace, modelConfig);

  const RooArgSet *nuisances = modelConfig.GetNuisanceParameters();
  std::unique_ptr<RooArgSet> constraintCopy;
  const RooArgSet *constraintPtr = nullptr;
  if (nuisances && nuisances->getSize() > 0) {
    constraintCopy = std::make_unique<RooArgSet>(*nuisances);
    constraintPtr = constraintCopy.get();
  }

  Combine::setNllBackend(backend);
  std::unique_ptr<RooAbsReal> nll;
  try {
    nll = combineCreateNLL(pdf, data, constraintPtr, /*offset=*/true);
  } catch (const std::exception &ex) {
    if (error) *error = ex.what();
    return std::nullopt;
  }
  if (!nll) {
    if (error) *error = "combineCreateNLL returned null pointer";
    return std::nullopt;
  }

  auto *primaryPoi = dynamic_cast<RooRealVar *>(poiSet.first());
  CascadeMinimizer minim(*nll, CascadeMinimizer::Constrained, primaryPoi);
  try {
    if (!minim.minimize(/*verbose=*/0)) {
      if (error) *error = "minimization failed";
      return std::nullopt;
    }
  } catch (const std::exception &ex) {
    if (error) *error = ex.what();
    return std::nullopt;
  }

  minim.hesse(/*verbose=*/0);
  auto fitResult = std::unique_ptr<RooFitResult>(minim.save());

  FitSummary summary;
  summary.minNll = nll->getVal();
  summary.fitResult = std::move(fitResult);

  for (RooAbsArg *arg : poiSet) {
    auto *poiVar = dynamic_cast<RooRealVar *>(arg);
    if (!poiVar)
      continue;
    summary.poiValues[poiVar->GetName()] = poiVar->getVal();
  }

  if (summary.fitResult) {
    const RooArgList &floats = summary.fitResult->floatParsFinal();
    for (int i = 0; i < floats.getSize(); ++i) {
      auto *var = dynamic_cast<RooRealVar *>(floats.at(i));
      if (!var)
        continue;
      const std::string name = var->GetName();
      const double err = var->getError();
      summary.parameterErrors[name] = err;
      const double val = var->getVal();
      std::cout << "  " << name << " = " << val;
      if (!std::isnan(err)) {
        std::cout << " +/-" << err;
      }
      std::cout << '\n';
    }
  }

  return summary;
}

void compareFits(const FitSummary &baseline,
                 const FitSummary &candidate,
                 const std::string &label,
                 double relTolerance) {
  for (const auto &entry : baseline.parameterErrors) {
    const std::string &name = entry.first;
    if (name.rfind("prop_", 0) != 0)
      continue;
    const auto it = candidate.parameterErrors.find(name);
    ASSERT_TRUE(it != candidate.parameterErrors.end())
        << "Parameter '" << name << "' missing in backend '" << label << "'";

    const double err1 = entry.second;
    const double err2 = it->second;
    const double diff = std::abs(err1 - err2);
    const double scale = std::max({1.0, std::abs(err1), std::abs(err2)});

    EXPECT_TRUE(diff / scale <= relTolerance)
        << "Parameter '" << name << "' error mismatch: combine=" << err1
        << ", " << label << "=" << err2;
  }
}

}  // namespace

struct TestConfig {
  BackendComparison backend;
  bool analyticBarlowBeeston;
  bool expectFailure;
};

class CreateNLLTest : public ::testing::TestWithParam<TestConfig> {
protected:
  void SetUp() override {
    oldAnalyticFlag_ = runtimedef::get("MINIMIZER_no_analytic");
    const bool enableAnalytic = GetParam().analyticBarlowBeeston;
    runtimedef::set("MINIMIZER_no_analytic", enableAnalytic ? 0 : 1);
  }

  void TearDown() override {
    runtimedef::set("MINIMIZER_no_analytic", oldAnalyticFlag_);
  }

  void runComparison() {
  CascadeMinimizer::initOptions();
  // Matches https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/bin/combine.cpp#L300-L317
  runtimedef::set("OPTIMIZE_BOUNDS", 1);
  runtimedef::set("ADDNLL_RECURSIVE", 1);
  runtimedef::set("ADDNLL_GAUSSNLL", 1);
  runtimedef::set("ADDNLL_HISTNLL", 1);
  runtimedef::set("ADDNLL_CBNLL", 1);
  runtimedef::set("TMCSO_AdaptivePseudoAsimov", 1);
  runtimedef::set("MINIMIZER_optimizeConst", 2);
  runtimedef::set("MINIMIZER_rooFitOffset", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_FACTOR", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_NONORM", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_BASICINT", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_KEEPZEROS", 1);
  runtimedef::set("ADDNLL_PRODNLL", 1);
  runtimedef::set("ADDNLL_HFNLL", 1);
  runtimedef::set("ADDNLL_HISTFUNCNLL", 1);
  runtimedef::set("ADDNLL_ROOREALSUM_CHEAPPROD", 1);

  const std::string fileName = "template-analysis_shape_autoMCStats.root";
  const std::string workspaceName = "w";
  const std::string modelConfigName = "ModelConfig";
  const std::string dataName = "data_obs";

  std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
  EXPECT_FALSE(!file || file->IsZombie()) << "ERROR: failed to open input file '" << fileName << "'.\n";

  auto *workspace = dynamic_cast<RooWorkspace *>(file->Get(workspaceName.c_str()));
  if (!workspace) {
    std::cerr << "ERROR: workspace '" << workspaceName << "' not found in '" << fileName << "'.\n";
    file->ls();
    //return 2;
  }

  // Matches https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/src/Combine.cc#L505-L603
  utils::check_inf_parameters(workspace->allVars(), /*verbosity=*/0);
  loadSnapshotIfExists(*workspace, "clean");

  auto *modelConfig =
      dynamic_cast<RooStats::ModelConfig *>(workspace->genobj(modelConfigName.c_str()));
  if (!modelConfig) {
    std::cerr << "ERROR: ModelConfig '" << modelConfigName << "' not found in workspace '"
              << workspaceName << "'.\n";
    //return 2;
  }

  RooAbsPdf *pdf = modelConfig->GetPdf();
  if (!pdf) {
    std::cerr << "ERROR: ModelConfig '" << modelConfigName << "' does not define a pdf.\n";
    //return 2;
  }

  RooAbsData *data = workspace->data(dataName.c_str());
  if (!data) {
    std::cerr << "ERROR: dataset '" << dataName << "' not found in workspace '" << workspaceName
              << "'.\n";
    //return 2;
  }

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CountErrors);

  const RooArgSet *poi = modelConfig->GetParametersOfInterest();
  if (!poi || poi->getSize() == 0) {
    std::cerr << "ERROR: ModelConfig '" << modelConfigName << "' has no parameters of interest.\n";
    //return 2;
  }
  for (RooAbsArg *arg : *poi) {
    auto *var = dynamic_cast<RooRealVar *>(arg);
    if (var)
      var->setConstant(false);
  }

  const TestConfig cfg = GetParam();

  std::string errMessage;
  auto combineFit = runFit(*workspace, *modelConfig, *pdf, *data, *poi, "combine", &errMessage);
  ASSERT_TRUE(combineFit) << "Failed to run fit with 'combine' backend: " << errMessage;

  errMessage.clear();
  const double relTolerance = 0.05; // 5%

  if (cfg.backend == BackendComparison::Cpu) {
    auto cpuFit = runFit(*workspace, *modelConfig, *pdf, *data, *poi, "cpu", &errMessage);
    ASSERT_TRUE(cpuFit) << "Failed to run fit with 'cpu' backend: " << errMessage;
    if (cfg.expectFailure) {
      ADD_FAILURE() << "Expected compare failure for cpu backend but fit succeeded";
    }
    compareFits(*combineFit, *cpuFit, "cpu", relTolerance);
  } else {
    auto codegenFit = runFit(*workspace, *modelConfig, *pdf, *data, *poi, "codegen", &errMessage);
    if (!codegenFit) {
      if (cfg.expectFailure) {
        GTEST_SKIP() << "Expected failure for codegen backend: " << errMessage;
      } else {
        ADD_FAILURE() << "Failed to run fit with 'codegen' backend: " << errMessage;
      }
      return;
    }
    if (cfg.expectFailure) {
      ADD_FAILURE() << "Expected compare failure for codegen backend but fit succeeded";
    }
    compareFits(*combineFit, *codegenFit, "codegen", relTolerance);
  }
  }

private:
  int oldAnalyticFlag_ = 0;
};

TEST_P(CreateNLLTest, BackendConsistency) {
  runComparison();
}

static std::string AnalyticParamName(const testing::TestParamInfo<TestConfig> &info) {
  const char *backend = info.param.backend == BackendComparison::Cpu ? "Cpu" : "Codegen";
  const char *analytic = info.param.analyticBarlowBeeston ? "AnalyticOn" : "AnalyticOff";
  const char *expect = info.param.expectFailure ? "ExpectFail" : "ExpectPass";
  return std::string(backend) + "_" + analytic + "_" + expect;
}

INSTANTIATE_TEST_SUITE_P(
  // Alternative evaluation backends are only supported from ROOT 6.30.
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 30, 0)
    DISABLED_BarlowBeestonMatrix,
#else
    BarlowBeestonMatrix,
#endif
    CreateNLLTest,
    ::testing::Values(TestConfig{BackendComparison::Cpu, true, false},
                      TestConfig{BackendComparison::Codegen, true, true},
                      TestConfig{BackendComparison::Cpu, false, false},
                      TestConfig{BackendComparison::Codegen, false, true}),
    AnalyticParamName);
