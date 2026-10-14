#include "../interface/CombineCodegenImpl.h"

#include <ROOT/RConfig.hxx>  // for ROOT_VERSION

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 36, 0)

#include "../interface/AsymPow.h"
#include "../interface/CMSHistBarlowBeestonGradientTerm.h"
#include "../interface/CMSHistErrorPropagator.h"
#include "../interface/CMSHistFunc.h"
#include "../interface/CMSHistSum.h"
#include "../interface/CombineMathFuncs.h"
#include "../interface/FastTemplate_Old.h"
#include "../interface/HGGRooPdfs.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/RooParametricHist.h"
#include "../interface/VerticalInterpHistPdf.h"
#include "../interface/VerticalInterpPdf.h"

#include <RooAbsRealLValue.h>
#include <RooArgList.h>
#include <RooUniformBinning.h>

#include <cctype>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void RooFit::Experimental::codegenImpl(AsymPow& arg, CodegenContext& ctx) {
  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::asymPow", arg.theta(), arg.kappaLow(), arg.kappaHigh()));
}

void RooFit::Experimental::codegenImpl(ProcessNormalization& arg, CodegenContext& ctx) {
  std::vector<double> logAsymmKappaLow;
  std::vector<double> logAsymmKappaHigh;
  logAsymmKappaLow.reserve(arg.logAsymmKappa().size());
  logAsymmKappaHigh.reserve(arg.logAsymmKappa().size());
  for (auto [lo, hi] : arg.logAsymmKappa()) {
    logAsymmKappaLow.push_back(lo);
    logAsymmKappaHigh.push_back(hi);
  }

  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::processNormalization",
                              arg.nominalValue(),
                              arg.thetaList().size(),
                              arg.asymmThetaList().size(),
                              arg.otherFactorList().size(),
                              arg.thetaList(),
                              arg.logKappa(),
                              arg.asymmThetaList(),
                              logAsymmKappaLow,
                              logAsymmKappaHigh,
                              arg.otherFactorList()));
}

void RooFit::Experimental::codegenImpl(FastVerticalInterpHistPdf2& arg, CodegenContext& ctx) {
  if (arg.smoothAlgo() < 0) {
    throw std::runtime_error("We only support smoothAlgo >= 0");
  }

  RooRealVar const& xVar = arg.x();

  int numBins = xVar.numBins();

  std::vector<double> nominalVec(numBins);
  std::vector<double> widthVec(numBins);
  std::vector<double> morphsVecSum;
  std::vector<double> morphsVecDiff;

  auto const& cacheNominal = arg.cacheNominal();

  for (int i = 0; i < numBins; ++i) {
    nominalVec[i] = cacheNominal.GetBinContent(i);
    widthVec[i] = cacheNominal.GetWidth(i);
  }

  std::size_t nCoefs = arg.coefList().size();

  morphsVecSum.reserve(numBins * nCoefs);
  morphsVecDiff.reserve(numBins * nCoefs);
  auto const& morphs = arg.morphs();
  for (unsigned int j = 0; j < nCoefs; ++j) {
    for (int i = 0; i < numBins; ++i) {
      morphsVecSum.push_back(morphs[j].sum[i]);
      morphsVecDiff.push_back(morphs[j].diff[i]);
    }
  }

  for (int i = 0; i < numBins; ++i) {
    nominalVec[i] = cacheNominal.GetBinContent(i);
  }

  // The bin index part - handle both uniform and non-uniform binning
  const RooAbsBinning& binning = xVar.getBinning();
  std::string binIdx;

  if (dynamic_cast<RooUniformBinning const*>(&binning)) {
    // UNIFORM BINNING - fast path using arithmetic
    double xLow = xVar.getMin();
    double xHigh = xVar.getMax();
    binIdx = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", xLow, xHigh, xVar, numBins, 1.);
  } else {
    // NON-UNIFORM BINNING - use binary search
    // Extract bin edges from binning object
    std::vector<double> binEdges(numBins + 1);
    for (int i = 0; i < numBins; ++i) {
      binEdges[i] = binning.binLow(i);
    }
    binEdges[numBins] = binning.binHigh(numBins - 1);

    // Pass vector to ctx
    binIdx = ctx.buildCall("RooFit::Detail::MathFuncs::rawBinNumber", xVar, binEdges, numBins + 1);
  }

  std::string arrName = ctx.getTmpVarName();

  std::stringstream code;
  code << "double " << arrName << "[" << numBins << "];\n";
  code << ctx.buildCall("RooFit::Detail::MathFuncs::fastVerticalInterpHistPdf2",
                        numBins,
                        nCoefs,
                        arg.coefList(),
                        nominalVec,
                        widthVec,
                        morphsVecSum,
                        morphsVecDiff,
                        arg.smoothRegion(),
                        arrName) +
              ";\n";

  ctx.addToCodeBody(code.str(), true);
  ctx.addResult(&arg, arrName + "[" + binIdx + "]");
}

void RooFit::Experimental::codegenImpl(FastVerticalInterpHistPdf2D2& arg, CodegenContext& ctx) {
  if (arg.smoothAlgo() < 0) {
    throw std::runtime_error("We only support smoothAlgo >= 0");
  }

  if (!arg.conditional()) {
    throw std::runtime_error("We only support conditional == true");
  }

  RooRealVar const& xVar = arg.x();
  RooRealVar const& yVar = arg.y();

  // Handle both uniform and non-uniform binning for X and Y
  const RooAbsBinning& binningX = xVar.getBinning();
  const RooAbsBinning& binningY = yVar.getBinning();

  auto const& cacheNominal = arg.cacheNominal();

  int numBinsX = cacheNominal.binX();
  int numBinsY = cacheNominal.binY();
  int numBins = numBinsY * numBinsY;

  std::vector<double> nominalVec(numBins);
  std::vector<double> widthVec(numBins);
  std::vector<double> morphsVecSum;
  std::vector<double> morphsVecDiff;

  for (int i = 0; i < numBins; ++i) {
    nominalVec[i] = cacheNominal.GetBinContent(i);
    widthVec[i] = cacheNominal.GetWidth(i);
  }

  std::size_t nCoefs = arg.coefList().size();

  morphsVecSum.reserve(numBins * nCoefs);
  morphsVecDiff.reserve(numBins * nCoefs);
  auto const& morphs = arg.morphs();
  for (unsigned int j = 0; j < nCoefs; ++j) {
    for (int i = 0; i < numBins; ++i) {
      morphsVecSum.push_back(morphs[j].sum[i]);
      morphsVecDiff.push_back(morphs[j].diff[i]);
    }
  }

  for (int i = 0; i < numBins; ++i) {
    nominalVec[i] = cacheNominal.GetBinContent(i);
  }

  // The bin index part - handle both uniform and non-uniform binning for X
  std::string binIdxX;
  if (dynamic_cast<RooUniformBinning const*>(&binningX)) {
    // UNIFORM BINNING for X - fast path
    double xLow = xVar.getMin();
    double xHigh = xVar.getMax();
    binIdxX = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", xLow, xHigh, arg.x(), numBinsX, 1.);
  } else {
    // NON-UNIFORM BINNING for X - use binary search
    std::vector<double> binEdgesX(numBinsX + 1);
    for (int i = 0; i < numBinsX; ++i) {
      binEdgesX[i] = binningX.binLow(i);
    }
    binEdgesX[numBinsX] = binningX.binHigh(numBinsX - 1);

    // Pass vector to ctx
    binIdxX = ctx.buildCall("RooFit::Detail::MathFuncs::rawBinNumber", arg.x(), binEdgesX, numBinsX + 1);
  }

  // Handle both uniform and non-uniform binning for Y
  std::string binIdxY;
  if (dynamic_cast<RooUniformBinning const*>(&binningY)) {
    // UNIFORM BINNING for Y - fast path
    double yLow = yVar.getMin();
    double yHigh = yVar.getMax();
    binIdxY = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", yLow, yHigh, arg.y(), numBinsY, 1.);
  } else {
    // NON-UNIFORM BINNING for Y - use binary search
    std::vector<double> binEdgesY(numBinsY + 1);
    for (int i = 0; i < numBinsY; ++i) {
      binEdgesY[i] = binningY.binLow(i);
    }
    binEdgesY[numBinsY] = binningY.binHigh(numBinsY - 1);

    // Pass vector to ctx
    binIdxY = ctx.buildCall("RooFit::Detail::MathFuncs::rawBinNumber", arg.y(), binEdgesY, numBinsY + 1);
  }

  std::stringstream binIdx;
  binIdx << "(" << binIdxY << " + " << yVar.numBins() << " * " << binIdxX << ")";

  std::string arrName = ctx.getTmpVarName();

  std::stringstream code;
  code << "double " << arrName << "[" << (numBinsX * numBinsY) << "];\n";
  code << ctx.buildCall("RooFit::Detail::MathFuncs::fastVerticalInterpHistPdf2D2",
                        numBinsX,
                        numBinsY,
                        nCoefs,
                        arg.coefList(),
                        nominalVec,
                        widthVec,
                        morphsVecSum,
                        morphsVecDiff,
                        arg.smoothRegion(),
                        arrName) +
              ";\n";

  ctx.addToCodeBody(code.str(), true);
  ctx.addResult(&arg, arrName + "[" + binIdx.str() + "]");
}

void RooFit::Experimental::codegenImpl(VerticalInterpPdf& arg, CodegenContext& ctx) {
  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::verticalInterpolate",
                              arg.coefList(),
                              arg.coefList().size(),
                              arg.funcList(),
                              arg.funcList().size(),
                              arg.pdfFloorVal(),
                              arg.quadraticRegion(),
                              arg.quadraticAlgo()));
}

void RooFit::Experimental::codegenImpl(RooParametricHist& arg, CodegenContext& ctx) {
  std::vector<double> diffs_flat;
  std::vector<double> sums_flat;
  if (arg.hasMorphs()) {
    arg.getFlattenedMorphs(diffs_flat, sums_flat);
  }

  std::string xName = ctx.getResult(arg.observable());
  std::stringstream bin_i;
  bin_i << ctx.buildCall("RooFit::Detail::MathFuncs::parametricHistFindBin", arg.getNBins(), arg.getBins(), xName);
  std::stringstream code;
  code << ctx.buildCall("RooFit::Detail::MathFuncs::parametricHistEvaluate",
                        bin_i.str(),
                        arg.getPars(),
                        arg.getBins(),
                        arg.getNBins(),
                        arg.getCoeffList(),
                        static_cast<int>(arg.getCoeffList().size()),
                        arg.hasMorphs() ? diffs_flat : std::vector<double>{},
                        arg.hasMorphs() ? sums_flat : std::vector<double>{},
                        arg.getWidths(),
                        arg.getSmoothRegion()) +
              ";\n";
  ctx.addToCodeBody(code.str());
  ctx.addResult(&arg, code.str());
}

void RooFit::Experimental::codegenImpl(RooPower& arg, CodegenContext& ctx) {
  ctx.addResult(&arg, ctx.buildCall("RooFit::Detail::MathFuncs::powerEvaluate", arg.base(), arg.exponent()));
}

std::string RooFit::Experimental::codegenIntegralImpl(VerticalInterpPdf& arg,
                                                      int code,
                                                      const char* rangeName,
                                                      CodegenContext& ctx) {
  return ctx.buildCall("RooFit::Detail::MathFuncs::verticalInterpPdfIntegral",
                       arg.coefList(),
                       arg.coefList().size(),
                       arg.funcIntListFromCache(),
                       arg.funcIntListFromCache().size(),
                       arg.pdfFloorVal(),
                       arg.integralFloorVal(),
                       arg.quadraticRegion(),
                       arg.quadraticAlgo());
}

std::string RooFit::Experimental::codegenIntegralImpl(RooParametricHist& arg,
                                                      int code,
                                                      const char* rangeName,
                                                      CodegenContext& ctx) {
  std::vector<double> diffs_flat;
  std::vector<double> sums_flat;
  if (arg.hasMorphs()) {
    arg.getFlattenedMorphs(diffs_flat, sums_flat);
  }
  double xMin = rangeName ? arg.getObs().getMin(rangeName) : -1;
  double xMax = rangeName ? arg.getObs().getMax(rangeName) : -1;  // if rangeName is null, these values won't be used

  return ctx.buildCall("RooFit::Detail::MathFuncs::parametricHistIntegral",
                       arg.getPars(),
                       arg.getBins(),
                       arg.getNBins(),
                       arg.getCoeffList(),
                       static_cast<int>(arg.getCoeffList().size()),
                       arg.hasMorphs() ? diffs_flat : std::vector<double>{},
                       arg.hasMorphs() ? sums_flat : std::vector<double>{},
                       arg.getWidths(),
                       arg.getSmoothRegion(),
                       rangeName ? std::string("\"") + rangeName + "\"" : std::string("nullptr"),
                       xMin,
                       xMax);
}

std::string RooFit::Experimental::codegenIntegralImpl(RooPower& arg,
                                                      int code,
                                                      const char* rangeName,
                                                      CodegenContext& ctx) {
  if (code != 1) {
    throw std::invalid_argument("RooPower: unsupported analytical integral code");
  }

  const auto& x = dynamic_cast<const RooAbsRealLValue&>(arg.base());
  return ctx.buildCall(
      "RooFit::Detail::MathFuncs::powerIntegral", x.getMin(rangeName), x.getMax(rangeName), arg.exponent());
}

// Code generation for the template histograms with bin-by-bin uncertainties
// (CMSHistSum, CMSHistErrorPropagator and their CMSHistFunc processes). All
// three are evaluated by RooFit::Detail::MathFuncs::cmsHistSum() from the
// same flattened description of their templates.

struct CMSHistTemplates {
  int nBins = 0;
  int nProcs = 0;
  int nMorphs = 0;
  std::vector<double> templates;     // slots of nBins values
  std::vector<double> nominalSlot;   // per process
  std::vector<double> morphSlot;     // per (process, morph): slot of the sum template, or -1
  std::vector<double> binErrors;     // per (process, bin)
  std::vector<double> vertType;      // per process: 0 quadratic-linear, 1 logarithmic
  std::vector<double> smoothRegion;  // per process
  std::vector<double> binType;       // per bin: 0 none, 1 total parameter, 2 per-process parameters
  std::vector<double> procBinType;   // per (bin, process): 2 Poisson, 3 Gaussian, 4 none
  std::vector<double> binParIndex;   // per (bin, process): index in the bin parameter list, or -1
  [[nodiscard]] int numSlots() const { return static_cast<int>(templates.size()) / nBins; }
};

static void appendTemplate(CMSHistTemplates &t, FastTemplate const &tmpl) {
  for (int j = 0; j < t.nBins; ++j) {
    t.templates.push_back(tmpl[j]);
  }
}

// The bin parameters are listed bin by bin, and within a bin process by
// process, skipping the processes without parameter (type 4), like in the
// initialize() functions of CMSHistSum and CMSHistErrorPropagator.
static void setBinTypes(CMSHistTemplates &t, std::vector<std::vector<unsigned>> const &binTypes) {
  t.binType.assign(t.nBins, 0.);
  t.procBinType.assign(static_cast<std::size_t>(t.nBins) * t.nProcs, 0.);
  t.binParIndex.assign(static_cast<std::size_t>(t.nBins) * t.nProcs, -1.);
  int next = 0;
  for (std::size_t j = 0; j < binTypes.size(); ++j) {
    std::vector<unsigned> const &types = binTypes[j];
    if (types[0] == 0U) {
      continue;
    }
    if (types[0] == 1U) {
      t.binType[j] = 1.;
      t.binParIndex[j * t.nProcs] = next++;
      continue;
    }
    t.binType[j] = 2.;
    for (std::size_t i = 0; i < types.size(); ++i) {
      t.procBinType[(j * t.nProcs) + i] = types[i];
      if (types[i] >= 1U && types[i] < 4U) {
        t.binParIndex[(j * t.nProcs) + i] = next++;
      }
    }
  }
}

// Appends the nominal template of a CMSHistFunc and the sum/diff templates of
// its vertical morphings, converted like in the constructor of CMSHistSum.
// morphList is the morphing parameter list of the generated code.
static void appendHistFunc(CMSHistTemplates &t, int ip, CMSHistFunc const &func, RooArgList const &morphList) {
  if (!func.hmorphList().empty() || func.isRebinned() || !func.externalMorphList().empty()) {
    throw std::runtime_error(std::string("codegen: CMSHistFunc ") + func.GetName() +
                             " uses horizontal morphing, rebinning or an external morph, which is not supported");
  }
  std::vector<FastTemplate> const &storage = func.storage();
  FastTemplate const &nominal = storage[func.getIdx(/*hindex=*/0, /*hpoint=*/0, /*vindex=*/0, /*vpoint=*/0)];
  t.nominalSlot[ip] = t.numSlots();
  appendTemplate(t, nominal);
  t.vertType[ip] = func.verticalType() == CMSHistFunc::VerticalSetting::LogQuadLinear ? 1. : 0.;
  t.smoothRegion[ip] = func.verticalSmoothRegion();
  for (std::size_t iv = 0; iv < func.vmorphList().size(); ++iv) {
    FastTemplate lo = storage[func.getIdx(/*hindex=*/0, /*hpoint=*/0, iv + 1, /*vpoint=*/0)];
    FastTemplate hi = storage[func.getIdx(/*hindex=*/0, /*hpoint=*/0, iv + 1, /*vpoint=*/1)];
    if (func.verticalType() == CMSHistFunc::VerticalSetting::LogQuadLinear) {
      hi.LogRatio(nominal);
      lo.LogRatio(nominal);
    } else {
      hi.Subtract(nominal);
      lo.Subtract(nominal);
    }
    FastTemplate sum = hi;
    FastTemplate diff = hi;
    FastTemplate::SumDiff(hi, lo, sum, diff);
    t.morphSlot[(ip * t.nMorphs) + morphList.index(func.vmorphList()[static_cast<int>(iv)].GetName())] = t.numSlots();
    appendTemplate(t, sum);
    appendTemplate(t, diff);
  }
}

// Names of the generated variables that hold the integral, the summed yields
// and the total errors of the histogram. They are derived from the node name,
// like RooFit does for its likelihood results, so that the translations of
// the integral and of the gradient term find them without shared state.
static std::string cmsHistArrayName(RooAbsArg const &arg, std::string const &what) {
  std::string name = arg.GetName();
  for (char &c : name) {
    if (std::isalnum(static_cast<unsigned char>(c)) == 0) {
      c = '_';
    }
  }
  return "cmsHist_" + name + "_" + what;
}

// The extended RooRealSumPdf of the autoMCStats channels needs the fix of the
// code placement of its integral in RooFit (realSumPdfTranslateImpl in
// roofit/codegen): before, the generated function doesn't compile.
static void checkRooFitVersionForRealSumPdf(RooAbsArg const &arg) {
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 41, 1)
  throw std::runtime_error(std::string("codegen: the translation of ") + arg.GetName() +
                           " needs ROOT >= 6.41.01 for the integral of an extended RooRealSumPdf");
#else
  (void)arg;
#endif
}

// Expression for the bin of the observable, like FastHisto::GetAt(). Unlike
// that function, which returns zero outside of the histogram range, the
// expression clamps to the first and last bin.
static std::string cmsHistBinIndex(RooFit::Experimental::CodegenContext &ctx,
                                   FastHisto const &histo,
                                   RooAbsReal const &xVar,
                                   int nBins) {
  std::vector<double> edges(nBins + 1);
  for (int i = 0; i <= nBins; ++i) {
    edges[i] = histo.GetEdge(i);
  }
  const double low = edges.front();
  const double high = edges.back();
  const double width = (high - low) / nBins;
  bool uniform = true;
  for (int i = 0; i < nBins; ++i) {
    uniform = uniform && std::abs(edges[i + 1] - edges[i] - width) < 1e-9 * (high - low);
  }
  if (uniform) {
    return ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", low, high, xVar, nBins, 1.);
  }
  return ctx.buildCall("RooFit::Detail::MathFuncs::rawBinNumber", xVar, edges, nBins + 1);
}

// Emits the evaluation of the histogram outside of the event loop and
// registers the value in the bin of the observable as the result of the node.
static void emitCMSHistSum(RooAbsArg const &arg,
                           RooFit::Experimental::CodegenContext &ctx,
                           CMSHistTemplates const &t,
                           std::string const &coefArg,
                           std::string const &morphArg,
                           std::string const &binParArg,
                           FastHisto const &histo,
                           RooAbsReal const &xVar,
                           bool poissonUsesCropped) {
  std::vector<double> widths(t.nBins);
  for (int j = 0; j < t.nBins; ++j) {
    widths[j] = histo.GetWidth(j);
  }
  std::string const binIdx = cmsHistBinIndex(ctx, histo, xVar, t.nBins);
  std::string const values = ctx.getTmpVarName();
  std::string const scratch = ctx.getTmpVarName();
  std::string const valsum = cmsHistArrayName(arg, "valsum");
  std::string const toterr = cmsHistArrayName(arg, "toterr");

  std::stringstream code;
  code << "double " << values << "[" << t.nBins << "];\n";
  code << "double " << scratch << "[" << 3 * t.nBins << "];\n";
  code << "double " << valsum << "[" << t.nBins << "];\n";
  code << "double " << toterr << "[" << t.nBins << "];\n";
  code << "const double " << cmsHistArrayName(arg, "integral") << " = "
       << ctx.buildCall("RooFit::Detail::MathFuncs::cmsHistSum",
                        t.nBins,
                        t.nProcs,
                        t.nMorphs,
                        coefArg,
                        morphArg,
                        binParArg,
                        t.templates,
                        t.nominalSlot,
                        t.morphSlot,
                        t.binErrors,
                        widths,
                        t.vertType,
                        t.smoothRegion,
                        t.binType,
                        t.procBinType,
                        t.binParIndex,
                        poissonUsesCropped ? 1 : 0,
                        scratch,
                        valsum,
                        toterr,
                        values)
       << ";\n";
  ctx.addToCodeBody(code.str(), /*isScopeIndep=*/true);
  ctx.addResult(&arg, values + "[" + binIdx + "]");
}

// The integral is computed together with the bin contents, see above. The
// node result is requested first, in case the integral is translated before
// the node itself.
static std::string cmsHistIntegral(RooAbsArg &arg, RooFit::Experimental::CodegenContext &ctx) {
  ctx.getResult(arg);
  return cmsHistArrayName(arg, "integral");
}

void RooFit::Experimental::codegenImpl(CMSHistSum &arg, CodegenContext &ctx) {
  checkRooFitVersionForRealSumPdf(arg);
  if (!arg.externalMorphList().empty()) {
    throw std::runtime_error(std::string("codegen: CMSHistSum ") + arg.GetName() +
                             " uses an external morph, which is not supported");
  }
  CMSHistTemplates t;
  t.nBins = static_cast<int>(arg.cache().size());
  t.nProcs = arg.numProcesses();
  t.nMorphs = arg.numMorphs();
  // The storage of CMSHistSum already has the layout of the generated code:
  // the nominal template of every process followed by the sum and diff
  // templates of its morphings.
  for (FastTemplate const &tmpl : arg.storage()) {
    appendTemplate(t, tmpl);
  }
  t.nominalSlot.assign(arg.processFields().begin(), arg.processFields().end());
  t.morphSlot.assign(arg.vmorphFields().begin(), arg.vmorphFields().end());
  for (int ip = 0; ip < t.nProcs; ++ip) {
    for (int j = 0; j < t.nBins; ++j) {
      t.binErrors.push_back(arg.binErrors()[ip][j]);
    }
    t.vertType.push_back(arg.verticalTypes()[ip] == CMSHistFunc::VerticalSetting::LogQuadLinear ? 1. : 0.);
    t.smoothRegion.push_back(arg.verticalSmoothRegions()[ip]);
  }
  setBinTypes(t, arg.binTypes());
  emitCMSHistSum(arg,
                 ctx,
                 t,
                 ctx.buildArg(arg.coefList()),
                 ctx.buildArg(arg.morphList()),
                 ctx.buildArg(arg.binParList()),
                 arg.cache(),
                 arg.getXVar(),
                 /*poissonUsesCropped=*/false);
}

// A CMSHistErrorPropagator over its processes, or a lone CMSHistFunc as a
// propagator with one process, a unit coefficient and no bin parameters.
static void translateHistFuncs(RooAbsArg &arg,
                               RooFit::Experimental::CodegenContext &ctx,
                               std::vector<CMSHistFunc const *> const &funcs,
                               std::string const &coefArg,
                               std::string const &binParArg,
                               std::vector<std::vector<unsigned>> const &binTypes,
                               RooAbsReal const &xVar) {
  // The union of the morphing parameters, ordered by name like in CMSHistSum.
  RooArgSet morphSet;
  for (CMSHistFunc const *func : funcs) {
    morphSet.add(func->vmorphList(), /*silent=*/true);
  }
  morphSet.sort();
  RooArgList const morphList(morphSet);

  CMSHistTemplates t;
  t.nBins = static_cast<int>(funcs[0]->cache().size());
  t.nProcs = static_cast<int>(funcs.size());
  t.nMorphs = static_cast<int>(morphList.size());
  t.nominalSlot.assign(t.nProcs, -1.);
  t.morphSlot.assign(static_cast<std::size_t>(t.nProcs) * t.nMorphs, -1.);
  t.vertType.assign(t.nProcs, 0.);
  t.smoothRegion.assign(t.nProcs, 1.);
  for (int ip = 0; ip < t.nProcs; ++ip) {
    appendHistFunc(t, ip, *funcs[ip], morphList);
    for (int j = 0; j < t.nBins; ++j) {
      t.binErrors.push_back(funcs[ip]->errors()[j]);
    }
  }
  setBinTypes(t, binTypes);
  emitCMSHistSum(
      arg, ctx, t, coefArg, ctx.buildArg(morphList), binParArg, funcs[0]->cache(), xVar, /*poissonUsesCropped=*/true);
}

void RooFit::Experimental::codegenImpl(CMSHistErrorPropagator &arg, CodegenContext &ctx) {
  checkRooFitVersionForRealSumPdf(arg);
  std::vector<CMSHistFunc const *> funcs;
  for (RooAbsArg const *func : arg.funcList()) {
    auto const *histFunc = dynamic_cast<CMSHistFunc const *>(func);
    if (!histFunc) {
      throw std::runtime_error(std::string("codegen: process ") + func->GetName() + " of " + arg.GetName() +
                               " is not a CMSHistFunc");
    }
    funcs.push_back(histFunc);
  }
  if (funcs.empty()) {
    throw std::runtime_error(std::string("codegen: CMSHistErrorPropagator ") + arg.GetName() + " has no processes");
  }
  translateHistFuncs(
      arg, ctx, funcs, ctx.buildArg(arg.coefList()), ctx.buildArg(arg.binParList()), arg.binTypes(), arg.getXVar());
}

void RooFit::Experimental::codegenImpl(CMSHistFunc &arg, CodegenContext &ctx) {
  std::vector<double> const unitCoef{1.};
  translateHistFuncs(arg, ctx, {&arg}, ctx.buildArg(unitCoef), "nullptr", {}, arg.getXVar());
}

void RooFit::Experimental::codegenImpl(CMSHistBarlowBeestonGradientTerm &arg, CodegenContext &ctx) {
  ctx.addResult(&arg, "0.0");
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistBarlowBeestonGradientTerm &arg,
                                                      int /*code*/,
                                                      const char * /*rangeName*/,
                                                      CodegenContext &ctx) {
  // The histogram comes first in the function list of the RooRealSumPdf that
  // holds the term, so its arrays have been emitted already (the name of the
  // original histogram is the name of its copy inside the likelihood). The
  // sum is written out as a loop: clad does not differentiate a call that
  // reads these arrays correctly.
  RooAbsReal const &hist = arg.histogram();
  std::string const valsum = cmsHistArrayName(hist, "valsum");
  std::string const toterr = cmsHistArrayName(hist, "toterr");
  std::string const slopes = ctx.buildArg(arg.slopes());
  std::string const ratios = ctx.buildArg(arg.ratios());
  std::string const sum = ctx.getTmpVarName();
  std::stringstream body;
  body << "double " << sum << " = 0.0;\n";
  for (std::size_t k = 0; k < arg.bins().size(); ++k) {
    std::string const j = std::to_string(arg.bins()[k]);
    body << sum << " += " << slopes << "[" << k << "] * ((" << toterr << "[" << j << "] > 0.0 ? " << valsum << "[" << j
         << "] / " << toterr << "[" << j << "] : 0.0) - " << ratios << "[" << k << "]);\n";
  }
  ctx.addToCodeBody(body.str(), /*isScopeIndep=*/true);
  return sum;
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistSum &arg,
                                                      int /*code*/,
                                                      const char * /*rangeName*/,
                                                      CodegenContext &ctx) {
  return cmsHistIntegral(arg, ctx);
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistErrorPropagator &arg,
                                                      int /*code*/,
                                                      const char * /*rangeName*/,
                                                      CodegenContext &ctx) {
  return cmsHistIntegral(arg, ctx);
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistFunc &arg,
                                                      int /*code*/,
                                                      const char * /*rangeName*/,
                                                      CodegenContext &ctx) {
  return cmsHistIntegral(arg, ctx);
}

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6,36,0)
