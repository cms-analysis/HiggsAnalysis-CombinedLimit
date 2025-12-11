#include "../interface/CombineCodegenImpl.h"

#include <ROOT/RConfig.hxx>  // for ROOT_VERSION

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 36, 0)

#include "../interface/AsymPow.h"
#include "../interface/CMSHistErrorPropagator.h"
#include "../interface/CMSHistFunc.h"
#include "../interface/CMSHistSum.h"
#include "../interface/CombineMathFuncs.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/VerticalInterpHistPdf.h"
#include "../interface/VerticalInterpPdf.h"

#include <RooUniformBinning.h>

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

void RooFit::Experimental::codegenImpl(CMSHistFunc& arg, CodegenContext& ctx) {
  arg.evaluate(); // trigger cache() creation
  std::vector<double> const& edges = arg.cache().GetBinEdges();

  // I don't know if these values are actually constant and we can take them
  // here to hardcode into the generated code...
  auto const& values = arg.cache().GetValues();

  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::cmsHistFunc",
                              arg.getXVar(),
                              edges.size() - 1,
                              edges,
                              values
                              ));
}

void RooFit::Experimental::codegenImpl(CMSHistErrorPropagator& arg, CodegenContext& ctx) {
  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::cmsHistErrorPropagator",
                              arg.getXVar(),
                              arg.coefList().size(),
                              arg.coefList(),
                              arg.funcList()));
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistErrorPropagator& arg,
                                                      int code,
                                                      const char* rangeName,
                                                      CodegenContext& ctx) {
  return "2.0";  // TODO: dummy for now
}

void RooFit::Experimental::codegenImpl(CMSHistSum& arg, CodegenContext& ctx) {
  arg.evaluate();
  std::vector<double> const& edges = arg.cache().GetBinEdges();
  std::size_t nBins = edges.size() - 1;
  RooArgList const& coefs = arg.coefList();
  std::size_t nSamples = coefs.size();
  std::vector<double> values(nBins * nSamples);
  for (std::size_t iSamples = 0; iSamples < nSamples; ++iSamples) {
    std::vector<double> sampleValues = arg.getFuncValList(iSamples);
    for (std::size_t iBin = 0; iBin < nBins; ++iBin)
      values[iBin + iSamples * nBins] = sampleValues[iBin];
  }

  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::cmsHistSum",
                              arg.getXVar(),
                              nBins,
                              nSamples,
                              coefs,
                              edges,
                              values
                              ));
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistSum& arg,
                                                      int code,
                                                      const char* rangeName,
                                                      CodegenContext& ctx) {
  return "3.0";  // TODO: dummy for now
}

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6,36,0)
