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

  // The bin index part
  // We also have to assert that x is uniformely binned!
  if (!dynamic_cast<RooUniformBinning const*>(&xVar.getBinning())) {
    throw std::runtime_error("We only support uniform binning!");
  }
  double xLow = xVar.getMin();
  double xHigh = xVar.getMax();
  std::string binIdx = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", xLow, xHigh, xVar, numBins, 1.);

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

  // We also have to assert that x and y are uniformely binned!
  if (!dynamic_cast<RooUniformBinning const*>(&xVar.getBinning())) {
    throw std::runtime_error("We only support uniform binning!");
  }
  if (!dynamic_cast<RooUniformBinning const*>(&yVar.getBinning())) {
    throw std::runtime_error("We only support uniform binning!");
  }

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

  // The bin index part
  double xLow = xVar.getMin();
  double xHigh = xVar.getMax();
  std::string binIdxX =
      ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", xLow, xHigh, arg.x(), numBinsX, 1.);
  double yLow = yVar.getMin();
  double yHigh = yVar.getMax();
  std::string binIdxY =
      ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", yLow, yHigh, arg.y(), numBinsY, 1.);

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
  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::cmsHistSum",
                              arg.getXVar(),
                              arg.coefList().size(),
                              arg.coefList(),
                              nullptr // arg.funcList() what should this be??
                              ));
}

std::string RooFit::Experimental::codegenIntegralImpl(CMSHistSum& arg,
                                                      int code,
                                                      const char* rangeName,
                                                      CodegenContext& ctx) {
  return "3.0";  // TODO: dummy for now
}

#endif  // ROOT_VERSION_CODE >= ROOT_VERSION(6,36,0)
