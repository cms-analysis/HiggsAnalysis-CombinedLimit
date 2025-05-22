#include "../interface/CombineCodegenImpl.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)

#include "../interface/AsymPow.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/VerticalInterpHistPdf.h"
#include "../interface/VerticalInterpPdf.h"
#include "../interface/CombineMathFuncs.h"

#include <RooUniformBinning.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
namespace RooFit {
namespace Experimental {
# define CODEGEN_IMPL(CLASS_NAME) void codegenImpl(CLASS_NAME &arg0, CodegenContext &ctx)
# define CODEGEN_INTEGRAL_IMPL(CLASS_NAME) std::string codegenIntegralImpl(CLASS_NAME &arg0, int code, const char *rangeName, CodegenContext &ctx)
# define ARG_VAR auto &arg = arg0;
#else
# define CODEGEN_IMPL(CLASS_NAME) void CLASS_NAME::translate(RooFit::Detail::CodeSquashContext &ctx) const
# define CODEGEN_INTEGRAL_IMPL(CLASS_NAME) std::string CLASS_NAME::buildCallToAnalyticIntegral(Int_t code, const char *rangeName, RooFit::Detail::CodeSquashContext &ctx) const
# define ARG_VAR auto &arg = *this;
#endif


CODEGEN_IMPL(AsymPow) {
  ARG_VAR;
  ctx.addResult(&arg,
                ctx.buildCall("RooFit::Detail::MathFuncs::asymPow", arg.theta(), arg.kappaLow(), arg.kappaHigh()));
}

CODEGEN_IMPL(ProcessNormalization) {
  ARG_VAR;

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

CODEGEN_IMPL(FastVerticalInterpHistPdf2) {
  ARG_VAR;

  if (arg.smoothAlgo() < 0) {
    throw std::runtime_error("We only support smoothAlgo >= 0");
  }

  RooRealVar const &xVar = arg.x();

  int numBins = xVar.numBins();

  std::vector<double> nominalVec(numBins);
  std::vector<double> widthVec(numBins);
  std::vector<double> morphsVecSum;
  std::vector<double> morphsVecDiff;

  auto const &cacheNominal = arg.cacheNominal();

  for (int i = 0; i < numBins; ++i) {
    nominalVec[i] = cacheNominal.GetBinContent(i);
    widthVec[i] = cacheNominal.GetWidth(i);
  }

  std::size_t nCoefs = arg.coefList().size();

  morphsVecSum.reserve(numBins * nCoefs);
  morphsVecDiff.reserve(numBins * nCoefs);
  auto const &morphs = arg.morphs();
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
  if (!dynamic_cast<RooUniformBinning const *>(&xVar.getBinning())) {
    throw std::runtime_error("We only support uniform binning!");
  }
  double xLow = xVar.getMin();
  double xHigh = xVar.getMax();
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
  std::string binIdx = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", xLow, xHigh, xVar, numBins, 1.);
#else
  std::string binIdx = ctx.buildCall("RooFit::Detail::MathFuncs::getUniformBinning", xLow, xHigh, xVar, numBins);
#endif

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

CODEGEN_IMPL(FastVerticalInterpHistPdf2D2) {
  ARG_VAR;

  if (arg.smoothAlgo() < 0) {
    throw std::runtime_error("We only support smoothAlgo >= 0");
  }

  if (!arg.conditional()) {
    throw std::runtime_error("We only support conditional == true");
  }

  RooRealVar const &xVar = arg.x();
  RooRealVar const &yVar = arg.y();

  // We also have to assert that x and y are uniformely binned!
  if (!dynamic_cast<RooUniformBinning const *>(&xVar.getBinning())) {
    throw std::runtime_error("We only support uniform binning!");
  }
  if (!dynamic_cast<RooUniformBinning const *>(&yVar.getBinning())) {
    throw std::runtime_error("We only support uniform binning!");
  }

  auto const &cacheNominal = arg.cacheNominal();

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
  auto const &morphs = arg.morphs();
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
  std::string binIdxX = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", xLow, xHigh, arg.x(), numBinsX, 1.);
#else
  std::string binIdxX = ctx.buildCall("RooFit::Detail::MathFuncs::getUniformBinning", xLow, xHigh, arg.x(), numBinsX);
#endif
  double yLow = yVar.getMin();
  double yHigh = yVar.getMax();
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
  std::string binIdxY = ctx.buildCall("RooFit::Detail::MathFuncs::uniformBinNumber", yLow, yHigh, arg.y(), numBinsY, 1.);
#else
  std::string binIdxY = ctx.buildCall("RooFit::Detail::MathFuncs::getUniformBinning", yLow, yHigh, arg.y(), numBinsY);
#endif

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

CODEGEN_IMPL(VerticalInterpPdf) {
  ARG_VAR;
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

CODEGEN_INTEGRAL_IMPL(VerticalInterpPdf) {
  ARG_VAR;
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

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
} // namespace RooFit
} // namespace Experimental
#endif

#endif // ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
