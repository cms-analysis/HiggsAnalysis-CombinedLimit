#ifndef CombineMathFuncs_h
#define CombineMathFuncs_h

#include <cmath>

namespace RooFit {
namespace Detail {
namespace MathFuncs {

inline double smoothStepFunc(double x, double smoothRegion)
{
   if (std::abs(x) >= smoothRegion)
      return x > 0 ? +1 : -1;
   double xnorm = x / smoothRegion;
   double xnorm2 = xnorm * xnorm;
   return 0.125 * xnorm * (xnorm2 * (3. * xnorm2 - 10.) + 15);
}

inline double
fastVerticalInterpHistPdf2(int nBins, int binIdx, int nCoefs, double const *coefs, double const *nominal, double const *binWidth,
                           double const *morphsSum, double const *morphsDiff, double smoothRegion, int offsetIdx)
{
   double out[nBins];

   coefs = coefs + offsetIdx;
   nominal = nominal + offsetIdx;
   binWidth = binWidth + offsetIdx;
   morphsSum = morphsSum + offsetIdx;
   morphsDiff = morphsDiff + offsetIdx;

   for (int iBin = 0; iBin < nBins; ++iBin) {
      out[iBin] = nominal[iBin];
   }

   double normSum = 0.0;

   for (int iBin = 0; iBin < nBins; ++iBin) {
      // apply all morphs one by one
      for (int iCoef = 0; iCoef < nCoefs; ++iCoef) {
         double const *sum = morphsSum + iCoef * nBins;
         double const *diff = morphsDiff + iCoef * nBins;
         double x = coefs[iCoef];
         double a = 0.5 * x;
         double b = smoothStepFunc(x, smoothRegion);

         out[iBin] += a * (diff[iBin] + b * sum[iBin]);
      }

      out[iBin] = std::max(1e-9, out[iBin]);
      normSum += out[iBin] * binWidth[iBin];
   }

   if (normSum > 0.0) {
      double normSumInv = 1. / normSum;
      for (int iBin = 0; iBin < nBins; ++iBin) {
         out[iBin] *= normSumInv;
      }
   }

   return out[binIdx];
}

inline double logKappaForX(double theta, double logKappaLow, double logKappaHigh)
{
   double logKappa = 0.0;

   if (std::abs(theta) >= 0.5) {
      logKappa = theta >= 0 ? logKappaHigh : -logKappaLow;
   } else {
      // interpolate between log(kappaHigh) and -log(kappaLow)
      //    logKappa(x) = avg + halfdiff * h(2x)
      // where h(x) is the 3th order polynomial
      //    h(x) = (3 x^5 - 10 x^3 + 15 x)/8;
      // chosen so that h(x) satisfies the following:
      //      h (+/-1) = +/-1
      //      h'(+/-1) = 0
      //      h"(+/-1) = 0
      double logKhi = logKappaHigh;
      double logKlo = -logKappaLow;
      double avg = 0.5 * (logKhi + logKlo);
      double halfdiff = 0.5 * (logKhi - logKlo);
      double twox = theta + theta;
      double twox2 = twox * twox;
      double alpha = 0.125 * twox * (twox2 * (3 * twox2 - 10.) + 15.);
      logKappa = avg + alpha * halfdiff;
   }

   return logKappa;
}

inline double asymPow(double theta, double kappaLow, double kappaHigh)
{
   return std::exp(logKappaForX(theta, std::log(kappaLow), std::log(kappaHigh)) * theta);
}

inline double processNormalization(double nominalValue, std::size_t nThetas, std::size_t nAsymmThetas,
                                   std::size_t nOtherFactors, double const *thetas, double const *logKappas,
                                   double const *asymmThetas, double const *asymmLogKappasLow,
                                   double const *asymmLogKappasHigh, double const *otherFactors)
{
   double logVal = 0.0;
   for (std::size_t i = 0; i < nThetas; i++) {
      logVal += thetas[i] * logKappas[i];
   }
   for (std::size_t i = 0; i < nAsymmThetas; i++) {
      double x = asymmThetas[i];
      logVal += x * logKappaForX(x, asymmLogKappasLow[i], asymmLogKappasHigh[i]);
   }
   double norm = nominalValue;
   norm *= std::exp(logVal);
   for (std::size_t i = 0; i < nOtherFactors; i++) {
      norm *= otherFactors[i];
   }
   return norm;
}

} // MathFuncs
} // Detail
} // RooFit

#endif
