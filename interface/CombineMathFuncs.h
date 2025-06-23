#ifndef CombineMathFuncs_h
#define CombineMathFuncs_h

#include <RooAbsReal.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooConstVar.h>
#include <RtypesCore.h>

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

inline void fastVerticalInterpHistPdf2(int nBins, int nCoefs, double const *coefs, double const *nominal,
                                       double const *binWidth, double const *morphsSum, double const *morphsDiff,
                                       double smoothRegion, double *out)
{
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
}

inline void fastVerticalInterpHistPdf2D2(int nBinsX, int nBinsY, int nCoefs, double const *coefs,
                                           double const *nominal, double const *binWidth, double const *morphsSum,
                                           double const *morphsDiff, double smoothRegion, double *out)
{
   int nBins = nBinsX * nBinsY;

   for (int iBin = 0; iBin < nBins; ++iBin) {
      out[iBin] = nominal[iBin];
   }

   for (int iBinX = 0; iBinX < nBinsX; ++iBinX) {

      double normSum = 0.0;

      for (int iBinY = 0; iBinY < nBinsY; ++iBinY) {
         int iBin = iBinY + nBinsY * iBinX;
         // apply all morphs one by one
         for (int iCoef = 0; iCoef < nCoefs; ++iCoef) {
            double const *sum = morphsSum + iCoef * nBinsY;
            double const *diff = morphsDiff + iCoef * nBinsY;
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
         for (int iBinY = 0; iBinY < nBinsY; ++iBinY) {
            int iBin = iBinY + nBinsY * iBinX;
            out[iBin] *= normSumInv;
         }
      }
   }
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

// Interpolation  (from VerticalInterpPdf)
inline Double_t interpolate(Double_t const coeff, Double_t const central, Double_t const fUp, 
                            Double_t const fDn, Double_t const quadraticRegion, Int_t const quadraticAlgo)
{
   if (quadraticAlgo == -1) {
      Double_t kappa = (coeff > 0 ? fUp/central : central/fDn);
      return pow(kappa, sqrt(pow(coeff, 2)));
   }

   if (fabs(coeff) >= quadraticRegion) {
      return coeff * (coeff > 0 ? fUp - central : central - fDn);
   }
   // quadratic interpolation coefficients between the three
   if (quadraticAlgo == 0) {
      // quadratic interpolation null at zero and continuous at boundaries, but not differentiable at boundaries
      // conditions:
      //   c_up (+quadraticRegion) = +quadraticRegion
      //   c_cen(+quadraticRegion) = -quadraticRegion
      //   c_dn (+quadraticRegion) = 0
      //   c_up (-quadraticRegion) = 0 
      //   c_cen(-quadraticRegion) = -quadraticRegion
      //   c_dn (-quadraticRegion) = +quadraticRegion
      //   c_up(0) = c_dn(0) = c_cen(0) = 0
      Double_t c_up  = + coeff * (quadraticRegion + coeff) / (2 * quadraticRegion);
      Double_t c_dn  = - coeff * (quadraticRegion - coeff) / (2 * quadraticRegion);
      Double_t c_cen = - coeff * coeff / quadraticRegion;
      return (c_up * fUp) + (c_dn * fDn) + (c_cen * central);
   } 
   if (quadraticAlgo == 1) { 
      // quadratic interpolation that is everywhere differentiable, but it's not null at zero
      // conditions on the function
      //   c_up (+quadraticRegion) = +quadraticRegion
      //   c_cen(+quadraticRegion) = -quadraticRegion
      //   c_dn (+quadraticRegion) = 0
      //   c_up (-quadraticRegion) = 0 
      //   c_cen(-quadraticRegion) = -quadraticRegion
      //   c_dn (-quadraticRegion) = +quadraticRegion
      // conditions on the derivatives
      //   c_up '(+quadraticRegion) = +1
      //   c_cen'(+quadraticRegion) = -1
      //   c_dn '(+quadraticRegion) = 0
      //   c_up '(-quadraticRegion) = 0
      //   c_cen'(-quadraticRegion) = +1
      //   c_dn '(-quadraticRegion) = -1
      Double_t c_up  = (quadraticRegion + coeff) * (quadraticRegion + coeff) / (4 * quadraticRegion);
      Double_t c_dn  = (quadraticRegion - coeff) * (quadraticRegion - coeff) / (4 * quadraticRegion);
      Double_t c_cen = - c_up - c_dn;
      return (c_up * fUp) + (c_dn * fDn) + (c_cen * central);
   }
   // P(6) interpolation that is everywhere differentiable and null at zero
   /* === how the algorithm works, in theory ===
   * let  dhi = h_hi - h_nominal
   *      dlo = h_lo - h_nominal
   * and x be the morphing parameter
   * we define alpha = x * 0.5 * ((dhi-dlo) + (dhi+dlo)*smoothStepFunc(x));
   * which satisfies:
   *     alpha(0) = 0
   *     alpha(+1) = dhi
   *     alpha(-1) = dlo
   *     alpha(x >= +1) = |x|*dhi
   *     alpha(x <= -1) = |x|*dlo
   *     alpha is continuous and has continuous first and second derivative, as smoothStepFunc has them
   * === and in practice ===
   * we already have computed the histogram for diff=(dhi-dlo) and sum=(dhi+dlo)
   * so we just do template += (0.5 * x) * (diff + smoothStepFunc(x) * sum)
   * ========================================== */
   Double_t cnorm = coeff/quadraticRegion;
   Double_t cnorm2 = pow(cnorm, 2);
   Double_t hi = fUp - central;
   Double_t lo = fDn - central;
   Double_t sum = hi+lo;
   Double_t diff = hi-lo;
   Double_t a = coeff/2.; // cnorm*quadraticRegion
   Double_t b = 0.125 * cnorm * (cnorm2 * (3.*cnorm2 - 10.) + 15.);
   Double_t result = a*(diff + b*sum);
   return result;
}

template <typename Operation>
inline Double_t opInterpolate(RooArgList const& coefList, RooArgList const& funcList, Double_t const pdfFloorVal,
                                    Double_t const quadraticRegion, Int_t const quadraticAlgo, const RooArgSet* normSet2=nullptr)
{
   // Do running sum of coef/func pairs, calculate lastCoef.
   RooAbsReal* func = &dynamic_cast<RooAbsReal&>(funcList[0]);
   Double_t central = func->getVal();
   Double_t value = central;

   Operation op;

   for (int iCoef = 0; iCoef < coefList.getSize(); ++iCoef) {
      Double_t coefVal = dynamic_cast<RooAbsReal&>(coefList[iCoef]).getVal(normSet2);
      RooAbsReal* funcUp = &dynamic_cast<RooAbsReal&>(funcList[(2 * iCoef) + 1]);
      RooAbsReal* funcDn = &dynamic_cast<RooAbsReal&>(funcList[(2 * iCoef) + 2]);
      value = op(value, interpolate(coefVal, central, funcUp->getVal(), funcDn->getVal(), quadraticRegion, quadraticAlgo));
   }
   return ( value > 0. ? value : pdfFloorVal);
}

inline Double_t additiveInterpolate(double const* coefList, std::size_t nCoeffs, double const* funcList, std::size_t nFuncs,
                                  Double_t const pdfFloorVal, Double_t const quadraticRegion, Int_t const quadraticAlgo)
{
   // Do running sum of coef/func pairs, calculate lastCoef.
   Double_t central = funcList[0];
   Double_t value = central;

   for (std::size_t iCoef = 0; iCoef < nCoeffs; ++iCoef) {
      double coefVal = coefList[iCoef];
      double funcUp = funcList[(2 * iCoef) + 1];
      double funcDn = funcList[(2 * iCoef) + 2];
      value += interpolate(coefVal, central, funcUp, funcDn, quadraticRegion, quadraticAlgo);
   }
   return ( value > 0. ? value : pdfFloorVal);
}

inline Double_t multiplicativeInterpolate(double const* coefList, std::size_t nCoeffs, double const* funcList, std::size_t nFuncs,
                                  Double_t const pdfFloorVal, Double_t const quadraticRegion, Int_t const quadraticAlgo)
{
   // Do running sum of coef/func pairs, calculate lastCoef.
   Double_t central = funcList[0];
   Double_t value = central;

   for (std::size_t iCoef = 0; iCoef < nCoeffs; ++iCoef) {
      double coefVal = coefList[iCoef];
      double funcUp = funcList[(2 * iCoef) + 1];
      double funcDn = funcList[(2 * iCoef) + 2];
      value *= interpolate(coefVal, central, funcUp, funcDn, quadraticRegion, quadraticAlgo);
   }
   return ( value > 0. ? value : pdfFloorVal);
}

inline Double_t verticalInterpolate(double const* coefList, std::size_t nCoeffs, double const* funcList, std::size_t nFuncs,
                                    double const pdfFloorVal, double const quadraticRegion, Int_t const quadraticAlgo)
{
   // Do running sum of coef/func pairs, calculate lastCoef.
   Double_t value = pdfFloorVal;
   if (quadraticAlgo >= 0) {
      value = RooFit::Detail::MathFuncs::additiveInterpolate(coefList, nCoeffs, funcList, nFuncs, pdfFloorVal, quadraticRegion, quadraticAlgo);
   } else {
      value = RooFit::Detail::MathFuncs::multiplicativeInterpolate(coefList, nCoeffs, funcList, nFuncs, pdfFloorVal, quadraticRegion, quadraticAlgo);
   }
   return value;
}

inline Double_t verticalInterpPdfIntegral(double const* coefList, std::size_t nCoeffs, double const* funcIntList, std::size_t nFuncs,
                                          double const pdfFloorVal, double const integralFloorVal,
                                          double const quadraticRegion, Int_t const quadraticAlgo)
{
   double value = RooFit::Detail::MathFuncs::additiveInterpolate(coefList, nCoeffs, funcIntList, nFuncs, 
                                                                 pdfFloorVal, quadraticRegion, quadraticAlgo);
   double normVal(1);
   double result = 0;
   if(normVal>0.) result = value / normVal;
   return result > 0. ? result : integralFloorVal;
}

} // namespace MathFuncs
} // namespace Detail
} // namespace RooFit

#endif
