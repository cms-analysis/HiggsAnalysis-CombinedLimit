#include "../interface/HGGRooPdfs.h"

#include "../interface/CombineMathFuncs.h"

#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"

#include <gtest/gtest.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <memory>
#include <set>

namespace {

  struct Reference {
    double exponent, valueAtTwo, integral, windowIntegral;
  };

  // c, x^c at x=2, integral from 1 to 4, integral from 2 to 3
  const std::array<Reference, 3> references{{
      {.exponent=-5., .valueAtTwo=1. / 32., .integral=255. / 1024., .windowIntegral=65. / 5184.},
      {.exponent=-1., .valueAtTwo=0.5, .integral=std::log(4.), .windowIntegral=std::log(1.5)},
      {.exponent=0., .valueAtTwo=1., .integral=3., .windowIntegral=1.},
  }};

  TEST(RooPowerTest, EvaluationAndIntegrals) {
    RooRealVar mass{/*name=*/"mass", /*title=*/"mass", /*value=*/2., /*minValue=*/1., /*maxValue=*/4.};
    RooRealVar exponent{/*name=*/"exponent", /*title=*/"exponent", /*value=*/-5., /*minValue=*/-6., /*maxValue=*/1.};
    RooPower pdf{/*name=*/"power", /*title=*/"power", mass, exponent};
    RooArgSet observables{mass};
    mass.setRange(/*name=*/"window", /*min=*/2., /*max=*/3.);
    std::unique_ptr<RooAbsReal> full{pdf.createIntegral(observables)};
    std::unique_ptr<RooAbsReal> window{pdf.createIntegral(observables, RooFit::Range(/*rangeName=*/"window"))};
    for (const auto& ref : references) {
      SCOPED_TRACE(ref.exponent);
      exponent.setVal(ref.exponent);
      EXPECT_NEAR(pdf.getVal(observables), ref.valueAtTwo / ref.integral, 1.e-12);
      EXPECT_NEAR(full->getVal(), ref.integral, 1.e-12);
      EXPECT_NEAR(window->getVal(), ref.windowIntegral, 1.e-12);
    }
  }

  TEST(RooPowerTest, IntegralAndDerivativeScan) {
    // Eight points per decade, including the exact limit and its adjacent floating-point values
    std::set<double> points{-1., 0., std::nextafter(-1., -2.), std::nextafter(-1., 0.)};
    for (int i = 0; i <= 120; ++i) {
      const double offset = std::pow(10., -1. - (i / 8.));
      points.insert(-1. - offset);
      points.insert(-1. + offset);
    }
    const std::array<std::array<double, 2>, 2> ranges{{{1., 4.}, {100., 180.}}};
    for (const auto& range : ranges) {
      SCOPED_TRACE(::testing::Message() << "range=[" << range[0] << ", " << range[1] << "]");
      // Evaluate the closed-form function with 100 decimal digits: the quotient
      // loses precision near c=-1 in double arithmetic.
      using Real = boost::multiprecision::number<boost::multiprecision::backends::cpp_dec_float<100>>;
      const Real logLow = log(Real(range[0]));
      const Real logHigh = log(Real(range[1]));
      for (double c : points) {
        SCOPED_TRACE(::testing::Message() << std::setprecision(17) << "c=" << c);
        const Real q = Real(c) + 1;
        Real integral;
        Real derivative;
        if (q == 0) {
          integral = logHigh - logLow;
          derivative = (logHigh * logHigh - logLow * logLow) / 2;
        } else {
          const Real lowPower = exp(q * logLow);
          const Real highPower = exp(q * logHigh);
          integral = (highPower - lowPower) / q;
          derivative = (q * (highPower * logHigh - lowPower * logLow) - (highPower - lowPower)) / (q * q);
        }
        const auto expectedIntegral = static_cast<double>(integral);
        const auto expectedDerivative = static_cast<double>(derivative);
        EXPECT_NEAR(RooFit::Detail::MathFuncs::powerIntegral(range[0], range[1], c),
                    expectedIntegral,
                    1.e-10 * std::abs(expectedIntegral));
        EXPECT_NEAR(RooFit::Detail::MathFuncs::powerIntegralDc(range[0], range[1], c),
                    expectedDerivative,
                    1.e-10 * std::abs(expectedDerivative));
      }
    }
  }

}  // namespace
