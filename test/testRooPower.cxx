#include "../interface/CachingNLL.h"
#include "../interface/CombineUtils.h"
#include "../interface/HGGRooPdfs.h"
#include "../interface/RooSimultaneousOpt.h"

#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>

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
      EXPECT_DOUBLE_EQ(pdf.getVal(), ref.valueAtTwo);
      EXPECT_NEAR(full->getVal(), ref.integral, 1.e-12);
      EXPECT_NEAR(window->getVal(), ref.windowIntegral, 1.e-12);
    }
  }

  TEST(RooPowerTest, CombineLikelihood) {
    RooRealVar mass{/*name=*/"mass", /*title=*/"mass", /*value=*/2., /*minValue=*/1., /*maxValue=*/4.};
    RooRealVar exponent{/*name=*/"exponent", /*title=*/"exponent", /*value=*/-5., /*minValue=*/-6., /*maxValue=*/1.};
    RooPower pdf{/*name=*/"power", /*title=*/"power", mass, exponent};
    RooConstVar yield{/*name=*/"yield", /*title=*/"yield", /*value=*/3.};
    RooAddPdf channelPdf{/*name=*/"channelPdf", /*title=*/"channelPdf", RooArgList{pdf}, RooArgList{yield}};
    RooCategory channel{/*name=*/"channel", /*title=*/"channel"};
    channel.defineType(/*label=*/"power");
    RooSimultaneousOpt model{/*name=*/"model", /*title=*/"model", channel};
    model.addPdf(channelPdf, /*catLabel=*/"power");
    RooArgSet observables{mass, channel};
    RooDataSet data{"data", "data", observables};
    for (double value : {1.5, 2., 3.5}) {
      mass.setVal(value);
      data.add(observables);
    }
    auto nll = combineCreateNLL(model, data);
    auto* cached = dynamic_cast<cacheutils::CachingSimNLL*>(nll.get());
    ASSERT_NE(cached, nullptr);
    cached->clearConstantZeroPoint();  // Removing Combine's offset term
    for (const auto& ref : references) {
      SCOPED_TRACE(ref.exponent);
      exponent.setVal(ref.exponent);
      // The extended term is zero at yield = observed count = 3.
      // -log(product_i(x_i^c / integral)); product_i(x_i) = 10.5.
      const double expected = (3. * std::log(ref.integral)) - (ref.exponent * std::log(10.5));
      EXPECT_NEAR(nll->getVal(), expected, 1.e-12);
    }
  }

}  // namespace
