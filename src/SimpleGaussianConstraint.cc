#include "../interface/SimpleGaussianConstraint.h"

#include <string>
#include <memory>
#include <stdexcept>

void SimpleGaussianConstraint::init() {
    if (!sigma.arg().InheritsFrom("RooRealVar") && !sigma.arg().InheritsFrom("RooConstVar")) {
        throw std::invalid_argument(std::string("SimpleGaussianConstraint created with non-RooRealVar non-RooConstVar sigma: ") + GetName());
    }
    if (!sigma.arg().isConstant()) {
        throw std::invalid_argument(std::string("SimpleGaussianConstraint created with non-constant sigma: ") + GetName());
    }
    Double_t sig = sigma;
    scale_ = -0.5/(sig*sig);
}

RooGaussian * SimpleGaussianConstraint::make(RooGaussian &c) {
    if (typeid(c) == typeid(SimpleGaussianConstraint)) {
        return &c;
    }
    try {
        std::unique_ptr<SimpleGaussianConstraint> opt(new SimpleGaussianConstraint(c));
        opt->SetNameTitle(c.GetName(),c.GetTitle());
        return opt.release();
    } catch (std::invalid_argument &ex) {
        return &c;
    }
}

ClassImp(SimpleGaussianConstraint)

#include <RooFitHS3/JSONIO.h>

namespace {

  auto RooFitHS3_wsexportkeys = R"({
    "SimpleGaussianConstraint": {
        "type": "gaussian_dist",
        "proxies": {
            "x": "x",
            "mean": "mean",
            "sigma": "sigma"
        }
    }
})";

auto RooFitHS3_wsfactoryexpression = R"({
    "gaussian_dist": {
        "class":"SimpleGaussianConstraint",
        "arguments":[
            "x",
            "mean",
            "sigma"
        ]
    }
})";
} // namespace

#include "static_execute.h"

STATIC_EXECUTE([]() {
  std::stringstream exportkeys;
  std::stringstream factoryexpressions;
  SimpleGaussianConstraint::Class();
  exportkeys << ::RooFitHS3_wsexportkeys;
  RooFit::JSONIO::loadExportKeys(exportkeys);
  factoryexpressions << ::RooFitHS3_wsfactoryexpression;
  RooFit::JSONIO::loadFactoryExpressions(factoryexpressions);
  return 0;
 });

