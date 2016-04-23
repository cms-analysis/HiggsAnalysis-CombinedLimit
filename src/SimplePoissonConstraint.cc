#include "HiggsAnalysis/CombinedLimit/interface/SimplePoissonConstraint.h"

#include <string>
#include <memory>
#include <stdexcept>
#include <TMath.h>

void SimplePoissonConstraint::init() {
    if (!x.arg().InheritsFrom("RooRealVar") && !x.arg().InheritsFrom("RooConstVar")) {
        throw std::invalid_argument(std::string("SimplePoissonConstraint created with non-RooRealVar non-RooConstVar sigma: ") + GetName());
    }
    if (!x.arg().isConstant()) {
        throw std::invalid_argument(std::string("SimplePoissonConstraint created with non-constant x: ") + GetName());
    }
    Double_t observed = x;
    logGamma_ = TMath::LnGamma(observed+1.);
}

RooPoisson * SimplePoissonConstraint::make(RooPoisson &c) {
    if (typeid(c) == typeid(SimplePoissonConstraint)) {
        return &c;
    }
    try {
        std::unique_ptr<SimplePoissonConstraint> opt(new SimplePoissonConstraint(c));
        opt->SetNameTitle(c.GetName(),c.GetTitle());
        return opt.release();
    } catch (std::invalid_argument &ex) {
        return &c;
    }
}

ClassImp(SimplePoissonConstraint)
