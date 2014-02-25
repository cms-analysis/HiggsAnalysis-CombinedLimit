#include "../interface/SimpleGaussianConstraint.h"

#include <string>
#include <stdexcept>

void SimpleGaussianConstraint::init() {
    if (!sigma.arg().isConstant()) {
        throw std::invalid_argument(std::string("SimpleGaussianConstraint created with non-constant sigma: ") + GetName());
    }
    Double_t sig = sigma;
    scale_ = -0.5/(sig*sig);
}

ClassImp(SimpleGaussianConstraint)
