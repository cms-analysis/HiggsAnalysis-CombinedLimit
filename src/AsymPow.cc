#include "../interface/AsymPow.h"
#include "../interface/logKappa.h"

#include <cmath>
#include <cassert>
#include <cstdio>

AsymPow::AsymPow(const char *name, const char *title, RooAbsReal &kappaLow, RooAbsReal &kappaHigh, RooAbsReal &theta) :
        RooAbsReal(name,title),
        kappaLow_("kappaLow","Base for theta < 0", this, kappaLow),
        kappaHigh_("kappaHigh","Base for theta > 0", this, kappaHigh),
        theta_("theta", "Exponent (unit gaussian)", this, theta)
        { }

AsymPow::AsymPow(const AsymPow &other, const char *newname) :
    RooAbsReal(other, newname),
    kappaLow_("kappaLow",this,other.kappaLow_),
    kappaHigh_("kappaHigh",this,other.kappaHigh_),
    theta_("theta",this,other.theta_)
    { }

AsymPow::~AsymPow() {}

TObject *AsymPow::clone(const char *newname) const
{
    return new AsymPow(*this,newname);
}

Double_t AsymPow::evaluate() const {
    Double_t x = theta_;
    return exp(logKappaForX(x) * x);
}

Double_t AsymPow::logKappaForX(Double_t x) const {
    const double logKhi =  log(kappaHigh_);
    const double logKlo = -log(kappaLow_);
#if 0
    // old version with discontinuous derivatives
    return (x >= 0 ? logKhi : logKlo);
#else
    return logKappa(x, logKhi, logKlo);
#endif
}

ClassImp(AsymPow)
