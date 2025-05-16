#include "../interface/AsymPow.h"

#include "../interface/CombineMathFuncs.h"

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

Double_t AsymPow::evaluate() const {
   return RooFit::Detail::MathFuncs::asymPow(theta_, kappaLow_, kappaHigh_);
}

ClassImp(AsymPow)
