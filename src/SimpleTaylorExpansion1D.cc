#include "HiggsAnalysis/CombinedLimit/interface/SimpleTaylorExpansion1D.h"

#include <cmath>
#include <cassert>
#include <cstdio>
#include <algorithm>

SimpleTaylorExpansion1D::SimpleTaylorExpansion1D(const char *name, const char *title, RooAbsReal & func, RooRealVar & xvar, double dx, int order) :
    RooAbsReal(name,title),
    x_("x", "Independent variable", this, xvar),
    x0_(xvar.getVal())
{ 
    std::fill_n(ci_, MaxOrder+1, 0.0);
    ci_[0] = func.getVal();
    if (order > 4) {
        std::cerr << "Only order 0..4 are implemented for the moment, out of lazyness" << std::endl;
        assert(false);
    }
    if (order > 0) {
        xvar.setVal(x0_ - dx);
        double ym = func.getVal();
        xvar.setVal(x0_ + dx);
        double yp = func.getVal();

        if (order <= 2) {
            ci_[1] = 0.5*(yp - ym)/dx;
            if (order == 2) {
                ci_[2] = 0.5*((yp - ci_[0]) + (ym - ci_[0]))/(dx*dx);
            }
        } else {
            xvar.setVal(x0_ - 2*dx);
            double y2m = func.getVal();
            xvar.setVal(x0_ + 2*dx);
            double y2p = func.getVal();
            
            double x11 = 0.5*(yp - ym), x12 = 0.5*(y2p - y2m);
            double x21 = ((yp  - ci_[0]) + ( ym - ci_[0]));
            double x22 = ((y2p - ci_[0]) + (y2m - ci_[0]));
            ci_[1] = (4*x11/3-x12/6)/dx;
            ci_[2] = (4*x21/3-x22/12)/(2*dx*dx);
            ci_[3] = (x12-2*x11)/(6*dx*dx*dx);
            if (order == 4) {
                ci_[4] = (x22-4*x21)/(24*dx*dx*dx*dx);
            }
        }
        xvar.setVal(x0_);
    }
}

SimpleTaylorExpansion1D::SimpleTaylorExpansion1D(const SimpleTaylorExpansion1D &other, const char *newname) :
    RooAbsReal(other, newname),
    x_("x",this,other.x_),
    x0_(other.x0_)
{ 
    std::copy_n(other.ci_, MaxOrder+1, ci_); 
}

SimpleTaylorExpansion1D::~SimpleTaylorExpansion1D() {}

TObject *SimpleTaylorExpansion1D::clone(const char *newname) const 
{
    return new SimpleTaylorExpansion1D(*this,newname);
}

Double_t SimpleTaylorExpansion1D::evaluate() const {
    Double_t dx = x_ - x0_; 
    Double_t ret = 0, dxn = 1;
    for (unsigned int i = 0; i <= MaxOrder; ++i) {
        ret += ci_[i] * dxn;
        dxn *= dx;
    }
    return ret;
}

ClassImp(SimpleTaylorExpansion1D)
