#include "HiggsAnalysis/CombinedLimit/interface/SimpleConstraintGroup.h"
#include "HiggsAnalysis/CombinedLimit/interface/Accumulators.h"

SimpleConstraintGroup::SimpleConstraintGroup() :
    RooAbsReal("unnamedGroup",""),
    _deps("deps", "deps", this)
{
}


SimpleConstraintGroup::SimpleConstraintGroup(const char *name, const char *title) :
    RooAbsReal(name,title),
    _deps("deps", "deps", this)
{
}

SimpleConstraintGroup::SimpleConstraintGroup(const SimpleConstraintGroup & other, const char *newname) :
    RooAbsReal(other, newname),
    _deps("deps", this, other._deps),
    _gaus(other._gaus),
    _pois(other._pois),
    _gaus0(other._gaus0),
    _pois0(other._pois0)
{
}

void SimpleConstraintGroup::add(const SimpleGaussianConstraint * gaus) {
    _deps.add(gaus->getX());
    _gaus.push_back(gaus);
    _gaus0.push_back(0);
}

void SimpleConstraintGroup::add(const SimplePoissonConstraint * pois) {
    std::cout << "Adding one poisson" << std::endl;
    _deps.add(pois->getMean());
    _pois.push_back(pois);
    _pois0.push_back(0);
}

void SimpleConstraintGroup::setZeroPoint() {
    auto ig0 = _gaus0.begin();
    for (auto ig = _gaus.begin(), eg = _gaus.end(); ig != eg; ++ig, ++ig0) {
        *ig0 = -(*ig)->getLogValFast();
    }
    auto ip0 = _pois0.begin();
    for (auto ip = _pois.begin(), ep = _pois.end(); ip != ep; ++ip, ++ip0) {
        *ip0 = -(*ip)->getLogValFast();
    }
    setValueDirty();
}

void SimpleConstraintGroup::clearZeroPoint() {
    std::fill(_gaus0.begin(), _gaus0.end(), 0.);
    std::fill(_pois0.begin(), _pois0.end(), 0.);
    setValueDirty();
}

Double_t SimpleConstraintGroup::evaluate() const {
    DefaultAccumulator<double> ret2 = 0;
    auto ig0 = _gaus0.begin();
    for (auto ig = _gaus.begin(), eg = _gaus.end(); ig != eg; ++ig, ++ig0) {
        ret2 += ((*ig)->getLogValFast() + *ig0);
    }
    auto ip0 = _pois0.begin();
    for (auto ip = _pois.begin(), ep = _pois.end(); ip != ep; ++ip, ++ip0) {
        ret2 += ((*ip)->getLogValFast() + *ip0);
    }
    return ret2.sum();
}



