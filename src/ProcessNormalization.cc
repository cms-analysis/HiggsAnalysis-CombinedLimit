#include "../interface/ProcessNormalization.h"

#include "../interface/CombineMathFuncs.h"

#include <cmath>
#include <cassert>
#include <cstdio>

ProcessNormalization::ProcessNormalization(const char *name, const char *title, double nominal) :
        RooAbsReal(name,title),
        nominalValue_(nominal),
        thetaList_("thetaList","List of nuisances for symmetric kappas", this), 
        asymmThetaList_("asymmThetaList","List of nuisances for asymmetric kappas", this), 
        otherFactorList_("otherFactorList","Other multiplicative terms", this)
{ 
}

ProcessNormalization::ProcessNormalization(const char *name, const char *title, RooAbsReal &nominal)
  : ProcessNormalization{name, title, 1.0}
{
   otherFactorList_.add(nominal);
}

ProcessNormalization::ProcessNormalization(const ProcessNormalization &other, const char *newname) :
        RooAbsReal(other, newname ? newname : other.GetName()),
        nominalValue_(other.nominalValue_),
        logKappa_(other.logKappa_),
        thetaList_("thetaList", this, other.thetaList_), 
        logAsymmKappa_(other.logAsymmKappa_),
        asymmThetaList_("asymmThetaList", this, other.asymmThetaList_), 
        otherFactorList_("otherFactorList", this, other.otherFactorList_)
{
}

void ProcessNormalization::addLogNormal(double kappa, RooAbsReal &theta) {
    if (kappa != 0.0 && kappa != 1.0) {
        logKappa_.push_back(std::log(kappa));
        thetaList_.add(theta);
    }
}

void ProcessNormalization::addAsymmLogNormal(double kappaLo, double kappaHi, RooAbsReal &theta) {
    if (fabs(kappaLo*kappaHi - 1) < 1e-5) {
        addLogNormal(kappaHi, theta);
    } else {
        logAsymmKappa_.push_back(std::make_pair(std::log(kappaLo), std::log(kappaHi)));
        asymmThetaList_.add(theta);
    }
}

void ProcessNormalization::addOtherFactor(RooAbsReal &factor) {
    otherFactorList_.add(factor);
}

void ProcessNormalization::fillAsymmKappaVecs() const
{
    if (logAsymmKappaLow_.size() != logAsymmKappa_.size()) {
       logAsymmKappaLow_.reserve(logAsymmKappa_.size());
       logAsymmKappaHigh_.reserve(logAsymmKappa_.size());
       for (auto [lo, hi] : logAsymmKappa_) {
          logAsymmKappaLow_.push_back(lo);
          logAsymmKappaHigh_.push_back(hi);
       }
    }
}

Double_t ProcessNormalization::evaluate() const
{
    thetaListVec_.resize(thetaList_.size());
    asymmThetaListVec_.resize(asymmThetaList_.size());
    otherFactorListVec_.resize(otherFactorList_.size());
    for (std::size_t i = 0; i < thetaList_.size(); ++i) {
        thetaListVec_[i] = static_cast<RooAbsReal const&>(thetaList_[i]).getVal();
    }
    for (std::size_t i = 0; i < asymmThetaList_.size(); ++i) {
        asymmThetaListVec_[i] = static_cast<RooAbsReal const&>(asymmThetaList_[i]).getVal();
    }
    for (std::size_t i = 0; i < otherFactorList_.size(); ++i) {
        otherFactorListVec_[i] = static_cast<RooAbsReal const&>(otherFactorList_[i]).getVal();
    }

    fillAsymmKappaVecs();
    return RooFit::Detail::MathFuncs::processNormalization(nominalValue_,
            thetaList_.size(), asymmThetaList_.size(), otherFactorList_.size(),
            thetaListVec_.data(), logKappa_.data(), asymmThetaListVec_.data(),
            logAsymmKappaLow_.data(), logAsymmKappaHigh_.data(),
            otherFactorListVec_.data());
}

void ProcessNormalization::dump() const {
    std::cout << "Dumping ProcessNormalization " << GetName() << " @ " << (void*)this << std::endl;
    std::cout << "\tnominal value: " << nominalValue_ << std::endl;
    std::cout << "\tlog-normals (" << logKappa_.size() << "):"  << std::endl;
    for (unsigned int i = 0; i < logKappa_.size(); ++i) {
        std::cout << "\t\t kappa = " << exp(logKappa_[i]) << ", logKappa = " << logKappa_[i] << 
                     ", theta = " << thetaList_.at(i)->GetName() << " = " << ((RooAbsReal*)thetaList_.at(i))->getVal() << std::endl;
    }
    std::cout << "\tasymm log-normals (" << logAsymmKappa_.size() << "):"  << std::endl;
    for (unsigned int i = 0; i < logAsymmKappa_.size(); ++i) {
        std::cout << "\t\t kappaLo = " << exp(logAsymmKappa_[i].first) << ", logKappaLo = " << logAsymmKappa_[i].first << 
                     ", kappaHi = " << exp(logAsymmKappa_[i].second) << ", logKappaHi = " << logAsymmKappa_[i].second << 
                     ", theta = " << asymmThetaList_.at(i)->GetName() << " = " << ((RooAbsReal*)asymmThetaList_.at(i))->getVal() << std::endl;
    }
    std::cout << "\tother terms (" << otherFactorList_.getSize() << "):"  << std::endl;
    for (int i = 0; i < otherFactorList_.getSize(); ++i) {  
        std::cout << "\t\t term " << otherFactorList_.at(i)->GetName() <<
                     " (class " << otherFactorList_.at(i)->ClassName() << 
                     "), value = " << ((RooAbsReal*)otherFactorList_.at(i))->getVal() << std::endl;
    }
    std::cout << std::endl;
}
ClassImp(ProcessNormalization)
