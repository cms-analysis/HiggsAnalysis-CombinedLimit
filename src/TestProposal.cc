#include "HiggsAnalysis/CombinedLimit/interface/TestProposal.h"
#include <RooArgSet.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <TIterator.h>
#include <RooRandom.h>
#include <RooStats/RooStatsUtils.h>

TestProposal::TestProposal(double divisor, const RooRealVar *alwaysStepMe) : 
    RooStats::ProposalFunction(),
    divisor_(1./divisor),
    poiDivisor_(divisor_),
    discreteModelIndicator_(0)
{
    alwaysStepMe_.add(*alwaysStepMe);
}
     
TestProposal::TestProposal(double divisor, const RooArgList &alwaysStepMe) : 
    RooStats::ProposalFunction(),
    divisor_(1./divisor),
    poiDivisor_(divisor_),
    alwaysStepMe_(alwaysStepMe),
    discreteModelIndicator_(0)
{
    if (alwaysStepMe.getSize() > 1) poiDivisor_ /= sqrt(double(alwaysStepMe.getSize()));
}
 

// Populate xPrime with a new proposed point
void TestProposal::Propose(RooArgSet& xPrime, RooArgSet& x )
{
   RooStats::SetParameters(&x, &xPrime);
   RooLinkedListIter it(xPrime.iterator());
   RooRealVar* var;
   int n = xPrime.getSize(), j = floor(RooRandom::uniform()*n);
   const RooRealVar *indicatorPrime = discreteModelIndicator_ ? (RooRealVar*)xPrime.find(discreteModelIndicator_->GetName()) : 0;
   for (int i = 0; (var = (RooRealVar*)it.Next()) != NULL; ++i) {
      if (i == j) {
        if (alwaysStepMe_.contains(*var)) break; // don't step twice
        if (discreteModelIndicator_ != 0) {
            if (var == indicatorPrime) {
                int k = floor(RooRandom::uniform()*discreteModelPoints_.size());
                xPrime.assignValueOnly(discreteModelPoints_[k]); 
                break;
            } else if (discreteModelVars_.contains(*var)) {
                break;
            } 
        } 
        double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
        val += RooRandom::gaussian() * len * divisor_;
        while (val > max) val -= len;
        while (val < min) val += len;
        var->setVal(val);
        break;
      }
   }
   it = alwaysStepMe_.iterator();
   for (RooRealVar *poi = (RooRealVar*)it.Next(); poi != NULL; poi = (RooRealVar*)it.Next()) {
        RooRealVar *var = (RooRealVar*) xPrime.find(poi->GetName());
        if (var == 0) {
            std::cout << "ERROR: missing POI " << poi->GetName() << " in xPrime" << std::endl;
            xPrime.Print("V");
            throw std::logic_error("Missing POI in ArgSet");
        }
        if (discreteModelIndicator_ != 0) {
            if (var == indicatorPrime) {
                int k = floor(RooRandom::uniform()*discreteModelPoints_.size());
                xPrime.assignValueOnly(discreteModelPoints_[k]); 
                continue;
            } else if (discreteModelVars_.contains(*var)) {
                continue;
            } 
        } 
        double val = var->getVal(), max = var->getMax(), min = var->getMin(), len = max - min;
        val += RooRandom::gaussian() * len * poiDivisor_;
        while (val > max) val -= len;
        while (val < min) val += len;
        var->setVal(val);
   }
}

Bool_t TestProposal::IsSymmetric(RooArgSet& x1, RooArgSet& x2) {
   return true;
}

// Return the probability of proposing the point x1 given the starting
// point x2
Double_t TestProposal::GetProposalDensity(RooArgSet& x1,
                                          RooArgSet& x2)
{
   return 1.0; // should not be needed
}

ClassImp(TestProposal)
