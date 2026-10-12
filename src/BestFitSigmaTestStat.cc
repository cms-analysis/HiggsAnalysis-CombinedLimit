#include "../interface/BestFitSigmaTestStat.h"
#include "../interface/ProfiledLikelihoodRatioTestStatExt.h"
#include "../interface/utils.h"

#include <RooRealVar.h>

#include <iostream>

//============================================================================================================================================
BestFitSigmaTestStat::BestFitSigmaTestStat(
    const RooArgSet & observables,
    RooAbsPdf &pdf, 
    const RooArgSet *nuisances, 
    const RooArgSet & params,
    int verbosity)
:
    pdf_(&pdf),
    verbosity_(verbosity)
{
    params.snapshot(snap_,false);
    ((RooRealVar*)snap_.find(params.first()->GetName()))->setConstant(false);
    poi_.add(snap_);
    if (nuisances) { nuisances_.add(*nuisances); snap_.addClone(*nuisances, /*silent=*/true); }
    params_.reset(pdf_->getParameters(observables));
}


Double_t BestFitSigmaTestStat::Evaluate(RooAbsData& data, RooArgSet& /*nullPOI*/)
{
    // Take snapshot of initial state, to restore it at the end 
    RooArgSet initialState; params_->snapshot(initialState);

   // Initialize parameters
    *params_ = snap_;

    // Initialize signal strength
    RooRealVar *rIn = (RooRealVar *) poi_.first();
    RooRealVar *r   = (RooRealVar *) params_->find(rIn->GetName());
    bool canKeepNLL = nllutils::createOrReuseNLL(*pdf_, data, nuisances_, nll_);

    double initialR = rIn->getVal();

    // Perform unconstrained minimization
    r->setVal(initialR == 0 ? 0.5 : 0.5*initialR); //best guess
    r->setConstant(false);

    std::cout << "Doing a fit for with " << r->GetName() << " in range [ " << r->getMin() << " , " << r->getMax() << "]" << std::endl;
    std::cout << "Starting point is " << r->GetName() << " = " << r->getVal() << std::endl;
    nllutils::minimizeNLL(*nll_, /*constrained=*/false, r, verbosity_);
    double bestFitR = r->getVal();
    std::cout << "Fit result was " << r->GetName() << " = " << r->getVal() << std::endl;

    //Restore initial state, to avoid issues with ToyMCSampler
    *params_ = initialState;

    if (!canKeepNLL) nll_.reset();

    return bestFitR;
}
