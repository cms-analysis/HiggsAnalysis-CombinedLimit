#include "HiggsAnalysis/CombinedLimit/interface/RooSimultaneousOpt.h"
#include "HiggsAnalysis/CombinedLimit/interface/CachingNLL.h"
#include <RooCmdConfig.h>

RooAbsReal* 
RooSimultaneousOpt::createNLL(RooAbsData& data, const RooLinkedList& cmdList) 
{
    RooCmdConfig pc(Form("RooSimultaneousOpt::createNLL(%s)",GetName())) ;
    pc.defineSet("cPars","Constrain",0,0);
    RooArgSet *cPars = pc.getSet("cPars");
    return new cacheutils::CachingSimNLL(this, &data, cPars);
}

RooSimultaneousOpt::~RooSimultaneousOpt()
{
}

void
RooSimultaneousOpt::addExtraConstraint(const RooAbsPdf &pdf)
{
    if (!_extraConstraints.contains(pdf)) _extraConstraints.add(pdf);
}

void
RooSimultaneousOpt::addExtraConstraints(const RooAbsCollection &pdfs)
{
    if (_extraConstraints.getSize() == 0) {
        // fast (hopefully)
        _extraConstraints.add(pdfs);
    } else {
        // slow
        RooLinkedListIter iter = pdfs.iterator();
        for (RooAbsArg *a = (RooAbsArg *) iter.Next(); a != 0; a = (RooAbsArg *) iter.Next()) {
            if (!_extraConstraints.contains(*a)) _extraConstraints.add(*a);
        }
    }
}
