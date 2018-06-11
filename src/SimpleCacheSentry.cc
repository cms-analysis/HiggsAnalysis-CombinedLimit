#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"

SimpleCacheSentry::SimpleCacheSentry() :  _deps("deps","deps",this)   {}

SimpleCacheSentry::SimpleCacheSentry(const char *name, const char *title)
    : RooAbsArg(name, title), _deps("deps", "deps", this) {}

SimpleCacheSentry::SimpleCacheSentry(const RooRealVar &var) :
    _deps("deps","deps",this)   
{
    addVar(var);
}

SimpleCacheSentry::SimpleCacheSentry(const RooAbsCollection &vars) :
    _deps("deps","deps",this)   
{
    addVars(vars);
}


SimpleCacheSentry::SimpleCacheSentry(const RooAbsArg &func, const RooArgSet *obs) :
    _deps("deps","deps",this)   
{
    addFunc(func,obs);
}

SimpleCacheSentry::SimpleCacheSentry(const SimpleCacheSentry &other, const char *newname) :
    _deps("deps",this,other._deps)   
{
}

void SimpleCacheSentry::addVars(const RooAbsCollection &vars) 
{
    TIterator *iter = vars.createIterator();
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
        if (_deps.containsInstance(*a)) continue;
        // RooRealVars can return true to isDerived() if the ranges or binning depend on
        // other parameters, so always add RooRealVars to the list
        if (dynamic_cast<RooRealVar*>(a)) {
          _deps.add(*a);
        } else if (a->isDerived()) {
          addFunc(*a);
        } else {
          _deps.add(*a);
        }
    }
    delete iter;
}

void SimpleCacheSentry::addFunc(const RooAbsArg &func, const RooArgSet *obs) 
{
    RooArgSet *deps = func.getParameters(obs,false);
    addVars(*deps);
    delete deps;
}

Bool_t SimpleCacheSentry::isIdentical(const RooAbsArg& other, 
            Bool_t /*assumeSameType*/) {
  bool ret = kFALSE;
  SimpleCacheSentry const& otherSentry = dynamic_cast<SimpleCacheSentry const&>(other);
  RooAbsCollection * common = _deps.selectCommon(otherSentry._deps);
  if ( (common->getSize() == otherSentry._deps.getSize()) && 
       (common->getSize() == _deps.getSize()) )
    ret = kTRUE;
  delete common;
  return ret;
}

