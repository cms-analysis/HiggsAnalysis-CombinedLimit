#ifndef ROO_SIM_OPT_PDF
#define ROO_SIM_OPT_PDF

/** Trick RooSimultaneous */

#include "RooSimultaneous.h"
#include "RooListProxy.h"

class RooSimultaneousOpt : public RooSimultaneous {
public:
	RooSimultaneousOpt(){};
  RooSimultaneousOpt(const char* name, const char* title, RooAbsCategoryLValue& indexCat) :
        RooSimultaneous(name, title, indexCat), _extraConstraints("extraConstraints", "Extra constraint terms independent on the category",this)  {}
  RooSimultaneousOpt(const RooSimultaneous &any, const char* name=0) :
        RooSimultaneous(any, name), _extraConstraints("extraConstraints", "Extra constraint terms independent on the category",this)  
        {
            const RooSimultaneousOpt * opt = dynamic_cast<const RooSimultaneousOpt *>(&any);
            if (opt) _extraConstraints.add(opt->extraConstraints(), true);
        }
  RooSimultaneousOpt(const RooSimultaneousOpt &other, const char* name=0) :
        RooSimultaneous(other, name), _extraConstraints("extraConstraints", this, other._extraConstraints)  {}

  virtual RooSimultaneousOpt *clone(const char* name=0) const { return new RooSimultaneousOpt(*this, name); }

  virtual ~RooSimultaneousOpt() ;

  virtual RooAbsReal* createNLL(RooAbsData& data, const RooLinkedList& cmdList) ;

  const RooArgList & extraConstraints() const { return _extraConstraints; }

  void addExtraConstraint(const RooAbsPdf &pdf) ;
  void addExtraConstraints(const RooAbsCollection &pdfs) ;
private:
  ClassDef(RooSimultaneousOpt,2) // Variant of RooSimultaneous that can put together binned and unbinned stuff 

  RooListProxy _extraConstraints;  //  List of extra constraint terms

};

#endif
