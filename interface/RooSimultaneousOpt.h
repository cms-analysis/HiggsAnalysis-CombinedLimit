#ifndef ROO_SIM_OPT_PDF
#define ROO_SIM_OPT_PDF

/** Trick RooSimultaneous */

#include "RooSimultaneous.h"
#include "RooListProxy.h"

class RooSimultaneousOpt : public RooSimultaneous {
public:
	RooSimultaneousOpt(){};
  RooSimultaneousOpt(const char* name, const char* title, RooAbsCategoryLValue& indexCat) :
        RooSimultaneous(name, title, indexCat),
        _extraConstraints("extraConstraints", "Extra constraint terms independent on the category",this),
        _channelMasks("channelMasks", "Terms or expressions that mask particular channels", this)  {}
  RooSimultaneousOpt(const RooSimultaneous &any, const char* name=0) :
        RooSimultaneous(any, name),
        _extraConstraints("extraConstraints", "Extra constraint terms independent on the category",this),
        _channelMasks("channelMasks", "Terms or expressions that mask particular channels", this)
        {
            const RooSimultaneousOpt * opt = dynamic_cast<const RooSimultaneousOpt *>(&any);
            if (opt) _extraConstraints.add(opt->extraConstraints(), true);
        }
  RooSimultaneousOpt(const RooSimultaneousOpt &other, const char* name=0) :
        RooSimultaneous(other, name),
        _extraConstraints("extraConstraints", this, other._extraConstraints),
        _channelMasks("channelMasks", this, other._channelMasks)  {}

  RooSimultaneousOpt *clone(const char* name=0) const override { return new RooSimultaneousOpt(*this, name); }

  ~RooSimultaneousOpt() override ;

#if ROOT_VERSION_CODE < ROOT_VERSION(6,30,0)
  RooAbsReal* createNLL(RooAbsData& data, const RooLinkedList& cmdList) override;
#else
  std::unique_ptr<RooAbsReal> createNLLImpl(RooAbsData& data, const RooLinkedList& cmdList) override;
#endif

  const RooArgList & extraConstraints() const { return _extraConstraints; }
  const RooArgList & channelMasks() const { return _channelMasks; }

  void addExtraConstraint(const RooAbsPdf &pdf) ;
  void addExtraConstraints(const RooAbsCollection &pdfs) ;
  void addChannelMasks(const RooArgList &args);
private:
  ClassDefOverride(RooSimultaneousOpt,3) // Variant of RooSimultaneous that can put together binned and unbinned stuff 

  RooListProxy _extraConstraints;  //  List of extra constraint terms
  RooListProxy _channelMasks;

};

#endif
