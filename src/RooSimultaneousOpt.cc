#include "../interface/RooSimultaneousOpt.h"
#include "../interface/CachingNLL.h"
#include "../interface/Combine.h"

#include <RooCmdConfig.h>

#if ROOT_VERSION_CODE < ROOT_VERSION(6,30,0)
RooAbsReal* 
RooSimultaneousOpt::createNLL(RooAbsData& data, const RooLinkedList& cmdList) 
#else
std::unique_ptr<RooAbsReal>
RooSimultaneousOpt::createNLLImpl(RooAbsData& data, const RooLinkedList& cmdList) 
#endif
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 36, 0)
  // Match the printout of https://github.com/root-project/root/blob/a6ec2bd88e82f46f131cece8b3cc6b357c756f02/roofit/roofitcore/src/FitHelpers.cxx#L594
  auto timingScope = std::make_unique<ROOT::Math::Util::TimingScope>(
      [this](std::string const& msg) { oocoutI(this, Fitting) << msg << std::endl; }, "Creation of NLL object took");
#endif
  RooCmdConfig pc(Form("RooSimultaneousOpt::createNLL(%s)", GetName()));
  pc.defineSet("cPars", "Constrain", 0, 0);
  RooArgSet* cPars = pc.getSet("cPars");
  auto nll = std::make_unique<cacheutils::CachingSimNLL>(this, &data, cPars);
  nll->setChannelMasks(this->channelMasks());
#if ROOT_VERSION_CODE < ROOT_VERSION(6,30,0)
    return nll.release();
#else
    return nll;
#endif
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
        for (RooAbsArg *a : pdfs) {
            if (!_extraConstraints.contains(*a)) _extraConstraints.add(*a);
        }
    }
}

void
RooSimultaneousOpt::addChannelMasks(const RooArgList &args)
{
    if (args.getSize() != indexCat().numTypes()) {
        std::cerr << "RooSimultaneousOpt: number of channel masks must equal number of channels\n";
        return;
    }
    _channelMasks.removeAll();
    _channelMasks.add(args);
}
