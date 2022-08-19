#include "HiggsAnalysis/CombinedLimit/interface/rVrFLikelihood.h"
#include <cstdio>

void
rVrFLikelihood::addChannel(const TH2* chi2, RooAbsReal &muV, RooAbsReal &muF)
{
    channels_.emplace_back(std::make_unique<Channel>(this,chi2,muV,muF));
}

rVrFLikelihood::rVrFLikelihood(const rVrFLikelihood& other, const char* name)
  : RooAbsReal{other, name}
{
    // never understood if RooFit actually cares of const-correctness or not.
    for (auto const& channel : other.channels_) {
        addChannel(channel->chi2, const_cast<RooAbsReal &>(channel->muV.arg()), const_cast<RooAbsReal &>(channel->muF.arg()));
    }
}

Double_t rVrFLikelihood::evaluate() const {
    double ret = 0;
    for (auto const& channel : channels_) {
        double x = channel->muF;
        double y = channel->muV;
        if (!(x > channel->chi2->GetXaxis()->GetXmin())) return 9999;
        if (!(x < channel->chi2->GetXaxis()->GetXmax())) return 9999;
        if (!(y > channel->chi2->GetYaxis()->GetXmin())) return 9999;
        if (!(y < channel->chi2->GetYaxis()->GetXmax())) return 9999;
        /*printf("looked up %s at (%g,%g), x [%g,%g], y [%g,%g]\n",
                    channel->chi2->GetName(), x, y,
                    channel->chi2->GetXaxis()->GetXmin(), channel->chi2->GetXaxis()->GetXmax(),
                    channel->chi2->GetYaxis()->GetXmin(), channel->chi2->GetYaxis()->GetXmax());*/
        ret += (const_cast<TH2 *>(channel->chi2))->Interpolate(x,y); // WTF Interpolate is not const??
    }
    return 0.5*ret; // NLL
}

ClassImp(rVrFLikelihood)
