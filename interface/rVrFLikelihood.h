#ifndef HiggsAnalysis_CombinedLimit_rVrFLikelihood_h
#define HiggsAnalysis_CombinedLimit_rVrFLikelihood_h

#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <TH2.h>
#include <vector>

class rVrFLikelihood : public RooAbsReal {

    public:
        rVrFLikelihood() {}
        rVrFLikelihood(const char *name, const char *title) : RooAbsReal{name, title} {}
        rVrFLikelihood(const rVrFLikelihood& other, const char* name=nullptr);

        void addChannel(const TH2* chi2, RooAbsReal &muV, RooAbsReal &muF);

        TObject* clone(const char* newname) const override {
            return new rVrFLikelihood(*this,newname);
        }

        // Must be public to get dictionaries to compile properly
        struct Channel {
            Channel() {}
            Channel(rVrFLikelihood *parent, const TH2* chi2_, RooAbsReal &muV_, RooAbsReal &muF_) :
                chi2(chi2_),
                muV("muV","signal strength modifier for qqH,VH",  parent, muV_),
                muF("muF","signal strength modifier for ggH,ttH", parent, muF_) {}
            const TH2* chi2;
            RooRealProxy muV, muF;
        };

    protected:
        Double_t evaluate() const override;

    private:
        std::vector<std::unique_ptr<Channel>> channels_;

        ClassDefOverride(rVrFLikelihood,2) // Asymmetric power
};

#endif
