#ifndef HiggsAnalysis_CombinedLimit_ProfiledLikelihoodRatioTestStatExt_h
#define HiggsAnalysis_CombinedLimit_ProfiledLikelihoodRatioTestStatExt_h

#include <memory>
#include <vector>

class RooMinimizer;
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooArgSet.h>
#include <RooStats/TestStatistic.h>
#include "RooSimultaneousOpt.h"
#include "CachingNLL.h"

namespace nllutils {
    bool robustMinimize(RooAbsReal &nll, RooMinimizer &minimizer, int verbosity=0, bool zeroPoint=false);
}

class ProfiledLikelihoodRatioTestStatOpt : public RooStats::TestStatistic {
    public:
        ProfiledLikelihoodRatioTestStatOpt(const RooArgSet &obs, RooAbsPdf &pdfNull, RooAbsPdf &pdfAlt, 
                const RooArgSet *nuisances, const RooArgSet & paramsNull = RooArgSet(), const RooArgSet & paramsAlt = RooArgSet(),
                int verbosity=0) ;
 
        Double_t Evaluate(RooAbsData& data, RooArgSet& nullPOI) override ;

        const TString GetVarName() const override {
            return TString::Format("-log(%s/%s)", pdfNull_->GetName(), pdfAlt_->GetName()); 
        }

        // Verbosity (default: 0)
        void setPrintLevel(Int_t level) { verbosity_ = level; }

    private:
        RooAbsPdf *pdfNull_, *pdfAlt_;
        RooArgSet snapNull_, snapAlt_; 
        RooArgSet nuisances_; 
        std::unique_ptr<RooArgSet> paramsNull_, paramsAlt_;
        std::unique_ptr<RooAbsReal> nllNull_, nllAlt_;
        Int_t verbosity_;

        // create NLL. if returns true, it can be kept, if false it should be deleted at the end of Evaluate
        bool createNLL(RooAbsPdf &pdf, RooAbsData &data, std::unique_ptr<RooAbsReal> &nll) ;

        double minNLL(std::unique_ptr<RooAbsReal> &nll) ;
}; // TestSimpleStatistics


class ProfiledLikelihoodTestStatOpt : public RooStats::TestStatistic {
    public:
        enum OneSidedness { twoSidedDef = 0, oneSidedDef = 1, signFlipDef = 2 };

        ProfiledLikelihoodTestStatOpt(const RooArgSet & observables,
                RooAbsPdf &pdf, 
                const RooArgSet *nuisances, 
                const RooArgSet & params, const RooArgSet & poi, const RooArgList &gobsParams, const RooArgList &gobs, int verbosity=0, OneSidedness oneSided = oneSidedDef) ; 

        Double_t Evaluate(RooAbsData& data, RooArgSet& nullPOI) override ;
        virtual std::vector<Double_t> Evaluate(RooAbsData& data, RooArgSet& nullPOI, const std::vector<Double_t> &rVals) ;

        const TString GetVarName() const override { return "- log (#lambda)"; }

        // Verbosity (default: 0)
        void setPrintLevel(Int_t level) { verbosity_ = level; }

        void SetOneSided(OneSidedness oneSided) { oneSided_ = oneSided; }
    private:

        RooAbsPdf *pdf_;
        RooArgSet snap_, poi_; // snapshot of parameters, and of the subset which are POI
        std::unique_ptr<RooArgSet>  params_;
        RooArgSet                 nuisances_; // subset of params which are nuisances (not a snapshot)
        RooArgSet                 poiParams_; // subset of params which are POI (not a snapshot)
        std::unique_ptr<RooAbsReal> nll_;
        RooArgList gobsParams_, gobs_;
        Int_t verbosity_;
        OneSidedness oneSided_;

        // create NLL. if returns true, it can be kept, if false it should be deleted at the end of Evaluate
        bool createNLL(RooAbsPdf &pdf, RooAbsData &data) ;
        double minNLL(bool constrained, RooRealVar *r=0) ;
}; // TestSimpleStatistics


#endif
