#ifndef HiggsAnalysis_CombinedLimit_BestFitSigmaTestStat
#define HiggsAnalysis_CombinedLimit_BestFitSigmaTestStat

#include <memory>

#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooArgSet.h>
#include <RooStats/TestStatistic.h>

class BestFitSigmaTestStat : public RooStats::TestStatistic {
    public:
        BestFitSigmaTestStat(const RooArgSet & observables,
                RooAbsPdf &pdf, 
                const RooArgSet *nuisances, 
                const RooArgSet & params, int verbosity=0) ; 

        Double_t Evaluate(RooAbsData& data, RooArgSet& nullPOI) override ;

        const TString GetVarName() const override { return "mu-hat`"; }

        // Verbosity (default: 0)
        void setPrintLevel(Int_t level) { verbosity_ = level; }
    private:

        RooAbsPdf *pdf_;
        RooArgSet snap_, poi_, nuisances_; 
        std::unique_ptr<RooArgSet> params_;
        std::unique_ptr<RooAbsReal> nll_;
        Int_t verbosity_;
}; // TestSimpleStatistics


#endif
