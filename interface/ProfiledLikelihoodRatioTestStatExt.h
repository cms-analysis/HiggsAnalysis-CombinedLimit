#ifndef HiggsAnalysis_CombinedLimit_ProfiledLikelihoodRatioTestStatExt_h
#define HiggsAnalysis_CombinedLimit_ProfiledLikelihoodRatioTestStatExt_h

#include <memory>
#include <vector>

class RooMinimizer;
class RooRealVar;
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooArgSet.h>
#include <RooStats/TestStatistic.h>

namespace nllutils {
    bool robustMinimize(RooAbsReal &nll, RooMinimizer &minimizer, int verbosity=0, bool zeroPoint=false);

    /// Create the NLL of `pdf` for `data` into `nll`, constraining `nuisances`.
    /// For a RooSimultaneousOpt, an existing `nll` is reused by just swapping
    /// the dataset. Returns true if `nll` can be kept for later evaluations,
    /// false if the caller should drop it once done with this dataset.
    bool createOrReuseNLL(RooAbsPdf &pdf, RooAbsData &data, const RooArgSet &nuisances, std::unique_ptr<RooAbsReal> &nll);

    /// Minimize `nll` with a CascadeMinimizer and return the value at the
    /// minimum. `poi` is only used in the unconstrained mode.
    double minimizeNLL(RooAbsReal &nll, bool constrained, RooRealVar *poi, int verbosity);
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

        double minNLL(RooAbsReal &nll) const ;
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

        /// Common start of both Evaluate() overloads: save the current
        /// parameter values into `initialState`, load the snapshot, fix the
        /// globally constrained nuisances to their global observables and set
        /// up `nll_` for `data`. Returns the POI to profile, and whether `nll_`
        /// can be kept after this evaluation.
        RooRealVar *setupEvaluation(RooAbsData &data, RooArgSet &initialState, bool &canKeepNLL) ;
        double minNLL(bool constrained, RooRealVar *r) ;
}; // TestSimpleStatistics


#endif
