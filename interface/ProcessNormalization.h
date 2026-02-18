#ifndef HiggsAnalysis_CombinedLimit_ProcessNormalization_h
#define HiggsAnalysis_CombinedLimit_ProcessNormalization_h

#include <RooAbsReal.h>
#include "RooListProxy.h"

//_________________________________________________
/*
BEGIN_HTML
<p>
ProcessNormalization is helper class for implementing process normalizations
</p>
END_HTML
*/
//
class ProcessNormalization : public RooAbsReal {
   public:
      ProcessNormalization() : nominalValue_(1) {}
      ProcessNormalization(const char *name, const char *title, double nominal=1) ;
      ProcessNormalization(const char *name, const char *title, RooAbsReal &nominal) ;
      ProcessNormalization(const ProcessNormalization &other, const char *newname = 0) ;

      TObject * clone(const char *newname) const override { return new ProcessNormalization(*this, newname); }

      void setNominalValue(double nominal) { nominalValue_ = nominal; }
      void addLogNormal(double kappa, RooAbsReal &theta) ;
      void addAsymmLogNormal(double kappaLo, double kappaHi, RooAbsReal &theta) ;
      void addOtherFactor(RooAbsReal &factor) ;
      void dump() const ;

      double nominalValue() const { return nominalValue_; }
      std::vector<double> const &logKappa() const { return logKappa_; }
      RooArgList const &thetaList() const { return thetaList_; }
      std::vector<std::pair<double, double>> const &logAsymmKappa() const { return logAsymmKappa_; }
      RooArgList const &asymmThetaList() const { return asymmThetaList_; }
      RooArgList const &otherFactorList() const { return otherFactorList_; }

    protected:
        Double_t evaluate() const override;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
        void doEval(RooFit::EvalContext &) const override;
#endif

    private:
        void fillAsymmKappaVecs() const;

        // ---- PERSISTENT ----
        double nominalValue_;                         
        std::vector<double> logKappa_; // Logarithm of symmetric kappas
        RooListProxy thetaList_;        // List of nuisances for symmetric kappas
        std::vector<std::pair<double,double> > logAsymmKappa_; // Logarithm of asymmetric kappas (low, high)
        RooListProxy asymmThetaList_;                           // List of nuisances for asymmetric kappas
        RooListProxy otherFactorList_;     // Other multiplicative terms 
        mutable std::vector<double> thetaListVec_; //! Don't serialize me
        mutable std::vector<double> asymmThetaListVec_; //! Don't serialize me
        mutable std::vector<double> otherFactorListVec_; //! Don't serialize me
        mutable std::vector<double> logAsymmKappaLow_; //! Don't serialize me
        mutable std::vector<double> logAsymmKappaHigh_; //! Don't serialize me

  ClassDefOverride(ProcessNormalization,1) // Process normalization interpolator 
};

#endif
