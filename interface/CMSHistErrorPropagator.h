#ifndef CMSHistErrorPropagator_h
#define CMSHistErrorPropagator_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "HiggsAnalysis/CombinedLimit/interface/FastTemplate.h"
#include "HiggsAnalysis/CombinedLimit/interface/Logging.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"

class CMSHistErrorPropagatorV;


class CMSHistErrorPropagator : public RooAbsReal {
public:
  CMSHistErrorPropagator();

  // Possibility to add this on a per-bin basis
  // Different types of bin errors to handle
  // [0] = Nothing
  // [1] = Barlow-Beeston
  // [2] = Mixed
  CMSHistErrorPropagator(const char* name, const char* title, RooRealVar& x,
                         RooArgList const& funcs, RooArgList const& coeffs);

  CMSHistErrorPropagator(CMSHistErrorPropagator const& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new CMSHistErrorPropagator(*this, newname);
  }

  virtual ~CMSHistErrorPropagator() {;}


  // std::vector<double> getErrorMultipliers(unsigned idx);
  // std::vector<double> const& getErrorAdditions(unsigned idx);
  void applyErrorShifts(unsigned idx, FastHisto const& nominal, FastHisto & result);

  Double_t evaluate() const;

  RooArgList * setupBinPars();

  std::unique_ptr<RooArgSet> getSentryArgs() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  friend class CMSHistErrorPropagatorV;

 protected:
  RooRealProxy x_;
  RooListProxy funcs_;
  RooListProxy coeffs_;
  RooListProxy binpars_;
  mutable std::vector<CMSHistFunc const*> vfuncs_; //!
  mutable std::vector<RooAbsReal const*> vcoeffs_; //!
  mutable std::vector<RooAbsReal *> vbinpars_; //!
  std::vector<unsigned> bintypes_;

  mutable std::vector<double> coeffvals_; //!
  mutable std::vector<std::vector<double>> valvec_; //!
  mutable std::vector<std::vector<double>> err2vec_; //!
  mutable FastHisto valsum_; //!
  mutable FastHisto scaledvalsum_; //!
  mutable std::vector<double> err2sum_; //!
  mutable std::vector<double> toterr_; //!
  mutable std::vector<std::vector<double>> binmods_; //!
  mutable std::vector<std::vector<double>> scaledbinmods_; //!
  mutable SimpleCacheSentry sentry_; //!
  mutable SimpleCacheSentry binsentry_; //!
  mutable std::vector<double> data_; //!

  int v;
  mutable bool initialized_; //! not to be serialized


  void initialize() const;
  void fillSumAndErr(int eval = 0) const;


 private:
  ClassDef(CMSHistErrorPropagator,1)
};

class CMSHistErrorPropagatorV {
 public:
  CMSHistErrorPropagatorV(const CMSHistErrorPropagator &,
                              const RooAbsData& data,
                              bool includeZeroWeights = false);
  void fill(std::vector<Double_t>& out) const;

 private:
  const CMSHistErrorPropagator& hpdf_;
  int begin_, end_, nbins_;
  struct Block {
    int index, begin, end;
    Block(int i, int begin_, int end_) : index(i), begin(begin_), end(end_) {}
  };
  std::vector<Block> blocks_;
  std::vector<int> bins_;
};

#endif
