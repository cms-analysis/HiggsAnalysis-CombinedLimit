#ifndef CMSHistSum_h
#define CMSHistSum_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "FastTemplate_Old.h"
#include "SimpleCacheSentry.h"
#include "CMSHistFunc.h"
#include "CMSHistV.h"

class CMSHistSum : public RooAbsReal {
private:
  struct BarlowBeeston {
    bool init = false;
    std::vector<unsigned> use;
    std::vector<double> dat;
    std::vector<double> valsum;
    std::vector<double> toterr;
    std::vector<double> err;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> tmp;
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> res;
    std::vector<double> gobs;
    std::set<RooAbsArg*> dirty_prop;
    std::vector<RooRealVar*> push_res;
  };
public:

  CMSHistSum();

  CMSHistSum(const char* name, const char* title, RooRealVar& x,
                         RooArgList const& funcs, RooArgList const& coeffs);

  CMSHistSum(CMSHistSum const& other, const char* name = 0);

  TObject* clone(const char* newname) const override {
    return new CMSHistSum(*this, newname);
  }

  ~CMSHistSum() override {;}

  Double_t evaluate() const override;

  std::map<std::string, Double_t> getProcessNorms() const;

  RooArgList * setupBinPars(double poissonThreshold);

  std::unique_ptr<RooArgSet> getSentryArgs() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const override;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const override;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const override;

  void setData(RooAbsData const& data) const;

  void setAnalyticBarlowBeeston(bool flag) const;

  inline FastHisto const& cache() const { return cache_; }

  RooArgList const& coefList() const { return coeffpars_; }
  // RooArgList const& funcList() const { return funcs_; }

  RooAbsReal const& getXVar() const { return x_.arg(); }

  static void EnableFastVertical();
  friend class CMSHistV<CMSHistSum>;

  void injectExternalMorph(int idx, CMSExternalMorph& morph);

 protected:
  RooRealProxy x_;

  RooListProxy morphpars_;
  RooListProxy coeffpars_;
  RooListProxy binpars_;

  int n_procs_;
  int n_morphs_;

  std::vector<FastTemplate> storage_;  // All nominal and vmorph templates
  std::vector<int> process_fields_; // Indicies for process templates in storage_
  std::vector<int> vmorph_fields_; // Indicies for vmorph templates in storage_

  std::vector<FastTemplate> binerrors_; // Bin errors for each process

  std::vector<CMSHistFunc::VerticalSetting> vtype_; // Vertical morphing type for each process
  std::vector<double> vsmooth_par_; // Vertical morphing smooth region for each process

  mutable std::vector<CMSHistFunc const*> vfuncstmp_; //!
  mutable std::vector<RooAbsReal const*> vcoeffpars_; //!
  mutable std::vector<RooAbsReal const*> vmorphpars_; //!
  mutable std::vector<std::vector<RooAbsReal *>> vbinpars_; //!
  std::vector<std::vector<unsigned>> bintypes_;

  mutable std::vector<double> coeffvals_; //!

  mutable std::vector<FastHisto> compcache_; //!
  mutable FastHisto staging_; //!
  mutable FastHisto valsum_; //!
  mutable FastHisto cache_;

  mutable std::vector<double> err2sum_; //!
  mutable std::vector<double> toterr_; //!
  mutable std::vector<std::vector<double>> binmods_; //!
  mutable std::vector<std::vector<double>> scaledbinmods_; //!

  mutable SimpleCacheSentry sentry_; //!
  mutable SimpleCacheSentry binsentry_; //!

  mutable std::vector<double> data_; //!

  mutable BarlowBeeston bb_; //!

  mutable bool initialized_; //! not to be serialized

  mutable bool analytic_bb_; //! not to be serialized

  mutable std::vector<double> vertical_prev_vals_; //! not to be serialized
  mutable int fast_mode_; //! not to be serialized
  static bool enable_fast_vertical_; //! not to be serialized

  RooListProxy external_morphs_;
  std::vector<int> external_morph_indices_;

  inline int& morphField(int const& ip, int const& iv) {
    return vmorph_fields_[ip * n_morphs_ + iv];
  }

  void initialize() const;
  void updateCache() const;
  inline double smoothStepFunc(double x, int const& ip) const;


  void runBarlowBeeston() const;

  void updateMorphs() const;


 private:
  ClassDefOverride(CMSHistSum,2)
};

#endif
