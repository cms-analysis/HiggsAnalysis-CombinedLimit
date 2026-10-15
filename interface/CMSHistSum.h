#ifndef CMSHistSum_h
#define CMSHistSum_h
#include <ostream>
#include <vector>
#include <memory>
#if __has_include(<ROOT/RVersion.hxx>)
#include <ROOT/RVersion.hxx>
#else
#include <RVersion.h>
#endif
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

namespace RooFit {
class EvalContext;
}

class CMSHistSum : public RooAbsReal {
private:
  struct BarlowBeeston {
    bool init = false;
    std::vector<unsigned> use;
    std::vector<double> dat;
    std::vector<double> valsum;
    std::vector<double> toterr;
    std::vector<double> res;
    std::vector<double> gobs;
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 32, 0)
  void doEval(RooFit::EvalContext &ctx) const override;
#endif

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

  /// See cmsHistBarlowBeestonGradientTerms() in CMSHistHelpers.h.
  void barlowBeestonGradientTerms(std::vector<double>& slopes, std::vector<double>& ratios) const;

  inline FastHisto const& cache() const { return cache_; }

  RooArgList const& coefList() const { return coeffpars_; }
  // RooArgList const& funcList() const { return funcs_; }

  // Per-process views (CMSHistSumWrapper), in the same order as coefList().
  // Uses the wrappers stored in the workspace if there are any, otherwise
  // they are created on first call and owned by this object.
  RooArgList wrapperList() const;

  // Fill result with the current shape of process idx, including its share
  // of the bin-wise shifts
  void fillProcessShape(unsigned idx, FastHisto& result) const;

  RooAbsReal const& getXVar() const { return x_.arg(); }

  // Read access to the persistent template data, for the code generation of
  // the "codegen" evaluation backend.
  int numProcesses() const { return n_procs_; }
  int numMorphs() const { return n_morphs_; }
  RooArgList const& morphList() const { return morphpars_; }
  RooArgList const& binParList() const { return binpars_; }
  RooArgList const& externalMorphList() const { return external_morphs_; }
  std::vector<FastTemplate> const& storage() const { return storage_; }
  std::vector<int> const& processFields() const { return process_fields_; }
  std::vector<int> const& vmorphFields() const { return vmorph_fields_; }
  std::vector<FastTemplate> const& binErrors() const { return binerrors_; }
  std::vector<CMSHistFunc::VerticalSetting> const& verticalTypes() const { return vtype_; }
  std::vector<double> const& verticalSmoothRegions() const { return vsmooth_par_; }
  std::vector<std::vector<unsigned>> const& binTypes() const { return bintypes_; }

  static void EnableFastVertical();
  friend class CMSHistV<CMSHistSum>;

  void injectExternalMorph(int idx, CMSExternalMorph& morph);

  void runBarlowBeeston() const;

protected:
  RooRealProxy x_;

  RooListProxy morphpars_;
  RooListProxy coeffpars_;
  RooListProxy binpars_;

  int n_procs_;
  int n_morphs_;

  std::vector<FastTemplate> storage_;  // All nominal and vmorph templates
  std::vector<int> process_fields_; // indices for process templates in storage_
  std::vector<int> vmorph_fields_; // indices for vmorph templates in storage_

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
  mutable std::vector<double> toterr_;  //!
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

  void updateMorphs() const;


 private:
  ClassDefOverride(CMSHistSum,2)
};

#endif
