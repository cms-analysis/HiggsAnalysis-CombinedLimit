#ifndef CMSHistFuncFactory_h
#define CMSHistFuncFactory_h
#include <algorithm>
#include <vector>
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooBinning.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "HiggsAnalysis/CombinedLimit/interface/CMSHistFunc.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/AsymPow.h"
#pragma GCC diagnostic pop

namespace ch {

class CMSHistFuncFactory {
public:
  void Run(CombineHarvester& cb, RooWorkspace& ws, std::map<std::string, std::string> process_vs_norm_postfix_map);
  void Run(CombineHarvester& cb, RooWorkspace& ws);
  void SetHorizontalMorphingVariable(std::map<std::string, RooAbsReal*> &hvar) { mass_var = hvar; }
  CMSHistFuncFactory();
private:
  std::string norm_postfix_ = "norm";
  unsigned v_;
  // RooAbsReal *mass_var;
  std::map<std::string, RooAbsReal*> mass_var;
  void RunSingleProc(CombineHarvester& cb, RooWorkspace& ws, std::string bin, std::string process);
  std::map<std::string, RooRealVar> obs_;
  unsigned hist_mode_;
  bool rebin_;

  TH1F AsTH1F(TH1 const* hist) {
    TH1F res;
    TH1F const* test_f = dynamic_cast<TH1F const*>(hist);
    TH1D const* test_d = dynamic_cast<TH1D const*>(hist);
    if (test_f) {
      test_f->Copy(res);
    } else if (test_d) {
      test_d->Copy(res);
    } else {
      throw std::runtime_error("TH1 is not a TH1F or a TH1D");
    }
    return res;
  }

  RooRealVar VarFromHist(TString name, TString title, TH1 const& hist) {
    RooBinning * binning;
    if (hist.GetXaxis()->GetXbins() && hist.GetXaxis()->GetXbins()->GetSize()) {
      binning = new RooBinning(hist.GetNbinsX(), hist.GetXaxis()->GetXbins()->GetArray());
    } else {
      binning = new RooBinning(hist.GetNbinsX(), hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
    }
    RooRealVar x(name, title, hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
    // x.getBinning().Print();
    x.setBinning(*binning);
    delete binning;
    return x;
  }
};

}

#endif
