#ifndef SimpleProdPdf_h
#define SimpleProdPdf_h

#include <iostream>
#include <list>
#include <memory>
#include <vector>
#include "Rtypes.h"
#include "RtypesCore.h"
#include "RooAbsArg.h"
#include "RooAbsCacheElement.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooObjCacheManager.h"
#include "RooProdPdf.h"

class SimpleProdPdf : public RooAbsPdf {
 public:
  SimpleProdPdf();
  SimpleProdPdf(const char* name, const char* title, RooArgList const& pdfList);
  SimpleProdPdf(const char* name, const char* title, RooProdPdf const& prodPdf);
  SimpleProdPdf(const char* name, const char* title, RooArgList& pdfList,
                std::vector<int> const& pdfSettings);
  SimpleProdPdf(const SimpleProdPdf& other, const char* name = 0);
  TObject* clone(const char* newname) const override { return new SimpleProdPdf(*this, newname); }
  inline ~SimpleProdPdf() override {}

  // The integrals over the observables are products of the integrals of the
  // components, built with the RooFit hooks for analytical integration so that
  // RooFit always gets a regular RooRealIntegral from createIntegral().
  bool forceAnalyticalInt(const RooAbsArg& dep) const override;
  Int_t getAnalyticalIntegralWN(RooArgSet& allVars,
                                RooArgSet& analVars,
                                const RooArgSet* normSet,
                                const char* rangeName = nullptr) const override;
  double analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName = nullptr) const override;

  const RooArgList& pdfList() const { return _pdfList; }
  void printMetaArgs(std::ostream& os) const override;

  // Copied from RooProdPdf - needed to get the correct toy/Asimov binning
  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override;

 protected:
  RooListProxy _pdfList;
  std::vector<int> _pdfSettings;  // 0 = normal, 1 = conditional
  Double_t evaluate() const override;

 private:
  class CacheElem : public RooAbsCacheElement {
   public:
    std::unique_ptr<RooAbsReal> integral;
    RooArgList containedArgs(Action) override { return RooArgList(*integral); }
  };
  // Named like the other data members of this class rather than with the m_ prefix
  mutable RooObjCacheManager _cacheMgr;  //! The cache manager for the integrals NOLINT(readability-identifier-naming)

  std::unique_ptr<RooAbsReal> makeIntegral(const RooArgSet& iset, const RooArgSet* nset, const char* rangeName) const;
  Int_t getIntegralCode(const RooArgSet& iset, const RooArgSet* nset, const char* rangeName) const;

  ClassDefOverride(SimpleProdPdf, 1)
};

#endif
