#ifndef SimpleProdPdf_h
#define SimpleProdPdf_h

#include <iostream>
#include <list>
#include "RooAbsPdf.h"
#include "RooListProxy.h"
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

  RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet* nset,
                             const RooNumIntConfig* cfg, const char* rangeName) const override;

  const RooArgList& pdfList() const { return _pdfList; }
  void printMetaArgs(std::ostream& os) const override;

  // Copied from RooProdPdf - needed to get the correct toy/Asimov binning
  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override;

 protected:
  RooListProxy _pdfList;
  std::vector<int> _pdfSettings;  // 0 = normal, 1 = conditional
  Double_t evaluate() const override;

 private:
  ClassDefOverride(SimpleProdPdf, 1)
};

#endif
