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
  virtual TObject* clone(const char* newname) const { return new SimpleProdPdf(*this, newname); }
  inline virtual ~SimpleProdPdf() {}

  RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet* nset,
                             const RooNumIntConfig* cfg, const char* rangeName) const;

  const RooArgList& pdfList() const { return _pdfList; }
  void printMetaArgs(std::ostream& os) const;

  // Copied from RooProdPdf - needed to get the correct toy/Asimov binning
  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const;

 protected:
  RooListProxy _pdfList;
  std::vector<int> _pdfSettings;  // 0 = normal, 1 = conditional
  virtual Double_t evaluate() const;

 private:
  ClassDef(SimpleProdPdf, 1)
};

#endif
