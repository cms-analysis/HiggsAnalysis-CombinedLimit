#ifndef ROO_TAYLOREXPANSION
#define ROO_TAYLOREXPANSION

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "TVectorD.h"

class RooRealVar;
// class RooFitResult;

#include <map>

class RooTaylorExpansion : public RooAbsReal {
 public:
  RooTaylorExpansion(){};
  RooTaylorExpansion(const char* name, const char* title, const RooArgList& x,
                     const RooArgList& x0,
                     std::vector<std::vector<int>> const& trackers,
                     const std::vector<double>& terms);

  RooTaylorExpansion(const RooTaylorExpansion& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const {
    return new RooTaylorExpansion(*this, newname);
  }
  inline virtual ~RooTaylorExpansion() {}

  void printMultiline(std::ostream& os, Int_t contents,
                                   Bool_t verbose, TString indent) const;

 protected:

  RooListProxy _x;
  RooListProxy _x0;
  std::vector<std::vector<int>> _trackers;
  std::vector<double> _terms;

  Double_t evaluate() const;

 private:
  ClassDef(RooTaylorExpansion,
           1)  // Multivariate Gaussian PDF with correlations
};

#endif
