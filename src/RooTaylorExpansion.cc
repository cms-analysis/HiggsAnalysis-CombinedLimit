#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "../interface/RooTaylorExpansion.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"
#include "RooConstVar.h"
#include "RooListProxy.h"
#include "TArrayI.h"

using namespace std;

ClassImp(RooTaylorExpansion);

//_____________________________________________________________________________
RooTaylorExpansion::RooTaylorExpansion(
    const char* name, const char* title, const RooArgList& x,
    const RooArgList& x0, std::vector<std::vector<int>> const& trackers,
    const std::vector<double>& terms)
    : RooAbsReal(name, title),
      _x("x", "Observables", this, kTRUE, kFALSE),
      _x0("x0", "initial point", this, kTRUE, kFALSE),
      _trackers(trackers),
      _terms(terms) {
  _x.add(x);
  _x0.add(x0);
  assert((_trackers.size()) == (_terms.size()));
}

//_____________________________________________________________________________
RooTaylorExpansion::RooTaylorExpansion(const RooTaylorExpansion& other,
                                       const char* name)
    : RooAbsReal(other, name),
      _x("x", this, other._x),
      _x0("x0", this, other._x0),
      _trackers(other._trackers),
      _terms(other._terms) {}

//_____________________________________________________________________________
Double_t RooTaylorExpansion::evaluate() const {
  double fullsum = 0.;

  unsigned nx = _x.getSize();
  std::vector<double> dx(nx);

  for (unsigned i = 0; i < nx; ++i) {
    dx[i] = ((RooAbsReal*)(_x.at(i)))->getVal() -
            ((RooAbsReal*)(_x0.at(i)))->getVal();
  }

  for (unsigned i = 0; i < _trackers.size(); ++i) {
    unsigned n = _trackers[i].size();
    double term = _terms[i];
    for (unsigned t = 0; t < n; ++t) {
      term *= dx[_trackers[i][t]];
    }
    fullsum += term;
  }
  return fullsum;
}


void RooTaylorExpansion::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  for (unsigned i = 0; i < _trackers.size(); ++i) {
    std::cout << " - [";
    for (unsigned j = 0; j < _trackers[i].size(); ++j) {
      std::cout << _trackers[i][j];
      if (j != (_trackers[i].size() - 1)) {
        std::cout << ", ";
      } else {
        std::cout << "] = ";
      }
    }
    std::cout << _terms[i] << "\n";
  }
}
