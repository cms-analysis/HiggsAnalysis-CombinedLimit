#ifndef ROOPIECEWISEPOLYNOMIAL_H
#define ROOPIECEWISEPOLYNOMIAL_H
#include <vector>
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

// Piecewise polynomial that looks like
// -- . -- ... -- . --
class RooPiecewisePolynomial : public RooAbsReal{
protected:
  RooRealProxy xvar;

  // First [0,...,nnodes-1] parameters are nodes
  // There should be 2*ndof_endfcn+(nfcn-2)*ndof_middlefcn more parameters for the free dofs in the functions
  // where ndof_endfcn=polyndof-1 and ndof_middlefcn=polyndof-2
  RooListProxy parList;

  const int nfcn; // Number of piecewise functions
  const int polyndof; // Ndof of the polynomial (e.g. 4: cubic)
  const int nnodes; // How many nodes are in between
  const int ndof_endfcn; // Number of degrees of freedom in fcns in the middle of the nodes
  const int ndof_middlefcn; // Number of degrees of freedom in fcns outside the nodes

  double eval(double x, std::vector<double> const& par)const;

public:
  RooPiecewisePolynomial();
  RooPiecewisePolynomial(const char* name, const char* title);
  RooPiecewisePolynomial(const char* name, const char* title, RooAbsReal& xvar_, RooArgList const& parList_, const int nfcn_, const int polyndof_);
  RooPiecewisePolynomial(RooPiecewisePolynomial const& other, const char* name=0);
  virtual TObject* clone(const char* newname)const{ return new RooPiecewisePolynomial(*this, newname); }
  inline virtual ~RooPiecewisePolynomial(){}

  double evaluate()const;

  ClassDef(RooPiecewisePolynomial, 1)

};

#endif
