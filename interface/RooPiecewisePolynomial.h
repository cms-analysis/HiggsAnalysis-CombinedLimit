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
  // For nfcn=4, polyndof=4, the order of parameters is
  // node 0, node 1, node 2 (nnodes=nfcn-1), and then
  // a0 + a1*x +        a2*x^2 + not_a_par*x^3 (= fcn 0)
  // b0 + b1*x + not_a_par*x^2 + not_a_par*x^3 (= fcn 1)
  // c0 + c1*x + not_a_par*x^2 + not_a_par*x^3 (= fcn 2)
  // d0 + d1*x + not_a_par*x^2 +        d3*x^3 (= fcn 3)
  // Continuity between the fcns is guaranteed, but notice that smoothness is not.
  RooListProxy parList;

  const int nfcn; // Number of piecewise functions
  const int polyndof; // Ndof of the polynomial (e.g. 4: cubic)
  const int nnodes; // How many nodes are in between (= nfcn-1)
  const int ndof_endfcn; // Number of degrees of freedom in fcns in the middle of the nodes (= polyndof-1)
  const int ndof_middlefcn; // Number of degrees of freedom in fcns outside the nodes (= polyndof-2)

  double eval(double x, std::vector<double> const& par)const;

public:
  RooPiecewisePolynomial(const int nfcn_=1, const int polyndof_=1);
  RooPiecewisePolynomial(const char* name, const char* title, const int nfcn_=1, const int polyndof_=1);
  RooPiecewisePolynomial(const char* name, const char* title, RooAbsReal& xvar_, RooArgList const& parList_, const int nfcn_, const int polyndof_);
  RooPiecewisePolynomial(RooPiecewisePolynomial const& other, const char* name=0);
  virtual TObject* clone(const char* newname)const{ return new RooPiecewisePolynomial(*this, newname); }
  inline virtual ~RooPiecewisePolynomial(){}

  double evaluate()const;
  double evaluate(double* x, double* p)const; // For calling in a TF1 object

  ClassDef(RooPiecewisePolynomial, 1)

};

#endif
