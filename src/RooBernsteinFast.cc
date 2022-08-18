//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// RooBernsteinFast implements a polynomial p.d.f of the form
// <pre>
// f(x) = sum_i a_i * x^i
//</pre>
// By default coefficient a_0 is chosen to be 1, as polynomial
// probability density functions have one degree of freedome
// less than polynomial functions due to the normalization condition
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include "TMath.h"

#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

using namespace std;

templateClassImp(RooBernsteinCoeffs)
templateClassImp(RooBernsteinFast)
