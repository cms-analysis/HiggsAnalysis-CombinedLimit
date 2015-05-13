/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBernsteinFast.cxx 44507 2012-06-04 12:30:41Z axel $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

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
