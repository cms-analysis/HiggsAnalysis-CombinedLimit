/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooPower.cxx 25185 2008-08-20 14:00:42Z wouter $
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
// Power function p.d.f
// END_HTML
//

#include "RooAbsPdf.h"
#include "RooFit.h"

#include <math.h>
#include <cassert>

#include "../interface/HGGRooPdfs.h"
#include "../interface/CombineMathFuncs.h"
#include "RooRealVar.h"
#include "Rtypes.h"
#include "RtypesCore.h"

ClassImp(RooPower)

    //_____________________________________________________________________________
    RooPower::RooPower(const char* name, const char* title, RooAbsReal& x, RooAbsReal& c)
    : RooAbsPdf(name, title),
      x(/*theName=*/"x", /*desc=*/"Dependent", this, x),
      c(/*theName=*/"c", /*desc=*/"Power", this, c) {}

//_____________________________________________________________________________
RooPower::RooPower(const RooPower& other, const char* name)
    : RooAbsPdf(other, name), x(/*theName=*/"x", this, other.x), c(/*theName=*/"c", this, other.c) {}

//_____________________________________________________________________________
Double_t RooPower::evaluate() const {
  //cout << "pow(x=" << x << ",c=" << c << ")=" << pow(x,c) << endl ;
  return RooFit::Detail::MathFuncs::powerEvaluate(x, c);
}

//_____________________________________________________________________________
Int_t RooPower::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x))
    return 1;
  return 0;
}

//_____________________________________________________________________________
Double_t RooPower::analyticalIntegral(Int_t code, const char* rangeName) const {
  if (code == 1) {
    return RooFit::Detail::MathFuncs::powerIntegral(x.min(rangeName), x.max(rangeName), c);
  }
  assert(0);
  return 0;
}
