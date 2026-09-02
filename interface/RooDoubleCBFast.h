#ifndef ROODOUBLECB
#define ROODOUBLECB

#include <ROOT/RConfig.hxx>  // for ROOT_VERSION

// RooCrystalBall was declared "final" until ROOT 6.38.02. The keyword was
// removed on master by dea293f14e05 (2026-02-18, "[RF] Remove `final` keyword
// from RooFit pdf classes") and cherry-picked to v6-38-00-patches as
// 802f9fd5d5d2 on the same day. Note that 6.38.00 was tagged earlier
// (2025-11-27) and therefore still has `final`; 6.38.02 is the first release
// without it. No backport was made to 6.36 or older -- v6-36-14, v6-34-10 and
// v6-32-24 all still declare it final.
//
// Deriving from it is a hard compile error on those versions:
//   error: cannot derive from 'final' base 'RooCrystalBall'
//          in derived type 'RooDoubleCBFast'
// so there we keep the original standalone implementation.
//
// The two variants are deliberately *different* schema versions:
//   v2 - derives from RooCrystalBall, no data members of its own
//   v1 - derives from RooAbsPdf, owns seven RooRealProxy members
// The read rule in src/classes_def.xml migrates v1 payloads into a v2 object
// and is therefore only active when this macro is 1.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 38, 2)
#define COMBINE_ROODOUBLECBFAST_FROM_ROOCRYSTALBALL 1
#else
#define COMBINE_ROODOUBLECBFAST_FROM_ROOCRYSTALBALL 0
#endif

#if COMBINE_ROODOUBLECBFAST_FROM_ROOCRYSTALBALL

#include "RooCrystalBall.h"
#include "RooAbsReal.h"

class RooDoubleCBFast : public RooCrystalBall {
public:
  RooDoubleCBFast() = default;
  RooDoubleCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2
           );

  // In-place migration from v1 (RooAbsPdf base) to v2 (RooCrystalBall base).
  // Uses TClass reflection to set RooCrystalBall's private proxy members
  // without calling the destructor, which would corrupt ROOT's TProcessID state.
  void migrateFromV1(RooAbsReal& x, RooAbsReal& mean, RooAbsReal& width,
                     RooAbsReal& alpha1, RooAbsReal& n1,
                     RooAbsReal& alpha2, RooAbsReal& n2);

  ClassDefOverride(RooDoubleCBFast, 2)
};

#else  // ROOT < 6.38.02: RooCrystalBall is final, keep the standalone class.

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooDoubleCBFast : public RooAbsPdf {
public:
  RooDoubleCBFast();
  RooDoubleCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2
           );
  RooDoubleCBFast(const RooDoubleCBFast& other, const char* name=0) ;
  TObject* clone(const char* newname) const override { return new RooDoubleCBFast(*this,newname); }
  inline ~RooDoubleCBFast() override { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const override ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const override ;

private:

  ClassDefOverride(RooDoubleCBFast,1)
};

#endif  // COMBINE_ROODOUBLECBFAST_FROM_ROOCRYSTALBALL
#endif
