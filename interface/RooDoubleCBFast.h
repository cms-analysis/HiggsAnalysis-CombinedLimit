#ifndef ROODOUBLECB
#define ROODOUBLECB

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

  ClassDefOverride(RooDoubleCBFast, 2)
};
#endif
