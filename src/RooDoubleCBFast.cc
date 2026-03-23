#include "../interface/RooDoubleCBFast.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"

using namespace RooFit;

 ClassImp(RooDoubleCBFast) 

 RooDoubleCBFast::RooDoubleCBFast(const char *name, const char *title, 
                    RooAbsReal& _x,
                    RooAbsReal& _mean,
                    RooAbsReal& _width,
                    RooAbsReal& _alpha1,
                    RooAbsReal& _n1,
                    RooAbsReal& _alpha2,
                    RooAbsReal& _n2
                    ) :
   RooCrystalBall(name, title, _x, _mean, _width, _alpha1, _n1, _alpha2, _n2)
 { 
 }

 
