#ifndef CombineTools_SOverBTools_h
#define CombineTools_SOverBTools_h
#include "TH1.h"

namespace ch {
struct SOverBInfo {
  double s;
  double b;
  double x_lo;
  double x_hi;
  SOverBInfo() { ; }
  SOverBInfo(TH1F const* sig, TH1F const* bkg, unsigned steps, double frac);
};

double IntegrateFloatRange(TH1F const* hist, double xmin, double xmax);
}

#endif
