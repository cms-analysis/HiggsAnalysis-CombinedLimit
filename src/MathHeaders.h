#ifndef MathHeaders_h
#define MathHeaders_h

#ifndef COMBINE_NO_VDT
#include "vdt/vdtMath.h"
#endif

#include <cmath>

#ifdef COMBINE_NO_VDT

#define my_exp std::exp
#define my_log std::log
inline double my_inv(double x) { return 1. / x; }

#else

#define my_exp vdt::fast_exp
#define my_log vdt::fast_log
#define my_inv vdt::fast_inv

#endif

namespace gbrmath {

  inline double fast_pow(double base, double exponent) {
    if (base == 0. && exponent > 0.)
      return 0.;
#ifdef COMBINE_NO_VDT
    else if (base > 0.)
      return std::exp(exponent * std::log(base));
#else
    else if (base > 0.)
      return vdt::fast_exp(exponent * vdt::fast_log(base));
#endif
    else
      return std::nan("");
  }

  //   inline float fast_powf(float base, float exponent) {
  //     if (base==0. && exponent>0.) return 0.;
  //     else if (base>0.) return vdt::fast_expf(exponent*vdt::fast_logf(base));
  //     else return std::nanf("");
  //   }

}  // namespace gbrmath

#endif
