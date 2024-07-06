#ifndef logKappa_h
#define logKappa_h

template<typename T>
T logKappa(const T x, const double &kappaHigh, const double &kappaLow) {
    if (fabs(x) >= 0.5) return (x >= 0 ? kappaHigh : -kappaLow);
    // interpolate between log(kappaHigh) and -log(kappaLow)
    //    logKappa(x) = avg + halfdiff * h(2x)
    // where h(x) is the 3th order polynomial
    //    h(x) = (3 x^5 - 10 x^3 + 15 x)/8;
    // chosen so that h(x) satisfies the following:
    //      h (+/-1) = +/-1
    //      h'(+/-1) = 0
    //      h"(+/-1) = 0
    const double logKhi =  kappaHigh;
    const double logKlo = -kappaLow;
    const double avg = 0.5*(logKhi + logKlo), halfdiff = 0.5*(logKhi - logKlo);
    const double twox = x+x, twox2 = twox*twox;
    const double alpha = 0.125 * twox * (twox2 * (3*twox2 - 10.) + 15.);
    return avg + alpha*halfdiff;
}

#endif
