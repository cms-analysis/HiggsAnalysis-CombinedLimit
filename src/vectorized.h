#include "vdt/vdtMath.h"

namespace vectorized {
    // oarray += coeff * iarray
    void mul_add(const uint32_t size, double coeff, double const * __restrict__ iarray, double* __restrict__ oarray) ;

    // nll_reduce = sum ( weights * log(pdfvals/sumCoeff) )
    double nll_reduce(const uint32_t size, double* __restrict__ pdfvals, double const * __restrict__ weights, double sumcoeff, double *  __restrict__ workingArea) ;

    // gaussians
    void gaussians(const uint32_t size, double mean, double sigma, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea, double * __restrict__ workingArea2) ;
}

