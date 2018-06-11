#include "vdt/vdtMath.h"

namespace vectorized {
    // oarray += coeff * iarray
    void mul_add(const uint32_t size, double coeff, double const * __restrict__ iarray, double* __restrict__ oarray) ;

    // oarray += (coeff * iarray)^2
    void mul_add_sqr(const uint32_t size, double coeff, double const * __restrict__ iarray, double* __restrict__ oarray) ;

    // oarray += sqrt(iarray)
    void sqrt(const uint32_t size, double const * __restrict__ iarray, double* __restrict__ oarray) ;

    // oarray *= iarray
    void mul_inplace(const uint32_t size, double const * __restrict__ iarray, double* __restrict__ oarray) ;

    // nll_reduce = sum ( weights * log(pdfvals/sumCoeff) )
    double nll_reduce(const uint32_t size, double* __restrict__ pdfvals, double const * __restrict__ weights, double sumcoeff, double *  __restrict__ workingArea) ;

    // gaussians
    void gaussians(const uint32_t size, double mean, double sigma, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea, double * __restrict__ workingArea2) ;

    // exponentials
    void exponentials(const uint32_t size, double lambda, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea) ;

    // powers
    void powers(const uint32_t size, double lambda, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea) ;

    // dot product of two vectors 
    double dot_product(const uint32_t size, double const * __restrict__ iarray, double const * __restrict__ iarray2) ;
}

