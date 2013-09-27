#include "vectorized.h"

void vectorized::mul_add(const uint32_t size, double coeff, double const * __restrict__ iarray, double* __restrict__ oarray) {
    for (uint32_t i = 0; i < size; ++i) {
        oarray[i] += coeff * iarray[i];
    } 
}
double vectorized::nll_reduce(const uint32_t size, double* __restrict__ pdfvals, double const * __restrict__ weights, double sumcoeff,  double *  __restrict__ workingArea) {
    double invsum = 1.0/sumcoeff;
    for (uint32_t i = 0; i < size; ++i) {
        pdfvals[i] *= invsum;
    }

    vdt::fast_logv(size, pdfvals, workingArea);

    for (uint32_t i = 0; i < size; ++i) {
        pdfvals[i] = weights[i] * workingArea[i];
    }


    double ret = 0;
    for (uint32_t i = 0; i < size; ++i) {
        ret += pdfvals[i];
    }

    return ret;
}
