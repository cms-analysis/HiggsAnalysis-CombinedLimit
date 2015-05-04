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

void vectorized::gaussians(const uint32_t size, double mean, double sigma, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea, double * __restrict__ workingArea2)
{
    double xscale = -0.5/(sigma*sigma);
    for (uint32_t i = 0; i < size; ++i) {
        workingArea[i] = xscale * std::pow(xvals[i] - mean, 2);
    }
    vdt::fast_expv(size, workingArea, workingArea2);
    double inorm = 1.0/norm;
    for (uint32_t i = 0; i < size; ++i) {
        out[i] = inorm*workingArea2[i];
    }
}

void vectorized::exponentials(const uint32_t size, double lambda, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea)
{
    //out[i] = std::exp(xvals[i]*lambda) * nfact; nfact = 1.0/norm
    double lognfact = -std::log(norm);
    for (uint32_t i = 0; i < size; ++i) {
        workingArea[i] = xvals[i] * lambda + lognfact;
    }
    vdt::fast_expv(size, workingArea, out);
}

void vectorized::powers(const uint32_t size, double exponent, double norm, const double* __restrict__ xvals, double * __restrict__ out, double * __restrict__ workingArea)
{
    //out[i] = std::pow(xvals[i],exponent) * nfact; // nfact = 1.0/norm
    double lognfact = -std::log(norm);
    vdt::fast_logv(size, xvals, workingArea);
    for (uint32_t i = 0; i < size; ++i) {
        workingArea[i] = workingArea[i]*exponent + lognfact;
    }
    vdt::fast_expv(size, workingArea, out);
}

double vectorized::dot_product(const uint32_t size, double const * __restrict__ vec1, double const *  __restrict__ vec2) {
    double ret = 0;
    for (uint32_t i = 0; i < size; ++i) {
        ret += vec1[i]*vec2[i];
    }
    return ret;
}


