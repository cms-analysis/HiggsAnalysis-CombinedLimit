#include "vectorized.h"

void vectorized::mul_add(const uint32_t size, double coeff, double const * __restrict__ iarray, double* __restrict__ oarray) {
    for (uint32_t i = 0; i < size; ++i) {
        oarray[i] += coeff * iarray[i];
    } 
}
