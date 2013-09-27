#include "vdt/vdtMath.h"

namespace vectorized {
    // oarray += coeff * iarray
    void mul_add(const uint32_t size, double coeff, double const * __restrict__ iarray, double* __restrict__ oarray) ;
}

