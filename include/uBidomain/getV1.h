#ifndef NBEM_GET_V1_H_
#define NBEM_GET_V1_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#include <fftw3.h>
#include <math.h>
#include <string.h>

extern "C" {

    int getV1(double *V1, int n);
}
#endif
