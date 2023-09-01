#ifndef MYFFTW_WRAPPER_H_
#define MYFFTW_WRAPPER_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 
#endif

extern "C" {

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* evaluates a trigonometrical expressions of the form
 * s1 = \sum_i=0^n-1 Fcoeff[i] cos(2*pi*i*xi)
 * s2 = \sum_i=1^n-1 Fcoeff[i] sin(2*pi*i*xi)
 * by the algorithm of Reinsch from Stoer 8. edition page 92f
 */
int ReinschEval(double *s1, double *s2, double *FCoeff, int n, double xi);
int fft(double *FVal, double *FCoeff, int n);
int ifft(double *FVal, double *FCoeff, int nCoeff, int n);

}
#endif
