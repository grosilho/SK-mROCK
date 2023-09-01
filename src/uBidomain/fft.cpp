#include "fft.h"
#include <iostream>
extern "C" {

/* evaluates a trigonometrical expressions of the form
 * s1 = \sum_i=0^n-1 Fcoeff[i] cos(2*pi*i*xi)
 * s2 = \sum_i=1^n-1 Fcoeff[i] sin(2*pi*i*xi)
 * by the algorithm of Reinsch from Stoer 8. edition page 92f
 */
int ReinschEval(double *s1, double *s2, double *FCoeff, int n, double xi) {
  int i = 0;
  long double lambda = 0;
  long double delta = 0;
  long double U = 0;

  *s1 = *s2 = 0;

  /* transform the angle into the interval [0,1) */
  xi -= (int)xi;

  /* case of cos(2*pi*xi) <= 0 */
  if (xi >= 0.25 && xi <= 0.75) {
    lambda = cos(M_PI * xi);
    lambda *= 4 * lambda;
    U = delta = 0;
    for (i = n - 1; i >= 0; --i) {
      U = delta - U;
      delta = lambda * U - delta + FCoeff[i];
    }
  }
  /* case of cos(2*pi*xi) > 0 */
  else {
    lambda = sin(M_PI * xi);
    lambda *= -4 * lambda;
    U = delta = 0;
    for (i = n - 1; i >= 0; --i) {
      U = delta + U;
      delta = lambda * U + delta + FCoeff[i];
    }
  }
  *s1 = delta - 0.5 * U * lambda;
  *s2 = U * sin(2 * M_PI * xi);
  return 0;
}

int fft(double *FVal, double *FCoeff, int n) {
  fftw_complex *beta = NULL;
  fftw_complex *val = NULL;
  fftw_plan plan;
  double *pFC = NULL;
  int M = 0;
  int i = 0;

  /* initialize fftw stuff and execute the fft                       */
  beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
  val = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
  plan = fftw_plan_dft_1d(n, val, beta, FFTW_FORWARD, FFTW_ESTIMATE);

  memset(val, 0, sizeof(fftw_complex) * n);
  for (i = 0; i < n; ++i) val[i][0] = FVal[i];

  fftw_execute(plan);

  memset(FCoeff, 0, sizeof(double) * n);

  M = n / 2;
  pFC = FCoeff + M;

  pFC[0] = 2. * beta[0][0] / n;
  pFC[-M] = 2. * beta[M][0] / n;
  for (i = 1; i < M; ++i) {
    pFC[-i] = 2. * beta[i][0] / n;
    pFC[i] = -2. * beta[i][1] / n;
  }

  fftw_destroy_plan(plan);
  fftw_free(beta);
  fftw_free(val);
  return 0;
}

int ifft(double *FVal, double *FCoeff, int nCoeff, int n) {
  fftw_complex *beta = NULL;
  fftw_complex *val = NULL;
  fftw_plan plan;
  double *pFC = NULL;
  int M = 0;
  int i = 0;

  /* initialize fftw stuff and execute the fft                       */
  beta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
  val = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
  plan = fftw_plan_dft_1d(n, beta, val, FFTW_BACKWARD, FFTW_ESTIMATE);

  M = nCoeff / 2;
  pFC = FCoeff + M;
  memset(beta, 0, sizeof(fftw_complex) * n); 
  for (i = 1; i < M; ++i) {
    beta[i][0] = 0.5 * pFC[-i];
    beta[i][1] = -0.5 * pFC[i];
    beta[n - i][0] = 0.5 * pFC[-i];
    beta[n - i][1] = 0.5 * pFC[i];
  }
  // handle beta[0] and beta[M] and beta[n-M] differently
  beta[0][0] = 0.5 * pFC[0];
  beta[0][1] = 0;

  beta[M][0] = 0.5 * pFC[-M];
  beta[M][1] = 0;

  beta[n - M][0] = 0.5 * pFC[-M];
  beta[n - M][1] = 0;

  fftw_execute(plan);
  for (i = 0; i < n; ++i) FVal[i] = val[i][0];

  fftw_destroy_plan(plan); 
  fftw_free(beta); 
  fftw_free(val); 
  return 0;
}
}