#include "getV1.h"

extern "C" {

    int getV1(double *V1, int n) {
        fftw_complex *beta = NULL;
        fftw_complex *val = NULL;
        fftw_plan plan;
        int M = 0;
        int i = 0;
        int j = 0;

        /* initialize fftw stuff and execute the fft                       */
        beta = (fftw_complex *) fftw_malloc(sizeof (fftw_complex) * n);
        val = (fftw_complex *) fftw_malloc(sizeof (fftw_complex) * n);
        plan = fftw_plan_dft_1d(n, beta, val, FFTW_BACKWARD, FFTW_ESTIMATE);

        memset(beta, 0, sizeof (fftw_complex) * n);

        M = n / 2;

        for (i = 1; i < M; ++i) {
            beta[i][0] = -1. / (double) i;
            beta[n - i][0] = beta[i][0];
        }

        beta[0][0] = 0;
        beta[M][0] = -0.5 * 1. / (double) M;

        fftw_execute(plan);

        for (i = 0; i < n; ++i)
            for (j = 0; j < n; ++j)
                V1[i * n + j] = val[abs(i - j)][0] / (-4 * M_PI * n);

        fftw_destroy_plan(plan);
        fftw_free(beta);
        fftw_free(val);

        return 0;
    }
}
