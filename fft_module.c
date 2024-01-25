/** \file fft_module.c
  * \brief Functions to compute FFT
  *
  * NOTES: 
  *
  * CALLED BY:  
  *
  */

#include <math.h>				// fabs
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>

/**
  * \brief Compute radix-2 FFT for real data and return trigonometric series
  *
  * This function computes a radix-2 FFT of length n on the real array y.
  * Moreover, it returns the coefficients \f$ A_k, B_k \f$ of the corresponding
  * trigonometric series,
  *
  * \f[ \sum_{k=0}^N B_k\sin{k\phi} + \sum_{k=0}^N A_k\cos{k\phi}. \f]
  *
  * @param[in]      n				Length of array y. Must be a power of 2.
  * @param[inout]   y               On entry, it contains the real array data.
  * On exit, it is replaced by the half-complex FFT sequence (see GSL
  * documentation).
  * @param[out]     A               Coefs of cosine series.
  * @param[out]     B               Coefs of sine series.
  */

void real_fft(int n, double y[n], double A[n/2+1], double B[n/2+1])
{
	/* auxiliary variables */
	int i;

	/* Compute the frequency coefficients of the real sequence y */
	/* For k < n/2, the real part of the k-th term is stored in location k, 
	   and the corresponding imaginary part is stored in location n-k.
	   The terms for k=0 and k=n/2 are both purely real, and are stored 
	   in locations 0 and n/2. */
  gsl_fft_real_radix2_transform(y, 1, n);	// stride == 1
 
    /* Write coefficients to stdout
  for (i = 0; i < n; i++)
    {
	    fprintf (stdout, "%e\n", y[i]);
    } */

 for (i = 0; i <= n/2; i++)
    {
		double re;
		double im;

		/* FFT storage Convention for first and last coefs is special */
		if(i == 0 || i == n/2) 
		{
			re = y[i];
			im = 0;
		}
		else
		{
			re = y[i];
			im = y[n-i];
		}

		/* Store coefficients of trigonometric series */ 
		if(i == 0 || i == n/2) 
		{
			A[i] = re/n;	// A0
			B[i] = 0.0;		// B0
		}
		else
		{
			A[i] = 2*re/n;
			B[i] = -2*im/n;
		}

		/* Write norm to stderr */
	    //fprintf (stderr, "%d %e %e\n", i, fabs(A[i]), fabs(B[i]));
    }
}

/**
 * \brief Given \f$ x \f$, evaluate trigonometric series at \f$ x \f$
 *
 * @param[in] N     Degree of Fourier series
 * @param[in] A     A = (A_0, A_1, ..., A_N), Fourier coeffs 
 * @param[in] B     B = (B_0, B_1, ..., B_N), Fourier coeffs
 * @param[in] x		Argument of Fourier series
 *
 * @returns		Value of the trig series evaluated at x
 */
double eval_trig_series (size_t N, double A[N+1], double B[N+1], double x)
{
    double res = 0.0;
    for(int j=0; j<=N; j++)
    {
        res = res + (A[j]*cos(j*x) + B[j]*sin(j*x));
    }
    return res;
}

