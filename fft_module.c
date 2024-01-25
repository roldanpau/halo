/** \file fft_module.c
  * \brief Functions to compute FFT
  *
  * NOTES: 
  *
  * CALLED BY:  
  *
  */

#include <math.h>				// fabs, M_PI
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>

/** \brief Period of nominal halo orbit */
static const double T = 0.3059226605957322E+01;

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
  * \brief Write coefs of trigonometric series to stdout
  *
  * @param[in]      N				Number of coefs in trig series.
  * @param[in]		A               For each coordinate, coefs of cosine series A_1, ..., A_N
  * @param[in]		B               For each coordinate, coefs of sine series B_1, ..., B_N
  */

void write_coefs(int N, double A[6][N], double B[6][N])
{
	/* Write coefs to stdout */
	for(int d=0; d<6; d++)
	{
		for (int i = 0; i < N; i++)
		{
			fprintf (stdout, "%le %le ", A[d][i], B[d][i]);
		}
		fprintf(stdout, "\n");
	}
}

/**
  * \brief Read coefs of trigonometric series from stdin
  *
  * @param[in]      N				Number of coefs in trig series.
  * @param[out]		A               For each coordinate, coefs of cosine series A_1, ..., A_N
  * @param[out]		B               For each coordinate, coefs of sine series B_1, ..., B_N
  */

void read_coefs(int N, double A[6][N], double B[6][N])
{
	/* Read coefs from stdin */
	for(int d=0; d<6; d++)
	{
		for (int i = 0; i < N; i++)
		{
			fscanf (stdin, "%le %le ", &(A[d][i]), &(B[d][i]));
		}
		fscanf(stdin, "\n");
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

/**
 * \brief Given time \f$ t \f$, evaluate orbit at \f$ t \f$
 *
 * This function returns an approximation to the point orbit(t), by evaluating
 * the trigonometric series of each component (x,px,y,py,z,pz).
 *
 * @param[in] N     Number of total coefs of Fourier series
 * @param[in] A     Fourier coeffs of every component
 * @param[in] B     Fourier coeffs of every component
 * @param[in] deg   Desired degree we want to use when evaluating series
 * @param[in] t     Argument of Fourier series
 * @param[out] p    Point orbit(t)
 */
void eval_orbit (int N, double A[6][N], double B[6][N], int deg, double t,
		double p[6]) 
{
	for(int i=0; i<6; i++)
		p[i] = eval_trig_series(deg, A[i], B[i], 2*M_PI*t/T);
}
