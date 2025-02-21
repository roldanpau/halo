/** \file eval_nominal_orbit.c
  * \brief Integrate nominal orbit from time 0 up to a given time T.
  *
  * Integrate nominal halo orbit from time 0 up to a given time T, to find the
  * final condition p(T).
  *
  * NOTES: 
  *		The nominal halo orbit is not actually integrated, we use its Fourier
  *		representation for efficiency. The Fourier coefficients must be given
  *		in stdin.
  *	
  * USAGE:	./eval_nominal_orbit < coefs.dat
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include "fft_module.h"

/** \brief Period of nominal halo orbit */
static const double T = 0.3059226605957322E+01;

int
main (int argc, char *argv[])
{
	const int n = 128;
	const int N = n/2+1;	// num of Fourier coefs
	const int deg = 15;		// degree of Fourier series

	// Fourier coefs:
	// A[0][:] contains coefs A_0, A_1, ..., A_N of function x(t)
	// A[1][:] contains coefs A_0, A_1, ..., A_N of function p_x(t)
	// ...
	// A[4][:] contains coefs A_0, A_1, ..., A_N of function z(t)
	// A[5][:] contains coefs A_0, A_1, ..., A_N of function p_z(t)
	double A[6][N]; 
	double B[6][N];							

	// p(T) is the Fourier-recovered point in the nominal orbit
	double p[6];

	/* Read Fourier coefs from stdin */
	read_coefs(N, A, B);

	/* Integrate nominal orbit, using Fourier */
	eval_orbit(N, A, B, deg, T/2, p);

	fprintf(stderr, "Final condition (nominal orbit): %e %e %e %e %e %e\n", p[0],
			p[1], p[2], p[3], p[4], p[5]);
}
