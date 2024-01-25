/** \file reward.c
  * \brief Compute RL reward of an initial condition close to nominal orbit
  *
  * Given an initial condition close to the nominal one, integrate both the
  * perturbed and nominal orbits, and compute the RL (negative) reward (or
  * loss) as the distance between both final conditions.
  *
  * NOTES: 
  *		The nominal halo orbit is not actually integrated, we use its Fourier
  *		representation for efficiency.
  *	
  * USAGE:	./reward <coefs.dat 
  *
  * CALLED BY:	
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>				// sqrt
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "int_rtbp.h"
#include "fft_module.h"

/** \brief Period of nominal halo orbit */
static const double T = 0.3059226605957322E+01;

/**
 * \brief Given vector \f$ v \f$, compute its Euclidean norm
 *
 * @param[in] n     dimension of vector
 * @param[in] v     vector
 *
 * @returns			Euclidean norm of vector
 */
double norm(int n, double v[n])
{
	double res = 0;
	for(int i=0; i<n; i++)
		res += (v[i]*v[i]);
	return sqrt(res);
}

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

	// q(T) is the point in the perturbed orbit
	double q[6];
	
	// v = p-q
	double v[6];

	// Auxiliary variables
	int i;

	// Masde passed me the i.c. as x,y,z,xd,yd,zd
	double q_Masde[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00};

	/* Read Fourier coefs from stdin */
	read_coefs(N, A, B);

	/* Use a random number generator to produce variates from a Normal
	 * distribution */
	const gsl_rng_type *rng_type;
	gsl_rng *r;

	/* create a generator chosen by the
    environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();

	rng_type = gsl_rng_default;
	r = gsl_rng_alloc (rng_type);

	/* My vectorfield expects the order x,px,y,py,z,pz. 
	   Recall that xd=px+y; yd=py-x; zd=pz;
       Hence px=xd-y; py=yd+x; pz=zd; 
	   */
	q[0] = q_Masde[0];
	q[1] = q_Masde[3]-q_Masde[1]; 
	q[2] = q_Masde[1]; 
	q[3] = q_Masde[4]+q_Masde[0]; 
	q[4] = q_Masde[2]; 
	q[5] + q_Masde[5];

	fprintf(stderr, "Initial condition (nominal orbit): %e %e %e %e %e %e\n", q[0],
			q[1], q[2], q[3], q[4], q[5]);

	/* Initial condition of perturbed orbit: we just add some noise to nominal
	 * i.c. */
	for(i=0; i<6; i++)
		q[i] += gsl_ran_gaussian(r, q[i]*1.e-5);

	fprintf(stderr, "Initial condition (perturbed orbit): %e %e %e %e %e %e\n", q[0],
			q[1], q[2], q[3], q[4], q[5]);

	/* Integrate nominal orbit, using Fourier */
	eval_orbit(N, A, B, deg, T, p);

	/* Integrate perturbed orbit, with params tol=1.e-14, hmin=hmax=1.e-2,
	 * bPrint=0 */
    int_rtbp(T, q, 1.e-14, 1.e-2, 1.e-2, 0);

	fprintf(stderr, "Final condition (nominal orbit): %e %e %e %e %e %e\n", p[0],
			p[1], p[2], p[3], p[4], p[5]);
	fprintf(stderr, "Final condition (perturbed orbit): %e %e %e %e %e %e\n", q[0],
			q[1], q[2], q[3], q[4], q[5]);

	/* RL negative reward of this orbit is the distance between final
	 * conditions */
	for(i=0; i<6; i++)
		v[i] = p[i]-q[i];

	fprintf(stderr, "Negative reward = |p(T)-q(T)| = %e\n", norm(6,v));

	gsl_rng_free(r);
}
