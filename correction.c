/** \file correction.c
  * \brief Given an IC on the halo, find a small correction after 90d that ensures it remains inside the LPO region.
  *
  * Given an initial condition that is very close to the nominal halo periodic
  * orbit, find a (tiny) correction 90 days from now that ensures we stay on
  * the LPO region up to 450 more days, effectively shadowing the halo.
  *
  * The correction can be found in several ways. We start using a bisection
  * method: Find an interval (slowly varying the velocity's module \f$|v|\f$)
  * such that at the endpoints we leave from the left/right of the region.
  * Starting with this interval, we perform a bisection procedure until we find 
  * a correction satisfying the 450d condition.
  *
  * NOTES: 
  *	
  * USAGE:	./correction
  *
  * CALLED BY:	
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>				// fmin
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "int_rtbp.h"			// DIM, int_rtbp_exit
#include "correction_module.h"	// correction
#include "utils_module.h"		// dblprint

/** \brief Period of nominal halo orbit. 
 *
 * The period is close to \f$\pi\f$, i.e. approximately 180 days. 
 */
static const double T = 0.3059226605957322E+01;

int
main (int argc, char *argv[])
{
	// Masde passed me the i.c. as x,y,z,xd,yd,zd
	double q_Masde[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00};

	double dv;				/* Modulus of the maneuver */
	double q90[DIM];		/* State after 90 days */
	double q90_new[DIM];	/* Corrected IC after 90 days */

	// Auxiliary variables
	int i;

	/* Use a random number generator to produce variates from a Normal
	 * distribution */
	const gsl_rng_type *rng_type;
	gsl_rng *r;

	/* create a generator chosen by the
    environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();

	rng_type = gsl_rng_default;
	r = gsl_rng_alloc (rng_type);

	/* Initial condition of perturbed orbit: we just add some noise to nominal
	 * i.c. */
	for(i=0; i<6; i++)
		q_Masde[i] += gsl_ran_gaussian(r, q_Masde[i]*1.e-7);

	gsl_rng_free(r);

	dv = correction(q_Masde, 0.5*T, 2.5*T, CORRECTION_ST, q90, q90_new);

	printf("Accepted maneuver:\n");
    printf("dv: %e\n", dv);
    printf("New IC after 90 days: \n");
	dblprint(q90_new, DIM);
    printf("\n");
}
