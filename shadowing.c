/** \file shadowing.c
  * \brief Construct a trajectory that ``shadows'' an LPO orbit.
  *
  * Construct a trajectory that ``shadows'' an LPO orbit for 20 years by making
  * tiny adjustments in velocity every half period (equivalent to roughly 90
  * days).
  *
  * We extend the time interval within the LPO region by applying a procedure
  * similar to that described in [Masdemont et al., Global analysis of direct
  * transfers..., 2021], Section 3.3. This correction procedure is implemented
  * in correction_module.c
  *
  * USAGE:	./shadowing
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>				// M_PI

#include "int_rtbp.h"			// DIM
#include "correction_module.h"	// correction, correction_opt
#include "utils_module.h"		// dblprint

/** \brief Period of nominal halo orbit. 
 *
 * The period is close to \f$\pi\f$, i.e. approximately 180 days. 
 */
static const double T = 0.3059226605957322E+01;

static const double twentyYrs = 20*2*M_PI;	///< 20 yrs (in normalized units)

int
main (int argc, char *argv[])
{
	// I.c. close to the nominal halo orbit, given as x,y,z,xd,yd,zd
	double X1[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00};

	double dv;				/* Modulus of the maneuver */
	double q90[DIM];		/* State after 90 days */
	double X2[DIM];			/* Corrected IC after 90 days */

	printf("IC: \n");
	dblprint(X1, DIM);
	printf("\n");

	double CORREC_TIME = T;
	double SHADOW_TIME = 3*T;

	double time=0.0;		/* Total time of extended orbit (in LPO region) */
	while(time < twentyYrs)
	{
		printf("Years: %f\n", time/(2*M_PI));
		/*
		dv = correction(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_VEL, q90, X2);

		fprintf(stderr, "CORRECTION_VEL Accepted maneuver dv: %e\n", dv);

		dv = correction(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_ST, q90, X2);

		fprintf(stderr, "CORRECTION_ST  Accepted maneuver dv: %e\n\n", dv);
		*/

		dv = int_correction_opt(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_ST,
				q90, X2);
		if(dv==0)
		{
			fprintf(stderr, "Correction procedure failed! Exiting.");
			exit(EXIT_FAILURE);
		}

		fprintf(stderr, "OPTIMAL CORRECTION_ST  Accepted maneuver dv: %e\n\n",
				dv);

		printf("New IC after 90 days: \n");
		dblprint(X2, DIM);
		printf("\n\n");

		dblcpy(X1, X2, DIM);

		time += CORREC_TIME;
	}
}
