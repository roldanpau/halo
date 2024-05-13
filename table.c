/** \file table.c
  * \brief Generate table containing training examples for supervised learning
  *
  * Generate a table of labelled examples of the form
  * x y z dx dy dz dv
  * associating to each state (t,x,y,z,dx,dy,dz) the correction dv (magnitude
  * applied in the direction of current velocity).
  * 
  * We generate these examples by constructing a trajectory that ``shadows'' an
  * LPO orbit by making tiny adjustments in velocity every GOLDEN_FRACT*T. The
  * frequency GOLDEN_FRACT is chosen so that the generated examples cover the
  * whole periodic orbit.
  *
  * USAGE:	./table > maneuvers.csv
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>				// M_PI

#include "int_rtbp.h"			// DIM
#include "correction_module.h"	// correction
#include "utils_module.h"		// dblprint

/** \brief Period of nominal halo orbit. 
 *
 * The period is close to \f$\pi\f$, i.e. approximately 180 days. 
 */
static const double T = 0.3059226605957322E+01;

static const double CORREC_TIME = 2*GOLDEN_FRACT*T;
static const double SHADOW_TIME = 4.25*T;

//static const double twentyYrs = 20*2*M_PI;	///< 20 yrs (in normalized units)
static const double twentyYrs = 600*2*M_PI;

int
main (int argc, char *argv[])
{
	// I.c. close to the nominal halo orbit, given as x,y,z,xd,yd,zd
	double X1[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00};

	double dv_vel;			/* Modulus of the CORRECTION_VEL maneuver */
	double dv_st;			/* Modulus of the CORRECTION_ST maneuver */
	double q90[DIM];		/* State after 90 days */
	double X2[DIM];			/* Corrected IC after 90 days */

	printf("t,x,y,z,dx,dy,dz,dv_st\n");	/* Table header */

	double time=0.0;		/* Total time of extended orbit (in LPO region) */

	while(time < 3*CORREC_TIME)
	// The first 3 corrections are not printed, since we are initially on the
	// halo and they are not significant.
	{	
		//dv_vel = correction(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_VEL, q90, X2);
		dv_st = correction(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_ST, q90, X2);

		dblcpy(X1, X2, DIM);

		time += CORREC_TIME;
	}
	while(time < twentyYrs)
	{
		//dv_vel = correction(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_VEL, q90, X2);
		dv_st = correction(X1, CORREC_TIME, SHADOW_TIME, CORRECTION_ST, q90, X2);

		printf("%e,%e,%e,%e,%e,%e,%e,%e\n", time, q90[0], q90[1], q90[2],
				q90[3], q90[4], q90[5], dv_st);

		dblcpy(X1, X2, DIM);

		time += CORREC_TIME;
	}
}
