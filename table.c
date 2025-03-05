/** \file table.c
  * \brief Generate table containing training examples for supervised learning
  *
  * Generate a table of labelled examples of the form
  * x y z dx dy dz dv
  * associating to each state (t,x,y,z,dx,dy,dz) the correction dv (magnitude
  * applied in the direction of current velocity).
  * 
  * We generate these examples by constructing a trajectory that ``shadows'' an
  * LPO orbit by making tiny adjustments in velocity every time the orbit
  * passes through Poincare section {y=0, vy>0}.
  *
  * USAGE:	./table > maneuvers.csv
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>				// M_PI

#include "int_rtbp.h"			// DIM
#include "section.h"			// section_t
#include "correction_module.h"	// int_sec_correction_opt
#include "utils_module.h"		// dblprint

/** \brief Period of nominal halo orbit. 
 *
 * The period is close to \f$\pi\f$, i.e. approximately 180 days. 
 */
static const double T = 0.3059226605957322E+01;

// SHADOW_TIME can be up to 3*T (for I.C. at t=0), but needs to be reduced to
// 2*T for I.C. at t=T/4.
// static const double SHADOW_TIME = 3*T;
static const double SHADOW_TIME = 2*T;
 
//static const double twentyYrs = 20*2*M_PI;	///< 20 yrs (in normalized units)
static const double twentyYrs = 1000*2*M_PI;

static const section_t SEC1 = {{-0.9916647163367744E+00,  0.0000000000000000E+00,
	0.8983543483564242E-03},{0.0,1.0,0.0}}; 

int
main (int argc, char *argv[])
{
	// I.c. on the nominal halo orbit, given as x,y,z,xd,yd,zd
	double X1[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
	0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
	0.0000000000000000E+00};

	// I.c. at QUARTER PERIOD (T/4) of the nominal halo orbit
	// (x,px,y,py,z,pz) = -9.905881e-01 -1.456353e-03 4.517090e-03 -9.909699e-01 -2.167137e-04 -2.007546e-03
	// In (x,y,z, dx=px+y, dy=py-x, dz) coords:
	//double X1[6] = {-9.905881e-01, 4.517090e-03, -2.167137e-04,
	//	-1.456353e-03+4.517090e-03, -9.909699e-01-(-9.905881e-01),
	//	-2.007546e-03};
 
	// I.c. at T/8 of the nominal halo orbit
	// (x,px,y,py,z,pz) = -9.914575e-01 -2.013749e-03 3.327175e-03 -9.850071e-01 5.350106e-04 -1.703465e-03
	// In (x,y,z, dx=px+y, dy=py-x, dz) coords:
	//double X1[6] = {-9.914575e-01, 3.327175e-03, 5.350106e-04,
	//	-2.013749e-03+3.327175e-03, -9.850071e-01-(-9.914575e-01),
	//	-1.703465e-03};
 
	section_t sec = SEC1;	/* Integrate to Poincare section 1 */

	double dv_vel;			/* Modulus of the CORRECTION_VEL maneuver */
	double dv_st;			/* Modulus of the CORRECTION_ST maneuver */
	double q90[DIM];		/* State after 90 days */
	double X2[DIM];			/* Corrected IC after 90 days */

	printf("t,x,y,z,dx,dy,dz,dv_st\n");	/* Table header */

	double time=0.0;		/* Total time of extended orbit (in LPO region) */
	double ti;				/* integration time to section */

	for(int i=0; i<3; i++)
	// The first 3 corrections are not printed, since we are initially on the
	// halo and they are not significant.
	{	
		dv_st = int_sec_correction_opt(X1, sec, SHADOW_TIME, CORRECTION_ST,
				q90, X2, &ti);
		if(dv_st==0)
		{
			fprintf(stderr, "Correction procedure failed! Exiting.");
			exit(EXIT_FAILURE);
		}

		dblcpy(X1, X2, DIM);

		time += ti;
	}
	while(time < twentyYrs)
	{
		dv_st = int_sec_correction_opt(X1, sec, SHADOW_TIME, CORRECTION_ST,
				q90, X2, &ti);
		if(dv_st==0)
		{
			fprintf(stderr, "Correction procedure failed! Exiting.");
			exit(EXIT_FAILURE);
		}

		printf("%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n", time, q90[0], q90[1], q90[2],
				q90[3], q90[4], q90[5], dv_st);

		dblcpy(X1, X2, DIM);

		time += ti;
	}
}
