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
#include "int_rtbp.h"			// int_rtbp_exit
#include "utils_module.h"		// norm, dblcpy
#include "cv_module.h"			// posvel_to_posmom

/** \brief Dimension of the RTBP. */
static const size_t DIM = 6;

/** \brief Period of nominal halo orbit. 
 *
 * The period is close to \f$\pi\f$, i.e. approximately 180 days. 
 */
static const double T = 0.3059226605957322E+01;

struct lag_params 
{
    double q90[6];  // state after 90 days (pos-vel)
};

/**
  * \brief Measure the difference (`lag') between 450 days (=2.5T) and our exit
  * time tout.
  *
  * Given a maneuver of magnitude dv, this function returns 
  * \f\[ \min\{ 2.5T - tout, 0} \f\]
  * multiplied by -1 (if we exit through the left) or by 1 (if we exit through
  * the right. 
  *
  * Zeros of this function correspond to maneuvers dv that lead to an exit time
  * tout greater than 2.5T. This function is used in the bisection procedure
  * (see correction.c).
  *
  * @param[in]      dv	    Magnitude of the maneuver (may be negative)
  * @param[in]      params	Parameters needed for this function (q90).
  *
  * @returns        lag \f$ \min\{ 2.5T - tout, 0} \f$
  */
double lag(double dv, void *params)
{
    struct lag_params *p = (struct lag_params *)params;

	double q[DIM];		// q is a point in the orbit (pos-mom)
	double q90[DIM];	// q90: state after 90 days (pos-vel)
	double q90_new[DIM];	// q90 after maneuver
	double v_mod;		// Modulus of velocity vector
	
	double tout;	/* exit time */
	int bLeft;	/* Do we leave through the left or right? */

	double Dv[3];	/* Maneuver (in the direction of current velocity) */

    dblcpy(q90, p->q90, DIM);

	v_mod = norm(3,&q90[3]);	// modulus of velocity vector

    dblcpy(q90_new, q90, DIM);

    /* Set maneuver to dv*(vx, vy, vz), where (vx, vy, vz) is a unitary
     * vector tangent to current velocity */
    Dv[0] = dv*q90[3]/v_mod;
    Dv[1] = dv*q90[4]/v_mod;
    Dv[2] = dv*q90[5]/v_mod;

    /* Apply maneuver to q90 */
    q90_new[3] += Dv[0];
    q90_new[4] += Dv[1];
    q90_new[5] += Dv[2];

    /* My vectorfield expects the order x,px,y,py,z,pz. */
    posvel_to_posmom(q90_new, q);

    /* Integrate orbit until exiting the LPO region, with params
     * tol=1.e-14, hmin=1.e-5, hmax=1.e-2, bPrint=0 */
    int_rtbp_exit(q, 1.e-14, 1.e-5, 1.e-1, 0, &tout, &bLeft);

    if(tout > 2.5*T)
    {
        return(0);  /* We found a zero! */
    }
    else 
    {
        if(bLeft) 
            return(-(2*T - tout));
        else
            return(2*T - tout);
    }
}

int
main (int argc, char *argv[])
{
	
	double q[DIM];		// q is a point in the orbit (pos-mom)
	double q90[DIM];	// q90: state after 90 days (pos-vel)
	double q90_new[DIM];	// q90 after maneuver
	double v_mod;		// Modulus of velocity vector
	
	double tout;	/* exit time */
	int bLeftA;	/* Do we leave through the left or right? (endpoint A) */
	int bLeftB;	/* Do we leave through the left or right? (endpoint B) */

	double Dv[3];	/* Maneuver (in the direction of current velocity) */
	double dv;		/* Modulus of the maneuver */

	// Auxiliary variables
	int i;

	// Masde passed me the i.c. as x,y,z,xd,yd,zd
	double q_Masde[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00};

	/* Use a random number generator to produce variates from a Normal
	 * distribution */
	const gsl_rng_type *rng_type;
	gsl_rng *r;

	/* create a generator chosen by the
    environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();

	rng_type = gsl_rng_default;
	r = gsl_rng_alloc (rng_type);

	/* My vectorfield expects the order x,px,y,py,z,pz. */
	posvel_to_posmom(q_Masde, q);

	/* Initial condition of perturbed orbit: we just add some noise to nominal
	 * i.c. */
	for(i=0; i<6; i++)
		q[i] += gsl_ran_gaussian(r, q[i]*1.e-10);

	/* Integrate orbit for 90 days, or approximately half the period */
	int_rtbp(0.5*T, q, 1.e-15, 1.e-5, 1.e-1, 0);

	posmom_to_posvel(q, q90);

	fprintf(stderr, "q90: %e %e %e %e %e %e\n", 
			q90[0], q90[1], q90[2], q90[3], q90[4], q90[5]);

	/* Integrate orbit until exiting the LPO region, with params tol=1.e-14,
	 * hmin=1.e-5, hmax=1.e-2, bPrint=1 */
    int_rtbp_exit(q, 1.e-14, 1.e-5, 1.e-1, 0, &tout, &bLeftA);

	fprintf(stderr, "Exit time: %e\n", tout);
	if(bLeftA)
		fprintf(stderr, "Exiting through the left\n");
	else
		fprintf(stderr, "Exiting through the right\n");

    if(tout > 2.5*T)
    {
        fprintf(stderr, "We stay in LPO region longer than 450 days! OK!\n");
        exit(0);
    }

    /* If we reach this part of the code, it means that a correction maneuver
     * is needed */

	v_mod = norm(3,&q90[3]);	// modulus of velocity vector

    /* Try corrections with positive modulus */
	for(dv=1.e-10; dv<1.e-5; dv+=1.e-10)
	{
		dblcpy(q90_new, q90, DIM);

		/* Set maneuver to dv*(vx, vy, vz), where (vx, vy, vz) is a unitary
		 * vector tangent to current velocity */
		Dv[0] = dv*q90[3]/v_mod;
		Dv[1] = dv*q90[4]/v_mod;
		Dv[2] = dv*q90[5]/v_mod;

		/* Apply maneuver to q90 */
		q90_new[3] += Dv[0];
		q90_new[4] += Dv[1];
		q90_new[5] += Dv[2];

		/* My vectorfield expects the order x,px,y,py,z,pz. */
		posvel_to_posmom(q90_new, q);

		/* Integrate orbit until exiting the LPO region, with params
		 * tol=1.e-14, hmin=1.e-5, hmax=1.e-2, bPrint=1 */
		int_rtbp_exit(q, 1.e-14, 1.e-5, 1.e-1, 0, &tout, &bLeftB);

		if(bLeftB != bLeftA)
			break;
	}

    if(bLeftB == bLeftA)   
    /* Could not find endpoint exiting through the other side */
    {
        /* Try corrections with negative modulus */
        for(dv=-1.e-10; dv>-1.e-5; dv-=1.e-10)
        {
            dblcpy(q90_new, q90, DIM);

            /* Set maneuver to dv*(vx, vy, vz), where (vx, vy, vz) is a unitary
             * vector tangent to current velocity */
            Dv[0] = dv*q90[3]/v_mod;
            Dv[1] = dv*q90[4]/v_mod;
            Dv[2] = dv*q90[5]/v_mod;

            /* Apply maneuver to q90 */
            q90_new[3] += Dv[0];
            q90_new[4] += Dv[1];
            q90_new[5] += Dv[2];

            /* My vectorfield expects the order x,px,y,py,z,pz. */
            posvel_to_posmom(q90_new, q);

            /* Integrate orbit until exiting the LPO region, with params
             * tol=1.e-14, hmin=1.e-5, hmax=1.e-2, bPrint=1 */
            int_rtbp_exit(q, 1.e-14, 1.e-5, 1.e-1, 0, &tout, &bLeftB);

            if(bLeftB != bLeftA)
                break;
        }

    }

	fprintf(stderr, "dv: %e\n", dv);
	fprintf(stderr, "Dv: %e %e %e\n", Dv[0], Dv[1], Dv[2]);

	fprintf(stderr, "Exit time: %e\n", tout);
	if(bLeftB)
		fprintf(stderr, "Exiting through the left\n");
	else
		fprintf(stderr, "Exiting through the right\n");

	gsl_rng_free(r);
}
