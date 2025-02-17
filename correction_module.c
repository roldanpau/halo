/** \file correction_module.c
  * \brief Given an IC close to the halo, find a small correction after
  * CORREC_TIME (usually 90 days) that ensures it remains inside the LPO region.
  *
  * Given an initial condition that is very close to the nominal halo periodic
  * orbit, find a (tiny) correction CORREC_TIME from now that ensures we stay
  * on the LPO region up to SHADOW_TIME more time, effectively shadowing the
  * halo.
  *
  * The correction can be found in several ways, depending on the flag
  * correction_t: 
  *
  * CORRECTION_VEL: apply a maneuver in the direction of the VELocity vector
  * (vx,vy,vz).
  * 
  * CORRECTION_ST: apply a maneuver in the direction of the STable direction
  * \f$ E^s \f$.
  *
  * USED BY: correction.c	
  *
  */

#include <stdlib.h>				// exit
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "int_rtbp.h"			// DIM, int_rtbp_exit
#include "cv_module.h"			// posvel_to_posmom
#include "utils_module.h"		// norm, dblcpy
#include "correction_module.h"	// correction_t

/**
  * \brief Apply correction dv to q90 and integrate trajectory until exiting LPO region.
  *
  * This function integrates the initial condition \f$
  * q_{90}=(x,y,z,v_x,v_y,v_z) \f$ after applying the maneuver until exiting
  * the LPO region, and returns the integration time and whether it exits
  * trough the left or the right of the region.
  *
  * The type of correction maneuver is specified through flag `correction'.
  *
  * \remark For convenience, the IC after applying the maneuver, \f$ q90 +
  * \Delta v\f$, is returned in q90_new.
  *
  * @param[in]		q90         I.c.  (pos-vel coordinates)
  * @param[in]      dv			Magnitude of the maneuver (may be negative)
  * @param[in]      correction	Type of correction maneuver
  * @param[out]     tout        Exit time
  * @param[out]     bLeft       Flag to specify if orbit leaves through
  * left region (bLeft = true) or right region (bLeft = false).
  * @param[out]		q90_new     q90 after maneuver (pos-vel coordinates)
  *
  * \return     Currently, simply returns 0 (no error control).
  */
int apply_correction(double q90[DIM], double dv, correction_t correction,
		double *tout, int *bLeft, double q90_new[DIM])
{
	double Dv[3];	/* Maneuver (in the direction of current velocity or stable
					   direction, depending on correction flag) */

	/* Auxiliary variables */
	double q[DIM];		// q is a point in the orbit (pos-mom)

	dblcpy(q90_new, q90, DIM);

	if(correction == CORRECTION_VEL)
	{
		double v_mod = norm(3,&q90[3]);		// Modulus of velocity vector
		
		/* Set maneuver to dv*(vx, vy, vz), where (vx, vy, vz) is a unitary
		 * vector tangent to current velocity */
		Dv[0] = dv*q90[3]/v_mod;
		Dv[1] = dv*q90[4]/v_mod;
		Dv[2] = dv*q90[5]/v_mod;

		/* Apply maneuver to q90 */
		q90_new[3] += Dv[0];
		q90_new[4] += Dv[1];
		q90_new[5] += Dv[2];
	}
	else if(correction == CORRECTION_ST)
	{
		/* w = stable direction */
		/* Masde says: Take only first 3 (position) components */
		double w[3] = {0.323988273471917, 0.173290594261497, 0};

		double w_mod = norm(3,w);		// Modulus of stable direction
		
		/* Set maneuver to dv*w, where w is the unitary stable dir */
		Dv[0] = dv*w[0]/w_mod;
		Dv[1] = dv*w[1]/w_mod;
		Dv[2] = dv*w[2]/w_mod;

		/* Apply maneuver to q90 */
		q90_new[3] += Dv[0];
		q90_new[4] += Dv[1];
		q90_new[5] += Dv[2];
	}

    /* My vectorfield expects the order x,px,y,py,z,pz. */
    posvel_to_posmom(q90_new, q);

    /* Integrate orbit until exiting the LPO region, with params
     * tol=1.e-14, hmin=1.e-5, hmax=1.e-1, bPrint=0 */
    int_rtbp_exit(q, 1.e-14, 1.e-5, 1.e-1, 0, tout, bLeft);
}

/** \struct lag_params
  * \brief Parameters to function "lag"
  */
struct lag_params 
{
    double q90[6];  ///< state after 90 days (pos-vel)
	double CORREC_TIME;		///< maneuver is applied after CORREC_TIME
	double SHADOW_TIME;		///< maneuver keeps orbit shadowing halo for SHADOW_TIME extra time
	correction_t correction;	///< type of correction maneuver
};

/**
  * \brief Measure the difference (`lag') between SHADOW_TIME and our exit
  * time tout.
  *
  * Given a maneuver of magnitude dv, this function returns 
  * \f\[ \min\{ SHADOW_TIME - tout, 0\} \f\]
  * multiplied by -1 (if we exit through the left) or by 1 (if we exit through
  * the right. 
  *
  * Zeros of this function correspond to maneuvers dv that lead to an exit time
  * tout greater than SHADOW_TIME. This function is used in the bisection
  * procedure (see function correction()).
  *
  * @param[in]      dv	    Magnitude of the maneuver (may be negative)
  * @param[in]      params	Parameters needed for this function.
  *
  * @returns        lag \f$ \min\{ SHADOW_TIME - tout, 0\} \f$
  */

double lag(double dv, void *params)
{
    struct lag_params *p = (struct lag_params *)params;

	double CORREC_TIME;
	double SHADOW_TIME;
	double tout;	/* exit time */
	int bLeft;		/* Do we leave through the left or right? */

	double lag;		/* return value */

	// auxiliary variables
	double q90_new[DIM];

	CORREC_TIME = p->CORREC_TIME;
	SHADOW_TIME = p->SHADOW_TIME;

    apply_correction(p->q90, dv, p->correction, &tout, &bLeft, q90_new);
    if(tout > SHADOW_TIME)
		lag=0;  /* We found a zero! */
    else 
    {
        if(bLeft) 
			lag=-(SHADOW_TIME - tout);
        else
			lag=SHADOW_TIME - tout;
    }
	return(lag);
}

/**
  * \brief Given an IC close to the halo, find a small correction after
  * CORREC_TIME (usually 90 days) that ensures it remains inside the LPO
  * region.
  *
  * Given an initial condition that is very close to the nominal halo periodic
  * orbit, find a (tiny) correction CORREC_TIME from now that ensures we stay
  * on the LPO region up to SHADOW_TIME more time, effectively shadowing the
  * halo.
  *
  * The correction can be found in several ways. We start using a bisection
  * method: Find an interval (slowly varying the velocity's module \f$|v|\f$)
  * such that at the endpoints we leave from the left/right of the region.
  * Starting with this interval, we perform a bisection procedure until we find 
  * a correction satisfying the SHADOW_TIME condition.
  *
  * \remark For convenience, the state after applying the maneuver CORREC_TIME
  * from now, \f$ q90 + \Delta v\f$, is returned in q90_new.
  *
  * @param[in]		q_Masde		Initial condition (pos-vel coordinates)
  *
  * @param[in]		CORREC_TIME		Apply correction after CORREC_TIME time
  * @param[in]		SHADOW_TIME		Extra time (after CORREC_TIME) we must
  * ensure we stay shadowing the halo.
  *
  * @param[in]      corr   		Type of correction maneuver
  * @param[out]		q90         State after CORREC_TIME (pos-vel)
  * @param[out]		q90_new     q90 after maneuver (pos-vel coordinates)
  *
  * @returns	dv	Modulus of correction maneuver
  */

double correction(double q_Masde[DIM], double CORREC_TIME, double SHADOW_TIME,
		correction_t corr, double q90[DIM], double q90_new[DIM])
{
	double q[DIM];		// q is a point in the orbit (pos-mom)
	double v_mod;		// Modulus of velocity vector
	
	double tout;	/* exit time */
	int bLeftA;	/* Do we leave through the left or right? (endpoint A) */
	int bLeftB;	/* Do we leave through the left or right? (endpoint B) */

	double dv;		/* Modulus of the maneuver */

	/* My vectorfield expects the order x,px,y,py,z,pz. */
	posvel_to_posmom(q_Masde, q);

	/* Integrate orbit for GOLDEN_FRACT*T time, approximately 69 days */
	int_rtbp(CORREC_TIME, q, 1.e-15, 1.e-5, 1.e-1, 0);

	posmom_to_posvel(q, q90);

	/*
	fprintf(stderr, "q90: %e %e %e %e %e %e\n", 
			q90[0], q90[1], q90[2], q90[3], q90[4], q90[5]); 
			*/

	/* Find out which side (left/right) we are exiting when using maneuver dv=0
	 * */
	apply_correction(q90, 0.0, corr, &tout, &bLeftA, q90_new);

	/*
	fprintf(stderr, "Exit time: %e\n", tout);
	if(bLeftA)
		fprintf(stderr, "Exiting through the left\n");
	else
		fprintf(stderr, "Exiting through the right\n");
		*/

    if(tout > SHADOW_TIME)
    {
        //fprintf(stderr, "We stay in LPO region longer than 450 days! OK!\n");
		return(0);
    }

    /* If we reach this part of the code, it means that a correction maneuver
     * is needed */

    /* Try corrections with positive modulus */
	for(dv=1.e-7; dv<1.e-3; dv+=1.e-7)
	{
		apply_correction(q90, dv, corr, &tout, &bLeftB, q90_new);

		if(bLeftB != bLeftA)
			break;
	}

    if(bLeftB == bLeftA)   
    /* Could not find endpoint exiting through the other side */
    {
        /* Try corrections with negative modulus */
        for(dv=-1.e-7; dv>-1.e-3; dv-=1.e-7)
        {
			apply_correction(q90, dv, corr, &tout, &bLeftB, q90_new);

            if(bLeftB != bLeftA)
                break;
        }

    }

	/*
	fprintf(stderr, "dv: %e\n", dv);
	fprintf(stderr, "Exit time: %e\n", tout);
	if(bLeftB)
		fprintf(stderr, "Exiting through the left\n");
	else
		fprintf(stderr, "Exiting through the right\n");
		*/

    if(tout > SHADOW_TIME)	/* No need to apply bisection */
    {
        //fprintf(stderr, "We stay in LPO region longer than 450 days! OK!\n");
		return(dv);
    }

	if(bLeftB == bLeftA)	/* Unrecoverable error */
	{
		fprintf(stderr, "Could not find two orbits that exit through "
				"sides of LPO region?? Exiting.\n");
		exit(EXIT_FAILURE);
	}

	/* Apply bisection procedure */

    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *Type;
    gsl_root_fsolver *s;
    double root = 0.0;
    double x_lo = fmin(0.0, dv);
    double x_hi = fmax(0.0, dv);
    gsl_function F;
    struct lag_params params;
    double lag_r;   /* lag for the current root */

    dblcpy(params.q90, q90, DIM);
	params.CORREC_TIME = CORREC_TIME;
	params.SHADOW_TIME = SHADOW_TIME;
	params.correction = corr;

    F.function = &lag;
    F.params = &params;

    Type = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (Type);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	/*
    printf ("using %s method\n", gsl_root_fsolver_name (s));

    printf ("%5s [%13s, %13s] %13s %13s\n", "iter", "lower", "upper", "root",
            "err(est)");
			*/

    do 
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      root = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      lag_r = lag(root, &params);
      status = gsl_root_test_residual (lag_r, 0.1);

	  /*
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.11f, %.11f] %.11f %.11f\n",
              iter, x_lo, x_hi, root, lag_r);
			  */
    }
    while (status == GSL_CONTINUE && iter < max_iter);

	if(iter == max_iter)	/* Unrecoverable error */
	{
		fprintf(stderr, "Maximum number of iterations reached "
				"in bisection procedure? Exiting.\n");
		exit(EXIT_FAILURE);
	}

	/* Correction maneuver has been found */
	dv = root;

	/* Set q90_new before returning */
	apply_correction(q90, dv, corr, &tout, &bLeftB, q90_new);
	return(dv);
}

/** \struct tout_params
  * \brief Parameters to function "tout"
  */
struct tout_params 
{
    double q[6];  ///< initial condition (pos-vel)
	correction_t correction;	///< type of correction maneuver
};

/**
  * \brief Measure (the inverse of) our exit time tout.
  *
  * Given a maneuver of magnitude dv, this function returns 
  * \f\[ \frac{1}{tout} \f\]
  * multiplied by -1 (if we exit through the left) or by 1 (if we exit through
  * the right. 
  *
  * Zeros of this function correspond to maneuvers dv that lead to an infinite
  * exit time tout. This function is used in the optimal bisection procedure
  * (see function correction_opt()).
  *
  * @param[in]      dv	    Magnitude of the maneuver (may be negative)
  * @param[in]      params	Parameters needed for this function.
  *
  * @returns        \f$ \frac{1}{tout} \f$
  */

double tout_f(double dv, void *params)
{
    struct tout_params*p = (struct tout_params *)params;

	double tout;	/* exit time */
	int bLeft;		/* Do we leave through the left or right? */

	double ret;		/* return value */

	// auxiliary variables
	double q_new[DIM];

    apply_correction(p->q, dv, p->correction, &tout, &bLeft, q_new);

	if(bLeft) 
		ret = -1.0/tout;
	else
		ret = 1.0/tout;

	return(ret);
}

/**
  * \brief Given an IC close to the halo, find a small correction that ensures it remains inside the LPO region.
  *
  * Given an initial condition that is very close to the nominal halo periodic
  * orbit, find a (tiny) correction that ensures we stay on the LPO region for
  * as long as possible, effectively shadowing the halo.
  *
  * We use a bisection method: Find an interval (slowly varying the velocity's
  * module \f$|v|\f$) such that at the endpoints we leave from the left/right
  * of the region.
  * Starting with this interval, we perform a bisection procedure until we find
  * the optimal correction that stays in the LPO region for as long as
  * possible.
  *
  * In practice, to find the optimal one, we iterate the bisection procedure
  * until the interval \f$ [dv_-, dv_+] \f$ is practically zero length.
  *
  * \remark For convenience, the state after applying the maneuver, \f$ q +
  * \Delta v\f$, is returned in q_new.
  *
  * \remark If the maneuver is not able to guarantee that orbit remains inside
  * the LPO region for at least SHADOW_TIME, the procedure is considered
  * unsuccessful, and we return dv=0. We also leave the orbit unmodified, i.e.
  * we return q_new=q.
  *
  * @param[in]		q			Initial condition (pos-vel coordinates)
  * @param[in]		SHADOW_TIME	Time we must ensure we stay shadowing the halo.
  * @param[in]      corr   		Type of correction maneuver
  * @param[out]		q_new		q after maneuver (pos-vel coordinates)
  *
  * @returns	dv	Modulus of correction maneuver (or 0, if procedure fails).
  */

double correction_opt(double q[DIM], double SHADOW_TIME, correction_t corr,
		double q_new[DIM])
{
	double v_mod;		// Modulus of velocity vector
	
	double tout;	/* exit time */
	int bLeftA;	/* Do we leave through the left or right? (endpoint A) */
	int bLeftB;	/* Do we leave through the left or right? (endpoint B) */

	double dv;		/* Modulus of the maneuver */

	/*
	fprintf(stderr, "q: %e %e %e %e %e %e\n", 
			q[0], q[1], q[2], q[3], q[4], q[5]); 
			*/

	/* Find out which side (left/right) we are exiting when using maneuver dv=0
	 * */
	apply_correction(q, 0.0, corr, &tout, &bLeftA, q_new);

	/*
	fprintf(stderr, "Exit time: %e\n", tout);
	if(bLeftA)
		fprintf(stderr, "Exiting through the left\n");
	else
		fprintf(stderr, "Exiting through the right\n");
		*/

	/* A correction maneuver is always applied */

    /* Try corrections with positive modulus */
	for(dv=1.e-7; dv<1.e-3; dv+=1.e-7)
	{
		apply_correction(q, dv, corr, &tout, &bLeftB, q_new);

		if(bLeftB != bLeftA)
			break;
	}

    if(bLeftB == bLeftA)   
    /* Could not find endpoint exiting through the other side */
    {
        /* Try corrections with negative modulus */
        for(dv=-1.e-7; dv>-1.e-3; dv-=1.e-7)
        {
			apply_correction(q, dv, corr, &tout, &bLeftB, q_new);

            if(bLeftB != bLeftA)
                break;
        }

    }

	/*
	fprintf(stderr, "dv: %e\n", dv);
	fprintf(stderr, "Exit time: %e\n", tout);
	if(bLeftB)
		fprintf(stderr, "Exiting through the left\n");
	else
		fprintf(stderr, "Exiting through the right\n");
		*/

	if(bLeftB == bLeftA)	/* Unrecoverable error */
	{
		fprintf(stderr, "Could not find two orbits that exit through "
				"sides of LPO region?? Exiting.\n");
		exit(EXIT_FAILURE);
	}

	/* Apply bisection procedure */

    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *Type;
    gsl_root_fsolver *s;
    double root = 0.0;
    double x_lo = fmin(0.0, dv);
    double x_hi = fmax(0.0, dv);
    gsl_function F;
    struct tout_params params;
    double tout_r;   /* 1/tout for the current root */

    dblcpy(params.q, q, DIM);
	params.correction = corr;

    F.function = &tout_f;
    F.params = &params;

    Type = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (Type);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	/*
    printf ("using %s method\n", gsl_root_fsolver_name (s));

    printf ("%5s [%13s, %13s] %13s %13s\n", "iter", "lower", "upper", "root",
            "err(est)");
			*/

    do 
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      root = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      tout_r = tout_f(root, &params);
      status = gsl_root_test_interval (x_lo, x_hi, 1.e-16, 1.e-5);

	  /*
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.11f, %.11f] %.11f %.11f\n",
              iter, x_lo, x_hi, root, lag_r);
			  */
    }
    while (status == GSL_CONTINUE && iter < max_iter);

	if(iter == max_iter)	/* We did our best */
	{
		fprintf(stderr, "Maximum number of iterations reached "
				"in bisection procedure. Returning last maneuver.\n");
	}

	/* Correction maneuver has been found */
	dv = root;

	/* Set q_new before returning */
	apply_correction(q, dv, corr, &tout, &bLeftB, q_new);

	if(tout < SHADOW_TIME)
	{
		fprintf(stderr, "tout = %f is smaller than SHADOW_TIME = %f!."
				"Returning dv=0\n", tout, SHADOW_TIME);
		dv=0;
		dblcpy(q_new, q, DIM);
	}

	return(dv);
}

/**
  * \brief Given an IC close to the halo, integrate it for CORREC_TIME and find a small correction that ensures it remains inside the LPO region.
  *
  * Given an initial condition that is very close to the nominal halo periodic
  * orbit, integrate it for time CORREC_TIME (usually around 90 days), and find
  * a (tiny) correction that ensures we stay on the LPO region for as long as
  * possible, effectively shadowing the halo.
  *
  * We use a bisection method: Find an interval (slowly varying the velocity's
  * module \f$|v|\f$) such that at the endpoints we leave from the left/right
  * of the region.
  * Starting with this interval, we perform a bisection procedure until we find
  * the optimal correction that stays in the LPO region for as long as
  * possible.
  *
  * In practice, to find the optimal one, we iterate the bisection procedure
  * until the interval \f$ [dv_-, dv_+] \f$ is practically zero length.
  *
  * \remark For convenience, the state after applying the maneuver CORREC_TIME
  * from now, \f$ q90 + \Delta v\f$, is returned in q90_new.
  *
  * \remark If the maneuver is not able to guarantee that orbit remains inside
  * the LPO region for at least SHADOW_TIME, the procedure is considered
  * unsuccessful, and we return dv=0. We also leave the orbit unmodified, i.e.
  * we return q90_new=q90=q_Masde.
  *
  * @param[in]		q_Masde		Initial condition (pos-vel coordinates)
  *
  * @param[in]		CORREC_TIME		Apply correction after CORREC_TIME time
  * @param[in]		SHADOW_TIME		Extra time (after CORREC_TIME) we must
  * ensure we stay shadowing the halo.
  *
  * @param[in]      corr   		Type of correction maneuver
  * @param[out]		q90         State after CORREC_TIME (pos-vel)
  * @param[out]		q90_new     q90 after maneuver (pos-vel coordinates)
  *
  * @returns	dv	Modulus of correction maneuver (or 0, if procedure fails).
  */

double int_correction_opt(double q_Masde[DIM], double CORREC_TIME, double
		SHADOW_TIME, correction_t corr, double q90[DIM], double q90_new[DIM])
{
	double q[DIM];		// q is a point in the orbit (pos-mom)
	double v_mod;		// Modulus of velocity vector
	
	double tout;	/* exit time */
	int bLeftA;	/* Do we leave through the left or right? (endpoint A) */
	int bLeftB;	/* Do we leave through the left or right? (endpoint B) */

	double dv;		/* Modulus of the maneuver */

	/* My vectorfield expects the order x,px,y,py,z,pz. */
	posvel_to_posmom(q_Masde, q);

	/* Integrate orbit for CORREC_TIME time */
	int_rtbp(CORREC_TIME, q, 1.e-15, 1.e-5, 1.e-1, 0);

	posmom_to_posvel(q, q90);

	/* Obtain optimal correction to q90 */
	dv = correction_opt(q90, SHADOW_TIME, corr, q90_new);

	return(dv);
}
