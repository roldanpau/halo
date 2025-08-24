/** \file bandit_module.c
  * \brief Function bandit(a) takes an action and returns corresponding reward.
  *
  * USED BY: bandit.c	
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>		// assert
#include <math.h>		// fabs

#include "int_rtbp.h"			// DIM
#include "section.h"			// section_t
#include "utils_module.h"		// dbldif, norm
#include "correction_module.h"	// int_sec_correction_opt
#include "cv_module.h"			// posvel_to_posmom, posmom_to_posvel
#include "prtbp_module.h"		// prtbp

/** \brief Period of nominal halo orbit.
 *
 * The period is close to \f$\pi\f$, i.e. approximately 180 days.
 */
static const double T = 0.3059226605957322E+01;

// SHADOW_TIME can be up to 3*T (for I.C. at t=0), but needs to be reduced to
// 2*T for I.C. at t=T/4.
// static const double SHADOW_TIME = 3*T;
static const double SHADOW_TIME = 2*T;

static const double w1=0;	//< Weight of 'fuel' cost
static const double w2=1;	//< Weight of 'closeness' cost

static const section_t SEC1 = {{-0.9916647163367744E+00,  0.0000000000000000E+00,
	0.8983543483564242E-03},{0.0,1.0,0.0}}; 
static const section_t SEC2 = {{-9.888427653607404e-01, 0.000000000000000e+00,
		-1.121019578563241e-03},{0.0,-1.0,0.0}}; 

// Section at T/4 (aprox)
static const section_t SEC3 = {{-9.905881e-01, 4.517090e-03, -2.167137e-04},
    {1.0,0.0,0.0}};

// Section at 3T/4 (aprox)
static const section_t SEC4 = {{-9.905881e-01, -4.517090e-03, -2.167137e-04},
    {-1.0,0.0,0.0}};

// Forward declarations
static int bandit1(double q[DIM], double *reward);
static int bandit2(double q[DIM], double *reward);
static int bandit3(double q[DIM], double *reward);
static int bandit4(double q[DIM], double *reward);

/**
 * \brief Take an action and return corresponding reward.
 *
 * Given the current state \f$ q \f$ in the orbit, this function takes action
 * \f$ a \f$, which is one of 4 possible actions: 
 *    - Action k means to apply correction maneuver at section SEC_k, and
 *    integrate q for one period.
 *
 * The reward associated to each action is defined as -cost, where cost is:
 *    - cost = w_1*dv + w_2*|q(T) - q(0)|
 *
 * @param[in,out]	q			Current state  (pos-vel coordinates)
 * @param[in]		sec			Current section (where q is)
 * @param[in]		a			Action (a=k, k=1,2,3,4)
 * @param[out]		reward      Reward
 *
 * \return      Error code 0 (success) or 1 (error performing action a).
 */
int bandit_2(double q[DIM], const section_t sec, int a, double *reward)
{
	double q_posmom[DIM];
	double ti;
	int err;

	fprintf(stderr, "Performing action %d\n", a);
	switch(a) {
		case 1: 
			/* Integrate to SEC1, if necessary */
			if(sec.n[1] != 1.0) {	/* sec != SEC1 */
				fprintf(stderr, "Integrating to SEC1\n");
				posvel_to_posmom(q, q_posmom);
				prtbp(mu, SEC1, 1, q_posmom, &ti);
				posmom_to_posvel(q_posmom, q);
			}
			else fprintf(stderr, "q is already in SEC1\n");

			err = bandit1(q,reward);
			if(err) {
				fprintf(stderr, "Error in function bandit1 !\n");
			}
			else {
				fprintf(stderr, "Reward obtained: %e\n", *reward);
			}
			return err;
			break;
		case 2: 
			/* Integrate to SEC2, if necessary */
			if(sec.n[1] != -1.0) {	/* sec != SEC2 */
				fprintf(stderr, "Integrating to SEC2\n");
				posvel_to_posmom(q, q_posmom);
				prtbp(mu, SEC2, 1, q_posmom, &ti);
				posmom_to_posvel(q_posmom, q);
			}
			else fprintf(stderr, "q is already in SEC2\n");

			err = bandit2(q,reward);
			if(err) {
				fprintf(stderr, "Error in function bandit2 !\n");
			}
			else {
				fprintf(stderr, "Reward obtained: %e\n", *reward);
			}
			return err;
			break;
		case 3: 
			/* Integrate to SEC3, if necessary */
			if(sec.n[0] != 1.0 || sec.n[1] != 0.0 || sec.n[2] != 0.0) {	
				/* sec != SEC3 */
				fprintf(stderr, "Integrating to SEC3\n");
				posvel_to_posmom(q, q_posmom);
				prtbp(mu, SEC3, 1, q_posmom, &ti);
				posmom_to_posvel(q_posmom, q);
			}
			else fprintf(stderr, "q is already in SEC3\n");

			err = bandit3(q,reward);
			if(err) {
				fprintf(stderr, "Error in function bandit3 !\n");
			}
			else {
				fprintf(stderr, "Reward obtained: %e\n", *reward);
			}
			return err;
			break;
		case 4: 
			/* Integrate to SEC4, if necessary */
			if(sec.n[0] != -1.0 || sec.n[1] != 0.0 || sec.n[2] != 0.0) {	
				/* sec != SEC4 */
				fprintf(stderr, "Integrating to SEC4\n");
				posvel_to_posmom(q, q_posmom);
				prtbp(mu, SEC4, 1, q_posmom, &ti);
				posmom_to_posvel(q_posmom, q);
			}
			else fprintf(stderr, "q is already in SEC4\n");

			err = bandit4(q,reward);
			if(err) {
				fprintf(stderr, "Error in function bandit4 !\n");
			}
			else {
				fprintf(stderr, "Reward obtained: %e\n", *reward);
			}
			return err;
			break;
		default: 
			fprintf(stderr, "Error: received unknown action a=%d\n", a);
			exit(EXIT_FAILURE);
	}
}

int bandit1(double q[DIM], double *reward)
{
	double dv;		/* Modulus of the maneuver */
	double q_new[DIM];
	double q_new_posmom[DIM];
	double ti;
	double dif[DIM];	/* |q(T) - halo(T) */
	double cost;

	// auxiliary variables
	double q_bak[DIM];	/* backup of q */
	dblcpy(q_bak, q, DIM);

	/* Obtain optimal correction to q */
	dv = correction_opt(q, SHADOW_TIME, CORRECTION_ST, q_new);
	if(dv == 0)
		return 1;	// Error performing action 1

	//printf("dv=%e\n", dv);

	/* Integrate orbit aprox one period */
	posvel_to_posmom(q_new, q_new_posmom);
	prtbp(mu, SEC1, 1, q_new_posmom, &ti);
	posmom_to_posvel(q_new_posmom, q);

	dbldif(DIM,q,q_bak,dif);
	cost = w1*fabs(dv) + w2*norm(DIM,dif);
	*reward = -cost;
	return 0;
}

int bandit2(double q[DIM], double *reward)
{
	double dv;		/* Modulus of the maneuver */
	double q_new[DIM];
	double q_new_posmom[DIM];
	double ti;
	double dif[DIM];	/* |q(T) - halo(T) */
	double cost;

	// auxiliary variables
	double q_bak[DIM];	/* backup of q */
	dblcpy(q_bak, q, DIM);

	/* Obtain optimal correction to q */
	dv = correction_opt(q, SHADOW_TIME, CORRECTION_ST, q_new);
	if(dv == 0)
		return 1;	// Error performing action 2

	//printf("dv=%e\n", dv);

	/* Integrate orbit aprox one period */
	posvel_to_posmom(q_new, q_new_posmom);
	prtbp(mu, SEC2, 1, q_new_posmom, &ti);
	posmom_to_posvel(q_new_posmom, q);

	dbldif(DIM,q,q_bak,dif);
	cost = w1*fabs(dv) + w2*norm(DIM,dif);
	*reward = -cost;
	return 0;
}

int bandit3(double q[DIM], double *reward)
{
	double dv;		/* Modulus of the maneuver */
	double q_new[DIM];
	double q_new_posmom[DIM];
	double ti;
	double dif[DIM];	/* |q(T) - halo(T) */
	double cost;

	// auxiliary variables
	double q_bak[DIM];	/* backup of q */
	dblcpy(q_bak, q, DIM);

	/* Obtain optimal correction to q */
	dv = correction_opt(q, SHADOW_TIME, CORRECTION_ST, q_new);
	if(dv == 0)
		return 1;	// Error performing action 2

	//printf("dv=%e\n", dv);

	/* Integrate orbit aprox one period */
	posvel_to_posmom(q_new, q_new_posmom);
	prtbp(mu, SEC3, 1, q_new_posmom, &ti);
	posmom_to_posvel(q_new_posmom, q);

	dbldif(DIM,q,q_bak,dif);
	cost = w1*fabs(dv) + w2*norm(DIM,dif);
	*reward = -cost;
	return 0;
}

int bandit4(double q[DIM], double *reward)
{
	double dv;		/* Modulus of the maneuver */
	double q_new[DIM];
	double q_new_posmom[DIM];
	double ti;
	double dif[DIM];	/* |q(T) - halo(T) */
	double cost;

	// auxiliary variables
	double q_bak[DIM];	/* backup of q */
	dblcpy(q_bak, q, DIM);

	/* Obtain optimal correction to q */
	dv = correction_opt(q, SHADOW_TIME, CORRECTION_ST, q_new);
	if(dv == 0)
		return 1;	// Error performing action 2

	//printf("dv=%e\n", dv);

	/* Integrate orbit aprox one period */
	posvel_to_posmom(q_new, q_new_posmom);
	prtbp(mu, SEC4, 1, q_new_posmom, &ti);
	posmom_to_posvel(q_new_posmom, q);

	dbldif(DIM,q,q_bak,dif);
	cost = w1*fabs(dv) + w2*norm(DIM,dif);
	*reward = -cost;
	return 0;
}

