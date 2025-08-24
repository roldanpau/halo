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

static const double w1=1;	//< Weight of 'fuel' cost
static const double w2=0;	//< Weight of 'closeness' cost

static const section_t SEC1 = {{-0.9916647163367744E+00,  0.0000000000000000E+00,
	0.8983543483564242E-03},{0.0,1.0,0.0}}; 
static const section_t SEC2 = {{-9.888427653607404e-01, 0.000000000000000e+00,
		-1.121019578563241e-03},{0.0,-1.0,0.0}}; 

// Forward declarations
static int bandit1(double q[DIM], double *reward);
static int bandit2(double q[DIM], double *reward);

/**
 * \brief Take an action and return corresponding reward.
 *
 * Given the current state \f$ q \f$ in the orbit, which must be in section
 * SEC1, this function takes action \f$ a \f$, which is one of two possible
 * actions: 
 *	- a=1 means apply one correction maneuver per period (at t=0)
 *	- a=2 means apply two correction maneuvers per period (at t=0 and t=T/2)
 *
 * The reward associated to each action is defined as -cost, where cost is:
 *	- When a=1, cost = |dv| + |q(T) - halo(T)|.
 *	- When a=2, cost = |dv_1 + dv_2| + 1/2*(|q(T/2) - halo(T/2)| + |q(T) -
 *	  halo(T)|).
 *
 * @param[in,out]	q			I.c.  (pos-vel coordinates)
 * @param[in]		a			Action (a=1 or a=2)
 * @param[out]		reward      Reward
 *
 * \return      Error code 0 (success) or 1 (error performing action a).
 */
int bandit(double q[DIM], int a, double *reward)
{
	// Make sure q is in SEC1
	double y=q[1];
	double vy=q[4];
	assert(y==0 && vy>0);	

	//printf("Performing action %d\n", a);
	switch(a) {
		case 1: 
			return bandit1(q,reward);
			break;
		case 2: 
			return bandit2(q,reward);
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

	/* Obtain optimal correction to q */
	dblcpy(q_bak, q, DIM);
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
	double dif[DIM];	/* |q(T/2) - halo(T/2) */
	double cost=0;

	// auxiliary variables
	double q_bak[DIM];	/* backup of q */

	/* Obtain optimal correction to q */
	dblcpy(q_bak, q, DIM);
	dv = correction_opt(q, SHADOW_TIME, CORRECTION_ST, q_new);
	if(dv == 0)
		return 1;	// Error performing action 1

	//printf("dv1=%e, ", dv);

	/* Integrate orbit aprox half period */
	posvel_to_posmom(q_new, q_new_posmom);
	prtbp(mu, SEC2, 1, q_new_posmom, &ti);
	posmom_to_posvel(q_new_posmom, q);

	//dbldif(DIM,q,q_bak,dif);
	//cost += (w1*fabs(dv) + w2*0.5*norm(DIM,dif));
	cost += w1*fabs(dv);

	/* Obtain optimal correction to q */
	dv = correction_opt(q, SHADOW_TIME, CORRECTION_ST, q_new);
	if(dv == 0)
		return 1;	// Error performing action 1

	//printf("dv2=%e\n", dv);

	posvel_to_posmom(q_new, q_new_posmom);
	prtbp(mu, SEC1, 1, q_new_posmom, &ti);
	posmom_to_posvel(q_new_posmom, q);

	dbldif(DIM,q,q_bak,dif);
	cost += (w1*fabs(dv) + w2*norm(DIM,dif));

	*reward = -cost;
	return 0;
}

