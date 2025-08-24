/** \file bandit2.c
  * \brief A simple bandit algorithm.
  *
  * Consider the following learning problem. You are faced repeatedly with a
  * choice among 4 different options, or actions: 
  *		- Action k means to apply correction maneuver at section SEC_k.
  *
  * After each choice you receive a numerical reward=-cost chosen from a
  * stationary distribution that depends on the action you selected:
  *		- cost = w_1*dv + w_2*|q(T) - q(0)|
  * where w_1, w_2 are weigths in [0,1].
  *
  * Your objective is to maximize the expected total reward over some time
  * period, for example, over NSEL action selections
  *
  * USAGE:	./bandit2 > fuel_value2
  * and then plot with bandit2.gpl.
  *
  */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>	// gsl_ran_gaussian

#include "section.h"
#include "bandit2_module.h"

const int SEED=2;
const int NSEL = 5000;	//< Number of action selections

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

int
main (int argc, char *argv[])
{
	const int NACT = 4;	//< Number of actions (2)

	// eps = how infrequently we want to explore random actions
	// 1-eps = how greedyly we want to exploit valuable actions
	const double EPS = 0.1;	

	// Halo at t=0, given as x,y,z,xd,yd,zd
	section_t sec = SEC1;
	double q[6] = {-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00};

	int a;	// action
	int n;	// number of actions taken so far

	double Q[NACT];		// Estimated value of each action
	int N[NACT];		// Number of times each action has been selected

    /* Use a random number generator to produce variates from a Normal
     * distribution */
	const gsl_rng_type *T;
	gsl_rng *r;

	int err;			// Return code of bandit function
	double reward;		// Reward after action
	double randnum;		// Random number between 0 and 1

	// Initialize estimated action values
	for(a=0; a<NACT; a++) {
		Q[a] = 0;
		N[a] = 0;
	}

    /* create a generator chosen by the
    environment variable GSL_RNG_TYPE */
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	// Set seed. Ch1nge this to any value > 1 to generate a different stream
	// of random numbers. Default seed = 0.
	gsl_rng_set(r, SEED);

	// Loop for NSEL action selections
	for(n=0; n<NSEL; n++) {
		double u = gsl_rng_uniform(r);	// produce uniform rand num in [0,1)

		// With prob 1-eps, do this
		if(u>EPS) {
			// Action a=argmax Q(a)
			if(Q[0]>Q[1] && Q[0]>Q[2] && Q[0]>Q[3]) a=0;
			else if(Q[1]>Q[0] && Q[1]>Q[2] && Q[1]>Q[3]) a=1;
			else if(Q[2]>Q[0] && Q[2]>Q[1] && Q[2]>Q[3]) a=2;
			else if(Q[3]>Q[0] && Q[3]>Q[1] && Q[3]>Q[2]) a=3;
			else {	
				// At least two actions with same value. Break ties
				// randomly. This would be very unlikely...
				randnum=gsl_rng_uniform(r);
				if(randnum<0.25) a=0;
				else if(0.25<randnum && randnum<0.5) a=1;
				else if(0.5<randnum && randnum<0.75) a=2;
				else a=3;
			}
			fprintf(stderr, "Exploiting action %d...\n", a+1);
		}
		// With prob eps, do this
		else	// Take a random action
		{
			fprintf(stderr, "Exploration: Taking random action...\n");
			randnum=gsl_rng_uniform(r);
			if(randnum<0.25) a=0;
			else if(0.25<randnum && randnum<0.5) a=1;
			else if(0.5<randnum && randnum<0.75) a=2;
			else a=3;
		}

		/* To simulate uncertain observations, we try to add some noise to
		 * nominal i.c.
		for(int i=0; i<6; i++)
			q[i] += gsl_ran_gaussian(r, q[i]*1.e-5); */
		
		// Perform action and collect reward
		err = bandit_2(q, sec, a+1, &reward);
		if(err) {
			fprintf(stderr, "Error in bandit function!\n");
			exit(EXIT_FAILURE);
		}
		// update sec
		switch(a) {
			case 0: 
				sec=SEC1;
				break;
			case 1: 
				sec=SEC2;
				break;
			case 2:
				sec=SEC3;
				break;
			case 3:
				sec=SEC4;
				break;
			default:
				fprintf(stderr, "Error updating sec!\n");
		}

		// Skip first 5 iterates, because initial condition is practically on
		// the halo orbit
		if(n<5) continue;

		// Update action values
		N[a] = N[a] + 1;
		Q[a] = Q[a] + 1/((double)N[a])*(reward-Q[a]);

		// Print action values
		/*
		fprintf(stderr, "Estimated action values: \n");
		fprintf(stderr, "Q[1] = %e, ", Q[0]);
		fprintf(stderr, "Q[2] = %e, ", Q[1]);
		fprintf(stderr, "Q[3] = %e, ", Q[2]);
		fprintf(stderr, "Q[4] = %e\n\n", Q[3]);
		*/

		// Output: besides action values Q, we also print updated states q.
		// These states are used in nominal_vs_corrected_orbits.gpl to check
		// that the corrected orbit does not drift away from nominal orbit too
		// much.
		printf("%e %e %e %e %e %e %e\n", Q[0], Q[1], Q[2], Q[3], q[0], q[1],
				q[2]);
	}

	gsl_rng_free(r);
}
