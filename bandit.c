/** \file bandit.c
  * \brief A simple bandit algorithm.
  *
  * Consider the following learning problem. You are faced repeatedly with a
  * choice among 2 different options, or actions: 
  *		- Action 1 means to apply correction maneuvers once per period.
  *		- Action 2 means to apply correction maneuvers twice per period.
  *
  * After each choice you receive a numerical reward=-cost chosen from a
  * stationary distribution that depends on the action you selected:
  *		- When action=1,  cost = dv + |q(T) - halo(T)|
  *		- When action=2,  cost = dv_1 + dv_2 + 1/2*(|q(T/2) - halo(T/2)| +
  *		|q(T) - halo(T)|)
  *
  * Your objective is to maximize the expected total reward over some time
  * period, for example, over 1000 action selections
  *
  * USAGE:	./bandit > fuel_value
  * and then plot with bandit.gpl.
  *
  */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>	// gsl_ran_gaussian

#include "bandit_module.h"

int
main (int argc, char *argv[])
{
	const int NACT = 2;	//< Number of actions (2)

	// eps = how infrequently we want to explore random actions
	// 1-eps = how greedyly we want to exploit valuable actions
	const double EPS = 0.2;	

	// Halo at t=0, given as x,y,z,xd,yd,zd
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

	// Loop for 1000 action selections
	for(n=0; n<1000; n++) {
		double u = gsl_rng_uniform(r);	// produce uniform rand num in [0,1)

		// With prob 1-eps, do this
		if(u>EPS) {
			// Action a=argmax Q(a)
			if(Q[0]>Q[1]) a=0;
			else if(Q[1]>Q[0]) a=1;
			else 	// break ties randomly
				a = (gsl_rng_uniform(r)<0.5 ? 0:1);
		}
		// With prob eps, do this
		else	// Take a random action
		{
			//printf("Taking random action... ");
			a = (gsl_rng_uniform(r)<0.5 ? 0:1);
		}

		/* To simulate uncertain observations, we try to add some noise to
		 * nominal i.c.
		for(int i=0; i<6; i++)
			q[i] += gsl_ran_gaussian(r, q[i]*1.e-5); */
		
		// Perform action and collect reward
		err = bandit(q, a+1, &reward);
		if(err) {
			fprintf(stderr, "Error in bandit function!\n");
			exit(EXIT_FAILURE);
		}
		//printf("Reward = %e\n", reward);

		// Update action values
		N[a] = N[a] + 1;
		Q[a] = Q[a] + 1/((double)N[a])*(reward-Q[a]);

		// Print action values
		//printf("Estimated action values: \n");
		//printf("Q[1] = %e, ", Q[0]);
		//printf("Q[2] = %e\n", Q[1]);
		printf("%e %e\n", Q[0], Q[1]);
	}

	gsl_rng_free(r);
}
