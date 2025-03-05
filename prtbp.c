// =================================================
// Poincare map of the Restricted Three Body Problem
// =================================================
// FILE:          $RCSfile: prtbp_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:39:05 $
//
// PURPOSE
// =======
// Let "sec" be a Poincare section in the RTBP, where
//      - sec = SEC1 means section {y=0, p_y>0}
//      - sec = SEC2 means section {y=0, p_y<0}.
// This program computes the n-th iterate of the Poincare map, $P^n(x)$.
//
// NOTES
// =====
// We do not check (and thus we do not impose) that the initial point "x" is
// on the Poincare section.
//
// OVERALL METHOD:
//
// 1. Input mass parameter, type of section "sec", number of iterates "n" and
// initial point "x" from stdin (in pos-vel).
// 2. Compute n-th iterate of the Poincare map, $P^n(x)$.
// 3. Output final point $P^n(x)$ (in pos-vel) and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // strcmp

#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off

#include "int_rtbp.h"	// DIM
#include "section.h"	// section_t
#include "prtbp_module.h"	// prtbp
#include "cv_module.h"	// posvel_to_posmom, posmom_to_posvel

static const section_t SEC1 = {{-0.9916647163367744E+00,  0.0000000000000000E+00,
	0.8983543483564242E-03},{0.0,1.0,0.0}}; 
static const section_t SEC2 = {{-9.888427653607404e-01, 0.000000000000000e+00,
		-1.121019578563241e-03},{0.0,-1.0,0.0}}; 

int main( )
{
   double mu;
   section_t sec;
   double ti;		// integration time
   double x[DIM];
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc
   double x_posmom[DIM];		// x in position-momentum

   // Input mass parameter, section, number of iterates and initial condition from stdin.
   if(scanf("%le %s %d %le %le %le %le %le %le", &mu, section_str, &n, x, x+1, x+2,
			   x+3, x+4, x+5) < 9)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   posvel_to_posmom(x,x_posmom);

   // Compute n-th iterate of Poincare map, $P^n(x)$.
   status=prtbp(mu,sec,n,x_posmom,&ti);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   posmom_to_posvel(x_posmom, x);

   // Output final point and integration time to stdout.
   status = printf("%.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", \
	 x[0], x[1], x[2], x[3], x[4], x[5], ti);
   if(status<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
