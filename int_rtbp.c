/** \file int_rtbp.c
  * \brief Functions to integrate trajectory of the RTBP.
  *
  * NOTES: 
  *
  * CALLED BY:	
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>	//sqrt

#include "rk78.h"

/* routine rtbp below must "see" this value */
static double mu = 3.040357143000000E-06; 

/** \brief Dimension of the RTBP. */
static const size_t DIM = 6;

/** 
  * \brief Integrate trajectory of the RTBP.
  *
  * Numerically integrate the RTBP using a Runge-Kutta 7-8 method with
  * automatic step-size control.
  *
  * @param[in]		total_time		Integration time.
  * @param[inout]   x				On entry, it contains i.c. On exit, f.c.
  * @param[in]		tol				Tolerance to integration error.
  * @param[in]		hmin			Minimum step size.
  * @param[in]		hmax			Maximum step size.
  * @param[in]		bPrint			Flag to print orbit (if bPrint = true)
  *
  * \return		Currently, simply returns 0 (no error control).
  */
int int_rtbp(double total_time, double *x, double tol, double hmin, double
		hmax, int bPrint)
{
   void rtbp(double t, double *x, int n, double *y);
   double t,h;

   ini_rk78(DIM);

   /* integrate initial condition up to time "total_time" */
   t=0.e0;
   h=hmax;

   while(t<total_time)
   {
	   if(bPrint)
	   {
		   printf("%24.16e", t);
		   for(int i=0; i<DIM; i++) printf("%24.16e", x[i]);
		   printf("\n");
	   }
	   if(t+h > total_time)
	   {
		   h = total_time - t;
		   hmin = hmax = h;
	   }
      rk78(&t,x,&h,tol,hmin,hmax,DIM,rtbp);
   }

   fflush(stdout);

   end_rk78(DIM);

   return(0);
}

/** 
  * \brief Integrate trajectory of the RTBP until exiting LPO region.
  *
  * Numerically integrate the RTBP using a Runge-Kutta 7-8 method with
  * automatic step-size control until trajectory exits the \f$L_1\f$ Libration
  * Point Orbit region of the Sun-Earth RTBP.
  *
  * @param[inout]   x				On entry, it contains i.c. On exit, f.c.
  * @param[in]		tol				Tolerance to integration error.
  * @param[in]		hmin			Minimum step size.
  * @param[in]		hmax			Maximum step size.
  * @param[in]		bPrint			Flag to print orbit (if bPrint = true)
  * @param[out]		tout			Exit time
  * @param[out]		bLeft			Flag to specify if orbit leaves through
  * left region (bLeft = true) or right region (bLeft = false).
  *
  * \return		Currently, simply returns 0 (no error control).
  */
int int_rtbp_exit(double *x, double tol, double hmin, double hmax, int bPrint,
		double *tout, int *bLeft)
{
   void rtbp(double t, double *x, int n, double *y);
   double t,h;

	/* x coordinate of L_1 */
   double x1 = -0.99003;

	/* Width of the LPO region (in AU), equivalent to 150*10^6 km */
   double LPO_width = 1.1e+06/150e+06;	

   ini_rk78(DIM);

   /* integrate initial condition up until exiting LPO region */
   t=0.e0;
   h=hmax;

   while(x[0] > x1-LPO_width/2 && x[0] < x1+LPO_width/2)	/* Inside LPO */
   {
	   if(bPrint)
	   {
		   printf("%24.16e", t);
		   for(int i=0; i<DIM; i++) printf("%24.16e", x[i]);
		   printf("\n");
	   }
      rk78(&t,x,&h,tol,hmin,hmax,DIM,rtbp);
   }

   fflush(stdout);

   end_rk78(DIM);

   *tout = t;	/* Return exit time */
   *bLeft = (x[0] <= x1-LPO_width/2 ? 1 : 0);	/* Did we leave through the
												   left? */
   return(0);
}

/*
   Compute the potential \Omega(x,y,z) associated to position (x,y,z).

parameters:
x,y,z: 	position coordinates.
mu: 	mass ratio.

Notes:
x,y,z are assumed to be synodic coordinates (origin at center of mass).
*/
double potential(double x, double y, double z, double mu)
{
   double r1, r2;
   r1 = sqrt((x-mu)*(x-mu) + y*y + z*z);
   r2 = sqrt((x+1-mu)*(x+1-mu) + y*y + z*z);
   return 0.5*(x*x + y*y) +
      (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);
}
void rtbp(double t, double *x, int n, double *y)
/*
this is the vectorfield of the restricted three body problem.

parameters:
t: adimensional time.
x: point in phase space (6 coordinates).
n: dimension of phase space (here it must be 6).
y: vectorfield at (t,x).
*/
{
   double r1,r13,r2,r23,aux3;
   r1=sqrt((x[0]-mu)*(x[0]-mu)+x[2]*x[2]+x[4]*x[4]);
   r13=r1*r1*r1;
   r2=sqrt((x[0]-mu+1)*(x[0]-mu+1)+x[2]*x[2]+x[4]*x[4]);
   r23=r2*r2*r2;
   aux3=(1-mu)/r13+mu/r23;
   y[0]=x[1]+x[2];
   y[1]=x[3]-(1-mu)*(x[0]-mu)/r13-mu*(x[0]-mu+1)/r23;
   y[2]=-x[0]+x[3];
   y[3]=-x[1]-aux3*x[2];
   y[4]=x[5];
   y[5]=-aux3*x[4];
   return;
}
