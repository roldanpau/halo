#ifndef INT_RTBP_H_INCLUDED
#define INT_RTBP_H_INCLUDED

#include <stdio.h>	/* size_t */

/* routine rtbp below must "see" this value */
extern const double mu; 

/** \brief Dimension of the RTBP. */
extern const size_t DIM;

int int_rtbp(double total_time, double *x, double tol, double hmin, double
		hmax, int bPrint);
int int_rtbp_exit(double *x, double tol, double hmin, double hmax, int bPrint,
		double *tout, int *bLeft);

#endif // INT_RTBP_H_INCLUDED
