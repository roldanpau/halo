#include <stdio.h>	/* size_t */

/* routine rtbp below must "see" this value */
static const double mu = 3.040357143000000E-06; 

/** \brief Dimension of the RTBP. */
static const size_t DIM = 6;

int int_rtbp(double total_time, double *x, double tol, double hmin, double
		hmax, int bPrint);
int int_rtbp_exit(double *x, double tol, double hmin, double hmax, int bPrint,
		double *tout, int *bLeft);
