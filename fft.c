/** \file fft.c
  * \brief Compute Fourier coefficients of \f$f(x)\f$ (REAL Fourier series)
  *
  * Read a curve from file, compute its Fourier coeffs \f$A_n, B_n\f$, and
  * print them to stdout.
  *
  * NOTES: 
  *		We assume that the number of points in curve is 128 (this comes from
  *		halo.py). Thus, the number of real Fourier coeffs is 65.
  *	
  * USAGE:	./fft orbit.dat > coefs.dat
  *
  * CALLED BY:	
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>			// M_PI
#include "fft_module.h"

/** \brief Period of nominal halo orbit */
static const double T = 0.3059226605957322E+01;

int
main (int argc, char *argv[])
{
	FILE *fp;
	const int n = 128;

	double t[n], x[n], px[n], y[n], py[n], z[n], pz[n];	// input data
	double A[n/2+1], B[n/2+1];							// output coefs

	// x(s) is the x-coord of an arbitrary point in the orbit.
	double s, xs;	

	/* auxiliary variables */
	double *tp, *xp, *pxp, *yp, *pyp, *zp, *pzp;
	int i;

	/* Read filename with interpolated curve */
	char *filename;

	if (argc != 2) {
		printf("Usage: %s datafile\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	filename = argv[1];

	/* Read input curve */
    tp = t; 
	xp = x; pxp = px;
	yp = y; pyp = py;
	zp = z; pzp = pz;
    fp = fopen(filename, "r");
	while(fscanf(fp,"%le %le %le %le %le %le %le", tp, xp, pxp, yp, pyp, zp,
				pzp) != EOF) 
	{
        /* advance pointers */
        tp++; 
		xp++; pxp++; 
		yp++; pyp++; 
		zp++; pzp++; 
    }
	fclose(fp);

	/* Backup an arbitrary point in the orbit for later */
	s = t[127];
	xs = x[127];

	/* Perform FFT */
	real_fft(n,x,A,B);

	/* Write coefs to stdout */
	for (i = 0; i <= n/2; i++)
    {
		fprintf (stdout, "%d %e %e\n", i, A[i], B[i]);
	}

	/* Check that we can recover function value x(s) from trig series */
	fprintf(stderr, "x-coord of an arbitrary point in the orbit: %e\n", xs);
	for(i=0; i<10; i++)
		fprintf(stderr, "x-coord recovered by trig series: %e\n",
				eval_trig_series(i, A, B, 2*M_PI*s/T));

	return 0;
}

