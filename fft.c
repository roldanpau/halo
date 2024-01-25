/** \file fft.c
  * \brief Compute Fourier coefficients of \f$orbit(t)\f$ (REAL Fourier series)
  *
  * Read orbit from file. For each coordinate (\f$ x, p_x, y, p_y, z, p_z)
  * \f$), compute its Fourier coeffs \f$A_n, B_n\f$, and print them to stdout.
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
#include "fft_module.h"

int
main (int argc, char *argv[])
{
	FILE *fp;
	const int n = 128;

	double t[n];		// times where orbit is evaluated

	// orbit data: 
	// orbit[0][:] contains x(t_0), x(t_1), ..., x(t_n)
	// orbit[1][:] contains p_x(t_0), p_x(t_1), ..., p_x(t_n)
	// ...
	// orbit[4][:] contains z(t_0), z(t_1), ..., z(t_n)
	// orbit[5][:] contains p_z(t_0), p_z(t_1), ..., p_z(t_n)
	double orbit[6][n];	

	// Fourier coefs:
	// A[0][:] contains coefs A_0, A_1, ..., A_N of function x(t)
	// A[1][:] contains coefs A_0, A_1, ..., A_N of function p_x(t)
	// ...
	// A[4][:] contains coefs A_0, A_1, ..., A_N of function z(t)
	// A[5][:] contains coefs A_0, A_1, ..., A_N of function p_z(t)
	double A[6][n/2+1]; 
	double B[6][n/2+1];							

	// p1(s) is the original point in the orbit
	double s; 
	double p1[6];	

	// p2(s) is the Fourier-recovered point in the orbit
	double p2[6];

	/* auxiliary variables */
	double *tp, *xp, *pxp, *yp, *pyp, *zp, *pzp;
	int i, d, deg;

	/* Read filename with interpolated curve */
	char *filename;

	if (argc != 2) {
		printf("Usage: %s datafile\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	filename = argv[1];

	/* Read input curve */
    tp = t; 
	xp = &orbit[0][0]; pxp = &orbit[1][0];
	yp = &orbit[2][0]; pyp = &orbit[3][0];
	zp = &orbit[4][0]; pzp = &orbit[5][0];
    fp = fopen(filename, "r");
	while(fscanf(fp,"%le %le %le %le %le %le %le", tp, xp, pxp, yp, pyp, zp,
				pzp) != EOF) 
	{
        /* advance pointers */
		tp++; xp++; pxp++; yp++, pyp++, zp++, pzp++;
    }
	fclose(fp);

	/* Backup an arbitrary point in the orbit for later */
	s = t[127];
	p1[0] = orbit[0][127];
	p1[1] = orbit[1][127];
	p1[2] = orbit[2][127];
	p1[3] = orbit[3][127];
	p1[4] = orbit[4][127];
	p1[5] = orbit[5][127];

	/* Perform FFT */
	real_fft(n,orbit[0],A[0],B[0]);
	real_fft(n,orbit[1],A[1],B[1]);
	real_fft(n,orbit[2],A[2],B[2]);
	real_fft(n,orbit[3],A[3],B[3]);
	real_fft(n,orbit[4],A[4],B[4]);
	real_fft(n,orbit[5],A[5],B[5]);

	/* Write coefs to stdout */
	write_coefs(n/2+1, A, B);

	/* Check that we can recover point p(s) from trig series */
	fprintf(stderr, "Arbitrary point in the orbit: %e %e %e %e %e %e\n", p1[0],
			p1[1], p1[2], p1[3], p1[4], p1[5]);
	for(deg=0; deg<10; deg++)
	{
		eval_orbit(n/2+1, A, B, deg, s, p2);
		fprintf(stderr, "Recovered point in the orbit: %e %e %e %e %e %e\n", p2[0],
				p2[1], p2[2], p2[3], p2[4], p2[5]);
	}

	return 0;
}

