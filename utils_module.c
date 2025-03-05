/*! 
  \file utils_module.c
  \brief Utility functions.
  */

#include <stdlib.h>	// malloc
#include <string.h>	// memcpy
#include <stdio.h>	// printf
#include <math.h>	// M_PI, floor

const double TWOPI = 2*M_PI;

double *dblcpy(double * dst, double const * src, size_t len)
{
   memcpy(dst, src, len * sizeof(double));
   return dst;
}

void dblprint(double const *x, size_t len)
{
   int i;
   for(i=0; i<len; i++)
      printf("%.15le ", x[i]);
}

// Floating-point modulo
// The result (the remainder) has same sign as the divisor.
// Similar to matlab's mod(); 
// Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
double Mod(double x, double y)
{
    if(y == 0)
        return x;

    double m= x - y * floor(x/y);
    return m;
}

// wrap [rad] angle to [-pi,pi)
double WrapPosNegPI(double fAng)
{
    return Mod(fAng + M_PI, TWOPI) - M_PI;
}

// wrap [rad] angle to [0..2\pi)
double WrapTwoPI(double fAng)
{
    return Mod(fAng, TWOPI);
}

/**
 * \brief Given vectors \f$ v1, v2 \f$, compute their difference v3 = v1-v2.
 *
 * @param[in] n			dimension of vectors v1,v2,v3
 * @param[in] v1		vector1
 * @param[in] v2		vector2
 * @param[out] v3		vector3
 */
void dbldif(int n, double v1[const], double v2[const], double v3[])
{
	for(int i=0; i<n; i++)
		v3[i] = v1[i] - v2[i];
}

/**
 * \brief Given vector \f$ v \f$, compute its Euclidean norm
 *
 * @param[in] n     dimension of vector
 * @param[in] v     vector
 *
 * @returns			Euclidean norm of vector
 */
double norm(int n, double v[const])
{
	double res = 0;
	for(int i=0; i<n; i++)
		res += (v[i]*v[i]);
	return sqrt(res);
}

/**
 * \brief Given vectors \f$ v_1, v_2 \f$, compute their dot product.
 *
 * @param[in] n     dimension of vectors
 * @param[in] v1     vector
 * @param[in] v2     vector
 *
 * @returns			Dot product v1*v2.
 */
double dot(int n, double v1[const], double v2[const])
{
	double res = 0;
	for(int i=0; i<n; i++)
		res += (v1[i]*v2[i]);
	return res;
}
