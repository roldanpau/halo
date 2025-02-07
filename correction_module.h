static const double GOLDEN_FRACT = 0.381966;	///< Fraction of a circle occupied by the golden angle
	
/** \enum correction_t
  * \brief Types of correction maneuver.
  *
  * The correction can be found in several ways, depending on the flag
  * correction_t: 
  *
  * CORRECTION_VEL: apply a maneuver in the direction of the velocity vector
  * (vx,vy,vz).
  * 
  * CORRECTION_ST: apply a maneuver in the direction of the stable direction
  * \f$ E^s \f$.
  */
typedef enum {CORRECTION_VEL, CORRECTION_ST} correction_t;

double lag(double dv, void *params);
double correction(double q_Masde[DIM], double CORREC_TIME, double SHADOW_TIME,
		correction_t corr, double q90[DIM], double q90_new[DIM]);

double tout(double dv, void *params);
double correction_opt(double q_Masde[DIM], double CORREC_TIME, double
		SHADOW_TIME, correction_t corr, int bNow, double q90[DIM], double q90_new[DIM]);
