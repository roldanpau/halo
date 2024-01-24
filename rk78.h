void ini_rk78(int n);
void end_rk78(int n);
double rk78(double *at, double x[], double *ah, double tol,
            double hmin, double hmax, int n,
            void (*deriv)(double, double *, int, double *));
