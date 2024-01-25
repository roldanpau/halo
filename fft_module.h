void real_fft(int n, double y[n], double A[n/2+1], double B[n/2+1]);
void write_coefs(int N, double A[6][N], double B[6][N]);
void read_coefs(int N, double A[6][N], double B[6][N]);
double eval_trig_series (size_t N, double A[N+1], double B[N+1], double x);
void eval_orbit (int n, double A[6][n], double B[6][n], int deg, double t,
		double p[6]);
