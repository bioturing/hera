#include "utility.h"
#include <math.h>
#include <stdio.h>

double pow_table[MAX_PAIRED_READ_LENGTH + 1];

void model_global_init()
{
	int i;
	pow_table[0] = 1;
	for (i = 1; i <= MAX_PAIRED_READ_LENGTH; ++i)
		pow_table[i] = pow_table[i-1] / 4;
}


static inline double calc_gamma2(double r_norm2, double v_norm2)
{
	if (v_norm2 <= EPSILON)
		return 1;
	
	double gamma2 = r_norm2 / v_norm2;
	return gamma2;
}

void accelerate(const double *a, const double *b, double *c, int n)
{
	#ifdef ACCELERATE
	int i;

	double r_norm2 = 0;
	double v_norm2 = 0;

	for (i = 0; i < n; i++) {
		double r = b[i] - a[i];
		double v = c[i] - b[i] - r;
		r_norm2 += r * r;
		v_norm2 += v * v;
	}

	double gamma2 = calc_gamma2(r_norm2, v_norm2);
	double gamma = -sqrt(gamma2);

	int zeroing = 0;

	for (i = 0; i < n; i++) {
		double r = b[i] - a[i];
		double v = c[i] - b[i] - r;

		double new_val = a[i] - 2 * gamma * r + gamma2 * v;
		if (new_val <= EPSILON) {
			zeroing += (c[i] > EPSILON);
			new_val = c[i] / ZERO_DIV;
		}

		c[i] = new_val;
	}
	//normalize(c, n);
	#ifdef VERBOSE
	fprintf(stderr, BLU "ZEROING %d\n" RESET, zeroing);
	#endif
	#endif
}

void lock_accelerate(const double *a, const double *b, double *c, const char *lock, int n)
{
	#ifdef ACCELERATE
	int i;

	double r_norm2 = 0;
	double v_norm2 = 0;

	for (i = 0; i < n; i++)
		if (!lock[i]) {
			double r = b[i] - a[i];
			double v = c[i] - b[i] - r;
			r_norm2 += r * r;
			v_norm2 += v * v;
		}

	double gamma2 = calc_gamma2(r_norm2, v_norm2);
	double gamma = -sqrt(gamma2);

	#ifdef VERBOSE
	fprintf(stderr, RED "GAMMA = %lf\n" RESET, gamma);
	#endif

	int zeroing = 0;

	for (i = 0; i < n; i++)
		if (!lock[i]) {
			double r = b[i] - a[i];
			double v = c[i] - b[i] - r;

			double new_val = a[i] - 2 * gamma * r + gamma2 * v;
			if (new_val <= EPSILON) {
				zeroing += (c[i] > EPSILON);
				new_val = c[i] / ZERO_DIV;
			}

			c[i] = new_val;
		}
	//lock_normalize(c, lock, n);
	#ifdef VERBOSE
	fprintf(stderr, BLU "ZEROING %d\n" RESET, zeroing);
	#endif
	#endif
}

const struct diff_info arr_diff(const double *a1, const double *a2, int n, double t1, double t2)
{
	struct diff_info info;
	info.r_change = 0;
	info.count1 = info.count2 = 0;
	int i;
	for (i = 0; i < n; i++) {
		double diff = fabs(a1[i] - a2[i]);
		
		if (a1[i] >= t2) {
			diff /= a1[i];
			info.count1 += diff >= t1;
			info.count2 += diff >= t2;

			if (diff > info.r_change)
				info.r_change = diff;
		}
	}
	return info;
}

void print_diff(const struct diff_info info)
{
	fprintf(stderr, "%16.5lg%10d%10d", info.r_change, info.count1, info.count2);
}