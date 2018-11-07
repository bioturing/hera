#ifndef _UTILITY_
#define _UTILITY_

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "../utils.h"

#define THRESHOLD1 1e-3
#define THRESHOLD2 1e-7

#define LOCK_THRESHOLD 1e-5
#define FILTERING_THRESHOLD 1e-6

//#define VERBOSE
#define PRINT_ROUND

#define LOOP_BEFORE_OFF 10 //mean x3 round em

#define CHECKPOINT_STEP 1

#define EPSILON 1e-300
#define RSPD_BUCKET 64
#define RSPD_LENGTH_THRESHOLD 500

#define MIN_EFFECTIVE_LENGTH 1
//#define EFFECTIVE_LENGTH_ROUND_UP

#define MAX_PATH_LENGTH 1000

#define MAX_READ_LENGTH 300
#define MAX_PAIRED_READ_LENGTH (2*(MAX_READ_LENGTH))

#define MIN_FRAGMENT_LENGTH 0
#define MAX_FRAGMENT_LENGTH 1000

#define ACCELERATE

//#define PRINT_PROGRESS
//#define SAVE_MODEL

#define ZERO_DIV 10

#define max(a,b) ((a) < (b)? (b) : (a))
#define min(a,b) ((a) < (b)? (a) : (b))

#define LEN_CORRECTION
#define RSP_CORRECTION

#define MAX_ALIGNMENT 1000

#define USE_RSPD //1s
#define USE_FRAG // 0.5s
#define USE_MIS //0.3-0.4s
#define USE_ORI //not significant <0.05
#define USE_NOISE 1 //not use noise

//#define ANALYSIS

//base time 1.56s

#if defined(USE_RSPD) && \
    defined(USE_FRAG) && \
    defined(USE_MIS)  && \
    defined(USE_ORI)
#define USE_ALL
#endif

#define get_trans(m, id) ((m)->trans->trans[(id)])
#define trans_cnt(m) ((m)->trans->cnt)

#define __s_check__(s, v) if ((s) != (v)) {fprintf(stderr, "\nRead error\n"); \
						exit(1);}

extern double pow_table[MAX_PAIRED_READ_LENGTH + 1];


extern void accelerate(const double *a, const double *b, double *c, int n);

extern void lock_accelerate(const double *a, const double *b, double *c,
			    const char *lock, int n);

extern void model_global_init();

static inline void normalize(double *a, int n)
{
	int  i;

	double sum = 0;

	for (i = 0; i < n; i++) 
		sum += a[i];

	#ifdef VERBOSE
	fprintf(stderr, "SUM %lf\n", sum);
	#endif

	for (i = 0; i < n; i++)
		a[i] /= sum;
}

static inline void lock_normalize(double *a, const char *lock, int n)
{
	int  i;

	double sum = 0;

	double l_sum = 0;

	for (i = 0; i < n; i++)
		if (!lock[i])
			sum += a[i];
		else
			l_sum += a[i];

	sum = 1/ sum * (1 - l_sum);

	for (i = 0; i < n; i++)
		if (!lock[i])
			a[i] *= sum;
}

static inline void norm_accumulate(double *a, double *b, int n)
{
	int i;
	double sum = 0;

	for (i = 0; i < n; i++)
		sum += a[i];

	assert(!isnan(sum));

	if (sum < EPSILON)
		sum = 1; //all zero

	

	double prev = 0;
	for (i = 0; i < n; i++) {
		a[i] /= sum;
		b[i] = prev += a[i];
		assert(!isnan(a[i]));
	}
}

static inline void accumulate(const double *a, double *b, int n)
{
	int i;
	double prev = 0;
	for (i = 0; i < n; i++)
		b[i] = prev += a[i];
}

static inline void array_add(double *a, const double *b, int n)
{
	int i;
	for (i = 0; i < n; i++)
		a[i] += b[i];
}

struct diff_info {
	double r_change;
	int count1, count2;
};

const struct diff_info arr_diff(const double *a1, const double *a2, int n,
				double threshold1, double threshold2);

void print_diff(const struct diff_info info);
#endif
