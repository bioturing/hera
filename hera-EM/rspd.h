#ifndef _RSPD_
#define _RSPD_

#include "utility.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

/*read start position distribution*/
struct rsp_dis
{
	double pdf[RSPD_BUCKET + 1], cdf[RSPD_BUCKET + 1];
	double x[RSPD_BUCKET + 1], y[RSPD_BUCKET + 1];
};

static inline double rsp_cdf(const struct rsp_dis *dis, int x, double inv) {
	double frac = x * inv;

	int b = frac;

	double over = frac - b;
	double prob =  dis->cdf[b] + over * dis->pdf[b + 1];

	return prob;
}

static inline double rsp_prob(const struct rsp_dis *dis, int pos, int f_len, int total_len)
{
	double inv = (double) RSPD_BUCKET / total_len;
	double pdf;

	/* divided for number of position in a bucket */
	pdf = dis->pdf[pos] * inv;

	double cdf = rsp_cdf(dis, total_len - f_len + 1, inv);
	return pdf / cdf;
}

void rsp_init(struct rsp_dis *dis);
void rsp_uniform_init(struct rsp_dis *dis);

static inline void rsp_update(struct rsp_dis *dis, struct rsp_dis *old_dis, int x, int f_len, int total_len, double prob)
{
	dis->x[x] += prob;

	#ifdef RSP_CORRECTION
	double inv = (double) RSPD_BUCKET / total_len;
	double meo = rsp_cdf(old_dis, total_len - f_len + 1, inv);
	int k = (total_len - f_len) * inv + 1;
	dis->y[k] += prob / meo;
	#endif
}

void rsp_join(struct rsp_dis *dis, const struct rsp_dis *part);
void rsp_finish(struct rsp_dis *dis);

void rsp_copy(struct rsp_dis *dis, const struct rsp_dis *source);

void rsp_release(struct rsp_dis *dis);

void rsp_save(const struct rsp_dis *dis, FILE *f);
void rsp_load(struct rsp_dis *dis, FILE *f);

void rsp_accelerate(const struct rsp_dis *d1, const struct rsp_dis *d2, struct rsp_dis *d3);

void rspd_diff(const struct rsp_dis *d1, const struct rsp_dis *d2);

#endif