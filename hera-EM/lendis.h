#ifndef _LENDIS_
#define _LENDIS_

#include "utility.h"
#include <stdio.h>
#include <assert.h>

#define len_span(d) ((d)->max_l - (d)->min_l + 1)

struct len_dis
{
	int min_l;
	int max_l;

	double *pdf, *cdf;

	double *x;
	double *y;
};

static inline double len_prob(const struct len_dis *dis, int len, int trans_len)
{
	double denom = dis->cdf[min(dis->max_l, trans_len)];

	if (denom <= EPSILON)
		return 0;
		
	return dis->pdf[len] / denom;
}

void len_init(struct len_dis *dis, int min_l, int max_l);
void len_uniform_init(struct len_dis *dis, int min_l, int max_l);
void len_normal_init(struct len_dis *dis, int min_l, int max_l, double mean, double sd);


static inline void len_update(struct len_dis *dis, int len, int trans_len, double prob)
{
	dis->x[len] += prob;

	#ifdef LEN_CORRECTION
	if (trans_len <= dis->max_l)
		dis->y[trans_len] += prob;
	#endif
}

void len_join(struct len_dis *dis, const struct len_dis *part);
void len_finish(struct len_dis *dis);
void len_release(struct len_dis *dis);

void len_save(const struct len_dis *dis, FILE *f);
void len_load(struct len_dis *dis, FILE *f);

void len_copy(struct len_dis *dis, const struct len_dis *source);

void len_accelerate(const struct len_dis *d1, const struct len_dis *d2, struct len_dis *d3);

void len_diff(const struct len_dis *d1, const struct len_dis *d2);

void len_move(struct len_dis *d1, struct len_dis *d2);
#endif