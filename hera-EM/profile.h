#ifndef _PROFILE_
#define _PROFILE_

#include "utility.h"
#include <stdio.h>

struct mis_dis {
	double q;
	double total;
	double p_pow[MAX_PAIRED_READ_LENGTH + 1];
	double q_pow[MAX_PAIRED_READ_LENGTH + 1];
};

static inline double mis_prob(const struct mis_dis *dis, int mismatch, int read_len)
{
	return dis->p_pow[read_len - mismatch] * dis->q_pow[mismatch];
}

static inline void mis_update(struct mis_dis *dis, int mismatch, int length, double prob)
{
	dis->q += prob * mismatch;
	dis->total += prob * length;
}

void mis_uniform_init(struct mis_dis *dis);

void mis_join(struct mis_dis *dis, const struct mis_dis *part);
void mis_finish(struct mis_dis *dis);
void mis_release(struct mis_dis *dis);
void mis_init(struct mis_dis *dis);

void mis_save(const struct mis_dis *dis, FILE *f);
void mis_load(struct mis_dis *dis, FILE *f);

void mis_copy(struct mis_dis *dis, const struct mis_dis *source);

void mis_accelerate(const struct mis_dis *d1, const struct mis_dis *d2, struct mis_dis *d3);
void mis_diff(const struct mis_dis *d1, const struct mis_dis *d2);

#endif