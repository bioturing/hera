#ifndef _ORIDIS_
#define _ORIDIS_

#include <stdio.h>
#include <assert.h>

struct ori_dis
{
	double p[2]; /* = {0.5, 0.5}; */
};

static inline double ori_pdf(const struct ori_dis *dis, int dir)
{
	return dis->p[dir];
}

void ori_init(struct ori_dis *dis);
void ori_uniform_init(struct ori_dis *dis);

static inline void ori_update(struct ori_dis *dis, int dir, double prob)
{
	dis->p[dir] += prob;
}

void ori_join(struct ori_dis *dis, const struct ori_dis *part);
void ori_finish(struct ori_dis *dis);
void ori_release(struct ori_dis *dis);

void ori_save(const struct ori_dis *dis, FILE *f);
void ori_load(struct ori_dis *dis, FILE *f);

void ori_copy(struct ori_dis *dis, const struct ori_dis *source);

void ori_accelerate(const struct ori_dis *d1, const struct ori_dis *d2, struct ori_dis *d3);


#endif