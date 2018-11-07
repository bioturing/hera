#ifndef _TRANSDIS_
#define _TRANSDIS_

#include <stdio.h>
#include <assert.h>

struct trans_dis {
	int trans_cnt;
	double *p; /*use p[-1] as noise*/
	char *lock; //is transcript locked
};

static inline double trans_pdf(const struct trans_dis *dis, int trans_id)
{
	return dis->p[trans_id];
}

void trans_init(struct trans_dis *dis, struct trans_dis *dis0);
void trans_master_init(struct trans_dis *dis, struct trans_dis *dis0);

void trans_uniform_init(struct trans_dis *dis, int trans_cnt);

static inline void trans_update(struct trans_dis *dis, int trans_id, double prob)
{
	dis->p[trans_id] += prob;
}

char *trans_locking(struct trans_dis *dis);

void trans_join(struct trans_dis *dis, const struct trans_dis *part);
void trans_finish(struct trans_dis *dis);
void trans_release(struct trans_dis* dis);
void trans_save(const struct trans_dis *dis, FILE *f);
void trans_load(struct trans_dis *dis, FILE *f);

void trans_copy(struct trans_dis *dis, const struct trans_dis *source);

static inline void ntrans_update(struct trans_dis * dis, double prob)
{
	dis->p[-1] += prob;
}

static inline double ntrans_pdf(const struct trans_dis * dis)
{
	return dis->p[-1];
}

void trans_accelerate(const struct trans_dis *d1, const struct trans_dis *d2, struct trans_dis *d3);

void trans_diff(const struct trans_dis *d1, const struct trans_dis *d2);
int is_stop(const struct trans_dis *d1, const struct trans_dis *d2);

void trans_move(struct trans_dis *d1, struct trans_dis *d2);

#endif