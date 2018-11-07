#include "profile.h"
#include <string.h>
#include <assert.h>
#include <math.h>

void mis_init(struct mis_dis *dis)
{
	memset(dis, 0, sizeof(struct mis_dis));
}

void to_prob(struct mis_dis *dis, double q)
{
	dis->q = q;
	int i;

	double p = 1 - q;
	dis->p_pow[0] = dis->q_pow[0] = 1;

	for (i = 1; i <= MAX_PAIRED_READ_LENGTH; i++) {
		dis->p_pow[i] = dis->p_pow[i - 1] * p;                                
		dis->q_pow[i] = dis->q_pow[i - 1] * q;
	}
}

void mis_uniform_init(struct mis_dis *dis)
{
	/* assume error rate is 50% */
	to_prob(dis, 0.5);
	dis->total = 0;
}

void mis_join(struct mis_dis *dis, const struct mis_dis *part)
{
	dis->q += part->q;
	dis->total += part->total;
}

void mis_finish(struct mis_dis *dis)
{
	dis->q /= dis->total;
	to_prob(dis, dis->q);
}

void mis_copy(struct mis_dis *dis, const struct mis_dis *source)
{
	memcpy(dis, source, sizeof(struct mis_dis));
}

void mis_release(struct mis_dis *dis)
{
}

void mis_save(const struct mis_dis *dis, FILE *f)
{
	fprintf(f, "%lg\n", dis->q);
}

void mis_load(struct mis_dis *dis, FILE *f)
{
	double q;
	__s_check__(fscanf(f, "%lf", &q), 1);
	to_prob(dis, q);
}

void mis_accelerate(const struct mis_dis *d1, const struct mis_dis *d2, struct mis_dis *d3)
{
}

void mis_diff(const struct mis_dis *d1, const struct mis_dis *d2)
{
	fprintf(stderr, "Error rate: %lf (%+lg)", d2->q, d2->q - d1->q);
};