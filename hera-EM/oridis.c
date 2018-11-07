#include "oridis.h"
#include "utility.h"
#include <assert.h>


void ori_init(struct ori_dis *dis)
{
	dis->p[0] = dis->p[1] = 0;
}

void ori_uniform_init(struct ori_dis *dis)
{
	dis->p[0] = dis->p[1] = 0.5;
}

void ori_join(struct ori_dis *dis, const struct ori_dis *part)
{
	array_add(dis->p, part->p, 2);
}

void ori_finish(struct ori_dis *dis)
{
	double total = dis->p[0] + dis->p[1];

	dis->p[0] /= total;
	dis->p[1] /= total;
}

void ori_copy(struct ori_dis *dis, const struct ori_dis *source)
{
	dis->p[0] = source->p[0];
	dis->p[1] = source->p[1];
}

void ori_release(struct ori_dis *dis)
{
	
}

void ori_save(const struct ori_dis *dis, FILE *f)
{
	fprintf(f, "%lg\n", dis->p[0]);
		
}

void ori_load(struct ori_dis *dis, FILE *f)
{
	__s_check__(fscanf(f, "%lf", &dis->p[0]), 1);
	dis->p[1] = 1 - dis->p[0];
}

void ori_accelerate(const struct ori_dis *d1, const struct ori_dis *d2, struct ori_dis *d3)
{

}