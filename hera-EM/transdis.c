#include "transdis.h"
#include <stdlib.h>
#include "utility.h"
#include <math.h>
#include <string.h>

void trans_alloc(struct trans_dis *dis, int trans_cnt)
{
	dis->trans_cnt = trans_cnt;
	dis->p = ((double*)calloc((trans_cnt + 1), sizeof(double))) + 1;
	dis->lock = NULL;
}

void trans_init(struct trans_dis *dis, struct trans_dis *dis0)
{
	trans_alloc(dis, dis0->trans_cnt);
	dis->lock = dis0->lock; //mem leak :)) whatever :)*
}

void trans_master_init(struct trans_dis *dis, struct trans_dis *dis0)
{
	trans_init(dis, dis0);
	int i;

	if (dis->lock) {
		for (i = -1; i < dis->trans_cnt; i++)
			if (dis->lock[i]) {
				dis->p[i] = dis0->p[i];
			}
	}
}

void trans_uniform_init(struct trans_dis *dis, int trans_cnt)
{
	trans_alloc(dis, trans_cnt);
	double prob = 1.0 / (trans_cnt + 1);

	int i;
	for (i = -1; i < trans_cnt; i++)
		dis->p[i] = prob;
}

char *trans_locking(struct trans_dis *dis)
{
	return dis->lock = (char*)calloc((dis->trans_cnt + 1),
			   sizeof(char)) + 1;
}

void trans_join(struct trans_dis *dis, const struct trans_dis *part)
{
	array_add(dis->p - 1, part->p - 1, dis->trans_cnt + 1);
}

void trans_finish(struct trans_dis *dis)
{
	if (dis->lock)
		lock_normalize(dis->p - 1, dis->lock - 1, dis->trans_cnt + 1);
	else
		normalize(dis->p - 1, dis->trans_cnt + 1);
	
}

void trans_release(struct trans_dis* dis)
{
	if (dis->p)
		free(dis->p - 1);
}

void trans_save(const struct trans_dis *dis, FILE *f)
{
	fprintf(f, "%d\n", dis->trans_cnt);
	
	int i;
	for (i = -1; i < dis->trans_cnt; i++)
		fprintf(f, "%lg\n", dis->p[i]);
}

void trans_load(struct trans_dis *dis, FILE *f)
{
	int trans_cnt;
	__s_check__(fscanf(f, "%d", &trans_cnt), 1);

	trans_uniform_init(dis, trans_cnt); //lazy to code

	int i;
	for (i = -1; i < trans_cnt; i++)
		__s_check__(fscanf(f, "%lf\n", &dis->p[i]), 1);
}

void trans_copy(struct trans_dis *dis, const struct trans_dis *source)
{
	dis->trans_cnt = source->trans_cnt;
	memcpy(dis->p - 1, source->p - 1,
		sizeof(double) * (dis->trans_cnt + 1));
}

void trans_accelerate(const struct trans_dis *d1, const struct trans_dis *d2,
			struct trans_dis *d3)
{
	if (d3->lock) {
		lock_accelerate(d1->p - 1, d2->p - 1, d3->p - 1, d3->lock - 1,
				d1->trans_cnt + 1);
		
		int i;
		for (i = -1; i < d3->trans_cnt; i++) {
			if (!d3->lock[i]) {
				double diff = fabs(d3->p[i] - d1->p[i]);
				d3->lock[i] = d3->p[i] < EPSILON ||
					      diff < LOCK_THRESHOLD * d3->p[i];
			}
		}
	}
	else
		accelerate(d1->p - 1, d2->p - 1, d3->p - 1, d1->trans_cnt + 1);
}

void trans_diff(const struct trans_dis *d1, const struct trans_dis *d2)
{
	print_diff(arr_diff(d1->p - 1, d2->p - 1, d1->trans_cnt + 1,
		THRESHOLD1, THRESHOLD2));
}

int is_stop(const struct trans_dis *d1, const struct trans_dis *d2)
{
	return arr_diff(d1->p - 1,
			d2->p - 1,
			d1->trans_cnt + 1,
			THRESHOLD1, THRESHOLD2
		).count1 == 0;
}

void trans_move(struct trans_dis *d1, struct trans_dis *d2)
{
	*d1 = *d2;
	d2->p = NULL;
}
