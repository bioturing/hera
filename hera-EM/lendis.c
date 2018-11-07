#include "lendis.h"
//#include <malloc.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define create(dis) ((double*) (calloc(len_span((dis)), sizeof(double))) \
				- (dis)->min_l)

void len_init(struct len_dis *dis, int min_l, int max_l)
{
	dis->min_l = min_l;
	dis->max_l = max_l;

	dis->pdf = create(dis);
	dis->cdf = create(dis);

	dis->x = create(dis);
	dis->y = create(dis);
}

void len_uniform_init(struct len_dis *dis, int min_l, int max_l)
{
	len_init(dis, min_l, max_l);
	double prob = 1.0 / len_span(dis);

	int i;
	for (i = min_l; i <= max_l; i++) {
		dis->pdf[i] = prob;
		dis->cdf[i] = prob * (i - min_l + 1);
	}
}

void len_normal_init(struct len_dis *dis, int min_l, int max_l, double mean, double sd)
{
	len_init(dis, min_l, max_l);

	int i;
	for (i = min_l; i <= max_l; i++)
		dis->pdf[i] = exp(-pow((i - mean)/(sd),2)/2);

	norm_accumulate(dis->pdf + min_l, dis->cdf + min_l, len_span(dis));
}

void len_join(struct len_dis *dis, const struct len_dis *part)
{
	int min_l = dis->min_l;
	int span = len_span(dis);
	array_add(dis->x + min_l, part->x + min_l, span);
	array_add(dis->y + min_l, part->y + min_l, span);
}

void len_finish(struct len_dis *dis) 
{
	/*len_trim(dis);*/
	int min_l = dis->min_l;
	int span = len_span(dis);

	int i;

	double z = 1;

	#ifdef LEN_CORRECTION
	double sum = 0;
	#endif

	for (i = dis->min_l; i <= dis->max_l; i++) {
		dis->pdf[i] = dis->x[i] / z;

		#ifdef LEN_CORRECTION
		sum += dis->pdf[i];
		if (sum > EPSILON)
			z -= dis->y[i]/sum;
		#endif
	}

	norm_accumulate(dis->pdf + min_l, dis->cdf + min_l, span);
}

void len_release(struct len_dis *dis) {
	if (dis->cdf)
		free(dis->cdf + dis->min_l);

	if (dis->pdf)
		free(dis->pdf + dis->min_l);

	if (dis->x)
		free(dis->x + dis->min_l);

	if (dis->y)
		free(dis->y + dis->min_l);
}

void len_save(const struct len_dis *dis, FILE *f)
{
	fprintf(f, "%d %d\n", dis->min_l, dis->max_l);

	int i;
	for (i = dis->min_l; i <= dis->max_l; i++)
		fprintf(f, "%lg\n", dis->pdf[i]);
	fprintf(f, "\n");
}

void len_load(struct len_dis *dis, FILE *f)
{
	int min, max;
	__s_check__(fscanf(f, "%d %d\n", &min, &max), 2);

	len_init(dis, min, max);

	int i;
	dis->pdf[min - 1] = dis->cdf[min - 1] = 0;
	for (i = min; i <= max; i++)
		__s_check__(fscanf(f, "%lf\n", &dis->pdf[i]), 1);

	accumulate(dis->pdf + min, dis->cdf + min, max - min + 1);
	
}

void len_copy(struct len_dis *dis, const struct len_dis *source)
{
	len_init(dis, source->min_l, source->max_l);

	memcpy(dis->pdf, source->pdf, len_span(dis) * sizeof(double));
	memcpy(dis->cdf, source->cdf, len_span(dis) * sizeof(double));
}

void len_accelerate(const struct len_dis *d1, const struct len_dis *d2, struct len_dis *d3)
{
	int min_l = d1->min_l;
	int span = d1->max_l - d1->min_l + 1;
	accelerate(d1->pdf + min_l, d2->pdf + min_l, d3->pdf + min_l, span);
	len_finish(d3);
}


void len_diff(const struct len_dis *d1, const struct len_dis *d2)
{
	print_diff(arr_diff(d1->pdf + d1->min_l,
			    d2->pdf + d1->min_l,
			    d1->max_l - d1->min_l + 1, 
			    THRESHOLD1,
			    THRESHOLD2
		)
	);
}

void len_move(struct len_dis *d1, struct len_dis *d2) 
{
	*d1 = *d2;
}
