#include "rspd.h"
#include <assert.h>
#include "utility.h"
#include <string.h>

void rsp_init(struct rsp_dis *dis)
{
	memset(dis, 0, sizeof(struct rsp_dis));
}

void rsp_uniform_init(struct rsp_dis *dis)
{
	double prob = 1.0 / RSPD_BUCKET;
	int i;
	dis->pdf[0] = dis->cdf[0] = 0;
	for (i = 1; i <= RSPD_BUCKET; i++) {
		dis->pdf[i] = prob;
		dis->cdf[i] = prob * i;
	}
}

void rsp_join(struct rsp_dis *dis, const struct rsp_dis *part)
{
	array_add(dis->x, part->x, RSPD_BUCKET + 1);
	array_add(dis->y, part->y, RSPD_BUCKET + 1);
}

double newton(struct rsp_dis *dis, double *z0)
{
	double sum = 0;
	double grad = 0;
	int i;

	double z = *z0;

	for (i = 1; i <= RSPD_BUCKET; i++) {
		dis->pdf[i] = (z < EPSILON? 0 : dis->x[i] / z);
		#ifdef RSP_CORRECTION
		grad += (dis->x[i] / (z * z));
		sum += dis->pdf[i];
		if (sum > EPSILON)
			z -= dis->y[i];
		#endif
	}
	#ifdef RSP_CORRECTION
	*z0 += (sum - 1)/grad;
	return fabs(sum - 1);
	#else
	return 0;
	#endif
}

void rsp_finish(struct rsp_dis *dis)
{
	int i;

	#ifdef RSP_CORRECTION
	double z = 0;
	for (int i = 1; i <= RSPD_BUCKET; ++i) {
		z += dis->y[i];
	}
	#else
	double z = 1;
	#endif

	int iter = 0;
	double diff = newton(dis, &z);


	while (diff > 1e-9) {
		++iter;
		diff = newton(dis, &z);
	}

	norm_accumulate(dis->pdf, dis->cdf, RSPD_BUCKET + 1);
}

void rsp_copy(struct rsp_dis *dis, const struct rsp_dis *source)
{
	memcpy(dis, source, sizeof(struct rsp_dis));

}

void rsp_release(struct rsp_dis *dis)
{
}

void rsp_save(const struct rsp_dis *dis, FILE *f)
{
	int i;
	for (i = 1; i <= RSPD_BUCKET; i++)
		fprintf(f, "%lg\n", dis->pdf[i]);

	fprintf(f, "\n");
}

void rsp_load(struct rsp_dis *dis, FILE *f)
{
	int i;

	dis->cdf[0] = dis->pdf[0] = 0;
	for (i = 1; i <= RSPD_BUCKET; i++)
		__s_check__(fscanf(f, "%lf\n", &dis->pdf[i]), 1);

	accumulate(dis->pdf, dis->cdf, RSPD_BUCKET + 1);
}

void rsp_accelerate(const struct rsp_dis *d1, const struct rsp_dis *d2, struct rsp_dis *d3)
{
	accelerate(d1->pdf + 1, d2->pdf + 1, d3->pdf + 1, RSPD_BUCKET);
	rsp_finish(d3);
}

void rspd_diff(const struct rsp_dis *d1, const struct rsp_dis *d2)
{
	print_diff(arr_diff(d1->pdf + 1, d2->pdf + 1, RSPD_BUCKET, THRESHOLD1, THRESHOLD2));
}