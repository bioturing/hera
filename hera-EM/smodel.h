#ifndef _SMODEL_
#define _SMODEL_

#include "model.h"
#include "utility.h"
#include "transcript.h"

struct s_model {
	int read_cnt;

	const struct ref_info_t *trans;
	struct trans_dis g_dis;
};

static inline double s_joint_prob(const struct s_model *m, const struct s_read_t *read, const union ualign_t *align)
{
	const struct s_align_t *s_align = &align->simp;
	double g = m->g_dis.p[s_align->trans_id];
	return g * s_align->c_prob;
}

static inline double s_noise_prob(const struct s_model *m, const struct s_read_t *read)
{
	double g = m->g_dis.p[-1];
	return g * pow_table[read->len];
}

static inline int s_model_update(struct s_model *m, const struct s_read_t *r, const union ualign_t *align, double prob)
{
	if (prob < FILTERING_THRESHOLD)
		return 2; //remove alignment
	if (m->g_dis.lock[align->simp.trans_id])
		return 1; //lock probabilty

	trans_update(&m->g_dis, align->simp.trans_id, prob);
	return 0;
}

static inline int s_model_update_noise(struct s_model *m, const struct s_read_t *read, double prob)
{
	#if USE_NOISE == 1
	if (prob < FILTERING_THRESHOLD)
		return 2; //remove
	if (m->g_dis.lock && m->g_dis.lock[-1])
		return 1;

	ntrans_update(&m->g_dis, prob);
	#endif
	return 0;
}

void s_model_init(struct s_model *m, struct s_model *m0);
void s_model_master_init(struct s_model *m, struct s_model *m0);

void s_model_join(struct s_model *full_m, const struct s_model *part_m);
void s_model_finish(struct s_model *m);

void s_model_release(struct s_model *m);

void s_model_accelerate(const struct s_model *m1, const struct s_model *m2, struct s_model *m3);

int s_model_converge(const struct s_model *m1, const struct s_model *m2);

void s_model_diff(const struct s_model *m1, const struct s_model *m2);

void merge_to_f_model(struct model *fm, struct s_model *sm);

char *copy_to_s_model(struct s_model *sm, struct model *fm);

#endif