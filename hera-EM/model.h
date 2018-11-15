#ifndef _MODEL_
#define _MODEL_

#include "lendis.h"
#include "rspd.h"
#include "oridis.h"
#include "transdis.h"
#include "profile.h"
#include "utility.h"
#include "result.h"
#include "transcript.h"
#include "align.h"

#include <math.h>

/* pair_ended full model */

struct model
{
	int read_cnt;
	int use_frag;
	const struct ref_info_t *trans;

	struct rsp_dis os_dis;

	struct trans_dis g_dis;
	struct ori_dis o_dis;
	struct len_dis f_dis;
	struct rsp_dis s_dis;
	struct mis_dis r_dis;
};


static inline double f_prob(const struct model *m, int f_len, int trans_id)
{
	#ifdef USE_FRAG
	int len_a = get_trans_len(m->trans, trans_id);
	return len_prob(&m->f_dis, f_len, len_a);
	#else
	return 1;
	#endif
}

static inline double s_prob(const struct model *m, int pos, int f_len, int trans_id)
{
	int total_len = get_trans_len(m->trans, trans_id);
	#ifdef USE_RSPD
	return rsp_prob(&m->s_dis, pos, f_len, total_len);
	#else
	return 1.0 / (total_len - f_len + 1);
	#endif
}

static inline double o_prob(const struct model *m, int dir)
{
	#ifdef USE_ORI
	return ori_pdf(&m->o_dis, dir);
	#else
	return 1.0; /*or 0.5, doesn't matter (only affect log likelihood) */
	#endif
}

static inline double g_prob(const struct model *m, int trans_id)
{
	return trans_pdf(&m->g_dis, trans_id);
}

static inline double ng_prob(const struct model *m)
{
	#if USE_NOISE == 1
	return ntrans_pdf(&m->g_dis);
	#else
	return 0;
	#endif
}

static inline double r_prob(const struct model *m, int mismatch, int read_len)
{
	#ifdef USE_MIS
	return mis_prob(&m->r_dis, mismatch, read_len);
	#else
	return 1;
	#endif
}

static inline double nr_prob(const struct model *m, int length)
{
	#if (USE_NOISE == 1) && defined(USE_MIS)
	return pow_table[length];
	#else
	return 1;
	#endif
}

static inline double cond_prob(const struct model *m, const struct s_read_t *read, const union ualign_t *a)
{
	const struct f_align_t *align = &a->full;
	int f_len = align->f_len;
	
	double f = 1, s = 1, o = 1, x = 1;

	f = f_prob(m, f_len, align->trans_id);

	s = s_prob(m, align->pos, f_len, align->trans_id);

	int ori = align->ori;

	o = o_prob(m, ori);

	x = r_prob(m, align->dis, read->len);

	return f * s * o * x;
}

static inline double cond_no_frag_prob(const struct model *m, const struct s_read_t *read, const union ualign_t *a)
{
	const struct f_align_t *align = &a->full;
	
	double s = 1, o = 1, x = 1;

	s = s_prob(m, align->pos, read->len, align->trans_id);

	int ori = align->ori;

	o = o_prob(m, ori);

	x = r_prob(m, align->dis, read->len);

	return s * o * x;
}

static inline double joint_prob(const struct model *m, const struct s_read_t *read, const union ualign_t *align)
{
	double g = g_prob(m, align->full.trans_id);
	return g * cond_prob(m, read, align);
}

static inline double joint_no_frag_prob(const struct model *m, const struct s_read_t *read, const union ualign_t *align)
{
	double g = g_prob(m, align->full.trans_id);
	return g * cond_no_frag_prob(m, read, align);
}

static inline double cond_noise_prob(const struct model *m, const struct s_read_t *read)
{
	#if USE_NOISE == 0
	return 0;
	#endif
 
	return nr_prob(m, read->len);
}

static inline double noise_prob(const struct model *m, const struct s_read_t *read)
{
	#if USE_NOISE == 0
	return 0;
	#endif

	double g = ng_prob(m);
	return g * cond_noise_prob(m, read);
}

void model_init(struct model *m, struct model *m0);
void model_master_init(struct model *m, struct model *m0);

void model_uniform_init(struct model *m, const struct ref_info_t *trans, int read_cnt, int use_frag);

static inline int model_no_frag_update(struct model *m, const struct s_read_t *r, const union ualign_t *a, double prob)
{
	const struct f_align_t *align = &a->full;
	if (prob < FILTERING_THRESHOLD)
		return 2; //remove alignment

	int trans_len = get_trans_len(m->trans, align->trans_id);

	trans_update(&m->g_dis, align->trans_id, prob);

	#ifdef USE_RSPD
	rsp_update(&m->s_dis, &m->os_dis, align->pos, r->len, trans_len, prob);
	#endif
	
	#ifdef USE_ORI
	ori_update(&m->o_dis, align->ori, prob);
	#endif
	
	#ifdef USE_MIS
	mis_update(&m->r_dis, align->dis, r->len, prob);
	#endif

	return 0; //keep
}

static inline int model_update(struct model *m, const struct s_read_t *r, const union ualign_t *a, double prob)
{
	const struct f_align_t *align = &a->full;
	if (prob < FILTERING_THRESHOLD)
		return 2; //remove alignment

	int trans_len = get_trans_len(m->trans, align->trans_id);
	
	int f_len = align->f_len;

	trans_update(&m->g_dis, align->trans_id, prob);

	#ifdef USE_FRAG
	len_update(&m->f_dis, f_len, trans_len, prob);
	#endif

	#ifdef USE_RSPD
	rsp_update(&m->s_dis, &m->os_dis, align->pos, f_len, trans_len, prob);
	#endif
	
	#ifdef USE_ORI
	ori_update(&m->o_dis, align->ori, prob);
	#endif
	
	#ifdef USE_MIS
	mis_update(&m->r_dis, align->dis, r->len, prob);
	#endif

	return 0; //keep
}

static inline int model_update_noise(struct model *m, const struct s_read_t *read, double prob)
{
	#if USE_NOISE == 1
	if (prob < FILTERING_THRESHOLD)
		return 2;
	ntrans_update(&m->g_dis, prob);
	#endif
	return 0;
}

const struct result model_to_result(const struct model *m);

int model_converge(const struct model *m1, const struct model *m2);

void model_diff(const struct model *m1, const struct model *m2);

void model_join(struct model *full_m, const struct model *part_m);
void model_finish(struct model *m);
void model_release(struct model *m);

void model_save(const struct model *m, const char *filename);
void model_load(struct model *m, const struct ref_info_t *a, const char *filename);
void model_bootstrap(struct model *m, const struct ref_info_t *a, const char *filename);

void model_accelerate(const struct model *m1, const struct model *m2, struct model *m3);
void model_copy(struct model *dest, const struct model *source);

#endif
