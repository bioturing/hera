#include "smodel.h"

void s_model_init(struct s_model *m, struct s_model *m0)
{
	m->trans = m0->trans;
	trans_init(&m->g_dis, &m0->g_dis);
}

void s_model_master_init(struct s_model *m, struct s_model *m0)
{
	m->read_cnt = m0->read_cnt;
	m->trans = m0->trans;
	trans_master_init(&m->g_dis, &m0->g_dis);
}

void s_model_join(struct s_model *full_m, const struct s_model *part_m)
{
	trans_join(&full_m->g_dis, &part_m->g_dis);
}

void s_model_finish(struct s_model *m)
{
	trans_finish(&m->g_dis);
}

void s_model_release(struct s_model *m)
{
	trans_release(&m->g_dis);
}

void s_model_accelerate(const struct s_model *m1, const struct s_model *m2, struct s_model *m3)
{
	trans_accelerate(&m1->g_dis, &m2->g_dis, &m3->g_dis);
}

int s_model_converge(const struct s_model *m1, const struct s_model *m2)
{
	return is_stop(&m1->g_dis, &m2->g_dis);
}

void s_model_diff(const struct s_model *m1, const struct s_model *m2)
{
	trans_diff(&m1->g_dis, &m2->g_dis);
}


void merge_to_f_model(struct model *fm, struct s_model *sm)
{
	trans_release(&fm->g_dis);
	trans_move(&fm->g_dis, &sm->g_dis);
}

char *copy_to_s_model(struct s_model *sm, struct model *fm)
{
	sm->read_cnt = fm->read_cnt;
	sm->trans = fm->trans;
	trans_init(&sm->g_dis, &fm->g_dis);
	trans_copy(&sm->g_dis, &fm->g_dis);
	return trans_locking(&sm->g_dis);
}