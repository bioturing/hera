#include "model.h"
#include "math.h"
#include <assert.h>
//#include <malloc.h>

void model_init(struct model *m, struct model *m0)
{
	m->trans = m0->trans;

	trans_init(&m->g_dis, &m0->g_dis);

	#ifdef USE_FRAG
	m->use_frag = m0->use_frag;
	if (m->use_frag)
		len_init(&m->f_dis, MIN_FRAGMENT_LENGTH, MAX_FRAGMENT_LENGTH);
	#endif

	#ifdef USE_ORI
	ori_init(&m->o_dis);
	#endif

	#ifdef USE_RSPD
	rsp_init(&m->s_dis);
	rsp_copy(&m->os_dis, &m0->s_dis);
	#endif

	#ifdef USE_MIS
	mis_init(&m->r_dis);
	#endif
}

void model_master_init(struct model *m, struct model *m0)
{
	m->read_cnt = m0->read_cnt;
	m->trans = m0->trans;

	trans_master_init(&m->g_dis, &m0->g_dis);

	m->use_frag = m0->use_frag;
	#ifdef USE_FRAG
	if (m->use_frag)
		len_init(&m->f_dis, MIN_FRAGMENT_LENGTH, MAX_FRAGMENT_LENGTH);
	#endif
	#ifdef USE_ORI
	ori_init(&m->o_dis);
	#endif
	#ifdef USE_RSPD
	rsp_init(&m->s_dis);
	rsp_copy(&m->os_dis, &m0->s_dis);

	#endif
	#ifdef USE_MIS
	mis_init(&m->r_dis);
	#endif
}

void model_uniform_init(struct model *m, const struct ref_info_t *trans, int read_cnt, int use_frag)
{
	m->trans = trans;	
	m->read_cnt = read_cnt;

	trans_uniform_init(&m->g_dis, trans->nref);

	#ifdef USE_FRAG
	m->use_frag = use_frag;
	if (m->use_frag)
		len_uniform_init(&m->f_dis, MIN_FRAGMENT_LENGTH, MAX_FRAGMENT_LENGTH);
	#endif

	#ifdef USE_ORI
	ori_uniform_init(&m->o_dis);
	#endif
	
	#ifdef USE_RSPD
	rsp_uniform_init(&m->s_dis);
	rsp_uniform_init(&m->os_dis);
	#endif
	
	#ifdef USE_MIS
	mis_uniform_init(&m->r_dis);
	#endif	
}

void model_join(struct model *full_m, const struct model *part_m)
{
	trans_join(&full_m->g_dis, &part_m->g_dis);

	#ifdef USE_ORI
	ori_join(&full_m->o_dis, &part_m->o_dis);
	#endif

	#ifdef USE_FRAG
	if (full_m->use_frag)
		len_join(&full_m->f_dis, &part_m->f_dis);
	#endif

	#ifdef USE_RSPD
	rsp_join(&full_m->s_dis, &part_m->s_dis);
	#endif
	
	#ifdef USE_MIS
	mis_join(&full_m->r_dis, &part_m->r_dis);
	#endif
}

void model_finish(struct model *m)
{
	trans_finish(&m->g_dis);

	#ifdef USE_ORI
	ori_finish(&m->o_dis);
	#endif

	#ifdef USE_FRAG
	if (m->use_frag)
		len_finish(&m->f_dis);
	#endif

	#ifdef USE_RSPD
	rsp_finish(&m->s_dis);
	#endif
	
	#ifdef USE_MIS
	mis_finish(&m->r_dis);
	#endif
}

void model_release(struct model *m)
{
	trans_release(&m->g_dis);
	#ifdef USE_ORI		
	ori_release(&m->o_dis);
	#endif

	#ifdef USE_FRAG
	if (m->use_frag)
		len_release(&m->f_dis);
	#endif

	#ifdef USE_RSPD	
	rsp_release(&m->s_dis);
	#endif
	#ifdef USE_MIS	
	mis_release(&m->r_dis);
	#endif
}

void model_save(const struct model *m, const char *filename)
{
	char path[MAX_PATH_LENGTH];
	FILE *f;

	#ifdef USE_ORI
	sprintf(path, "%s.ori.txt", filename);
	f = xfopen(path, "w");
	if (f) {
		ori_save(&m->o_dis, f);
		xfclose(f);
	} else {
		__WARNING("Unable to open file %s to write\n", path);
	}
	#endif

	#ifdef USE_FRAG
	sprintf(path, "%s.frag.txt", filename);
	f = xfopen(path, "w");
	if (f) {
		len_save(&m->f_dis, f);
		xfclose(f);
	} else {
		__WARNING("Unable to open file %s to write\n", path);
	}
	#endif

	#ifdef USE_RSPD
	sprintf(path, "%s.rspd.txt", filename);
	f = xfopen(path, "w");
	if (f) {
		rsp_save(&m->s_dis, f);
		xfclose(f);
	} else {
		__WARNING("Unable to open file %s to write\n", path);
	}
	#endif

	#ifdef USE_ALL
	sprintf(path, "%s.all.txt", filename);
	f = xfopen(path, "w");
	if (f) {
		trans_save(&m->g_dis, f);
		ori_save(&m->o_dis, f);
		len_save(&m->f_dis, f);	
		rsp_save(&m->s_dis, f);
		mis_save(&m->r_dis, f);
		xfclose(f);
	} else {
		__WARNING("Unable to open file %s to write\n", path);
	}
	#endif
}

void model_load(struct model *m, const struct ref_info_t *a, const char *filename)
{
	char path[MAX_PATH_LENGTH];
	sprintf(path, "%s.all.txt", filename);	
	
	FILE *f = xfopen(path, "r");

	if (f) {
		m->trans = a;

		trans_load(&m->g_dis, f);
		ori_load(&m->o_dis, f);
		len_load(&m->f_dis, f);
		rsp_load(&m->s_dis, f);
		mis_load(&m->r_dis, f);
		
		xfclose(f);
	} else {
		__WARNING("Unable to open file %s to read\n", path);
	}
}

void model_bootstrap(struct model *m, const struct ref_info_t *a, const char *filename)
{	
	m->trans = a;	

	FILE *f = xfopen(filename, "r");
	if (f) {
		trans_load(&m->g_dis, f);
		xfclose(f);
	} else {
		__ERROR("Unable to open file %s to read\n", filename);
	}

	#ifdef USE_FRAG
	len_uniform_init(&m->f_dis, MIN_FRAGMENT_LENGTH, MAX_FRAGMENT_LENGTH);
	#endif

	#ifdef USE_ORI
	ori_uniform_init(&m->o_dis);
	#endif
	
	#ifdef USE_RSPD
	rsp_uniform_init(&m->s_dis);
	#endif
	
	#ifdef USE_MIS
	mis_uniform_init(&m->r_dis);
	#endif	
}

void model_accelerate(const struct model *m1, const struct model *m2, struct model *m3)
{
	trans_accelerate(&m1->g_dis, &m2->g_dis, &m3->g_dis);

	#ifdef USE_FRAG
	if (m3->use_frag && m2->use_frag && m1->use_frag)
		len_accelerate(&m1->f_dis, &m2->f_dis, &m3->f_dis);
	#endif

	#ifdef USE_MIS
	mis_accelerate(&m1->r_dis, &m2->r_dis, &m3->r_dis);
	#endif

	#ifdef USE_ORI
	ori_accelerate(&m1->o_dis, &m2->o_dis, &m3->o_dis);
	#endif

	#ifdef USE_RSPD
	rsp_accelerate(&m1->s_dis, &m2->s_dis, &m3->s_dis);
	#endif
}

void model_copy(struct model *dest, const struct model *source)
{
	dest->read_cnt = source->read_cnt;
	dest->trans = source->trans;

	trans_copy(&dest->g_dis, &source->g_dis);

	#ifdef USE_FRAG
	len_copy(&dest->f_dis, &source->f_dis);
	#endif

	#ifdef USE_ORI
	ori_copy(&dest->o_dis, &source->o_dis);
	#endif
	
	#ifdef USE_RSPD
	rsp_copy(&dest->s_dis, &source->s_dis);
	rsp_copy(&dest->os_dis, &source->os_dis);
	#endif

	#ifdef USE_MIS
	mis_copy(&dest->r_dis, &source->r_dis);
	#endif
}

const struct result model_to_result(const struct model *m)
{
	struct result res;
	
	res.p = m->g_dis.p;
	res.e_count = malloc(m->trans->nref * sizeof(double));
	res.eff_length = malloc(m->trans->nref * sizeof(double));
	res.tpm = malloc(m->trans->nref * sizeof(double));
	res.fpkm = malloc(m->trans->nref * sizeof(double));
	
	res.t = m->trans;
	
	calc_e_count(&m->g_dis, m->read_cnt, res.e_count);
	calc_eff_length(m->trans, &m->f_dis, res.eff_length);
	calc_tpm(&m->g_dis, res.eff_length, res.tpm);
	calc_fpkm(&m->g_dis, res.eff_length, res.fpkm);
	 
	return res;
}

int model_converge(const struct model *m1, const struct model *m2)
{
	return is_stop(&m1->g_dis, &m2->g_dis);
}

void model_diff(const struct model *m1, const struct model *m2)
{
	trans_diff(&m1->g_dis, &m2->g_dis);

	#ifdef USE_RSPD
	rspd_diff(&m1->s_dis, &m2->s_dis);
	#endif
	#ifdef USE_FRAG
	len_diff(&m1->f_dis, &m2->f_dis);
	#endif
	#ifdef USE_MIS
	mis_diff(&m1->r_dis, &m2->r_dis);
	#endif
}
