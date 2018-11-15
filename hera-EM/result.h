#ifndef _RESULT_
#define _RESULT_

#include "transdis.h"
#include "lendis.h"
#include "transcript.h"

struct result
{
	const struct ref_info_t *t;

	double *p;
	double *eff_length;
	double *e_count;
	double *tpm;
	double *fpkm;
};

void calc_eff_length(const struct ref_info_t *t, const struct len_dis *f_dis, double *el);

void calc_tpm(const struct trans_dis *g_dis, const double *el, double *tpm);
void calc_fpkm(const struct trans_dis *g_dis, const double *el, double *fpkm);

void calc_e_count(const struct trans_dis *g_dis, int read_cnt, double *e_count);

int print_full_result(const struct result *r, const char *filename);
int print_rounded_tpm(const struct result *r, const char *filename);
void write_abundance(const char *prefix, const struct result *r);


void result_release(const struct result *r);

#endif
