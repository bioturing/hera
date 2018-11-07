#include "result.h"
#include "../utils.h"
#include "../bin.h"

void calc_eff_length(const struct ref_info_t *t, const struct len_dis *f_dis, double *el)
{
	int i;
	
	#ifdef USE_FRAG	

	int min_l = f_dis->min_l;
	int max_l = f_dis->max_l;

	double *p = f_dis->pdf;
	double *c = f_dis->cdf;
	double *g = (double *)malloc(sizeof(double) * (max_l - min_l + 1)) - min_l;

	g[min_l] = p[min_l] * min_l;
	for (i = min_l + 1; i <= max_l; i++) 
		g[i] = g[i - 1] + p[i] * i;

	#endif

	for (i = 0; i < t->nref; i++) {
		int len = get_trans_len(t, i);
		#ifdef USE_FRAG
		double eff = 0;
		/*double total = 0;
		eff += p[j] * (len - j + 1) = p[j] * len - p[j] * j + p[j]
		=> eff = c[max] * (len + 1) - g[max] */
		int max_f = min(len, max_l);
		eff = (len >= min_l ? (c[max_f] * (len + 1) - g[max_f]) : 0);
		/*for (j = m->f_dis.minL; j <= min(len, m->f_dis.maxL); j++) {
			//total += m->f_dis.pdf[j];
			eff += m->f_dis.pdf[j] * (len - j + 1);
		}*/

		if (eff < MIN_EFFECTIVE_LENGTH)
		#ifdef EFFECTIVE_LENGTH_ROUND_UP
			eff = MIN_EFFECTIVE_LENGTH;
		#else
			eff = 0;
		#endif
		
		el[i] = eff;
		#else
		el[i] = len;
		#endif
	}

	#ifdef USE_FRAG
	free(g + min_l);
	#endif

}

void calc_fpkm(const struct trans_dis *g_dis, const double *el, double *fpkm)
{
	double total = 0;
	int  i;

	for (i = 0; i < g_dis->trans_cnt; i++)
		total += g_dis->p[i];

	for (i = 0; i < g_dis->trans_cnt; i++)
		if (el[i] > EPSILON)
			fpkm[i] = g_dis->p[i] / total / el[i] * 1e9;
		else
			fpkm[i] = 0;
}

void calc_tpm(const struct trans_dis *g_dis, const double *el, double *tpm)
{
	double total = 0;
	int i;
	for (i = 0; i < g_dis->trans_cnt; i++)
		if (el[i] > EPSILON)
			total += tpm[i] = g_dis->p[i] / el[i];
		else
			tpm[i] = 0;

	total /= 1e6;

	for (i = 0; i < g_dis->trans_cnt; i++)
		tpm[i] /= total;
}

void calc_e_count(const struct trans_dis *g_dis, int read_cnt, double *e_count)
{
	int i;
	for (i = 0; i < g_dis->trans_cnt; i++)
		e_count[i] = g_dis->p[i] * read_cnt;
}

#if defined(APP_BUILD)
int print_full_result(const struct result *r)
{
	int i, g_iter, cur_gene_id, prev_gene_id;
	float t_express, g_express, g_count, g_tpm[r->t->nref], temp[r->t->nref];
	memset(g_tpm, 0, r->t->nref * sizeof(float));

	/* trans count */
	for (i = 0; i < r->t->nref; ++i)
		temp[i] = (float)r->e_count[i];
	__BIN_WRITE(temp, sizeof(float), r->t->nref);
	/* trans tpm */
	for (i = 0; i < r->t->nref; ++i)
		temp[i] = (float)r->tpm[i];
	__BIN_WRITE(temp, sizeof(float), r->t->nref);

	t_express = g_express = 0.0;
	g_count = 0.0;

	for (i = g_iter = 0; i < r->t->nref; ++i) {
		cur_gene_id = r->t->gene[i];

		if (r->e_count[i] > 1)
			++t_express;

		if (i > 0 && cur_gene_id != prev_gene_id) {
			if (g_count > 1)
				++g_express;
			/* gene count */
			__BIN_WRITE(&g_count, sizeof(float), 1);
			g_count = 0;
			g_iter++;
		}

		g_tpm[g_iter] += (float)r->tpm[i];
		g_count += (float)r->e_count[i];
		prev_gene_id = cur_gene_id;
	}

	/* gene count */
	__BIN_WRITE(&g_count, sizeof(float), 1);
	if (g_count > 1)
		++g_express;

	/* gene tpm */
	__BIN_WRITE(g_tpm, sizeof(float), r->t->gene_map->n_gene);
	/* total trans express */
	__BIN_WRITE(&t_express, sizeof(float), 1);
	/* total gene express */
	__BIN_WRITE(&g_express, sizeof(float), 1);

	return 1;
}

int print_rounded_tpm(const struct result *r)
{
	// FIXME: implement me
}

void write_abundance(const struct result *r)
{
	print_full_result(r);

#ifdef WRITE_ONLY_TPM
	print_rounded_tmp(r);
#endif
}
#else
int print_full_result(const struct result *r, const char *filename)
{
	FILE *f = xfopen(filename, "w");

	if (!f)
		return 0;
	
	int name_length = r->t->gene_map->l_gene_name;
	int id_length = r->t->gene_map->l_gene_id;

	char *gene_name = calloc(name_length + 1, sizeof(char));
	char *gene_id = calloc(r->t->gene_map->l_gene_id + 1, sizeof(char));

	fprintf(f, "#trans_id\tstrand\tgene_id\tgene_name\ttrans_length\teffective_len\tread_count\ttpm\tfpkm\n");
	int i;
	for (i = 0; i < r->t->nref; i++) {
		int id = r->t->gene[i];

		char strand = r->t->gene_map->strand[id];

		memcpy(gene_name, r->t->gene_map->gene_name 
				+ id * name_length, name_length);

		memcpy(gene_id, r->t->gene_map->gene_id 
				+ id * id_length, id_length);

		fprintf(f, "%s\t%c\t%s\t%s\t%d\t%lf\t%lf\t%lf\t%lf\n",
			get_trans_id(r->t, i),
			strand,
			gene_id,
			gene_name,
			get_trans_len(r->t, i),
			r->eff_length[i],
			r->e_count[i],
			r->tpm[i],
			r->fpkm[i]
		);
	}

	free(gene_name);
	free(gene_id);
	xfclose(f);
	
	return 1;
}

int print_rounded_tpm(const struct result *r, const char *filename)
{
	FILE *f = xfopen(filename, "w");

	if (!f)
		return 0;
	
	int i;
	for (i = 0; i < r->t->nref; i++)
		fprintf(f, "%s\t%.2lf\n",
			get_trans_id(r->t, i),
			r->tpm[i]
		);

	xfclose(f);
	return 1;
}

void write_abundance(const char *prefix, const struct result *r)
{
	char *file_name = str_concate(prefix, ".tsv");

	if (!print_full_result(r, file_name))
		fprintf(stderr, RED "Cannot write result to file %s" RESET, file_name);
	else
		fprintf(stderr, GREEN "Result written to file %s\n" RESET, file_name);

	free(file_name);

#ifdef WRITE_ONLY_TPM
	file_name = str_concate(prefix, ".abundance.tsv");
	print_rounded_tpm(r, file_name);
	free(file_name);
#endif
}
#endif




void result_release(const struct result *r)
{
	free(r->e_count);
	free(r->eff_length);
	free(r->tpm);
	free(r->fpkm);
}
