#include "bow_parser.h"
#include "transcript.h"

static inline int compare_name(uint8_t *r_name, int name_len, bam1_t *aln)
{
	if (name_len != aln->core.l_qname)
		return 0;

	if (strncmp((char *)r_name, (char *)aln->data, name_len) != 0)
		return 0;
	return 1;
}

void add_alignment(bam1_t *aln, union ualign_t *read_align,
		int *n, bam_hdr_t *header)
{
	if ((aln->core.flag & BAM_FUNMAP) == BAM_FUNMAP ||
	    (aln->core.flag & BAM_FPROPER_PAIR) != BAM_FPROPER_PAIR)
		return;

	if (aln->core.isize < 0)
		return;


	int s1, s2, ori;

	s1 = bam_aux2i(bam_aux_get(aln, "NM"));
	s2 = -bam_aux2i(bam_aux_get(aln, "YS"));


	ori = (aln->core.flag & BAM_FREVERSE) == BAM_FREVERSE? 1: 0;

	read_align[*n].full.trans_id = aln->core.tid;
	read_align[*n].full.ori = ori;
	read_align[*n].full.pos = aln->core.pos * RSPD_BUCKET / 
					header->target_len[aln->core.tid] + 1;
	read_align[*n].full.f_len = aln->core.isize;
	read_align[*n].full.dis = s1 + s2;

	++*n;
}

void batch_extend(struct align_batch_t *a, int cnt)
{
	if (a->size + cnt <= a->cap)
		return;

	uint64_t new_cap = a->size + cnt + a->cap;

	a->align = realloc(a->align, new_cap * sizeof(*a->align));
	a->cap = new_cap;
}

void write_read(struct align_batch_t *a, const union ualign_t *alg, int l, int cnt)
{
	if (cnt > MAX_ALIGNMENT)
		return;

	batch_extend(a, cnt + 1);

	a->align[a->size].read.cnt = cnt + 1; //jump not count
	a->align[a->size].read.len = l;
	a->align[a->size].read.prob = 0;
	memcpy(a->align + a->size + 1, alg, cnt * sizeof(*alg));

	a->size += cnt + 1;
	++a->read_cnt;
}

struct align_batch_t *parse_bam(const char *path, int thread_cnt, struct ref_info_t *ref)
{
	samFile *f;
	bam_hdr_t *header; 
	bam1_t *aln;
	uint8_t *r_name;
	union ualign_t *read_align;
	struct align_batch_t * batch;
	int i, name_len, ret, alg_cnt, current, name_cap, skip, l1, l2;

	f = hts_open(path,"r");
	header = sam_hdr_read(f);
	aln = bam_init1();

	batch = calloc(thread_cnt, sizeof(struct align_batch_t));
	current = 0;

	read_align = calloc(MAX_ALIGNMENT, sizeof(union ualign_t));
	r_name = malloc(1);
	name_len = 0;
	name_cap = 1;
	ret = sam_read1(f, header, aln);

	//read ref_inf;

	ref->nref = header->n_targets;
	ref->name = calloc(ref->nref, sizeof(*ref->name));
	ref->len = calloc(ref->nref, sizeof(*ref->len));

	for (i = 0; i < header->n_targets; ++i) {
		ref->name[i] = strdup(header->target_name[i]);
		ref->len[i] = header->target_len[i];
	}

	i = 0;
	alg_cnt = 0;

	while(ret > 0) {
		name_len = aln->core.l_qname;
		if (name_len > name_cap) {
			name_cap = name_len;
			r_name = realloc(r_name, name_cap);
		}

		memcpy(r_name, aln->data, name_len);
		read_align[0].read.len = aln->core.l_qseq;

		l1 = l2 = 0;
		skip = false;
		while (ret > 0 && compare_name(r_name, name_len, aln)){
			if (alg_cnt >= MAX_ALIGNMENT)
				skip = true;
			else
				add_alignment(aln, read_align, &alg_cnt, header);

			if (aln->core.isize >= 0)
				l1 = aln->core.l_qseq;
			else 
				l2 = aln->core.l_qseq;
			
			ret = sam_read1(f, header, aln);
		}
		

		if (alg_cnt == 0) {
			++batch[current].noise; //unmapped read
			skip = true;
		}

		if (!skip) {
			write_read(&batch[current], read_align, l1 + l2, alg_cnt);
			++current;
			if (current == thread_cnt)
				current = 0;
		}

		alg_cnt = 0;

		++i;
		if (i % 1000000 == 0) {
			fprintf(stderr, "Done load %d reads\r", i);
		}
	}

	fprintf(stderr, "Done load %d reads\n", i);
	
	bam_destroy1(aln);
	sam_close(f);

	free(read_align);
	free(r_name);

	for(i = 0; i < thread_cnt; ++i)
		batch_resize(&batch[i]);

	return batch;
}
