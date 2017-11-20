#include "bam_write.h"

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

void contruct_gene_map(FILE *idx)
{
	unsigned int ret, i;
	GENE_MAP = malloc(sizeof(Gene_map));

	// Get chromosome name
	ret = fread(&GENE_MAP->n_chr, 1, sizeof(int), idx);
	ret = fread(&GENE_MAP->l_chr, 1, sizeof(int), idx);
	GENE_MAP->chr_name = malloc(GENE_MAP->l_chr*GENE_MAP->n_chr);
	GENE_MAP->chr_len = calloc(GENE_MAP->n_chr, sizeof(int));
	ret = fread(GENE_MAP->chr_name, GENE_MAP->l_chr, GENE_MAP->n_chr, idx);
	ret = fread(GENE_MAP->chr_len, sizeof(int), GENE_MAP->n_chr, idx);

	// Get gene name + id + strand
	ret = fread(&GENE_MAP->n_gene, sizeof(int), 1, idx);
	ret = fread(&GENE_MAP->l_gene, sizeof(int), 1, idx);
	GENE_MAP->gene_name = malloc(GENE_MAP->l_gene*GENE_MAP->n_gene);
	ret = fread(GENE_MAP->gene_name, GENE_MAP->l_gene,
						 GENE_MAP->n_gene, idx);

	// Get transcriptome -> genome map
	GENE_MAP->gene = calloc(REF_INF->n, sizeof(int));
	GENE_MAP->chr = calloc(REF_INF->n, sizeof(int));
	GENE_MAP->n_exon = calloc(REF_INF->n, sizeof(short));
	GENE_MAP->start = calloc(REF_INF->n, sizeof(int*));
	GENE_MAP->end = calloc(REF_INF->n, sizeof(int*));

	for (i = 0; i < REF_INF->n; ++i){
		ret = fread(&GENE_MAP->n_exon[i], 1, sizeof(int), idx);
		GENE_MAP->start[i] = calloc(GENE_MAP->n_exon[i], sizeof(int));
		GENE_MAP->end[i] = calloc(GENE_MAP->n_exon[i], sizeof(int));
		ret = fread(GENE_MAP->start[i], sizeof(int),
						  GENE_MAP->n_exon[i], idx);
		ret = fread(GENE_MAP->end[i], sizeof(int),
						  GENE_MAP->n_exon[i], idx);
	}
	ret = fread(GENE_MAP->chr, sizeof(int), REF_INF->n, idx);
	ret = fread(GENE_MAP->gene, sizeof(int), REF_INF->n, idx);
}

static inline int bam_reg2bin(unsigned int beg, unsigned int end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

unsigned int make_flag(unsigned int pair, unsigned int proper, unsigned int map,
			 unsigned int mmap, unsigned int rev, unsigned int mrev,
					  unsigned int first, unsigned int last)
{
	unsigned int flag = last;
	flag = (flag << 1) | first;
	flag = (flag << 1) | mrev;
	flag = (flag << 1) | rev;
	flag = (flag << 1) | mmap;
	flag = (flag << 1) | map;
	flag = (flag << 1) | proper;
	flag = (flag << 1) | pair;

	return flag;
}

void add_padding(unsigned int *cstring, unsigned int *n, unsigned int *l,
		unsigned int *coord, unsigned int trans, unsigned int len,
			   unsigned int type, int *size, unsigned int err)
{
	if ((type != MATCH && type != DIFF) || err == 1){
		cstring[(*n)++] = len << CHUNK | type;
		return;
	}

	*l += len;

	while (*l > GENE_MAP->end[trans][*coord] + 1){
		len -= *l - GENE_MAP->end[trans][*coord] - 1;
		cstring[(*n)++] = len << CHUNK | type;

		len = GENE_MAP->start[trans][*coord + 1] -
				GENE_MAP->end[trans][*coord] - 1;
		*size += len;
		cstring[(*n)++] = len << CHUNK | SKIP;

		len = *l - GENE_MAP->end[trans][*coord] - 1;
		*l += GENE_MAP->start[trans][*coord + 1] -
				GENE_MAP->end[trans][*coord] - 1;
		++*coord;
	}

	cstring[(*n)++] = len << CHUNK | type;
}

void convert_cigar(Bam_core *core, Candidate *r, unsigned int p,
			unsigned int gene, unsigned int convert)
{
	unsigned int i, k, j, match, n, l, trans;

	trans = r->pos[p][gene];
	if (convert == 0){
		core->ref = GENE_MAP->chr[trans];
		core->pos = gene_coord(trans, r->p[p][gene], &core->coord);
	} else {
		chr_coord(r->pos[p][gene], &core->ref, &core->pos);
	}

	for (i = match = k = 0; i < r->cigar[p]->n_cigar; ++i){
		if (r->cigar[p]->type[i] == MATCH){
			n = r->cigar[p]->count[i++];
			while (i < r->cigar[p]->n_cigar &&
				 r->cigar[p]->type[i] != DIFF &&
				 r->cigar[p]->type[i] != DELETE){
				if (r->cigar[p]->type[i] == MATCH)
					n += r->cigar[p]->count[i];
				++i;
			}
			--i;
			match += n;
			sprintf(core->miss + core->len,"%u%c", n, '\0');
		} else if (r->cigar[p]->type[i] == DIFF){
			sprintf(core->miss + core->len,"%c%c",
					  r->cigar[p]->miss[k++], '\0');
			for(j = 1; j < r->cigar[p]->count[i]; ++j)
				sprintf(core->miss + core->len + 2*j - 1,
				  "0%c%c", r->cigar[p]->miss[k++], '\0');
		} else if (r->cigar[p]->type[i] == DELETE){
			sprintf(core->miss + core->len,"^%c", '\0');
			for(j = 0; j < r->cigar[p]->count[i]; ++j)
				sprintf(core->miss + core->len + j + 1,
				 "%c%c", r->cigar[p]->miss[k++], '\0');
		}
		core->len = strlen(core->miss);
	}

	for (i = n = 0, l = core->pos; i < r->cigar[p]->n_cigar; ++i){
		k = 0;
		while (i < r->cigar[p]->n_cigar && (r->cigar[p]->type[i] == DIFF 
					       || r->cigar[p]->type[i] == MATCH))
			k += r->cigar[p]->count[i++];

		if (k > 0)			
			add_padding(core->cstring, &n, &l, &core->coord, trans,
			 		  k, MATCH, &core->align_len, convert);

		if (i < r->cigar[p]->n_cigar)
			add_padding(core->cstring, &n, &l, &core->coord, trans,
				   r->cigar[p]->count[i] ,r->cigar[p]->type[i],
						    &core->align_len, convert);
	}

	core->n_cigar = n;
	core->match = match;
}

inline unsigned int get_qual(unsigned int n)
{
	double error = 1 - (1.0/n);
	if (error == 0)
		error = 0.0000001;
	return -10*log10(error);
}

void compact_inf(struct read_inf *read, Bam_core *core, unsigned int n_miss,
						unsigned int name_len)
{
	unsigned int add = BLOCK + CHUNK;
	unsigned int i, len, l;

	l = name_len + CHUNK*core->n_cigar + 3*read->len + add + core->len;
	core->data = malloc(l);

	// Read name
	memcpy(core->data + add, read->name, name_len);
	len = add + name_len;

	// Cigar
	memcpy(core->data + len, core->cstring, CHUNK*core->n_cigar);
	len += CHUNK*core->n_cigar;

	// Sequence
	for (i = 0, core->data[len] = 0; i < read->len; ++i){
		core->data[len] = (core->data[len] << CHUNK) |
				  bam_nt16_table[read->seq[i]];
		if ((i&1) == 1)
			core->data[++len] = 0;
	}
	if ((i&1) == 1){
		core->data[len]  = core->data[len] << CHUNK;
		++len;
	}

	// Quality string
	for (i = 0; i < read->len; ++i)
		core->data[len + i] = read->qual[i] - 33;
	len += read->len;

	// Tag field
	if (core->match > 0){
		// Number of different base on reference
		memcpy(core->data + len, "NMi", 3);
		memcpy(core->data + len + 3, &n_miss, CHUNK);
		len += 3 + CHUNK;

		// Different position
		memcpy(core->data + len, "MDZ", 3);
		memcpy(core->data + len + 3, core->miss, core->len + 1);
		len += 3 + core->len + 1;
	}

	// Align score
	memcpy(core->data + len, "ASi", 3);
	memcpy(core->data + len + 3, &core->match, CHUNK);
	len += 3 + CHUNK;

	core->block_len = len;
}

Bam_core *init_bamCore(unsigned int len)
{
	Bam_core *core = malloc(sizeof(Bam_core));
	core->cstring = calloc(len, sizeof(int));
	core->miss = malloc(2*len);
	core->miss[0] = '\0';
	core->data = NULL;

	core->bin = core->qual = core->n_cigar = core->len = 0;
	core->block_len = core->match = core->align_len = 0;
	core->ref = core->pos = -1;

	return core;
}

void destroy_bamCore(Bam_core *core)
{
	if (core->data != NULL)
		free(core->data);
	free(core->miss);
	free(core->cstring);
	free(core);
}

void pack_bam_core(Bam_core *core, Candidate *r, unsigned int p,
		struct read_inf *read, unsigned int gene, unsigned int name_len)
{
	unsigned int n;
	unsigned int convert = r->pair & 1;

	if (r->n[p] > 0)
		convert_cigar(core, r, p, gene, convert);

	core->bin = bam_reg2bin(core->pos, core->pos + core->match);

	if (convert == 1)
		n = r->n[p];
	else
		n = check_quanlity(r, p);
	core->qual = get_qual(n);

	compact_inf(read, core, r->cigar[p]->n_miss, name_len);
}

void pack_bam_block(Bam_core *core, Bam_core *mate_core, unsigned int name_len,
			    unsigned int read_len, unsigned int flag, int size,
					      char *stream, unsigned int *slen)
{
	int x[8];

	x[0] = core->ref;
	x[1] = core->pos;
	x[2] = (unsigned int) core->bin << 16 | core->qual << 8 | name_len;
	x[3] = (unsigned int) flag << 16 | core->n_cigar;
	x[4] = read_len;
	x[5] = mate_core->ref;
	x[6] = mate_core->pos;
	x[7] = size;

	core->len = core->block_len - CHUNK;

	memcpy(core->data, &core->len, CHUNK); 
	memcpy(core->data + CHUNK, x, BLOCK); 

	if (*slen + core->block_len > MAX_COMPRESS){
		bam_write(BAM, stream, *slen);
		*slen = 0;
	}
	memcpy(stream + *slen, core->data, core->block_len);
	*slen += core->block_len;
}

void write_alignment_single(struct read_inf *read, unsigned int flag, Candidate *r,
			 		 unsigned int p, unsigned int gene,
					  char *stream, unsigned int *slen)
{
	Bam_core *core[2];
	core[0] = init_bamCore(read->len);
	core[1] = init_bamCore(read->len);

	unsigned int name_len;

	name_len = strlen(read->name);
	if (name_len > 2 && read->name[name_len - 2] == '/'){
		read->name[name_len - 2] = '\0';
		--name_len;
	} else {
		++name_len;
	}

	// Assign
	pack_bam_core(core[0], r, p, read, gene, name_len);
	pack_bam_block(core[0], core[1], name_len, read->len, flag,
					  	 0, stream, slen);
	destroy_bamCore(core[0]);
	destroy_bamCore(core[1]);
}

void write_alignment_pair(struct read_inf *read1, struct read_inf *read2,
	     unsigned int *flag, Candidate *r, unsigned int gene,
				char *stream, unsigned int *slen)
{
	Bam_core *core[2];
	core[0] = init_bamCore(read1->len);
	core[1] = init_bamCore(read2->len);

	unsigned int name_len;
	int size = 0;

	name_len = strlen(read1->name);
	if (name_len > 2 && read1->name[name_len - 2] == '/'){
		read1->name[name_len - 2] = '\0';
		read1->name[name_len - 2] = '\0';
		--name_len;
	} else {
		++name_len;
	}

	// Assign
	pack_bam_core(core[0], r, 0, read1, gene, name_len);
	pack_bam_core(core[1], r, 1, read2, gene, name_len);

	if (core[0]->ref >= 0 && core[0]->ref == core[1]->ref){
		if (core[0]->pos < core[1]->pos)
			size = core[1]->pos + read2->len + core[1]->align_len
							      - core[0]->pos;
		else
			size = core[0]->pos + read1->len + core[0]->align_len
							      - core[1]->pos;
	}

	pack_bam_block(core[0], core[1], name_len, read1->len, flag[0],
						  size, stream, slen);
	pack_bam_block(core[1], core[0], name_len, read2->len, flag[1],
						 -size, stream, slen);

	destroy_bamCore(core[0]);
	destroy_bamCore(core[1]);
}

void bam_write_pair(struct read_inf *read1, struct read_inf *read2, Candidate *r,
	unsigned int proper, unsigned int first, unsigned int rev,
	      unsigned int mrev, char *stream, unsigned int *slen)
{
	if (WRITE_BAM == 1)
		return;

	unsigned int flag[2];
	unsigned int i;

	// Make flag
	flag[0] = make_flag(1, proper, r->n[0] == 0? 1: 0, r->n[1] == 0? 1: 0,
							     rev, mrev, 1, 0);

	flag[1] = make_flag(1, proper, r->n[1] == 0? 1: 0, r->n[0] == 0? 1: 0,
							     mrev, rev, 0, 1);

	// Write alignment
	write_alignment_pair(read1, read2, flag, r, 0, stream, slen);
	flag[0] |= 1 << 8;
	flag[1] |= 1 << 8;

	if (proper == 0 || r->n[0] < 2)
		return;

	// Secondary aligment
	for (i = 1; i < r->n[0]; ++i){
		while (i < r->n[0] && (r->pair & 1) == 0 && 
			check_overlapGene(r->pos[0][i], r->pos[0][i - 1]) == 1)
			++i;

		if (i == r->n[0]) break;
		write_alignment_pair(read1, read2, flag, r, i, stream, slen);
	}
}

void bam_write_single(struct read_inf *read, Candidate *r, unsigned int rev,
	          unsigned int p, char *stream, unsigned int *slen)
{
	if (WRITE_BAM == 1)
		return;

	unsigned int flag;
	unsigned int i;

	// Make flag
	flag = make_flag(0, 0, r->n[p] == 0? 1: 0, 1, rev, 0, 1, 0);

	// Write alignment
	write_alignment_single(read, flag, r, p, 0, stream, slen);
	flag |= 1 << 8;

	if (r->n[p] < 2)
		return;

	// Secondary aligment
	for (i = 1; i < r->n[p]; ++i){
		while (i < r->n[p] && check_overlapGene(r->pos[p][i], 
					      r->pos[p][i - 1]) == 1)
			++i;

		if (i == r->n[p]) break;
		write_alignment_single(read, flag, r, p, i, stream, slen);
	}
}

void init_bam_header(char *idx_dir, int argc, char *argv[])
{
	char buf[5000], *ref;
	unsigned int i, len, l;
	FILE *idx;

	// Get reference information
	len = strlen(idx_dir);
	if (idx_dir[len-1] == '/')
		--len;
	memcpy(buf, idx_dir, len);
	memcpy(buf + len, "/reference.inf\0", 15);
	idx = fopen(buf, "rb");

	// Magic word
	memcpy(buf, "BAM\001", CHUNK);
	if (WRITE_BAM == 0)
		bam_write(BAM, buf, CHUNK);

	// Command line
	sprintf(buf, "@HD\tVN:1.0\tSO:unsorted\n@PG\tID:hera\tPN:hera\tVN:%s\tCL:%s%c", 
		HERA_VERSION, argv[0], '\0');
	len = strlen(buf);
	l = len - strlen(argv[0]);

	for (i = 1; i < argc; ++i){
		sprintf(buf + len, " %s%c", argv[i], '\0');
		len = strlen(buf);
	}
	buf[len++] = '\n';

	// Reference
	i = fread(&l, 1, sizeof(int), idx);
	ref = malloc(l);
	i = fread(ref, l, 1, idx);
	
	if (WRITE_BAM == 0){
		bam_write(BAM, &len, CHUNK);
		bam_write(BAM, buf, len);
		bam_write(BAM, ref, l);
	}

	// Contruct gene map
	contruct_gene_map(idx);

	free(ref);
	fclose(idx);
}