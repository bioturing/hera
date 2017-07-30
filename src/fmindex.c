#include "hash_align.h"

char *readGenome(char *input, unsigned int *n, unsigned int code){
	int l;
	unsigned int i;
	gzFile fp;
	kseq_t *seq;
	char *s;

	s = malloc(1);
	*n = 0;
	fp = gzopen(input, "r");
	if (!fp) {
		perror("Error opening Genome file\n");
		exit(EXIT_FAILURE);
	}

	seq= kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0){
		s = realloc(s, l + *n + 1);
		for (i = 0; i < l; ++i, ++*n)
			s[*n] = code == 0? seq->seq.s[i]:
				(char) (ascii_table[seq->seq.s[i]] & 3);
	}
	gzclose(fp);

	return s;
}

void FMIndex(char *seq, unsigned int n, saidx64_t *sa,
		unsigned int *C, unsigned int **Occ, unsigned int **tCount)
{
	unsigned int i, k, c, mod, block;

	for (i = mod = block = 0; i < n; ++i){
		c = sa[i] > 0? seq[sa[i] - 1]: N_CHAR;
		if (c < N_CHAR) ++C[c];

		if (mod == 0)
			for (k = 0; k < N_CHAR; ++k)
				Occ[k][block] = C[k];
		else if (c < N_CHAR)
			tCount[c][block] |= 1 << mod;

		++mod;
		if (mod == BLOCK){
			++block;
			mod = 0;
		}
	}

	for (i = 1; i < N_CHAR; ++i) C[i] += C[i-1];
}

void writeFile(char *filename, saidx64_t *sa, unsigned int n,
		unsigned int *C, unsigned int **Occ, unsigned int **tCount)
{
	unsigned int i, l, x;
	FILE *f = fopen(strcat(filename, ".genome"), "wb");
	if (!f){
		perror("Error in writting file\n");
		exit(EXIT_FAILURE);
	}

	fwrite(&n, sizeof(int), 1, f);
	for (i = 0, l = n/BLOCK + 1; i < N_CHAR; ++i){
		fwrite(Occ[i], sizeof(int), l, f);
		fwrite(tCount[i], sizeof(int), l, f);
	}
	fwrite(C, sizeof(int), N_CHAR, f);
	fclose(f);

	f = fopen(strcat(filename, ".sa"), "wb");
	if (!f){
		perror("Error in writting file\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < n; ++i){
		x = (unsigned int) sa[i];
		fwrite(&x, sizeof(int), 1, f);
	}
	fclose(f);
}

void indexGenome(char *input, char *output)
{
	unsigned char *s;
	unsigned int n, tmp, i, l;
	unsigned int C[N_CHAR], **Occ, **tCount;
	saidx64_t *sa;

	s = readGenome(input, &n, 1);
	sa = calloc(n + 1, sizeof(saidx64_t));

	tmp = divsufsort64(s, sa+1, n);
	sa[0] = n;
	++n;

	Occ = calloc(N_CHAR, sizeof(int*));
	tCount = calloc(N_CHAR, sizeof(int*));
	for (i = 0, l = n/BLOCK + 1; i < N_CHAR; ++i){
		Occ[i] = calloc(l, sizeof(int));
		tCount[i] = calloc(l, sizeof(int));
  	}
	memset(C, 0, N_CHAR*sizeof(int));

	FMIndex(s, n, sa, C, Occ, tCount);
	free(s);

	writeFile(output, sa, n, C, Occ, tCount);
}

unsigned int calOcc(unsigned int **Occ, unsigned int **tCount,
				unsigned int c, unsigned int pos)
{
	unsigned int i, mod, block, add;

	block = pos/BLOCK;
	mod = pos & (BLOCK - 1);

	if (mod == 0) return Occ[c][block];

	for (add = 0, i = 1; i <= mod; ++i)
		add += (tCount[c][block] >> i) & 1;

	return Occ[c][block] + add;
}

unsigned int query(char *seq, unsigned int len,
			FMindex *index, unsigned long *ret)
{
	unsigned int left, right, c, l, r, rank, k;
	int i;
	unsigned int *pos;

	i = len - 1;
	c = ascii_table[seq[i]] & 3;
	left = c == 0? 1: index->C[c-1] + 1;
	right = index->C[c];

	for (--i; i >= 0; --i, left = l, right = r) {
		c = ascii_table[seq[i]] & 3;
		rank = c == 0? 1: index->C[c-1] + 1;
		l = rank + calOcc(index->Occ, index->tCount, c, left-1);
		r = rank + calOcc(index->Occ, index->tCount, c, right) - 1;

		if (l > r )
			return len - i - 1;
		*ret = (long) left << BLOCK | right;
	}

	return len; 
}

FMindex *load_FMindex(char *index_file, unsigned int l, char *genome)
{
	unsigned int i, len, ret;
	FILE *index, *saidx;
	FMindex *fmindex = malloc(sizeof(FMindex));

	memcpy(index_file + l, ".genome\0", 8);
	index = fopen(index_file, "rb");

	memcpy(index_file + l + 7, ".sa\0", 4);
	saidx = fopen(index_file, "rb");

	if (!index || !saidx){
		printf("Can not open genome file or file does not exist. You should run index command with option --full-index\n");
		exit(EXIT_FAILURE);
	}

	ret = fread(&fmindex->n, sizeof(int), 1, index);
	fmindex->Occ = calloc(N_CHAR, sizeof(int*));
	fmindex->tCount = calloc(N_CHAR, sizeof(int*));
	fmindex->sa = calloc(fmindex->n, sizeof(int));

	for (i = 0, l = fmindex->n/BLOCK + 1; i < N_CHAR; ++i){
		fmindex->Occ[i] = calloc(l, sizeof(int));
		fmindex->tCount[i] = calloc(l, sizeof(int));
		ret = fread(fmindex->Occ[i], sizeof(int), l, index);
		ret = fread(fmindex->tCount[i], sizeof(int), l, index);
	}
	ret = fread(fmindex->C, sizeof(int), N_CHAR, index);
	ret = fread(fmindex->sa, sizeof(int), fmindex->n, saidx);
	fclose(index);
	fclose(saidx);

	fmindex->seq = readGenome(genome, &len, 0);

	return fmindex;
}
