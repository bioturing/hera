#ifndef FMINDEX_H_
#define FMINDEX_H_

#include <divsufsort64.h>//"../lib/libdivsufsort/include/divsufsort64.h"
#include "hash_align.h"

typedef struct {
	unsigned int n;
	unsigned int C[N_CHAR];
	unsigned int *sa;
	unsigned int **Occ;
	unsigned int **tCount;
	char *seq;
} FMindex;

void indexGenome(char *input, char *output);
FMindex *load_FMindex(char *index_file, unsigned int l, char *genome);
unsigned long query(char *seq, unsigned int len, FMindex *index);
#endif
