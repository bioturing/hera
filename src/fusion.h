#ifndef FUSION_H_
#define FUSION_H_

#include <stdio.h>
#include <stdlib.h>

#include "hash_align.h"

#define ABS(x) ((int)(x)<0 ? -(x) : (x))
#define CHUNK 4
#define TRUE_FUSION 7

typedef struct {
	unsigned int *read;
	unsigned int *left;
	unsigned int *right;
	unsigned int n[3];
	unsigned long id;
} Fusion_inf;

typedef struct {
	unsigned int n;
	unsigned int n_read;
	unsigned int *map;
	Fusion_inf *detail;
	char **read;
} Fusion;

void fusion_add(Candidate *r, Read_inf *read1, Read_inf *read2);

#endif
