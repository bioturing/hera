#ifndef EM_H_
#define EM_H_

#include "../lib/hdf5/include/hdf5.h"
#include "../lib/hdf5/include/hdf5_hl.h"

#include "hash_align.h"

#define COMPRESS 6

typedef struct {
	unsigned int *table;
	unsigned int n;
} Discrete_dist;

typedef struct {
	double *grama;
	double *overlap;
} EM_val;

typedef struct {
    unsigned int order;
    hid_t *bs;
} Thread_data3;

Discrete_dist *DIST;

EM_val *estimate_count(unsigned int *count, unsigned int order);
void write_result(char *out_dir, unsigned int n_bs,
				 EM_val *em_val, char *prefix);

#endif
