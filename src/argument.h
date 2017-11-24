#ifndef ARGUMENT_H_
#define ARGUMENT_H_

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "usage.h"

struct arg_quant_data {
	int32_t nLeft, nRight;
	char **left_file, **right_file;
	int32_t threads;
	int32_t bootstrap;
	char *idx_dir;
	char *out_dir;
	char *genome;
	char *prefix;
	bool bam;
	bool verbose;
	int32_t compress;
};

struct arg_index_data {
	char *transcript;
	char *genome;
	char *out_dir;
};

struct arg_quant_data *arg_get_quant(int32_t pos, int32_t argc, char *argv[]);

struct arg_index_data *arg_get_index(int32_t pos, int32_t argc, char *argv[]);

#endif /* ARGUMENT_H_ */