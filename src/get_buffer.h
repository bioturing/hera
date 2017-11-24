#ifndef GET_BUFFER_H_
#define GET_BUFFER_H_

#include <stdio.h>
#include <stdlib.h>
// #include <zlib.h> zlib of system
#include "../lib/zlib/zlib.h"
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "macros.h"

enum gb_file_type {
	TYPE_FASTQ,
	TYPE_FASTA
};

struct gb_file_inf {
	gzFile file;
	int32_t end;
	bool finish_read;
	char *buffer;
	char *file_name;
};

/* pair */

struct gb_pair_data {
	struct gb_file_inf *left, *right;
	bool finish_flag;
	bool warning_flag;
	enum gb_file_type format;
};

struct gb_pair_buf {
	char *left, *right;
};

struct gb_pair_data *gb_init_pair(char *file_name1, char *file_name2);

struct gb_pair_buf *gb_get_pair(struct gb_pair_data *gb_data);

void gb_destroy_pair_buf(struct gb_pair_buf *data);

void gb_destroy_pair_data(struct gb_pair_data *data);

/* single */

struct gb_single_data {
	struct gb_file_inf *single;
	bool finish_flag;
	enum gb_file_type format;
};

struct gb_single_buf {
	char *single;
};

struct gb_single_data *gb_init_single(char *file_name);

struct gb_single_buf *gb_get_single(struct gb_single_data *gb_data);

void gb_destroy_single_buf(struct gb_single_buf *data);

void gb_destroy_single_data(struct gb_single_data *data);

#endif /* GET_BUFFER_H_ */