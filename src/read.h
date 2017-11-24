#ifndef READ_H_
#define READ_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

enum read_exit_code {
	READ_SUCCESS,
	READ_END,
	READ_FAIL
};

struct read_inf {
	char *seq;
	char *name;
	char *qual;
	int32_t len;
};

void free_read_inf(struct read_inf *read);

enum read_exit_code get_read_from_fq_buffer(struct read_inf *read,
					    char *buffer,
					    int32_t *pos);

enum read_exit_code get_read_from_fa_buffer(struct read_inf *read,
					    char *buffer,
					    int32_t *pos);

#endif /* HREAD_H_ */