#ifndef _BOW_PARSER_
#define _BOW_PARSER_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "align.h"
#include "utility.h"
#include "transcript.h"

struct align_batch_t *parse_bam(const char *path, int thread_cnt, struct ref_info_t *ref);

#endif