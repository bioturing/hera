#ifndef _EM_
#define _EM_

#include "model.h"
#include "smodel.h"

#include "result.h"
#include <stdio.h>

void squarem(struct model *m, struct align_batch_t *a, int thread_cnt, int full_loop);

#endif