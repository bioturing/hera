#include "align.h"
#include "utility.h"
#include <stdlib.h>

void align_batch_release(struct align_batch_t *a)
{
	free(a->align);
}


void batch_resize(struct align_batch_t *a)
{
	int new_cap = a->size;
	a->align = realloc(a->align, new_cap * sizeof(*a->align));
	a->cap = new_cap;
}
