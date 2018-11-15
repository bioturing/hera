#ifndef _TRANS_
#define _TRANS_
#include <stdint.h>

#define get_trans_len(t, i) ((t)->len[(i)])
#define get_trans_id(t, i) ((t)->name[(i)])


struct ref_info_t {
	int nref;			// Number of transcripts
	char **name;			// Transcripts name
        int *len;
};

struct ref_info_t read_ref_info(char *filename);
void release_ref_info(struct ref_info_t *ref);


#endif