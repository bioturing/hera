#ifndef BAM_WRITE_H_
#define BAM_WRITE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hash_align.h"

/* Read flag */
#define MUL 0x1			// Reads have its mate 		
#define PROPER 0x2		// Read and its mate are proper pair 
#define UNMAP 0x4		// Read is unmapped
#define M_UNMAP 0x8		// Read's mate is unmapped
#define REVERSE 0x10		// Read is reversed
#define M_REVERSE 0x20		// Read's mate is reversed
#define R1 0x40			// Read is first read
#define R2 0x80			// Read is second read
#define SECOND 0x100		// Alginment is chimeric
#define DUP 0x400		// Read is pcr duplicated
#define SUP 0x800		// Read is supplement read

typedef struct {
	int bin, qual, ref, pos, n_cigar;
	unsigned int len, block_len, coord, match, align_len;
	char *miss, *data;
	unsigned int *cstring;
} Bam_core;

extern unsigned int BUFF_LEN;
extern char *BAM_BUF;

void init_output(char *idx_dir, char *out_dir, int argc, char *argv[]);
void bam_write_pair(Read_inf read1, Read_inf read2, Candidate *r,
	unsigned int proper, unsigned int first, unsigned int rev,
	     unsigned int mrev, char *stream, unsigned int *slen);
void bam_write_single(Read_inf read, Candidate *r, unsigned int rev,
	          unsigned int p, char *stream, unsigned int *slen);

#endif
