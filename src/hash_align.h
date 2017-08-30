#ifndef HASH_ALIGN_H_
#define HASH_ALIGN_H_

#if defined(USE_JEMALLOC)
#include "../lib/jemalloc/include/jemalloc/jemalloc.h"
#define malloc(size) je_malloc(size)
#define calloc(count,size) je_calloc(count,size)
#define realloc(ptr,size) je_realloc(ptr,size)
#define free(ptr) je_free(ptr)
#endif

#if defined(DEBUG)
#define DEBUG_PRINT(fmt, args...) fprintf(stdout, fmt, ##args)
#define ERROR_PRINT(fmt, args...) fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, __FILE__, __LINE__, __func__,##args)
#else
#define DEBUG_PRINT(fmt, args...)
#define ERROR_PRINT(fmt, args...)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>

#include "../lib/zlib/zlib.h"
#include "../lib/libdivsufsort/include/divsufsort64.h"
#include "bgzf.h"
#include "kseq.h"
#include "xxhash.h"

#define GROUP 10000
#define FORCE_EXIT true
#define STACK_SIZE 16777216

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY

KSEQ_INIT(gzFile, gzread)
#endif

#define VERSION "hera-v1.0.1\0"

#define ALPHABET "ACGT"
#define BLOCK 32
#define HBLOCK 16
#define CHUNK 4
#define N_CHAR 4
#define KMER 31
#define ERR 7
#define MIN_FRAG 10
#define MAX_FRAG 1000
#define MIL 33554432

#define POP_HEAD(id) ((id << (2 * (BLOCK - KMER + 1))) >> (2 * (BLOCK - KMER)))
#define MIN3(a, b, c) (a < b? (a < c? a: c): (a < c? b: (c < b? c: b)))
#define ABS(x) ((int)(x)<0 ? -(x) : (x))
#define MIN2(a, b) (a < b? a : b)
#define MAX2(a, b) (a > b? a : b)

/* Cigar operator */
#define CIGAR "MIDNSHP=XB"
#define MATCH	0
#define INSERT	1
#define DELETE	2
#define SKIP	3
#define SCLIP	4
#define HCLIP	5
#define PAD	6
#define EQUAL	7
#define DIFF	8
#define BACK	9

#define MAGIC  "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0"

#define MAX_FIND 1000000
#define MAX_SPLIT 2000000
#define MAX_COMPRESS 65535
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

#define GZIP_ID1 31
#define GZIP_ID2 139
#define CM_DEFLATE 8
#define FLG_FEXTRA 4
#define OS_UNKNOWN 255
#define BGZF_ID1 66
#define BGZF_ID2 67
#define BGZF_LEN 2
#define BGZF_XLEN 6
#define GZIP_WINDOW_BITS -15
#define Z_DEFAULT_MEM_LEVEL 8

#define bam_open(fn, mode) fopen(fn, mode)
#define bam_close(fp) hera_close(fp)
#define bam_write(fp, buf, size) hera_write(fp, buf, size)

typedef FILE *bamFile;

typedef struct {
	unsigned short n_cigar;		// Number of type
	unsigned short n_miss;		// Number of miss char
	unsigned short *count;		// Length of each type
	char *type; 			// Cigar type
	char *miss;
} Cigar;

typedef struct {
	int *p[2];			// Real position on transcript
	unsigned int *pos[2]; 		// Transcript id
	unsigned int n[2];		// Number of transcript/position
	int err[2];			// Score for alignment
	Cigar *cigar[2];		// Cigar string
	int pair;
} Candidate;

typedef struct {
	unsigned int n; 		// Number of transcript
	unsigned int *ref;		// Transcript id
	double count; 			// Number of read
	unsigned int id;
} Class;

typedef struct {
	unsigned int n;
	unsigned int *map;
	Class *cls;
} Gclass;

typedef struct {
	unsigned int n; 		// Number of transcript
	unsigned int name_len;		// Len of each reference name
	char *ref_name; 		// Reference name
	unsigned int *len; 		// Transcipt length
	char *seq;			// Transcipt sequence
	unsigned int *find;		// Where to find a transcipt sequence
	double *eff_len;		// Effective length
	unsigned int *count; 		// Read count
} Ref_inf;

typedef struct {
	unsigned int **start;		// Start position of each exon
	unsigned int **end;		// End position of  each exon
	unsigned int *gene;		// Gene that transcript belong
	unsigned int *chr;		// Chromosome that transcript belong
	unsigned short *n_exon; 	// Number of exon in each transcipt
	unsigned int n_chr;		// Number of chromosome
	unsigned int n_gene;		// Number of gene
	unsigned int l_gene;		// Gene name block len
	unsigned int l_chr;		// Chromosome name block len
	unsigned int *chr_len;		// Chromosome length
	char *chr_name;			// Chromosome name
	char *gene_name;		// Gene name	
} Gene_map;

typedef struct {
	unsigned int id;
	unsigned int start;
	unsigned int end;
} Bucket;

typedef struct {
	Bucket **bucket;
	unsigned int *n;
	unsigned int *pos;
	unsigned int *p;
	unsigned int n_pos;
} Kmer_hash;

typedef struct {
	char *name;
	char *qual;
	char *seq;
	unsigned short len;
	unsigned short name_len;
} Read_inf;

typedef struct {
	unsigned int n;
	unsigned int C[N_CHAR];
	unsigned int *sa;
	unsigned int **Occ;
	unsigned int **tCount;
	char *seq;
} FMindex;

typedef struct {
	unsigned int *read;
	unsigned int n;
	unsigned int trans[2];
} Fusion_inf;

typedef struct {
	unsigned int n;
	unsigned int n_read;
	Fusion_inf *detail;
	char **read;
} Fusion;

typedef struct {
        unsigned int id;
        unsigned int len;
        unsigned int cover;
        char *seq;
} Child;

typedef struct {
        unsigned int id;
        unsigned int cover;
        unsigned int len;
} Parent;

typedef struct {
        unsigned long id;
        Parent parent[4];
        Child child[4];
        unsigned short n_child, n_parent;
} Node;

typedef struct {
        unsigned int n;
        Node *node;
} Graph;

typedef struct {
	int trans;
	unsigned int start;
	unsigned int end;
	int count;
} Fusion_pair;

typedef struct {
	Read_inf r1[GROUP], r2[GROUP];
	unsigned int n;
} Thread_data;

typedef struct {
	Read_inf r[2*GROUP];
	unsigned int n;
} Thread_data2;

extern const unsigned int H;
extern double MEAN_FRAG, MEAN_LEN;
extern unsigned int N_READ, NTHREAD, PAIRED, MAPPED, WRITE_BAM, CLASS_MAX, seed;
extern int COMPRESS_LEVEL;
extern pthread_mutex_t LOCK;
extern pthread_rwlock_t LOCK_RW;
extern unsigned int ascii_table[256];

extern Ref_inf *REF_INF;
extern Gene_map *GENE_MAP;
extern Kmer_hash *KMER_HASH;
extern Gclass *CLASS;
extern Fusion *FUSION;
extern FMindex *FMINDEX;
extern bamFile BAM;
extern FILE *SUMMARY;
extern FILE *OUT_FUSION;
extern kseq_t *R1, *R2;

/* Initialize */
void init_hash();
void init_refInf();

/* Destroy */
void destroy_hash();
void destroy_refInf();

/* Hash align */
void get_ref_seq(char *file_path);
void make_class();
void get_alignment_pair(char *fq1, char *fq2);
void get_alignment(char *fq);
void get_position(Read_inf read, Candidate *r, unsigned short str, short space);

/* Genome align */
void indexGenome(char *input, char *output);
FMindex *load_FMindex(char *index_file, unsigned int l, char *genome);
unsigned int query(char *seq, unsigned int len,
			 FMindex *index, unsigned long *ret);
void genome_map(Read_inf r1, Read_inf r2, Candidate **r);

/* Fusion */
void fusion_add(Candidate *r, Read_inf read1, Read_inf read2);
void *assembly_fusion(void *idx);

/* Cigar */
Cigar *init_cigar();
void add_cigar(Cigar *cigar, char type, unsigned short add, char *c);
void concat_cigar(Cigar *des, Cigar *target);
void free_cigar(Cigar *cigar);

/* File I/O */
void write_to_file(char *out_file);
void load_index(char *idx_file, char *genome);

/* Some utility */
double Norm(double *arr, unsigned int n);
void merge_sort(unsigned int **arr1,  unsigned int n1,
		 unsigned int **arr2, unsigned int n2, unsigned int idx);
void chr_coord(unsigned int pos, unsigned int *chr, unsigned int *p);
unsigned int gene_coord(unsigned int trans, int pos, unsigned int *coord);
inline unsigned long get_index(unsigned int start, char *seq);
inline unsigned int hash_get(unsigned long id, unsigned int *ret,
						unsigned int add);
unsigned int check_quanlity(Candidate *r, unsigned int p);
unsigned int check_overlapGene(unsigned int trans1, unsigned int trans2);
Candidate *init_candidate();
void destroy_candidate(Candidate *read);
void reverse_str(char *seq, unsigned int len);
unsigned long reverse(unsigned long kmer);
char *value_to_key(unsigned long value);

Cigar *SW_align(char *ref, unsigned int rstart, char *str, unsigned int ref_len,
					unsigned int len, int *score, short err);

#endif
