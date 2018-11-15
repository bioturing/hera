#ifndef _UTILS_H_
#define _UTILS_H_

#if defined(_MSC_VER)
#pragma warning(disable:4996)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>
#if defined(_MSC_VER)
#include <time.h>
#include <windows.h>
#include <getopt.h>
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <semaphore.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif

/* color terminal */
#if defined(_MSC_VER)
#define RED			""
#define GREEN			""
#define YELLOW			""
#define BLUE			""
#define MAGENTA			""
#define CYAN			""
#define WHITE			""
#define RESET			""
#else
#define RED			"\x1B[31m"
#define GREEN			"\x1B[32m"
#define YELLOW			"\x1B[33m"
#define BLUE			"\x1B[34m"
#define MAGENTA			"\x1B[35m"
#define CYAN			"\x1B[36m"
#define WHITE			"\x1B[37m"
#define RESET			"\x1B[0m"
#endif

#define MAX_INT32		2147483647
#define MIN_INT32		-2147483648

#define MASK32			4294967295ULL

#define BUFSZ			4096

#define THREAD_STACK_SIZE	16777216

#define FORWARD			0
#define REVERSE			1
#define LEFT			0
#define RIGHT			1

//#define WRITE_TEXT
//#define RUN_AFTER_ANALYSIS

/*
 * Built in macros
 */

#define __abs(x) 		((x) < 0 ? -(x) : (x))

#define __min(a, b) 		((a) < (b) ? (a) : (b))

#define __max(a, b) 		((a) > (b) ? (a) : (b))

#define __min3(a, b, c)		__min(__min((a), (b)), (c))

#define __max3(a, b, c)		__max(__max((a), (b)), (c))

#define __round_up_32(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))

#define __is_sep(c)		((c) == ' ' || (c) == '\t')

#define normalize_mapq(s)	do {					       \
	if ((s) < 0) (s) = 0;						       \
	if ((s) > 60) (s) = 60;						       \
} while (0)

/*
 * Built-in macros function
 */

#define __ALLOC(ptr, sz)	(ptr) = xmalloc(sizeof(*(ptr)) * (sz))

#define __REALLOC(ptr, sz)	(ptr) = xrealloc((ptr), sizeof(*(ptr)) * (sz))

/* push back val to ptr, ptr has sz element, realloc + 1 */
#define __PUSH_BACK(ptr, sz, val) do {					       \
	assert((sz) >= 0);						       \
	__REALLOC((ptr), (sz) + 1);					       \
	(ptr)[(sz)++] = (val);						       \
} while(0)

#define __FREE_AND_NULL(ptr) do {					       \
	free(p);							       \
	(p) = NULL;							       \
} while (0)

#if defined(_MSC_VER)
#define __VERBOSE(fmt, ...) do {					       \
	fprintf(stderr, fmt, __VA_ARGS__);				       \
	fflush(stderr);							       \
} while (0)

#if defined(NDEBUG)
#define __DEBUG(fmt, ...) 0
#else
#define __DEBUG(fmt, ...) do {						       \
	fprintf(stderr, "DEBUG: " fmt, __VA_ARGS__);			       \
	fflush(stderr);							       \
} while (0)

#endif

#define __WARNING(fmt, ...) fprintf(stderr, "WARNING: " fmt, __VA_ARGS__)

#define __ERROR(fmt, ...) do {						       \
	fprintf(stderr, "ERROR: " fmt "\n", __VA_ARGS__);		       \
	exit(EXIT_FAILURE);						       \
} while(0)
#else
#define __VERBOSE(fmt, args...) fprintf(stderr, fmt, ##args)

#define __DEBUG(fmt, args...)	fprintf(stderr, "DEBUG: " fmt, ##args)

#define __WARNING(fmt, args...) fprintf(stderr, "WARNING: " fmt, ##args)

#define __ERROR(fmt, args...) do {					       \
	fprintf(stderr, "ERROR: " fmt "\n", ##args);			       \
	exit(EXIT_FAILURE);						       \
} while(0)

#endif




#define __PERROR(fmt) do {						       \
	fprintf(stderr, "ERROR: ");					       \
	perror(fmt);							       \
	exit(EXIT_FAILURE);						       \
} while(0)

#define __SWAP(x, y) do {						       \
	assert(sizeof(x) == sizeof(y));					       \
	int8_t temp[sizeof(x)];						       \
	memcpy(temp, &(y), sizeof(x));					       \
	memcpy(&(y), &(x), sizeof(x));					       \
	memcpy(&(x), temp, sizeof(x));					       \
} while (0)

/*
 * Built-in function
 */

/* support openning file with UNICODE path */
FILE *xfopen(const char *file_path, const char *mode);

/* fflush before close file */
void xfclose(FILE *f);

/* check fread function read enough nmemb */
size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* check fwrite function write enough nmemb */
size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* auto remove /n character if found */
ssize_t xgetline(char **str, size_t *size, FILE *stream);

/* get time */
double realtime();

/* make directory if is not exist */
void make_dir(const char *path);

/* reverse compelemnt */
char *get_rev_complement(const char *seq, int len);

/* reverse string */
char *get_rev(const char *seq, int len);

/* return new char* concate s1 and s2 */
char *str_concate(const char *s1, const char *s2);

/* remove redundant / character */
void normalize_dir(char *path);

/* get size of all file from file_path */
size_t fetch_size(char **file_path, int n_file);

/*
 * Global variable
 */

extern int8_t nt4_table[256];
extern char *nt4_char, *rev_nt4_char;

#endif /* _UTILS_H_ */
