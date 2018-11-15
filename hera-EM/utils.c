#include "utils.h"

int8_t nt4_table[256] = {
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4
};

char *nt4_char = "ACGTN";
char *rev_nt4_char = "TGCAN";

FILE *xfopen(const char *file_path, const char *mode) {
	FILE *fi = NULL;
#ifdef WIN32
	fi = fopen(file_path, mode);
#else
	fi = fopen(file_path, mode);
#endif
	return fi;
}

void xfclose(FILE *f)
{
	fflush(f);
	fclose(f);
}

size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__ERROR("fread, wrong file or file is corrupted");
	return ret;
}

size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__ERROR("fwrite, could not write data to file");
	return ret;
}

ssize_t xgetline(char **str, size_t *size, FILE *stream)
{
	ssize_t ret = -1;
#if defined(_MSC_VER)
	ret = _getline(str, size, stream);
#else
	ret = getline(str, size, stream);
#endif
	if (ret == 0 || ret == -1)
		return ret;
	if ((*str)[ret - 1] == '\n')
		(*str)[--ret] = '\0';
	return ret;
}

char *str_concate(const char *str1, const char *str2)
{
	size_t len1 = strlen(str1);
	size_t len2 = strlen(str2);
	char *str3 = malloc(len1 + len2 + 1);
	strcpy(str3, str1);
	strcpy(str3 + len1, str2);
	str3[len1 + len2] = '\0';
	return str3;
}

void normalize_dir(char *path)
{
	int len = strlen(path), i, j;
	for (i = 0; i < len - 1; ) {
		if (path[i] == '/' && path[i + 1] == '/') {
			for (j = i; j < len; ++j)
				path[j] = path[j + 1];
			--len;
		} else {
			++i;
		}
	}
}

void make_dir(const char *path)
{
	struct stat st = {0};
	if (stat(path, &st) == -1) {
		if (mkdir(path, 0700)) {
			perror("Could not make output directory");
			exit(EXIT_FAILURE);
		}
	}
}

char *get_rev(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	int i, k;
	char *ret = malloc(len + 1);
	for (i = 0, k = len - 1; i < len; ++i, --k)
		ret[i] = seq[k];
	ret[len] = '\0';
	return ret;
}

char *get_rev_complement(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	char *ret = malloc(len + 1);
	int i, k;
	for (i = 0, k = len - 1; i < len; ++i, --k) 
		ret[i] = rev_nt4_char[nt4_table[(int)seq[k]]];
	ret[len] = '\0';
	return ret;
}

double realtime()
{
	struct timeval tp;
#if defined(_MSC_VER)
	_gettimeofday(&tp,  NULL);
#else
	gettimeofday(&tp,  NULL);
#endif
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

size_t fetch_size(char **file_path, int n_file)
{
	size_t ret = 0;
	FILE *fid;
	int i;
	for (i = 0; i < n_file; ++i) {
		fid = xfopen(file_path[i], "rb");
		if (!fid)
			__ERROR("Unable to open file: %s", file_path[i]);
		fseek(fid, 0L, SEEK_END);
		ret += ftell(fid);
		xfclose(fid);
	}
	return ret;
}
