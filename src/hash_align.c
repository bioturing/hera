#include "hash_align.h"
#include "ssw.h"

unsigned int ascii_table[256] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

const unsigned int H = (int) 1 << 27;
unsigned int NTHREAD = 1;
int COMPRESS_LEVEL = -1;
double MEAN_FRAG = 0;
double MEAN_LEN = 0;
unsigned int N_READ = 0;
unsigned int PAIRED = 0;
unsigned int MAPPED = 0;
unsigned int WRITE_BAM = 0;
unsigned int seed = 0;

Gclass *CLASS;
pthread_mutex_t LOCK = PTHREAD_MUTEX_INITIALIZER;
pthread_rwlock_t LOCK_RW = PTHREAD_RWLOCK_INITIALIZER;

Ref_inf *REF_INF;
Gene_map *GENE_MAP;
Kmer_hash *KMER_HASH;
Fusion *FUSION;
FMindex *FMINDEX;
bamFile BAM;
FILE *SUMMARY;
FILE *OUT_FUSION;
kseq_t *R1, *R2;

void hera_write(bamFile bam_file, void *input_buf, unsigned int block_length)
{
	unsigned int crc, remaining, compressed_length, input_length, total;
	uint8_t *buffer;

	total = MAX_COMPRESS;
	buffer = malloc(total);

	// Init gzip header
	buffer[0] = GZIP_ID1;
	buffer[1] = GZIP_ID2;
	buffer[2] = CM_DEFLATE;
	buffer[3] = FLG_FEXTRA;
	buffer[4] = 0;
	buffer[5] = 0;
	buffer[6] = 0;
	buffer[7] = 0;
	buffer[8] = 0;
	buffer[9] = OS_UNKNOWN;
	buffer[10] = BGZF_XLEN;
	buffer[11] = 0;
	buffer[12] = BGZF_ID1;
	buffer[13] = BGZF_ID2;
	buffer[14] = BGZF_LEN;
	buffer[15] = 0;
	buffer[16] = 0; // placeholder for block length
	buffer[17] = 0;

	// loop to retry for blocks that do not compress enough
	input_length = block_length;
	compressed_length = 0;
	while (1) {
		z_stream zs;
		zs.zalloc = NULL;
		zs.zfree = NULL;
		zs.next_in = input_buf;
		zs.avail_in = input_length;
		zs.next_out = (void*)&buffer[BLOCK_HEADER_LENGTH];
		zs.avail_out = total - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
		int status = deflateInit2(&zs, COMPRESS_LEVEL, Z_DEFLATED,
				    GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL,
						      Z_DEFAULT_STRATEGY);

		if (status != Z_OK) {
			ERROR_PRINT("Deflate init failed\n");
			exit(EXIT_FAILURE);
		}

		status = deflate(&zs, Z_FINISH);
		if (status != Z_STREAM_END) {
			deflateEnd(&zs);
			if (status == Z_OK) {
				input_length -= 1024;
				if (input_length <= 0) {
					ERROR_PRINT("Input reduction failed\n");
					exit(EXIT_FAILURE);
				}
				continue;
			}
			ERROR_PRINT("Deflate failed\n");
			exit(EXIT_FAILURE);
		}

		status = deflateEnd(&zs);
		if (status != Z_OK) {
			ERROR_PRINT("Deflate end failed\n");
			exit(EXIT_FAILURE);
		}
		compressed_length = zs.total_out;
		compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
		break;
	}

	packInt16((uint8_t*)&buffer[16], compressed_length-1);
	crc = crc32(0L, NULL, 0L);
	crc = crc32(crc, input_buf, input_length);
	packInt32((uint8_t*)&buffer[compressed_length-8], crc);
	packInt32((uint8_t*)&buffer[compressed_length-4], input_length);

	remaining = block_length - input_length;
	if (remaining > 0) {
		if (remaining > input_length) {
			ERROR_PRINT("Remainder too large\n");
			exit(EXIT_FAILURE);
		}
		memcpy(input_buf, input_buf + input_length, remaining);
	}

	pthread_rwlock_wrlock(&LOCK_RW);
	fwrite(buffer, 1, compressed_length, bam_file);
	pthread_rwlock_unlock(&LOCK_RW);
	free(buffer);

	if (remaining > 0)
		hera_write(bam_file, input_buf, remaining);
}

void hera_close(bamFile bam_file)
{
	fwrite(MAGIC, 1, 28, bam_file);
	fclose(bam_file);
}

/******************************************************************************
 * ************************************************************************** *
 * *                               UTILITY                                  * *
 * ************************************************************************** * 
 ******************************************************************************/

double Norm(double *arr, unsigned int n)
{
	unsigned int i;
	double sum, total;

	for (i = 0, sum = 0; i < n; i++)
		sum += arr[i];

	sum /= n;
	for (i = 0, total = 0; i < n; i++)
		total += (arr[i] - sum)*(arr[i] - sum);
	total /= n - 1;

	return sqrt(total);
}

unsigned long key_to_value(char *key)
{
	unsigned long value = 0;
	unsigned short i;
	unsigned int x;

	for (i = 0; i < KMER; i++) {
		x = ascii_table[key[i]];
		if (x > 3)
			x = rand()&3;
		value = (value << 2) | (x&3);
	}

	return value;
}

char *value_to_key(unsigned long value)
{
	char *key = malloc(KMER + 1);
	unsigned short i;

	for (i = 0; i < KMER; i++, value >>= 2)
		key[KMER - 1 - i] = ALPHABET[value & 3];
	key[KMER] = '\0';
	return key;
}

unsigned long reverse(unsigned long kmer)
{
	unsigned long rev = 0;
	unsigned int i;

	for (i = 0; i < KMER; i++) {
		rev <<= 2;
		rev |= ~kmer & 3;
		kmer >>= 2;
	}
	return rev;
}

void reverse2(char *seq, unsigned int len)
{
        unsigned short i;
        char tmp;
        for (i = 0; i < (len>>1); ++i){
                tmp = seq[len-1-i];
                seq[len-1-i] = seq[i];
                seq[i] = tmp;
        }
}

void reverse_str(char *seq, unsigned int len)
{
	unsigned short i;
	char tmp;

	for (i = 0; i < (len>>1); ++i){
		tmp = ascii_table[seq[i]] > 3 ? seq[i]:
                                    ALPHABET[3 - ascii_table[seq[i]]];
		seq[i] = ascii_table[seq[len-1-i]] > 3 ? seq[len-1-i]:
                              ALPHABET[3 - ascii_table[seq[len-1-i]]];
		seq[len-1-i] = tmp;
	}

	if (len & 1 == 1 && ascii_table[seq[len >> 1]] <= 3)
		seq[len >> 1] = ALPHABET[3 - ascii_table[seq[len >> 1]]];
}


void reverse_qual(char *seq, unsigned int len)
{
	unsigned short i;
	char tmp;

	for (i = 0; i < (len>>1); ++i){
		tmp = seq[i];
		seq[i] = seq[len-1-i];
		seq[len-1-i] = tmp;
	}
}

void add_sort(unsigned int **arr1, unsigned int **arr2, unsigned int idx,
				        unsigned int *i, unsigned int *k)
{
	arr1[0][*i] = arr2[0][*k];
	arr1[1][*i] = arr2[1][*k];
	arr1[2][*i] = idx;
	++*i;
	++*k;
}

void merge_sort(unsigned int **arr1,  unsigned int n1,
		 unsigned int **arr2, unsigned int n2, unsigned int idx)
{
	unsigned int i, k, n, dis, *new_arr[3];

	new_arr[0] = calloc(n1 + n2, sizeof(int));
	new_arr[1] = calloc(n1 + n2, sizeof(int));
	new_arr[2] = calloc(n1 + n2, sizeof(int));

	for (i = k = n = 0; i < n1 && k < n2;){
		if (arr1[0][i] < arr2[0][k]){
			add_sort(new_arr, arr1, arr1[2][i], &n, &i);
		} else if (arr1[0][i] > arr2[0][k]){
			add_sort(new_arr, arr2, idx, &n, &k);
		} else {
			dis = ABS((arr2[1][k] - arr1[1][i]) -
					  (idx - arr1[2][i]));
			if (dis < 2*ERR){
				add_sort(new_arr, arr1, arr1[2][i], &n, &i);
			} else if (arr1[1][i] > arr2[1][k]){
				add_sort(new_arr, arr2, idx, &n, &k);
			} else {
				add_sort(new_arr, arr1, arr1[2][i], &n, &i);
			}
		}
	}

	while (i < n1)
		add_sort(new_arr, arr1, arr1[2][i], &n, &i);

	while (k < n2)
		add_sort(new_arr, arr2, idx, &n, &k);

	free(arr1[0]);
	free(arr1[1]);
	free(arr1[2]);

	arr1[0] = new_arr[0];
	arr1[1] = new_arr[1];
	arr1[2] = new_arr[2];
}

Cigar *init_cigar()
{
	Cigar *cigar = malloc(sizeof(Cigar));
	cigar->n_cigar = cigar->n_miss = 0;

	return cigar;
}

void chr_coord(unsigned int pos, unsigned int *chr, unsigned int *p)
{
	unsigned int sum, i;

	for (i = sum = 0; i < GENE_MAP->n_chr; ++i){
		if (sum + GENE_MAP->chr_len[i]> pos)
			break;
		sum += GENE_MAP->chr_len[i];
	}
	*chr = i;
	*p = pos - sum;
}

unsigned int gene_coord(unsigned int trans, int pos, unsigned int *coord)
{
	unsigned int len, i;

	if (pos <= 0){
		*coord = 0;
		return GENE_MAP->start[trans][*coord];
	}

	for (i = len = 0; len <= pos; ++i)
		len += GENE_MAP->end[trans][i] - GENE_MAP->start[trans][i] + 1;

	*coord = i - 1;
	return GENE_MAP->end[trans][*coord] - (len - pos) + 1;
}

void swap_trans(Candidate *r1, unsigned int p1, 
			Candidate *r2, unsigned int p2)
{
	unsigned int *tmp, n;

	tmp = r1->pos[p1];
	n = r1->n[p1];

	r1->pos[p1] = r2->pos[p2];
	r1->n[p1] = r2->n[p2];

	r2->pos[p2] = tmp;
	r2->n[p2] = n;
}

void check_mate(Candidate *r1, Read_inf read1, Read_inf read2,
	   		      unsigned int len, unsigned site)
{
	if (len < KMER)
		return;

	Candidate *r2;
	Read_inf anchor;

	r2 = init_candidate();
	anchor.len = len;
	anchor.seq = read1.seq + (site > 0? read1.len - len: 0);

	get_position(anchor, r2, 0, ERR);
	get_position(anchor, r2, 1, ERR);

	if (r2->err[0] < r2->err[1])
		swap_trans(r1, 0, r2, 0);
	else
		swap_trans(r1, 0, r2, 1);

	fusion_add(r1, read1, read2);

	if (r2->err[0] < r2->err[1])
		swap_trans(r1, 0, r2, 0);
	else
		swap_trans(r1, 0, r2, 1);

	destroy_candidate(r2);
}

void check_proper_pair(unsigned int i, Candidate *r,
				     Read_inf read1, Read_inf read2)
{
	unsigned int n;

	n = r->cigar[1]->n_cigar - 1;
	if (r->cigar[0]->type[0] == SCLIP){
		check_mate(r, read1, read2, r->cigar[0]->count[0], 0);
	} else if (r->cigar[1]->type[n] == SCLIP){
                reverse_str(read2.seq, read2.len);
		check_mate(r, read2, read1, r->cigar[1]->count[n], 1);
                reverse_str(read2.seq, read2.len);
	}	
}

unsigned int check_overlapGene(unsigned int trans1, unsigned int trans2)
{
	unsigned int start[2], end[2];
	if (trans1 == trans2 ||
	    GENE_MAP->gene[trans1] == GENE_MAP->gene[trans2])
		return 1;

	if (GENE_MAP->chr[trans1] != GENE_MAP->chr[trans2])
		return 0;

	start[0] = GENE_MAP->start[trans1][0];
	start[1] = GENE_MAP->start[trans2][0];

	end[0] = GENE_MAP->end[trans1][GENE_MAP->n_exon[trans1] - 1];
	end[1] = GENE_MAP->end[trans2][GENE_MAP->n_exon[trans2] - 1];

	if (start[1] > end[0] || start[0] > end[1])
		return 0;
	return 1;
}

unsigned int check_quanlity(Candidate *r, unsigned int p)
{
	unsigned int i, k;

	if (r->n[p] == 0)
		return 0;

	for (i = k = 1; i < r->n[p]; ++i)
		if (check_overlapGene(r->pos[p][i], r->pos[p][i-1]) == 0)
			++k;
	return k;
}

void store_fusion(Candidate *r, Read_inf read1, Read_inf read2)
{
	unsigned int i, k;
	if (r->pair > 1)
		check_proper_pair(i, r, read1, read2);
	else
		fusion_add(r, read1, read2);
}

unsigned int *intersect(unsigned int *arr1, unsigned int *arr2,
                                     unsigned int *size1, unsigned int *size2,
                                                    int **p, unsigned int add)
{
	unsigned int i, k, j, l, n, dis, has;
	unsigned int *ret = calloc(1, sizeof(int));
	int *pos = calloc(1, sizeof(int));

	i = k = n = 0;
	while (i < *size1 && k < *size2) {
		if (arr1[i] < arr2[k]){
			++i;
		} else if (arr1[i] > arr2[k]) {
			++k;
		} else {
			for (j = i, has = 0; j < *size1 &&
                                        arr1[j] == arr1[i] && has == 0; ++j) {
				for (l = k; l < *size2 &&
                                                    arr2[l] == arr2[k]; ++l) {
					dis = ABS(p[0][j] - p[1][l]);
					if (dis == add) {
						++n;

						ret = realloc(ret,
                                                            n * sizeof (int));
						ret[n - 1] = arr1[j];
						pos = realloc(pos,
                                                            n * sizeof (int));
						pos[n - 1] = p[0][j];
						has = 1;
						++l;
						break;
					} else if (dis > add &&
                                                            p[1][l] > p[0][j])
						break;
				}
			}
			i = j;
			k = l;
		}
	}

	if (n > 0) {
		*size1 = n;
		free(p[0]);
		free(arr1);
		p[0] = pos;
		return ret;
	}

	*size2 = 0;
	free(ret);
	free(pos);
	return arr1;
}

unsigned int intersect3(unsigned int *pos, unsigned int *p, unsigned int *len,
					     Candidate *r, unsigned short str)
{
	unsigned int i, k, j, m, n, dis;

	i = k = n = 0;
	while (i < *len && k < r->n[str]) {
		if (pos[i] < r->pos[str][k]){
			++i;
		} else if (pos[i] > r->pos[str][k]){
			++k;
		} else {
			++n;
			pos[n - 1] = pos[i];
			p[n - 1] = p[i];
			r->pos[str][n - 1] = r->pos[str][k];
			r->p[str][n - 1] = r->p[str][k];

			++i;
			++k;
		}
	}

	if (n > 0)
		*len = n;
	return n;
}

/******************************************************************************
 * ************************************************************************** *
 * *                            PAIRWISE ALIGN                              * *
 * ************************************************************************** * 
 ******************************************************************************/

/******************************************************************************
 * Append align result to cigar string                                        *
 * @ cigar: cigar struct to store cigar string                                *
 * @ type : type of variant (match, miss match, insert, delete...)            *
 * @ add  : length of variant                                                 *
 * @ c    : store the origin nucleotide of reference (NULL if match)          *
 ******************************************************************************/
void add_cigar(Cigar *cigar, char type, unsigned short add, char *c)
{
	if (add == 0)
		return;

	unsigned int n = cigar->n_cigar;

	if (n > 1 && cigar->type[n-1] == SCLIP){
		if (type != DELETE)
			cigar->count[n-1] += add;
		return;
	}

	if (c != NULL){
		if (cigar->n_miss == 0)
			cigar->miss = malloc(1);
		else
			cigar->miss = realloc(cigar->miss, cigar->n_miss + add);
		memcpy(cigar->miss + cigar->n_miss, c, add);
		cigar->n_miss += add;
	}


	if (n > 0 && cigar->type[n-1] == type){
		cigar->count[n-1] += add;
		return;
	}

	++cigar->n_cigar;
	if (cigar->n_cigar == 1){
		cigar->count = calloc(1, sizeof(short));
		cigar->type = malloc(1);
	} else {
		cigar->count = realloc(cigar->count,
                                     cigar->n_cigar*sizeof(short));
		cigar->type = realloc(cigar->type, cigar->n_cigar);
	}
	cigar->type[n] = type;
	cigar->count[n] = add;
}

void reverse_cigar(Cigar *cigar)
{
	unsigned int i, tmp;
	char c;

	for (i = 0; i < cigar->n_cigar/2; ++i){
		tmp = cigar->count[i];
		cigar->count[i] = cigar->count[cigar->n_cigar - 1 - i];
		cigar->count[cigar->n_cigar - 1 - i] = tmp;

		c = cigar->type[i];
		cigar->type[i] = cigar->type[cigar->n_cigar - 1 - i];
		cigar->type[cigar->n_cigar - 1 - i] = c;
	}

	for (i = 0; i < cigar->n_miss/2; ++i){
		c = cigar->miss[i];
		cigar->miss[i] = cigar->miss[cigar->n_miss - 1 - i];
		cigar->miss[cigar->n_miss - 1 - i] = c;
	}
}

void concat_cigar(Cigar *des, Cigar *target)
{
	unsigned int i, k;

	reverse_cigar(target);
	for (i = k = 0; i < target->n_cigar; ++i){
		if (target->type[i] == DELETE || target->type[i] == DIFF){
			add_cigar(des, target->type[i],
				      target->count[i], target->miss + k);
			k += target->count[i];
		} else {
			add_cigar(des, target->type[i], target->count[i], NULL);
		}
	}
}

void free_cigar(Cigar *cigar)
{
	if (cigar->n_cigar > 0){
		free(cigar->count);
		free(cigar->type);
	}

	if (cigar->n_miss > 0)
		free(cigar->miss);

	free(cigar);
}

/******************************************************************************
 * Smith Waterman global alignment with modify to allowed error               *
 * @ ref   : reference string (1 char = 2 bit)                                *
 * @ start : start position on reference to algin                             *
 * @ str   : query string                                                     *
 * @ len   : query string length                                              *
 * @ score : return alignment score                                           *
 * @ err   : maximum allowed error (err < 0: query string is read preffix)    *
 ******************************************************************************/
Cigar *SW_align(char *ref, unsigned int rstart, char *str, unsigned int ref_len,
					unsigned int len, int *score, short err)
{
	int i, k, j, ferr;
	unsigned int s, pen, min, xmin, add;
	unsigned short matrix[len + 1][ref_len + 1];
	unsigned short error = ABS(err) + 1;
	char r[ref_len], q[len];

	memcpy(r, ref + rstart, ref_len);
	memcpy(q, str, len);
	if (err < 0){
		reverse2(r, ref_len);
		reverse2(q, len);
	}

	s = *score = 0;
	Cigar *cigar = init_cigar();
	if (len == 0)
		return cigar;

	for (i = 0; i < len + 1; ++i)
		memset(matrix[i], 0, (ref_len + 1)*sizeof(short));

	for (i = 0; i <= error + 1; ++i)
		matrix[0][i] = i;

	for (k = 1, ferr = -1; k <= len; ++k){
		min = error;
		j = 0;
		matrix[k][s] = matrix[k-1][s] + 1;

		for (i = s + 1; i <= ref_len; ++i){
			if (matrix[k-1][i-1] == error + 1){
				matrix[k][i] = error + 1;
				break;
			}

			pen = 0;
			if (matrix[k-1][i-1] < error)
				pen = ascii_table[r[i-1]] ==
				      ascii_table[q[k - 1]] ? 0: 1;

			matrix[k][i] = MIN3(matrix[k-1][i-1] + pen,
                            matrix[k][i-1] + 1, matrix[k-1][i] + 1);

			if (matrix[k][i] == error && j == 0)
				j = i;
			if (matrix[k][i] < error)
				j = 0;
			if (matrix[k][i] <= min){
				min = matrix[k][i];
				xmin = i;
			}
			if (ferr == -1 && i == k && matrix[k][i] > 0)
				ferr = i;
		}

		if (min == error){
			++k;
			break;
		}

		while (matrix[k][s] >= error && s < i)
			++s;
		if (j > xmin && j < ref_len)
			matrix[k][j] = error + 1;
	}

	--k;
	if (min >= error){
		while (k >= 5 && (min > k/10 || min * 2 > k - ferr)){
			for (i = j = 0, --k; i < ref_len; ++i){
				if (matrix[k][i] == error + 1)
					break;
				if (matrix[k][i] > 0)
					j = error;
				if (matrix[k][i] <= min && j > 0 &&
				    ascii_table[r[i-1]] == ascii_table[q[k-1]] &&
				    matrix[k][i] == matrix[k - 1][i - 1]){
					min = matrix[k][i];
					xmin = i;
				}
			}
		}

		if (k < 5){
			add_cigar(cigar, SCLIP, len, NULL);
			*score = len;
			return cigar;
		}
		add_cigar(cigar, SCLIP, len - k, NULL);
		*score = len - k;
		i = k;
	} else {
		if (xmin < len && k == len){
			min += len - xmin;
			xmin = len;
		}
		i = len;
	}
	k = xmin;
	*score += min;

	if (WRITE_BAM == 1){
		add_cigar(cigar, MATCH, k, NULL);
		return cigar;
	}

	while (i > 0 && k > 0){
		pen = ascii_table[r[k - 1]] == ascii_table[q[i - 1]] ? 0: 1;
		if (matrix[i][k] == matrix[i-1][k-1] + pen){
			if (pen == 0)
				add_cigar(cigar, MATCH, 1, NULL);
			else 
				add_cigar(cigar, DIFF, 1, r + k - 1);
			--i;
			--k;
		} else if (matrix[i][k] == matrix[i - 1][k] + 1){
			add_cigar(cigar, INSERT, 1, NULL);
			--i;
		} else {
			add_cigar(cigar, DELETE, 1, r + k - 1);
			--k;
		}
	}

	if (i > 0)
		add_cigar(cigar, INSERT, i, NULL);

	for (i = 1; i <= k; ++i)
		add_cigar(cigar, DELETE, 1, r + k - i);
	return cigar;
}

unsigned int re_mapSW(Candidate *r, unsigned short str, Read_inf read,
						    	     int type)
{
	unsigned int i, k, new_n, ref_start, len, *pos;
	int score1, score2, range, start, *p, skip;
	Cigar *cigar, *new_cigar;

	cigar = new_cigar = NULL;
	range = type > 0? 0: read.len - MAX_FRAG;

	for (i = new_n = 0; i < r->n[1 - str]; ++i){
		start = r->p[1 - str][i] + range;
		len = MAX_FRAG;

		if (start < 0)
			start = 0;

		if (start + MAX_FRAG > REF_INF->len[r->pos[1 - str][i]])
			len = REF_INF->len[r->pos[1 - str][i]] - start;

		ref_start = r->pos[1 - str][i] == 0? 0:
			    REF_INF->find[r->pos[1 - str][i] - 1];

		cigar = ssw_align(REF_INF->seq + ref_start + start, len, read.seq,
			read.len, &score1, &skip);

		if (score1 > read.len/2 || start + skip < 0 ||
		    start + skip + read.len > REF_INF->len[r->pos[1 - str][i]]){
			free_cigar(cigar);
			continue;
		}

		if (new_cigar == NULL){
			new_cigar = cigar;
			score2 = score1;

			r->pos[1 - str][0] = r->pos[1 - str][i];
			r->p[1 - str][0] = r->p[1 - str][i];
			pos = calloc(r->n[1 - str], sizeof(int));
			p = calloc(r->n[1 - str], sizeof(int));
			pos[0] = r->pos[1 - str][i];
			p[0] = start + skip;
			++new_n;
		} else if (score1 == score2){
			r->pos[1 - str][new_n] = r->pos[1 - str][i];
			r->p[1 - str][new_n] = r->p[1 - str][i];
			pos[new_n] = r->pos[1 - str][i];
			p[new_n] = start + skip;
			++new_n;
			free_cigar(cigar);
		} else if (score1 < score2) {
			free_cigar(new_cigar);
			new_cigar = cigar;
			score2 = score1;

			r->pos[1 - str][0] = r->pos[1 - str][i];
			r->p[1 - str][0] = r->p[1 - str][i];
			pos[0] = r->pos[1 - str][i];
			p[0] = start + skip;			
			new_n = 1;
		} else {
			free_cigar(cigar);
		}
	}

	if (new_n == 0) return new_n;

	if (r->n[str] > 0){
		free(r->pos[str]);
		free(r->p[str]);
	}

	r->n[1 - str] = r->n[str] = new_n;
	r->pos[str] = pos;
	r->p[str] = p;
	r->err[str] = score2;

	if (type <= 0)
		reverse_cigar(new_cigar);

	free_cigar(r->cigar[str]);
	r->cigar[str] = new_cigar;

	return new_n;
}


/******************************************************************************
 * ************************************************************************** *
 * *                              HASH TABLE                                * *
 * ************************************************************************** * 
 ******************************************************************************/

void init_hash()
{
	KMER_HASH = malloc(sizeof (Kmer_hash));
	KMER_HASH->n = calloc(H + 1, sizeof (int));
	KMER_HASH->bucket = calloc(H + 1, sizeof (Bucket*));
}

inline unsigned int hash_get(unsigned long id, unsigned int *ret, unsigned int add)
{
	unsigned int block, mod;

	mod  = id & (H - 1);
	block = id / H;

	*ret = 0;
	if (KMER_HASH->n[mod] == 0)
		return 0;

	unsigned int start, end, mid, pos;
	start = 0;
	end = KMER_HASH->n[mod];

	while (start != end) {
		mid = (start + end) >> 1;
		if (KMER_HASH->bucket[mod][mid].id == block) {
			start = end = mid;
			break;
		}

		if (KMER_HASH->bucket[mod][mid].id > block)
			end = mid;
		else
			start = mid + 1;
	}

	if (start < KMER_HASH->n[mod] && KMER_HASH->bucket[mod][mid].id == block){
		KMER_HASH->bucket[mod][mid].end += add;
		*ret = 1;
	}

	return start;
}



inline void hash_put(unsigned long id, unsigned int ins, unsigned int order)
{
	unsigned int i, mod, block;

	mod = id & (H - 1);
	block = id / H;

	i = KMER_HASH->n[mod];
	++KMER_HASH->n[mod];

	if (i == 0){
		KMER_HASH->bucket[mod] = calloc(1, sizeof (Bucket));
	} else {
		KMER_HASH->bucket[mod] = realloc(KMER_HASH->bucket[mod],
						    (i + 1) * sizeof (Bucket));
	        memmove(KMER_HASH->bucket[mod] + ins + 1,
					 KMER_HASH->bucket[mod] + ins,
                                     	  	  (i - ins) * sizeof (Bucket));
	}

	KMER_HASH->bucket[mod][ins].id = block;
	KMER_HASH->bucket[mod][ins].end = 1;
}

void destroy_hash()
{
	unsigned int i, k;

	for (i = 0; i < H + 1; i++) 
		if (KMER_HASH->n[i] > 0)
			free(KMER_HASH->bucket[i]);

	free(KMER_HASH->bucket);
	free(KMER_HASH->n);
	free(KMER_HASH->pos);
	free(KMER_HASH->p);
	free(KMER_HASH);
}

void init_refInf()
{
	REF_INF = malloc(sizeof (Ref_inf));
	REF_INF->n = 0;
	REF_INF->name_len = 0;
	REF_INF->seq = malloc(1);
	REF_INF->len = calloc(1, sizeof (int));
	REF_INF->find = calloc(1, sizeof (int));
}

Candidate *init_candidate()
{
	Candidate *r = malloc(sizeof(Candidate));
	r->err[0] = r->err[1] = r->n[0] = r->n[1] = 0;
	r->cigar[0] = init_cigar();
	r->cigar[1] = init_cigar();
	r->pair = 0;

	return r;
}

void destroy_refInf()
{
	free(REF_INF->len);
	free(REF_INF->seq);
	free(REF_INF->find);
	free(REF_INF->count);
	free(REF_INF->eff_len);
	free(REF_INF->ref_name);
	free(REF_INF);
}

void destroy_candidate(Candidate *read)
{
	if (read->n[0] > 0) {
		free(read->pos[0]);
		free(read->p[0]);
	}

	if (read->n[1] > 0) {
		free(read->pos[1]);
		free(read->p[1]);
	}

	free_cigar(read->cigar[0]);
	free_cigar(read->cigar[1]);

	free(read);
}

void destroy_readInf(Read_inf read)
{
	free(read.seq);
	free(read.qual);
	free(read.name);
}

/******************************************************************************
 * ************************************************************************** *
 * *                                INDEX                                   * *
 * ************************************************************************** * 
 ******************************************************************************/

inline unsigned long get_index(unsigned int start, char *seq)
{
	unsigned short i, border;
	unsigned long idx = 0;

	border = start + KMER;
	for (i = start; i < border; ++i)
		idx = (idx << 2) | ascii_table[seq[i]];
	return idx;
}

void get_name(char *file_path)
{
	unsigned int n, total, sum;
	int len;
	unsigned short *name_len;
	gzFile *fa;
	kseq_t *ref;
	char *name;

	// Init
	fa = gzopen(file_path, "r");
	if (!fa) {
		ERROR_PRINT("Can not open file %s\nExit.\n", file_path);
		exit(EXIT_FAILURE);
	}
	ref = kseq_init(fa);
	n = total = sum = 0;
	name = malloc(1);
	name_len = calloc(1, sizeof(short)); 

	// Read sequence
	while ((len = kseq_read(ref)) >= 0){
		++n;
		REF_INF->len = realloc(REF_INF->len, n * sizeof (int));
		REF_INF->find = realloc(REF_INF->find, n * sizeof (int));

		name_len = realloc(name_len, n * sizeof(short));
		name_len[n-1] = ref->name.l + 1;
		name = realloc(name, name_len[n-1] + sum);
		memcpy(name + sum, ref->name.s, name_len[n-1]);
		sum += name_len[n-1];
		if (name_len[n-1] > REF_INF->name_len)
			REF_INF->name_len = name_len[n-1];

		REF_INF->len[n-1] = ref->seq.l;
		REF_INF->seq = realloc(REF_INF->seq, total + REF_INF->len[n-1]);
		memcpy(REF_INF->seq + total, ref->seq.s, REF_INF->len[n-1]);
		total += REF_INF->len[n-1];
		REF_INF->find[n-1] = total;
	}
	REF_INF->n = n;

	REF_INF->ref_name = malloc(REF_INF->n * REF_INF->name_len);
	for (len = sum = 0; len < n; sum += name_len[len], ++len)
		memcpy(REF_INF->ref_name + len * REF_INF->name_len,
					name + sum, name_len[len]);

	free(name);
	free(name_len);
	gzclose(fa);
	kseq_destroy(ref);
}

void construct_hash(unsigned int total_length)
{
	unsigned long idx;
	unsigned int i, k, p, ch, pos, start, count, mod, ret;

	KMER_HASH->n_pos = total_length;
	KMER_HASH->pos = calloc(total_length, sizeof(int));
	KMER_HASH->p = calloc(total_length, sizeof(int));

	for (i = count = 0; i < H +1; ++i){
		for (k = 0; k < KMER_HASH->n[i]; ++k){
			KMER_HASH->bucket[i][k].start = count; 
			count += KMER_HASH->bucket[i][k].end;
			KMER_HASH->bucket[i][k].end =
				 KMER_HASH->bucket[i][k].start;
		}
	}

	for (i = 0; i < REF_INF->n; ++i){
		p = i == 0? 0: REF_INF->find[i-1];
		for (k = idx = 0; k < REF_INF->len[i]; ++k){
			ch = ascii_table[REF_INF->seq[p + k]] & 3;

			idx = POP_HEAD(idx) | ch;

			if (k + 1 >= KMER) {
				pos = hash_get(idx, &ret, 0);
				mod = idx & (H - 1);
				start = KMER_HASH->bucket[mod][pos].end;
				KMER_HASH->pos[start] = i;
				KMER_HASH->p[start] = k;
				++KMER_HASH->bucket[mod][pos].end; 
			}
		}
	}
}

void get_ref_seq(char *file_path)
{
	unsigned long idx;
	unsigned int i, k, p, ch, pos, count, total, ret;

	get_name(file_path);

	for (i = count = total = 0; i < REF_INF->n; ++i){
		p = i == 0? 0: REF_INF->find[i-1];
		for (k = idx = 0; k < REF_INF->len[i]; ++k){
			ch = ascii_table[REF_INF->seq[p + k]] & 3;

			if (k + 1 > KMER)
				idx = (idx << (2 * (BLOCK - KMER + 1)))
                                                   >> (2 * (BLOCK - KMER + 1));

			idx = (idx << 2) | ch;

			if (k + 1 >= KMER) {
				++total;
				pos = hash_get(idx, &ret, 1);
				if (ret == 0){
					hash_put(idx, pos, count);
					++count;
				}
			}
		}
	}
	
	construct_hash(total);
}

void write_to_file(char *out_file)
{
	FILE *out = fopen(out_file, "wb");
	unsigned int i, k;
	int j = -1;

	// Write reference information
	fwrite(&REF_INF->n, sizeof (int), 1, out);
	fwrite(&REF_INF->name_len, sizeof (int), 1, out);
	fwrite(REF_INF->len, sizeof (int), REF_INF->n, out);
	fwrite(REF_INF->find, sizeof (int), REF_INF->n, out);
	fwrite(REF_INF->ref_name, REF_INF->name_len, REF_INF->n, out);
	fwrite(REF_INF->seq, sizeof (char), REF_INF->find[REF_INF->n-1], out);

	destroy_refInf();

	// Write hash table
	for (i = 0; i < H + 1; ++i) {
		if (KMER_HASH->n[i] == 0) continue;
		fwrite(&i, sizeof (int), 1, out);
		fwrite(&KMER_HASH->n[i], sizeof (int), 1, out);
		fwrite(KMER_HASH->bucket[i], sizeof (Bucket),
					 KMER_HASH->n[i], out);
	}
	fwrite(&j, sizeof (int), 1, out);

	fwrite(&KMER_HASH->n_pos, sizeof (int), 1, out);
	fwrite(KMER_HASH->pos, sizeof (int), KMER_HASH->n_pos, out);
	fwrite(KMER_HASH->p, sizeof (int), KMER_HASH->n_pos, out);

	destroy_hash();	
}

void init_fusion()
{
	FUSION = calloc(1, sizeof(Fusion));
	FUSION->n = FUSION->n_read = 0;
	FUSION->detail = calloc(1, sizeof(Fusion_inf));
	FUSION->read = calloc(1, sizeof(char*));
}

void init_class()
{
	CLASS = calloc(1, sizeof(Gclass));
	CLASS->n = 0;
	CLASS->cls = calloc(1, sizeof(Class));
	CLASS->map = calloc(H, sizeof(int));
}

void load_hash(FILE *input)
{
	int i, ret;

	ret = fread(&i, sizeof (int), 1, input);
	while (i >= 0) {
		ret = fread(&KMER_HASH->n[i], sizeof (int), 1, input);
		KMER_HASH->bucket[i] = calloc(KMER_HASH->n[i], sizeof (Bucket));
		ret = fread(KMER_HASH->bucket[i], sizeof (Bucket),
							 KMER_HASH->n[i], input);
		ret = fread(&i, sizeof (int), 1, input);
	}

	ret = fread(&KMER_HASH->n_pos, sizeof (int), 1, input);
	KMER_HASH->pos = calloc(KMER_HASH->n_pos, sizeof(int));
	KMER_HASH->p = calloc(KMER_HASH->n_pos, sizeof(int));
	ret = fread(KMER_HASH->pos, sizeof (int), KMER_HASH->n_pos, input);
	ret = fread(KMER_HASH->p, sizeof (int), KMER_HASH->n_pos, input);
}

void load_index(char *idx_dir, char *genome)
{
	unsigned int l = strlen(idx_dir);
	int ret;

	char file_path[l + 20];
	memcpy(file_path, idx_dir, l);
	memcpy(file_path + l, "/index\0", 7);
	FILE *input = fopen(file_path, "rb");

	if (!input) {
 		fprintf(stdout, "Can not open index file or file does not exist!!!\n");
 		exit(EXIT_FAILURE);
 	}

	// Read reference information
	ret = fread(&REF_INF->n, sizeof (int), 1, input);
	ret = fread(&REF_INF->name_len, sizeof (int), 1, input);
	REF_INF->ref_name = malloc(REF_INF->n * REF_INF->name_len);
	REF_INF->len = calloc(REF_INF->n, sizeof (int));
	REF_INF->find = calloc(REF_INF->n, sizeof (int));
	REF_INF->count = calloc(REF_INF->n, sizeof (int));
	REF_INF->eff_len = calloc(REF_INF->n, sizeof (double));

	ret = fread(REF_INF->len, sizeof (int), REF_INF->n, input);
	ret = fread(REF_INF->find, sizeof (int), REF_INF->n, input);
	ret = fread(REF_INF->ref_name, REF_INF->name_len, REF_INF->n, input);

	REF_INF->seq = malloc(REF_INF->find[REF_INF->n-1]);
	ret = fread(REF_INF->seq, 1, REF_INF->find[REF_INF->n-1], input);

	load_hash(input);
	fclose(input);

	init_class();
	init_fusion();

	FMINDEX = NULL;
	if (genome != NULL)
		FMINDEX = load_FMindex(file_path, l + 6, genome);
}

/******************************************************************************
 * ************************************************************************** *
 * *                              HASH ALIGN                                * *
 * ************************************************************************** * 
 ******************************************************************************/

void get_SWalign(Read_inf read, unsigned int **pos, unsigned int id, 
		unsigned int n, unsigned int *store, unsigned int *max,
					Candidate *r, unsigned int str)
{
	unsigned int i, ref_start, len, add;
	int start, score, skip, err, sub;
	Cigar *cigar, *tmp;

	if (n > 1 && pos[1][id + n -1] - pos[1][id] > read.len)
		return;

	start = pos[1][id] - (pos[2][id] + KMER - 1);

	if (start < 0 || start + read.len > REF_INF->len[pos[0][id]])
		return;

	add = start < ERR? start: ERR;
	ref_start = (pos[0][id] == 0? 0: REF_INF->find[pos[0][id] - 1]) + start;
	err = sub = 0;

	if (pos[2][id] > 0){
		cigar = SW_align(REF_INF->seq, ref_start - add, read.seq,
				    pos[2][id] + add, pos[2][id], &score,
				          pos[2][id] < 2*ERR? -3 : -ERR);
                err += score; 

		if (err > read.len/2 || err > *max){
			free_cigar(cigar);
			return;
		}
		for (i = len = 0; i < cigar->n_cigar; ++i){
			if (cigar->type[i] == INSERT || cigar->type[i] == SCLIP)
				sub += cigar->count[i];
			else if (cigar->type[i] == DELETE)
				sub -= cigar->count[i];
		}
	} else {
		cigar = init_cigar();
	}

	for (i = 0; i < n; ++i){
		start = i;
		while (i + 1 < n && pos[2][id + i] + KMER >= pos[2][id + i + 1])
			++i;
		add_cigar(cigar, MATCH, 
			  pos[2][id + i] - pos[2][id + start] + KMER, NULL);

		start = pos[2][id + i] + KMER;
		len = ((i + 1) < n? pos[2][id + i + 1] : read.len) - start;
		add = (i + 1) < n? 0: ERR;
		tmp = SW_align(REF_INF->seq, ref_start + start,
			    read.seq + start, len + add, len, &score,
				     	        len < 2*ERR? 3: ERR);
		err += score;

		if (err > read.len/2 || err > *max + 1){
			free_cigar(tmp);
			free_cigar(cigar);
			return;
		}

		concat_cigar(cigar, tmp);
		free_cigar(tmp);
	}

	if (*max == read.len/2 || err < *max){
		tmp = r->cigar[str];
		r->cigar[str] = cigar;
		cigar = tmp;
		*store = 0;
		*max = err;
	}

	pos[0][*store] = pos[0][id];
	pos[1][*store] = pos[1][id] - (pos[2][id] + KMER - 1) + sub;

	++*store;
	free_cigar(cigar);
}

void get_position(Read_inf read, Candidate *r, unsigned short str, short space)
{
	unsigned int i, k, pos, add, max, mod, ret;
	unsigned long idx;
	unsigned int *pos1[3], *pos2[3], n1, n2, n;

	if (str == 1)
		reverse_str(read.seq, read.len);

	for (i = n1 = n2 = 0; i + KMER <= read.len; i += ABS(space)){
		idx = get_index(i, read.seq);
		pos = hash_get(idx, &ret, 0);

		if (ret == 0)
			continue;

		mod = idx & (H - 1);
		n1 = KMER_HASH->bucket[mod][pos].end -
				KMER_HASH->bucket[mod][pos].start;
		pos1[0] = KMER_HASH->pos + KMER_HASH->bucket[mod][pos].start;
		pos1[1] = KMER_HASH->p + KMER_HASH->bucket[mod][pos].start;

		if (n2 == 0){
			n2 = n1;
			pos2[0] = calloc(n1, sizeof(int));
			pos2[1] = calloc(n1, sizeof(int));
			pos2[2] = calloc(n1, sizeof(int));
			memcpy(pos2[0], pos1[0], n1 * sizeof(int));
			memcpy(pos2[1], pos1[1], n1 * sizeof(int));
			for (k = 0; k < n1; ++k)
				pos2[2][k] = i;
			n1 = 0;
		}

		if (n1 == 0)
			continue;

		merge_sort(pos2, n2, pos1, n1, i);

		n2 += n1;
	}

	if (n2 == 0)
		goto exit;
	n1 = 0;
	if (r->n[1 - str] > 0 && space < KMER && space > 0){
		for (i = k =  n1 = 0, max = read.len/2; i < n2; ++i){
			if (pos2[0][i] < r->pos[1 - str][k])
				continue;
			while (pos2[0][i] > r->pos[1 - str][k] &&
					 	  k < r->n[1 - str])
				++k;
			if (k == r->n[1 - str])
				break;

			n = 1;
			while (i + 1 < n2 && pos2[0][i] == pos2[0][i + 1] &&
			    pos2[2][i + 1] > pos2[2][i] &&
			    pos2[1][i + 1] > pos2[1][i + 1 - n] &&
			    pos2[1][i + 1] - pos2[1][i + 1 - n] < read.len){
				++i;
				++n;
			}

			get_SWalign(read, pos2, i + 1 - n, n, &n1, &max, r, str);
			for (; n > 0; --n)
				pos2[2][i + 1 - n] = read.len;
		}
	}

	if (n1 == 0) {
		for (i = n1 = 0, max = read.len/2; i < n2; ++i){
			if (pos2[2][i] == read.len)
				continue;

			n = 1;
			while (i + 1 < n2 && pos2[0][i] == pos2[0][i + 1] &&
			    pos2[2][i + 1] > pos2[2][i] &&
			    pos2[1][i + 1] > pos2[1][i + 1 - n] &&
			    pos2[1][i + 1] - pos2[1][i + 1 - n] < read.len &&
			    pos2[2][i + 1] != read.len){
				++i;
				++n;
			}

			if (n < 2 && read.len / KMER > 2)
				continue;

			get_SWalign(read, pos2, i - n + 1, n, &n1, &max, r, str);
		}
	}

	if (n1 == 0 && space < KMER){
		for (i = n1 = 0, max = read.len/2; i < n2; ++i){
			if (pos2[2][i] == read.len)
				continue;

			n = 1;
			while (i + 1 < n2 && pos2[0][i] == pos2[0][i + 1] &&
			    pos2[2][i + 1] > pos2[2][i] &&
			    pos2[1][i + 1] > pos2[1][i + 1 - n] &&
			    pos2[1][i + 1] - pos2[1][i + 1 - n] < read.len &&
			    pos2[2][i + 1] != read.len){
				++i;
				++n;
			}

			if (n >= 2)
				continue;

			get_SWalign(read, pos2, i, n, &n1, &max, r, str);
		}
	}

	if (n1 == 0){
		free(pos2[0]);
		free(pos2[1]);
		free(pos2[2]);
		goto exit;
	}

	if (r->n[str] > 0){
		free(r->pos[str]);
		free(r->p[str]);
	}

	pos2[0] = realloc(pos2[0], n1 * sizeof(int));
	pos2[1] = realloc(pos2[1], n1 * sizeof(int));

	r->n[str] = n1;
	r->p[str] = pos2[1];
	r->pos[str] = pos2[0];
	r->err[str] = max;

	free(pos2[2]);
	if (str == 1)
		reverse_str(read.seq, read.len);
	return;

exit:
	r->n[str] = 0;
	r->err[str] = read.len;
	if (str == 1)
		reverse_str(read.seq, read.len);
	return;
}

unsigned int get_position2(Candidate *r, unsigned short str,
					 Read_inf read, int type)
{
	unsigned int i, k, n, ret;
	Cigar *cigar, *tmp;

	cigar = init_cigar();

	tmp = r->cigar[str];
	r->cigar[str] = cigar;
	n = r->n[str];

	get_position(read, r, str, 2*ERR);
	ret = r->n[str];
	if (ret == 0){
		if (type > 0)
			reverse_str(read.seq, read.len);

		ret = re_mapSW(r, str, read, type);

		if (type > 0)
			reverse_str(read.seq, read.len);
	}

	if (ret > 0){
		free_cigar(tmp);
	} else {
		r->n[str] = n;
		free_cigar(r->cigar[str]);
		r->cigar[str] = tmp;
	}

	return ret;
}

void re_map(Candidate *r, Read_inf read1, Read_inf read2)
{
	unsigned int ret, i;

	reverse_str(read2.seq, read2.len);

	i = 0;
	if (r->err[0] < r->err[1])
		i = 1;

	ret = re_mapSW(r, i, i == 0? read1: read2, i == 0? -1: 1);

	if (ret == 0)
		ret = re_mapSW(r, 1 - i, i == 0? read2: read1, i == 0? 1: -1);

	if (ret > 0)
		r->pair = 2;
	reverse_str(read2.seq, read2.len);
}

void intersect2(Candidate *r, Read_inf read1, Read_inf read2, unsigned int add)
{
	unsigned int i, k, j, m, n, dis;

	i = k = n = 0;
	while (i < r->n[0] && k < r->n[1]) {
		if (r->pos[0][i] < r->pos[1][k]){
			++i;
		} else if (r->pos[0][i] > r->pos[1][k]){
			++k;
		} else {
			dis = ABS(r->p[0][i] - r->p[1][k]);
			if (dis <= add || add == 0) {
				++n;
				r->pos[0][n - 1] = r->pos[0][i];
				r->p[0][n - 1] = r->p[0][i];
				r->pos[1][n - 1] = r->pos[1][k];
				r->p[1][n - 1] = r->p[1][k];

				++i;
				++k;
			} else if (r->p[0][i] < r->p[1][k]){
				++i;
			} else {
				++k;
			}
		}
	}

	if (n > 0){
		r->n[0] = r->n[1] = n;
		r->pair = 2;
		return;
	} else if (add > 0){
		return intersect2(r, read1, read2, 0);
	} else {
		return re_map(r, read1, read2);
	}
}

void add_class(Candidate *r, unsigned int pair)
{
	pthread_mutex_lock(&LOCK);
	++PAIRED;
	if (pair == 1)
		MEAN_FRAG += ABS(r->p[0][0] - r->p[1][0]);

	if (r->n[0] == 1) {
		++REF_INF->count[r->pos[0][0]];
	} else if (r->err[0] + r->err[1] < ERR){
		unsigned int id, i, tmp;

		id = XXH32((char*) r->pos[0], r->n[0] * sizeof(int), 0);
		tmp = id;
		i = CLASS->map[id & (H - 1)];
		while (CLASS->map[tmp & (H - 1)] > 0 &&
					   CLASS->cls[i - 1].id != id){
			++tmp;
			i = CLASS->map[tmp & (H - 1)];
		}

		if (i == 0){
			i = CLASS->n;
			++CLASS->n;
			CLASS->map[tmp & (H - 1)] = CLASS->n;
			CLASS->cls = realloc(CLASS->cls, CLASS->n*sizeof(Class));
			CLASS->cls[i].count = 1;
			CLASS->cls[i].ref = r->pos[0];
			CLASS->cls[i].n = r->n[0];
			CLASS->cls[i].id = id;

			free(r->p[0]);
			r->n[0] = 0;
		} else {
			++CLASS->cls[i - 1].count;
		}
	}
	pthread_mutex_unlock(&LOCK);
}

void concate_candidate(Candidate *r1, Candidate *r2)
{
	r1->pos[0] = realloc(r1->pos[0], (r1->n[0] + r2->n[0]) * sizeof (int));
	r1->p[0] = realloc(r1->p[0], (r1->n[0] + r2->n[0]) * sizeof (int));
	r1->pos[1] = realloc(r1->pos[1], (r1->n[1] + r2->n[1]) * sizeof (int));
	r1->p[1] = realloc(r1->p[1], (r1->n[1] + r2->n[1]) * sizeof (int));

	memcpy(r1->pos[0] + r1->n[0], r2->pos[0], r2->n[0] * sizeof (int));
	memcpy(r1->p[0] + r1->n[0], r2->p[0], r2->n[0] * sizeof (int));
	memcpy(r1->pos[1] + r1->n[1], r2->pos[1], r2->n[1] * sizeof (int));
	memcpy(r1->p[1] + r1->n[1], r2->p[1], r2->n[1] * sizeof (int));

	r1->n[0] += r2->n[0];
	r1->n[1] += r2->n[1];
}

void move_candidate(Candidate *r, unsigned int p1, unsigned int p2)
{
	Cigar *tmp;

	if (r->n[p1] > 0){
		free(r->pos[p1]);
		free(r->p[p1]);
	}

	r->n[p1] = r->n[p2];
	r->p[p1] = r->p[p2];
	r->pos[p1] = r->pos[p2];

	tmp = r->cigar[p1];
	r->cigar[p1] = r->cigar[p2];
	r->cigar[p2] = tmp;

	r->n[p2] = 0;
}

void swap_candidate(Candidate *r1, unsigned int p1,
			 Candidate *r2, unsigned int p2)
{
	Cigar *tmp_cigar;
	unsigned int tmp_n, *tmp_p;

	tmp_n = r1->n[p1];
	r1->n[p1] = r2->n[p2];
	r2->n[p2] = tmp_n;

	tmp_p = r1->pos[p1];
	r1->pos[p1] = r2->pos[p2];
	r2->pos[p2] = tmp_p;

	tmp_p = r1->p[p1];
	r1->p[p1] = r2->p[p2];
	r2->p[p2] = tmp_p;

	tmp_cigar = r1->cigar[p1];
	r1->cigar[p1] = r2->cigar[p2];
	r2->cigar[p2] = tmp_cigar;
}

void check_in_transcriptome(Read_inf r1, Read_inf r2, Candidate *r[])
{
	unsigned int ret;

	get_position(r1, r[0], 0, KMER);
	get_position(r2, r[0], 1, KMER);
	get_position(r1, r[1], 1, KMER);
	get_position(r2, r[1], 0, KMER);

	if (r[0]->n[0] == 0 && r[0]->n[1] == 0 &&
		 r[1]->err[0] + r[1]->err[1] > ERR){
		get_position(r1, r[0], 0, 2*ERR);
		get_position(r2, r[0], 1, 2*ERR);
	}
	if (r[1]->n[0] == 0 && r[1]->n[1] == 0 &&
		 r[0]->err[0] + r[0]->err[1] > ERR){
		get_position(r1, r[1], 1, 2*ERR);
		get_position(r2, r[1], 0, 2*ERR);
	}

	if (r[0]->n[0] == 0 && r[0]->n[1] > 0 && r[0]->err[1] < ERR && 
				    r[1]->err[0] + r[1]->err[1] > ERR)
		ret = get_position2(r[0], 0, r1, -1);
	if (r[0]->n[1] == 0 && r[0]->n[0] > 0 && r[0]->err[0] < ERR && 
				    r[1]->err[0] + r[1]->err[1] > ERR)
		ret = get_position2(r[0], 1, r2, 1);
	if (r[1]->n[0] == 0 && r[1]->n[1] > 0 && r[1]->err[1] < ERR && 
				    r[0]->err[0] + r[0]->err[1] > ERR)
		ret = get_position2(r[1], 0, r2, -1);
	if (r[1]->n[1] == 0 && r[1]->n[0] > 0 && r[1]->err[0] < ERR && 
				    r[0]->err[0] + r[0]->err[1] > ERR)
		ret = get_position2(r[1], 1, r1, 1);
}

void check_in_genome(Read_inf r1, Read_inf r2, Candidate *r[])
{
	if (FMINDEX == NULL)
		return;

	unsigned int n_map = 0;
	if (r[0]->n[0] > 0 && r[0]->err[0] < 2*ERR)
		++n_map;
	if (r[0]->n[1] > 0 && r[0]->err[1] < 2*ERR)
		++n_map;
	if (r[1]->n[0] > 0 && r[1]->err[0] < 2*ERR)
		++n_map;
	if (r[1]->n[1] > 0 && r[1]->err[1] < 2*ERR)
		++n_map;

	if (n_map < 2)
		genome_map(r1, r2, r);
}

void align_read(Read_inf r1, Read_inf r2, char *stream, unsigned int *slen)
{
	if (r1.name == NULL || r2.name == NULL ||
		 r1.seq == NULL || r2.seq == NULL)
		return;

	// Get position of each read
	Candidate *r[2];
	r[0] = init_candidate();
	r[1] = init_candidate();

	check_in_transcriptome(r1, r2, r);

	if (r[0]->n[0] > 0 && r[0]->n[1] > 0)
		intersect2(r[0], r1, r2, MAX_FRAG);
	if (r[1]->n[0] > 0 && r[1]->n[1] > 0)
		intersect2(r[1], r2, r1, MAX_FRAG);

	if (r[0]->err[0] + r[0]->err[1] < r[1]->err[0] + r[1]->err[1])
		r[1]->pair = -1;
	else if (r[0]->err[0] + r[0]->err[1] > r[1]->err[0] + r[1]->err[1])
		r[0]->pair = -1;

	check_in_genome(r1, r2, r);

	if ((r[0]->pair&1) == 0 && r[0]->n[0] > 0 && r[0]->n[1] > 0)
		store_fusion(r[0], r1, r2);
	if ((r[1]->pair&1) == 0 && r[1]->n[0] > 0 && r[1]->n[1] > 0)
		store_fusion(r[1], r2, r1);

	if (r[0]->pair > 1 && r[1]->pair > 1) {
		reverse_str(r2.seq, r2.len);
		reverse_qual(r2.qual, r2.len);
		bam_write_pair(r1, r2, r[0], 1, 1, 0, 1, stream, slen);
		reverse_str(r2.seq, r2.len);
		reverse_qual(r2.qual, r2.len);
		reverse_str(r1.seq, r1.len);
		reverse_qual(r1.qual, r1.len);
		bam_write_pair(r2, r1, r[1], 1, 0, 0, 1, stream, slen);
		if ((r[0]->pair & 1) + (r[1]->pair & 1) == 0){
			concate_candidate(r[0], r[1]);
			add_class(r[0], 1);
                }
	} else if (r[0]->pair > 1) {
		reverse_str(r2.seq, r2.len);
		reverse_qual(r2.qual, r2.len);
		bam_write_pair(r1, r2, r[0], 1, 1, 0, 1, stream, slen);
		if ((r[0]->pair & 1) == 0)
			add_class(r[0], 1);
	} else if (r[1]->pair > 1){
		reverse_str(r1.seq, r1.len);
		reverse_qual(r1.qual, r1.len);
		bam_write_pair(r2, r1, r[1], 1, 0, 0, 1, stream, slen);
		if ((r[1]->pair & 1) == 0)
			add_class(r[1], 1);
	} else if (r[0]->n[0] > 0 && r[0]->n[1] > 0 && 
		r[0]->err[0] + r[0]->err[1] < r[1]->err[0] + r[1]->err[1]){
		reverse_str(r2.seq, r2.len);
		reverse_qual(r2.qual, r2.len);
		bam_write_pair(r1, r2, r[0], 0, 1, 0, 1, stream, slen);
	} else if (r[1]->n[0] > 0 && r[1]->n[1] > 0){
		reverse_str(r1.seq, r1.len);
		reverse_qual(r1.qual, r1.len);
		bam_write_pair(r2, r1, r[1], 0, 0, 0, 1, stream, slen);
	} else if (r[0]->n[0] > 0 && r[1]->n[0] > 0 &&
                        (r[0]->pair & 1) == (r[1]->pair & 1)){
		swap_candidate(r[0], 1, r[1], 0);
		r[0]->pair = -10;
		if ((r[0]->pair & 1) == 0)
			store_fusion(r[0], r1, r2);
		bam_write_pair(r1, r2, r[0], 0, 1, 0, 0, stream, slen);
	} else if (r[0]->n[1] > 0 && r[1]->n[1] > 0 &&
                        (r[0]->pair & 1) == (r[1]->pair & 1)){
		reverse_str(r1.seq, r1.len);
		reverse_qual(r1.qual, r1.len);
		reverse_str(r2.seq, r2.len);
		reverse_qual(r2.qual, r2.len);
		swap_candidate(r[0], 0, r[1], 1);
		r[0]->pair = -20;
		if ((r[0]->pair & 1) == 0)
			store_fusion(r[0], r1, r2);
		bam_write_pair(r1, r2, r[0], 0, 0, 1, 1, stream, slen);
	} else if (r[0]->n[0] > 0){
		bam_write_pair(r1, r2, r[0], 0, 1, 0, 0, stream, slen);
	} else if (r[0]->n[1] > 0){
		reverse_str(r2.seq, r2.len);
		reverse_qual(r2.qual, r2.len);
		bam_write_pair(r1, r2, r[0], 0, 1, 0, 1, stream, slen);
	} else if (r[1]->n[0] > 0){
		bam_write_pair(r2, r1, r[1], 0, 0, 0, 0, stream, slen);
	} else {
		reverse_str(r1.seq, r1.len);
		reverse_qual(r1.qual, r1.len);
		bam_write_pair(r2, r1, r[1], 0, 0, 0, 1, stream, slen);
	}

	if (r[0]->n[0] > 0 || r[0]->n[1] > 0 || r[1]->n[0] > 0 ||r[1]->n[1] > 0){
		pthread_mutex_lock(&LOCK);
		++MAPPED;
		pthread_mutex_unlock(&LOCK);
	}

	destroy_candidate(r[0]);
	destroy_candidate(r[1]);
}

void align_read2(Read_inf read, char *stream, unsigned int *slen)
{
	if (read.name == NULL || read.seq == NULL)
		return;

	Candidate *r;
	r = init_candidate();

	get_position(read, r, 0, 2*ERR);
	get_position(read, r, 1, 2*ERR);

	if (r->n[0] == 0 && r->n[1] == 0) {
		bam_write_single(read, r, 0, 0, stream, slen);
		destroy_candidate(r);
		return;
	}

	if (r->n[0] > 0){
		bam_write_single(read, r, 0, 0, stream, slen);
	} else if (r->n[1] > 0){
		reverse_str(read.seq, read.len);
		reverse_qual(read.qual, read.len);
		bam_write_single(read, r, 1, 1, stream, slen);
	}

	if (r->err[0] > r->err[1]) {
		move_candidate(r, 0, 1);
	} else if (r->err[0] == r->err[1]) {
		r->pos[0] = realloc(r->pos[0],
					   (r->n[0] + r->n[1]) * sizeof (int));
		r->p[0] = realloc(r->p[0], (r->n[0] + r->n[1]) * sizeof (int));

		memcpy(r->pos[0] + r->n[0], r->pos[1], r->n[1] * sizeof (int));
		memcpy(r->p[0] + r->n[0], r->p[1], r->n[1] * sizeof (int));
		r->n[0] += r->n[1];
	}

	add_class(r, 0);
	destroy_candidate(r);
}

unsigned int get_reads(Read_inf *read1, Read_inf *read2)
{
	pthread_mutex_lock(&LOCK);
	int l1, l2;
	unsigned int n;
	
	n = 0;
	while (n < GROUP){
		l1 = kseq_read(R1);
		l2 = kseq_read(R2);
		while (l1 == -2 || l2 == -2){
			++N_READ;
			l1 = kseq_read(R1);
			l2 = kseq_read(R2);
		}
		if (l1 == -1 || l2 == -1)
			break;

		++N_READ;
		MEAN_LEN += l1 + l2;

		read1[n].seq = R1->seq.s;
		read1[n].len = R1->seq.l;
		read1[n].name = R1->name.s;
		read1[n].qual = R1->qual.s;
		read2[n].seq = R2->seq.s;
		read2[n].len = R2->seq.l;
		read2[n].qual = R2->qual.s;
		read2[n].name = R2->name.s;

		R1->seq.s = R2->seq.s = 0;
		R1->qual.s = R2->qual.s = R1->name.s = R2->name.s = 0;
		R1->qual.m = R2->qual.m = R1->name.m = R2->name.m = 0;
		++n;
	}

	printf("\t\rNumber of processed pairs\t: %u", N_READ);
	fflush(stdout);
	if ((l1 == -1 || l2 == -1) && l1 != l2)
		fprintf(stderr, "\n[WARNING] Number of reads in two files are not equal.\n");
	pthread_mutex_unlock(&LOCK);
	return n;
}

unsigned int get_reads_single(Read_inf *read)
{
	pthread_mutex_lock(&LOCK);
	int l;
	unsigned int n;
	
	n = 0;
	while (n < GROUP){
		l = kseq_read(R1);
		while (l == -2)
			l = kseq_read(R1);
		if (l == -1)
			break;

		++N_READ;
		MEAN_LEN += l;

		read[n].seq = R1->seq.s;
		read[n].len = R1->seq.l;
		read[n].name = R1->name.s;
		read[n].qual = R1->qual.s;

		R1->seq.s = R1->qual.s = R1->name.s = 0;
		R1->qual.m = R1->name.m = 0;
		++n;
	}

	printf("\t\rNumber of processed reads\t: %u", N_READ);
	fflush(stdout);
	pthread_mutex_unlock(&LOCK);
	return n;
}

void *align_pool(void *data)
{
	unsigned int i, n;
	Read_inf read1[GROUP];
	Read_inf read2[GROUP];
	char *stream = malloc(MAX_COMPRESS);
	unsigned int slen = 0;

	do {
		n = get_reads(read1, read2);
		for (i = 0; i < n; ++i){
			align_read(read1[i], read2[i], stream, &slen);
			destroy_readInf(read1[i]);
			destroy_readInf(read2[i]);
		}
	} while (n > 0);

	if (slen > 0 && WRITE_BAM == 0)
		bam_write(BAM, stream, slen);

	free(stream);
	pthread_exit(NULL);
	return NULL;
}

void *align_pool2(void *data)
{
	unsigned int i, n;
	Read_inf read[GROUP];
	char *stream = malloc(MAX_COMPRESS);
	unsigned int slen = 0;

	do {
		n = get_reads_single(read);
		for (i = 0; i < n; ++i){
			align_read2(read[i], stream, &slen);
			destroy_readInf(read[i]);
		}
	} while (n > 0);

	if (slen > 0 && WRITE_BAM == 0)
		bam_write(BAM, stream, slen);

	free(stream);
	pthread_exit(NULL);
	return NULL;
}

void get_alignment_pair(char *fq1, char *fq2)
{
	unsigned int i, k, rc;
	gzFile *input1, *input2;
	pthread_t thr[NTHREAD];
	pthread_attr_t attr;

	input1 = gzopen(fq1, "r");
	if (!input1) {
		ERROR_PRINT("Can not open file %s\nExit.\n", fq1);
		exit(EXIT_FAILURE);
	}

	input2 = gzopen(fq2, "r");
	if (!input2) {
		ERROR_PRINT("Can not open file %s\nExit.\n", fq2);
		exit(EXIT_FAILURE);
	}

	// Initialize
	R1 = kseq_init(input1);
	R2 = kseq_init(input2);
	memset(&thr, '\0', NTHREAD);
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (i = 0; i < NTHREAD; ++i) 
		rc = pthread_create(&thr[i], &attr, align_pool, NULL);
	for (; i > 0; --i)
		rc = pthread_join(thr[i - 1], NULL);
	pthread_attr_destroy(&attr);

	gzclose(input1);
	gzclose(input2);
	kseq_destroy(R1);
	kseq_destroy(R2);

	MEAN_LEN /= 2 * N_READ;
	MEAN_FRAG = MEAN_FRAG / PAIRED + MEAN_LEN;
	for (i = 0; i < REF_INF->n; ++i) {
		if (REF_INF->len[i] < MEAN_FRAG)
			REF_INF->eff_len[i] = 0;
		else
			REF_INF->eff_len[i] = REF_INF->len[i] - MEAN_FRAG + 1;
	}

	printf("\rNumber of aligned pairs\t: %u/%u\n", MAPPED, N_READ);
	printf("Mean read length\t: %f\n", MEAN_LEN);
	printf("Mean fragment length\t: %f\n", MEAN_FRAG);

	fprintf(SUMMARY, "\nMAPPING RESULT:\n");
	fprintf(SUMMARY, "\tTotal number of read pairs\t\t\t: %u\n", N_READ);
	fprintf(SUMMARY, "\t\t- Number of proper mapped pairs\t: %u (%f)\n",
					  PAIRED, (float) PAIRED/N_READ);
	fprintf(SUMMARY, "\t\t- Number of unproper mapped pairs\t: %u (%f)\n",
		      MAPPED -  PAIRED, (float)(MAPPED -  PAIRED)/N_READ);
	fprintf(SUMMARY, "\t\t- Number of unmapped pairs\t\t: %u (%f)\n",
			N_READ - MAPPED, (float) (N_READ - MAPPED)/N_READ);
	fprintf(SUMMARY, "\tMean read length\t\t\t: %f\n", MEAN_LEN);
	fprintf(SUMMARY, "\tMean fragment length\t\t\t: %f\n", MEAN_FRAG);
}

void get_alignment(char *fq)
{
	unsigned int i, k, rc;
	kseq_t *r;
	gzFile *input;
	pthread_t thr[NTHREAD];
	pthread_attr_t attr;

	input = gzopen(fq, "r");
	if (!input) {
		ERROR_PRINT("Can not open file %s\nExit.\n", fq);
		exit(EXIT_FAILURE);
	}

	// Initialize
	R1 = kseq_init(input);
	memset(&thr, '\0', NTHREAD);
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (i = 0; i < NTHREAD; ++i) 
		rc = pthread_create(&thr[i], &attr, align_pool2, NULL);
	for (; i > 0; --i)
		rc = pthread_join(thr[i - 1], NULL);
	pthread_attr_destroy(&attr);

	gzclose(input);
	kseq_destroy(R1);

	MEAN_LEN /= N_READ;
	MEAN_FRAG = MEAN_LEN;
	for (i = 0; i < REF_INF->n; ++i) {
		if (REF_INF->len[i] < MEAN_FRAG)
			REF_INF->eff_len[i] = 0;
		else
			REF_INF->eff_len[i] = REF_INF->len[i] - MEAN_FRAG + 1;
	}
	printf("\rNumber of aligned reads\t: %u/%u\n", PAIRED, N_READ);
	printf("Mean read lenght\t: %f\n", MEAN_LEN);

	fprintf(SUMMARY, "\nMAPPING RESULT:\n");
	fprintf(SUMMARY, "\tTotal number of reads\t\t\t: %u\n", N_READ);
	fprintf(SUMMARY, "\t\t- Number of mapped reads\t: %u (%f)\n",
					  PAIRED, (float) PAIRED/N_READ);
	fprintf(SUMMARY, "\t\t- Number of unmapped reads\t: %u (%f)\n",
			N_READ - PAIRED, (float) (N_READ - PAIRED)/N_READ);
	fprintf(SUMMARY, "\tMean read length\t\t\t: %f\n", MEAN_LEN);
}
