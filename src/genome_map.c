#include "hash_align.h"

void g_add_sort(unsigned int **arr, unsigned int value, unsigned int idx,
				        unsigned int *i, unsigned int *k)
{
	arr[0][*i] = value;
	arr[1][*i] = idx;
	++*i;
	++*k;
}

void g_merge_sort(unsigned int **arr1,  unsigned int n1,
		 unsigned int **arr2, unsigned int n2, unsigned int idx)
{
	unsigned int i, k, n, dis, *new_arr[2];

	new_arr[0] = calloc(n1 + n2, sizeof(int));
	new_arr[1] = calloc(n1 + n2, sizeof(int));

	for (i = k = n = 0; i < n1 && k < n2;){
		if (arr1[0][i] <= arr2[0][k] - 1){
			k += arr1[0][i] == arr2[0][k] - 1? 1: 0;
			g_add_sort(new_arr, arr1[0][i], arr1[1][i], &n, &i);
		} else{
			g_add_sort(new_arr, arr2[0][k] - 1, idx, &n, &k);
		}
	}

	while (i < n1)
		g_add_sort(new_arr, arr1[0][i], arr1[1][i], &n, &i);

	while (k < n2)
		g_add_sort(new_arr, arr2[0][k] - 1, idx, &n, &k);

	free(arr1[0]);
	free(arr1[1]);

	arr1[0] = new_arr[0];
	arr1[1] = new_arr[1];
}

void g_get_SWalign(Read_inf read, unsigned int **pos, unsigned int id, 
		unsigned int n, unsigned int *store, unsigned int *max,
					Candidate *r, unsigned int str)
{
	unsigned int i, len, add, p;
	int start, score, skip, err, sub, ref_start;
	Cigar *cigar, *tmp;

	if (n > 1 && pos[0][id + n -1] - pos[0][id] > read.len)
		return;

	p = pos[1][id] >> HBLOCK;
	ref_start = pos[0][id] - p;

	if (ref_start < 0 || ref_start + read.len > FMINDEX->n)
		return;

	add = ref_start < ERR? ref_start: ERR;
	err = sub = 0;

	if (p > 0){
		cigar = SW_align(FMINDEX->seq, ref_start - add, read.seq,
			p + add, p, &score, p < 2*ERR? -3 : -ERR, &skip, read);
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
		while (i + 1 < n && (pos[1][id + i] >> HBLOCK) +
		   (pos[1][id + i] & 0xffff) >= (pos[1][id + i + 1] >> HBLOCK))
			++i;
		add_cigar(cigar, MATCH,  (pos[1][id + i] >> HBLOCK) -
					   (pos[1][id + start] >> HBLOCK) + 
						(pos[1][id + i] & 0xffff), NULL);

		start = (pos[1][id + i] >> HBLOCK) + (pos[1][id + i] & 0xffff);
		if (start >= read.len)
			break;

		len = ((i + 1) < n? pos[1][id + i + 1] >> HBLOCK :
							    read.len) - start;
		add = (i + 1) < n? 0: ERR;
		tmp = SW_align(FMINDEX->seq, ref_start + start,
			        read.seq + start, len, len, &score,
			     	   len < 2*ERR? 3: ERR, &skip, read);
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

	pos[0][*store] = pos[0][id] - (pos[1][id] >> HBLOCK) + sub;

	++*store;
	free_cigar(cigar);
}

void g_intersect(Candidate *r, Read_inf read1, Read_inf read2, unsigned int add)
{
	unsigned int i, k, j, m, n, dis;

	i = k = n = 0;
	while (i < r->n[0] && k < r->n[1]) {
		dis = ABS(r->pos[0][i] - r->pos[1][k]);
		if (dis > add && r->pair > 1)
			dis -= r->pair - 1;
		if (dis <= add) {
			++n;
			r->pos[0][n - 1] = r->pos[0][i];
			r->pos[1][n - 1] = r->pos[1][k];

			++i;
			++k;
		} else if (r->pos[0][i] < r->pos[1][k]){
			++i;
		} else {
			++k;
		}
	}

	if (n > 0){
		r->n[0] = r->n[1] = n;
		r->pair = 3;
	} else {
		r->pair = 1;
	}
}

unsigned int g_intersect2(unsigned int *arr1, unsigned int *n1,
			unsigned int *arr2, unsigned int n2, int str)
{
	unsigned int i, n, min_dis, chr1, chr2, p1, p2, l[2];
	int dis;

	i = str > 0? 0: 1;
	l[0] = l[1] = n = 0;
	min_dis = MAX_SPLIT;
	while (l[i] < *n1 && l[1 - i] < n2) {
		dis = ((arr2[l[1 - i]] - 1) - arr1[l[i]]) * str;
		if (dis < 0) {
			++l[1];
		} else if (dis > MAX_SPLIT){
			++l[0];
		} else {
			chr_coord(arr1[l[i]], &chr1, &p1);
			chr_coord(arr2[l[1 - i]], &chr2, &p2);
			if (chr1 != chr2){
				++l[0];
				continue;
			}
			if (dis < min_dis){
				min_dis  = dis;
				n = 1;
				arr1[0] = str > 0? arr1[l[i]]:
					   arr2[l[1 - i]] - 1;		
			} else if (dis == min_dis){
				arr1[n++] = str > 0? arr1[l[i]]:
					     arr2[l[1 - i]] - 1;
			}
			++l[0];
			++l[1];
		}
	}

	if (n > 0){
		*n1 = n;
		return min_dis;
	}
	return 0;
}

void insert_padding(Cigar *cigar, unsigned int len, unsigned int pos)
{
	unsigned int n = cigar->n_cigar;
	cigar->type = realloc(cigar->type, n + 1);
	cigar->count = realloc(cigar->count, (n + 1)*sizeof(short));

	memmove(cigar->type + pos + 1, cigar->type + pos, n - pos);
	memmove(cigar->count + pos + 1, cigar->count + pos,
					  (n - pos)*sizeof(short));

	cigar->type[pos] = SKIP;
	cigar->count[pos] = len;
	++cigar->n_cigar;
}

unsigned int check_split(Read_inf read, Candidate *r, unsigned short str)
{
	unsigned int n = r->cigar[str]->n_cigar - 1;
	if (r->cigar[str]->type[0] != SCLIP &&
			r->cigar[str]->type[n] != SCLIP)
		return 1;

	unsigned int p, l, len, start, end, ret;
	unsigned long range;

	if (r->cigar[str]->type[0] == SCLIP){
		p = 0;
		len = r->cigar[str]->count[0];
	} else {
		len = r->cigar[str]->count[n];
		p = read.len - len;
	}

	if (len < KMER)
		return 1;

	l = query(read.seq + p, len, FMINDEX, &range);
	end = (unsigned int) range;
        start = (unsigned int) (range >> BLOCK);

	if (l != len || end - start + 1 > MAX_FIND)
		return 1;
	
	ret = g_intersect2(r->pos[str], &r->n[str], FMINDEX->sa + start,
					end - start + 1, p == 0? -1: 1);
	if (ret > 0){
		p = p == 0? 0: n;
		r->cigar[str]->type[p] = MATCH;
		r->cigar[str]->count[p] = len;
		insert_padding(r->cigar[str], ret, p == 0? 1: n);
	}
	return ret + 1;
}

void g_get_position(Read_inf read, Candidate *r, unsigned short str, short space)
{
	unsigned int i, k, len, add, max, start, end, pidx;
	unsigned long idx, pos;
	unsigned int *pos1[2], *pos2[2], n1, n2, n;

	if (str == 1)
		reverse_str(read.seq, read.len);

	for (i = n1 = n2 = len = 0; i + KMER <= read.len; i += space){
		len = query(read.seq, read.len - i, FMINDEX, &pos);

		end = (unsigned int) pos;
                start = (unsigned int) (pos >> BLOCK);

		n1 = end - start + 1;
		if (len < space || n1 > MAX_FIND)
			continue;

		pos1[0] = FMINDEX->sa + start;
		pidx = (read.len - i - len) << HBLOCK | len;

		if (n2 == 0){
			n2 = n1;
			pos2[0] = calloc(n1, sizeof(int));
			pos2[1] = calloc(n1, sizeof(int));
			memcpy(pos2[0], pos1[0], n1 * sizeof(int));
			for (k = 0; k < n1; ++k){
				--pos2[0][k];
				pos2[1][k] = pidx;
			}
			n1 = 0;
		}
		
                if (len == read.len)
                        break;

		if (n1 == 0)
			continue;

		g_merge_sort(pos2, n2, pos1, n1, pidx);

		n2 += n1;
	}

	n1 = 0;
	if (r->n[1 - str] > 0 && space < KMER){
		for (i = k = 0, max = read.len/2; i < n2; ++i){
			if (pos2[0][i] + 2*MAX_FRAG < r->pos[1 - str][k])
				continue;
			while (pos2[0][i] > r->pos[1 - str][k] + 2*MAX_FRAG &&
						 	  k < r->n[1 - str])
				++k;

			if (k == r->n[1 - str])
				break;

			n = 1;
			while (i + 1 < n2 && (pos2[1][i + 1] >> HBLOCK) >
						(pos2[1][i] >> HBLOCK) &&
                    	  pos2[0][i + 1] - pos2[0][i + 1 - n] < read.len){
				++i;
				++n;
			}

			g_get_SWalign(read, pos2, i + 1 - n, n, &n1,
						      &max, r, str);
		}
	} else if (len < read.len && n2 > 0){
                for (i = 0, max = read.len/2; i < n2; ++i){
                        n = 1;
			while (i + 1 < n2 && (pos2[1][i + 1] >> HBLOCK) >
						(pos2[1][i] >> HBLOCK) &&
                         pos2[0][i + 1] - pos2[0][i + 1 - n] < read.len){
                                ++i;
                                ++n;
                        }

			if (n < 2 && (pos2[1][i] & 0xffff) <= KMER)
                                continue;

                        g_get_SWalign(read, pos2, i - n + 1, n, &n1,
						      &max, r, str);
                }
        } else if (n2 > 0) {
                n1 = n2;
                r->cigar[str] = init_cigar();
                add_cigar(r->cigar[str], MATCH, read.len, NULL);
                max = 0;
        }

	if (n1 > 0){
                pos2[0] = realloc(pos2[0], n1 * sizeof(int));

                r->n[str] = n1;
                r->pos[str] = pos2[0];
                r->p[str] = calloc(1, sizeof(int));
                r->err[str] = max;
		r->pair = MAX2(r->pair, check_split(read, r, str));
	} else {
                r->err[str] = read.len;
		if (n2 > 0)
                	free(pos2[0]);
        }

        if (n2 > 0)
        	free(pos2[1]);

	if (str == 1)
		reverse_str(read.seq, read.len);
	return;
}

void mate_pair(Read_inf r1, Read_inf r2, Candidate *r)
{
	if (r->n[0] == 0 && r->n[1] > 0 && r->err[1] < ERR)
		g_get_position(r1, r, 0, KMER/2);
	if (r->n[1] == 0 && r->n[0] > 0 && r->err[0] < ERR)
		g_get_position(r2, r, 1, KMER/2);

	if (r->n[0] > 0 && r->n[1] > 0)
		g_intersect(r, r1, r2, MAX_FRAG);
}

void genome_map(Read_inf r1, Read_inf r2, Candidate *r[])
{
	Candidate *r_tmp[2];
	r_tmp[0] = init_candidate();
	r_tmp[1] = init_candidate();

    	g_get_position(r1, r_tmp[0], 0, KMER);
	g_get_position(r2, r_tmp[0], 1, KMER);
	g_get_position(r1, r_tmp[1], 1, KMER);
	g_get_position(r2, r_tmp[1], 0, KMER);

	if (r_tmp[0]->pair == 0 && r_tmp[1]->pair == 0){
		destroy_candidate(r_tmp[0]);
		destroy_candidate(r_tmp[1]);
		return;
	}

	if (r_tmp[0]->pair > 0)
		mate_pair(r1, r2, r_tmp[0]);

	if (r_tmp[1]->pair > 0)
		mate_pair(r2, r1, r_tmp[1]);

	if (r_tmp[0]->pair != 3 && r_tmp[1]->pair != 3){
		destroy_candidate(r_tmp[0]);
		destroy_candidate(r_tmp[1]);
		return;
	}

	destroy_candidate(r[0]);
	destroy_candidate(r[1]);
	r[0] = r_tmp[0];
	r[1] = r_tmp[1];
}