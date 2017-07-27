#include "hash_align.h"

unsigned int get_pos(unsigned int trans1, unsigned int trans2)
{
	unsigned int i, n, t1, t2;
	for (i = 0; i < FUSION->n; ++i){
		t1 = FUSION->detail[i].trans[0];
		t2 = FUSION->detail[i].trans[1];
		if (check_overlapGene(trans1, t1) == 1 ||
		    check_overlapGene(trans1, t2) == 1 ||
		    check_overlapGene(trans2, t1) == 1 ||
		    check_overlapGene(trans2, t2) == 1)
		    return i;
	}

	n = FUSION->n;
	++FUSION->n;

	FUSION->detail = realloc(FUSION->detail, (n + 1)*sizeof(Fusion_inf));
	FUSION->detail[n].trans[0] = trans1;
	FUSION->detail[n].trans[1] = trans2;
	FUSION->detail[n].n = 0;
	FUSION->detail[n].read = calloc(1, sizeof(int));

	return n;
}

unsigned int add_read(Read_inf read1, Read_inf read2)
{
	unsigned int n = FUSION->n_read;
	FUSION->n_read += 2;
	FUSION->read = realloc(FUSION->read, FUSION->n_read*sizeof(char*));

	FUSION->read[n] = malloc(read1.len + 1);
	memcpy(FUSION->read[n], read1.seq, read1.len + 1);

	FUSION->read[n + 1] = malloc(read2.len + 1);
	memcpy(FUSION->read[n + 1], read2.seq, read2.len + 1);

	return n;
}

void add_pair(unsigned int pos, unsigned int read_pos)
{
	unsigned int n = FUSION->detail[pos].n;
	if (n > 0 && FUSION->detail[pos].read[n - 1] == read_pos)
		return;
	++FUSION->detail[pos].n;
	FUSION->detail[pos].read = realloc(FUSION->detail[pos].read,
						(n + 1)*sizeof(int));
	FUSION->detail[pos].read[n] = read_pos;
}

void fusion_add(Candidate *r, Read_inf read1, Read_inf read2)
{
	unsigned int t1, t2;
	unsigned int i, k, p, pos;
	unsigned long id;

	pthread_mutex_lock(&LOCK);
	p = add_read(read1, read2);
	pthread_mutex_unlock(&LOCK);

	for (i = 0; i < r->n[0]; ++i){
		if (i > 0 && check_overlapGene(r->pos[0][i],
				r->pos[0][i - 1]) == 1)
			continue;
		for (k = 0; k < r->n[1]; ++k){
			if (check_overlapGene(r->pos[0][i], r->pos[1][k]) == 1 ||
					(k > 0 && check_overlapGene(r->pos[1][k],
							r->pos[1][k - 1]) == 1))
				continue;

			t1 = r->pos[0][i];
			t2 = r->pos[1][k];

			pthread_mutex_lock(&LOCK);
			pos = get_pos(t1, t2);
			add_pair(pos, p);
			pthread_mutex_unlock(&LOCK);
		}
	}
}