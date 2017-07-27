#include "hash_align.h"

#define RATIO 5
#define COMPAIR(x, y) (x > RATIO*y? 1: (y > RATIO*x? -1: 0))

inline float ratio(unsigned int cover, unsigned int len)
{
        return (float) cover / (len - KMER);
}

unsigned long check_node(unsigned long id, unsigned int *pos,
                                 Graph *graph, unsigned int *hash_map)
{
	unsigned long tmp = id;

	*pos = hash_map[id & (MIL - 1)];
	while (*pos > 0 && graph->node[*pos - 1].id != id){
		++tmp;
		*pos = hash_map[tmp & (MIL - 1)];
	}
	return tmp;
}

unsigned int new_node(Graph *graph, unsigned long id,
                         unsigned long tmp, unsigned int *hash_map)
{
        unsigned int i = graph->n;
        ++graph->n;
	hash_map[tmp & (MIL - 1)] = graph->n;
        graph->node = realloc(graph->node, graph->n*sizeof(Node));

        graph->node[i].id = id;
        memset(graph->node[i].child, 0, N_CHAR*sizeof(Child));
        memset(graph->node[i].parent, 0, N_CHAR*sizeof(Parent));
	graph->node[i].n_child = graph->node[i].n_parent = 0;

        return graph->n;
}

unsigned int add_node(unsigned long id, Graph *graph, unsigned int *hash_map)
{
        unsigned int pos;
        unsigned long tmp;

	tmp = check_node(id, &pos, graph, hash_map);
        if (pos == 0)
                pos = new_node(graph, id, tmp, hash_map);
        return pos - 1;
}

void add_edge(unsigned long u, unsigned long v,
                 Graph *graph, unsigned int *hash_map)
{
        unsigned int u_pos, v_pos, n, i;

        u_pos = add_node(u, graph, hash_map);
        v_pos = add_node(v, graph, hash_map);

        n = graph->node[u_pos].n_child;
        for (i = 0; i < n; ++i){
                if (graph->node[u_pos].child[i].id == v_pos){
                        ++graph->node[u_pos].child[i].cover;
                        return;
                }
        }

        ++graph->node[u_pos].n_child;
        graph->node[u_pos].child[i].id = v_pos;
        graph->node[u_pos].child[i].cover = 1;
        graph->node[u_pos].child[i].len = KMER + 1;
        graph->node[u_pos].child[i].seq = value_to_key(u);
        graph->node[u_pos].child[i].seq[KMER] = ALPHABET[v & 3];

        n = graph->node[v_pos].n_parent++;
        graph->node[v_pos].parent[n].id = u_pos;
}

unsigned int get_child(Graph *graph, Parent parent, unsigned int id)
{
        unsigned int i, p_id;

        p_id = parent.id;
        for (i = 0; i < graph->node[p_id].n_child; ++i)
                if (graph->node[p_id].child[i].id == id &&
                        graph->node[p_id].child[i].cover == parent.cover &&
                             graph->node[p_id].child[i].len == parent.len)
                        return i;
        return 0;
}

void remove_edge(Graph *graph, unsigned int i, unsigned int k)
{
        unsigned int j, w, l;
        unsigned int child = graph->node[i].child[k].id;

        w = graph->node[i].child[k].cover;
        l = graph->node[i].child[k].len;
        free(graph->node[i].child[k].seq);
        memmove(graph->node[i].child + k, graph->node[i].child + k + 1,
                        (graph->node[i].n_child - k - 1)*sizeof(Child));


        --graph->node[i].n_child;
        for (j = 0; j < graph->node[child].n_parent; ++j)
                if (graph->node[child].parent[j].id == i &&
                        graph->node[child].parent[j].cover == w &&
                             graph->node[child].parent[j].len == l)
                        break;
        memmove(graph->node[child].parent + j,
                        graph->node[child].parent + j + 1,
                        (graph->node[child].n_parent - j - 1)*sizeof(Parent));
        --graph->node[child].n_parent;
        
}

void split_read(char *read, Graph *graph, unsigned int *hash_map)
{
        unsigned int i, len, ch;
        unsigned long u, v;

        for (i = u = v = len = 0; read[i] != '\0'; ++i){
                ch = ascii_table[read[i]];
                if (ch >= N_CHAR){
                        len = v = 0;
                        continue;
                }
                u = v;
                v = POP_HEAD(v) | ch;
                ++len;
                if (len > KMER){
                        add_edge(u, v, graph, hash_map);
                        add_edge(reverse(v), reverse(u), graph, hash_map);
                }
        }
}

void edge_cover(Graph *graph, unsigned int p, Child child)
{
        unsigned int i;

        for (i = 0; i < graph->node[child.id].n_parent; ++i){
                if (graph->node[child.id].parent[i].id == p){
                        graph->node[child.id].parent[i].cover = child.cover;
                        graph->node[child.id].parent[i].len = child.len;
                        return;
                }
        }
}

void contruct_graph(unsigned int idx, Graph *graph)
{
        unsigned int i, k;
        unsigned int *hash_map = calloc(MIL, sizeof(int));

        for (i = 0; i < FUSION->detail[idx].n; ++i){
                k = FUSION->detail[idx].read[i];
                split_read(FUSION->read[k], graph, hash_map);
                split_read(FUSION->read[k + 1], graph, hash_map);
        }

        for (i = 0; i < graph->n; ++i)
                for (k = 0; k < graph->node[i].n_child; ++k)
                        edge_cover(graph, i, graph->node[i].child[k]);
        free(hash_map);
}

void merge_edge(Child *source, unsigned int i, Child *target)
{
        unsigned int len, add;

        len = source[i].len;
        add = target[0].len - KMER;

        source[i].seq = realloc(source[i].seq, len + add + 1);
        memcpy(source[i].seq + len, target[0].seq + KMER, add);
        free(target[0].seq);

        source[i].len += add;
        source[i].cover += target[0].cover;
}

void condense_graph(Graph *graph, unsigned int idx)
{
        
        unsigned int i, child, parent, c_id, w, l;
        int n_id;

        parent = graph->node[idx].parent[0].id;
        child = graph->node[idx].child[0].id;

        if (parent == child)
                return;

        c_id = get_child(graph, graph->node[idx].parent[0], idx);
        graph->node[parent].child[c_id].id = child;
        merge_edge(graph->node[parent].child, c_id, graph->node[idx].child);
        w = graph->node[parent].child[c_id].cover;
        l = graph->node[parent].child[c_id].len;

        for (i = 0; i < graph->node[child].n_parent; ++i){
                if (graph->node[child].parent[i].id == idx &&
                    graph->node[child].parent[i].cover == 
                        graph->node[idx].child[0].cover &&
                    graph->node[child].parent[i].len == 
                        graph->node[idx].child[0].len){
                        graph->node[child].parent[i].id = parent;
                        graph->node[child].parent[i].cover = w;
                        graph->node[child].parent[i].len = l;
                        break;
                }
        }
        graph->node[idx].child[0].len = 0;
        graph->node[idx].n_child = graph->node[idx].n_parent = 0;
}

void simplify_graph(Graph *graph)
{
        unsigned int i, k, n_edge;

        for (i = 0; i < graph->n; ++i){
                if (graph->node[i].n_child == 1 && 
                        graph->node[i].n_parent == 1)
                        condense_graph(graph, i);
        }
}

unsigned int check_child(Child child, Graph *graph)
{
        unsigned int c = child.id;

        if (graph->node[c].n_child == 0 && graph->node[c].n_parent == 1 &&
                 (ratio(child.cover, child.len) < 2 || child.len <= KMER))
                return 1;
        return 0;
}

unsigned int check_parent(unsigned int p, unsigned int id, Graph *graph)
{
        if (graph->node[p].n_parent == 0 && graph->node[p].n_child == 1 &&
           (ratio(graph->node[p].child[0].cover, graph->node[p].child[0].len) < 2 ||
                                        graph->node[p].child[0].len <= KMER))
                return 1;
        return 0;
}

unsigned int remove_odd_child(Graph *graph, unsigned int id)
{
        if (graph->node[id].n_child < 2)
                return 0;

        int i;
        float max_ratio, r;
        unsigned int count;

        for (i = 0, max_ratio = 0.0; i < graph->node[id].n_parent; ++i){
                r = ratio(graph->node[id].parent[i].cover, 
                                  graph->node[id].parent[i].len);
                max_ratio = MAX2(max_ratio, r);
        }

        for (i = count = 0; i < graph->node[id].n_child; ++i){
                r = ratio(graph->node[id].child[i].cover,
                                  graph->node[id].child[i].len);
                if (r*RATIO >= max_ratio)
                        ++count;
        }

        if (count > 0){
                count = 0;
                for (i = graph->node[id].n_child - 1; i >= 0; --i){
                        r = ratio(graph->node[id].child[i].cover,
                                        graph->node[id].child[i].len);
                        if (r*RATIO < max_ratio){
                                ++count;
                                remove_edge(graph, id, i);
                        }
                }
        }
        return count;
}

unsigned int remove_odd_parent(Graph *graph, unsigned int id)
{
        if (graph->node[id].n_parent < 2)
                return 0;

        int i;
        float max_ratio, r;
        unsigned int count, c_id;

        for (i = 0, max_ratio = 0.0; i < graph->node[id].n_child; ++i){
                r = ratio(graph->node[id].child[i].cover,
                                  graph->node[id].child[i].len);
                max_ratio = MAX2(max_ratio, r);
        }

        for (i = count = 0; i < graph->node[id].n_parent; ++i){
                r = ratio(graph->node[id].parent[i].cover,
                                  graph->node[id].parent[i].len);
                if (r*RATIO >= max_ratio)
                        ++count;
        }

        if (count > 0){
                count = 0;
                for (i = graph->node[id].n_parent - 1; i >= 0; --i){
                        r = ratio(graph->node[id].parent[i].cover,
                                        graph->node[id].parent[i].len);
                        if (r*RATIO < max_ratio){
                                ++count;
                                c_id = get_child(graph, graph->node[id].parent[i], id);
                                remove_edge(graph, graph->node[id].parent[i].id, c_id);
                        }
                }
        }
        return count;
}

unsigned int remove_child(Graph *graph, unsigned int id)
{
        int i;
        unsigned int count, max_w, max_l;

        for (i = count = max_w = max_l = 0; i < graph->node[id].n_child; ++i){
                count += check_child(graph->node[id].child[i], graph);
                max_w = MAX2(max_w, graph->node[id].child[i].cover);
                max_l = MAX2(max_l, graph->node[id].child[i].len);
        }

        if (count == graph->node[id].n_child){
                for (count = 0, i = graph->node[id].n_child - 1; i >= 0; --i){
                        if ((max_w == 1 && graph->node[id].child[i].len < max_l) ||
                            (max_w > 1 && graph->node[id].child[i].cover < max_w)){
                                remove_edge(graph, id, i);
                                ++count;
                        }
                }
        } else if (count > 0){
                for (i = graph->node[id].n_child - 1; i >= 0; --i)
                        if (check_child(graph->node[id].child[i], graph) == 1)
                                remove_edge(graph, id, i);
        } else {
                count = 0;
        }
        count += remove_odd_child(graph, id);
        return count;
}

unsigned int remove_parent(Graph *graph, unsigned int id)
{
        unsigned int count, max_w, max_l;
        int i;

        for (i = count = max_w = max_l = 0; i < graph->node[id].n_parent; ++i){
                count += check_parent(graph->node[id].parent[i].id, id, graph);
                max_w = MAX2(max_w, graph->node[id].parent[i].cover);
                max_l = MAX2(max_l, graph->node[id].parent[i].len);
        }

        if (count == graph->node[id].n_parent){
                for (count = 0, i = graph->node[id].n_parent - 1; i >= 0; --i){
                        if ((max_w == 1 && graph->node[id].parent[i].len < max_l) ||
                             (max_w > 1 && graph->node[id].parent[i].cover < max_w)){
                                remove_edge(graph, graph->node[id].parent[i].id, 0);
                                ++count;
                             }
                }
        } else if (count > 0){
                for (i = graph->node[id].n_parent - 1; i >= 0; --i)
                        if (check_parent(graph->node[id].parent[i].id, id, graph) == 1)
                                remove_edge(graph, graph->node[id].parent[i].id, 0);
        } else {
                count = 0;
        }
        count += remove_odd_parent(graph, id);
        return count;        
}

unsigned int remove_bubble(Graph *graph, unsigned int id)
{
        unsigned int i, child, c_id;
        float r, max_ratio;
        
        max_ratio = ratio(graph->node[id].child[0].cover,
                                graph->node[id].child[0].len);
        c_id = graph->node[id].child[0].id;
        child = 0;

        if (graph->node[c_id].n_parent != graph->node[id].n_child)
                return 0;

        for (i = 1; i < graph->node[id].n_child; ++i){
                if (graph->node[id].child[i].id != c_id)
                        return 0;
                
                r = ratio(graph->node[id].child[i].cover,
                                graph->node[id].child[i].len);
                if (r > max_ratio){
                        child = i;
                        max_ratio = r;
                }
        }
        
        for (i = graph->node[id].n_child; i > 0; --i)
                if (i - 1 != child)
                        remove_edge(graph, id, i - 1);
        return 1;
}

unsigned int remove_branch(Graph *graph)
{
        unsigned int i, count;
        for (i = count = 0; i < graph->n; ++i){
                if (graph->node[i].n_child > 1)
                        count += remove_child(graph, i);

                if (graph->node[i].n_parent > 1)
                        count += remove_parent(graph, i);

                if (graph->node[i].n_child > 1)
                        count += remove_bubble(graph, i);
        }

        return count;
}

unsigned int match_end(char *seq, unsigned int start, unsigned int end,
				unsigned int ref, unsigned int r_start)
{
	unsigned int i, k, count;
	unsigned int s = ref == 0? 0: REF_INF->find[ref - 1];
	char *ref_seq = REF_INF->seq + s;

	k = r_start + 1;
	i = start + KMER;
	for (count = 0; i < end; ++i, ++k, ++count){
		if (seq[i] != ref_seq[k])
			return count;
	}

	return count;
}

unsigned int match_start(char *seq, unsigned int start,
				unsigned int ref, unsigned int r_start)
{
	int i, k, count;
	unsigned int s = ref == 0? 0: REF_INF->find[ref - 1];
	char *ref_seq = REF_INF->seq + s;

	k = r_start - KMER;
	i = start - 1;

	for (count = 0; i >= 0; --i, --k, ++count){
		if (seq[i] != ref_seq[k])
			return count + KMER;
	}

	return count + KMER;
}

void get_pair(char *seq, unsigned int **pos, unsigned int n,
			 Fusion_pair *trans, unsigned int len)
{
	unsigned int i, k, p, count, add_start, add_end;
	int score;

	for (i = p = 0; i < n;){
		count = 1;
		score = 0;
		p = pos[2][i] > len/3? 1: 0;
		++i;
		while(i < n && pos[0][i] == pos[0][i - 1]){
			score += KMER - ABS((int)(pos[1][i] - pos[1][i - 1]) -
					    (int)(pos[2][i] - pos[2][i - 1]));
			++count;
			++i;
		}

		add_start = match_start(seq, pos[2][i - count], 
					  pos[0][i - 1], pos[1][i - count]);
		add_end = match_end(seq, pos[2][i - 1], len,
					   pos[0][i - 1], pos[1][i - 1]);
		score += add_start + add_end;
		if (score > trans[p].count){
			trans[p].trans = pos[0][i - 1];
			trans[p].end = pos[1][i - 1];
			trans[p].end += add_end;
			
			trans[p].start  = pos[1][i - count];
			trans[p].start -= add_start;
			trans[p].count = score;
		}
	}
}

void get_transcript(char *seq, unsigned int len, Fusion_pair *trans)
{
	unsigned int i, k, p, mod, ret;
	unsigned long idx;
	unsigned int *pos1[3], *pos2[3], n1, n2;

	for (i = n1 = n2 = 0; i + KMER <= len; i += KMER){
		idx = get_index(i, seq);
		p = hash_get(idx, &ret, 0);

		if (ret == 0)
			continue;

		mod = idx & (H - 1);
		n1 = KMER_HASH->bucket[mod][p].end -
				KMER_HASH->bucket[mod][p].start;
		pos1[0] = KMER_HASH->pos + KMER_HASH->bucket[mod][p].start;
		pos1[1] = KMER_HASH->p + KMER_HASH->bucket[mod][p].start;

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
	get_pair(seq, pos2, n2, trans, len);

	if (n2 > 0){
		free(pos2[0]);
		free(pos2[1]);
		free(pos2[2]);
	}
}

unsigned int get_exon_border(unsigned int pos, unsigned int trans,
				unsigned int exon, unsigned int type)
{
	unsigned int dis_start = ABS(pos - GENE_MAP->start[trans][exon]);
	unsigned int dis_end = ABS(pos - GENE_MAP->end[trans][exon]);
	if (type == 0){
		if (exon > 0 && dis_start < dis_end && dis_start < MIN_FRAG)
			return GENE_MAP->end[trans][exon - 1] + 1;
		else if (dis_end < MIN_FRAG)
			return GENE_MAP->end[trans][exon] + 1;
	} else if (type == 1){
		if (dis_start < dis_end && dis_start < MIN_FRAG)
			return GENE_MAP->start[trans][exon] + 1;
		else if (exon < GENE_MAP->n_exon[trans] && dis_end < MIN_FRAG)
			return GENE_MAP->start[trans][exon + 1] + 1;
	}
	return pos + 1;
}

void print_fusion(char *seq, unsigned int t1, unsigned int t2,
                             unsigned int g1, unsigned int g2,
			 unsigned int start, unsigned int end,
			   unsigned int len, unsigned int str)
{
	unsigned int exon, p1, p2;
	unsigned int chr1 = GENE_MAP->chr[t1];
	unsigned int chr2 = GENE_MAP->chr[t2];
	char name1[16], name2[16];

	p1 = gene_coord(t1, end, &exon);
	p1 = get_exon_border(p1, t1, exon, len == 0? 1 - str: str);

	p2 = gene_coord(t2, start, &exon);
	p2 = get_exon_border(p2, t2, exon, len == 1? str: 1 - str);

	memcpy(name1, GENE_MAP->gene_name + (g1*GENE_MAP->l_gene + 1), 15);
	memcpy(name2, GENE_MAP->gene_name + (g2*GENE_MAP->l_gene + 1), 15);
	name1[15] = '\0';
	name2[15] = '\0';


	pthread_mutex_lock(&LOCK);
	fprintf(OUT_FUSION, "%s\t%u\t%s-(%c)%s\t%s\t%u\t%s-(%c)%s\t%s\n",
	            GENE_MAP->chr_name + GENE_MAP->l_chr*chr1, p1, name1,
	  		        GENE_MAP->gene_name[g1*GENE_MAP->l_gene],
	  	        GENE_MAP->gene_name + (g1*GENE_MAP->l_gene + 16),
	            GENE_MAP->chr_name + GENE_MAP->l_chr*chr2, p2, name2,
	  		        GENE_MAP->gene_name[g2*GENE_MAP->l_gene], 
	  	        GENE_MAP->gene_name + (g2*GENE_MAP->l_gene + 16),
	 					     	            seq);
	
	pthread_mutex_unlock(&LOCK);	
}

void get_gene_pair(char *seq)
{
	Fusion_pair trans[2][2];
	unsigned int len, t;

	len = strlen(seq);
	trans[0][0].trans = trans[0][1].trans = -1;
	trans[1][0].trans = trans[1][1].trans = -1;
	trans[0][0].count = trans[0][1].count = -len*KMER;
	trans[1][0].count = trans[1][1].count = -len*KMER;

	get_transcript(seq, len, trans[0]);

	reverse_str(seq, len);
	get_transcript(seq, len, trans[1]);

	t = 0;
	if (trans[0][0].count + trans[0][1].count <
		trans[1][0].count + trans[1][1].count)
		t = 1;

	if (trans[t][0].count < trans[1 - t][1].count){
		trans[t][0] =  trans[1 - t][1];
		len = trans[t][0].start;
		trans[t][0].start = trans[t][0].end;
		trans[t][0].end = len;
		len = 0;
	}

	if (trans[t][1].count < trans[1 - t][0].count){
		trans[t][1] = trans[1 - t][0];
		len = trans[t][1].start;
		trans[t][1].start = trans[t][1].end;
		trans[t][1].end = len;
		len = 1;
	}

	if (trans[t][0].trans == -1 || trans[t][1].trans == -1)
		return;

	unsigned int t1 = trans[t][0].trans;
	unsigned int t2 = trans[t][1].trans;
	unsigned int g1 = GENE_MAP->gene[t1];
	unsigned int g2 = GENE_MAP->gene[t2];
	int str1 = GENE_MAP->gene_name[g1*GENE_MAP->l_gene] == '+'? 1: -1;
	int str2 = GENE_MAP->gene_name[g1*GENE_MAP->l_gene] == '+'? 1: -1;

	if (g1 == g2 || (str1 < 0 && len != 0 && str2 > 0 && len != 1))
		return;

	if (str1 > 0 || (str1 < 0 && len == 0))
		print_fusion(seq, t1, t2, g1, g2,
			trans[t][1].start, trans[t][0].end, len, 0);
	else
		print_fusion(seq, t2, t1, g2, g1,
			trans[t][0].end, trans[t][1].start, len, 1);
}

void *assembly_fusion(void *id)
{
        unsigned int process, i, k, n;
        long idx;
        char *seq;
        Graph *graph = malloc(sizeof(Graph));
        graph->node = calloc(1, sizeof(Node));
        graph->n = 0;
        idx = (long) id;

        contruct_graph(idx, graph);
        do {
                simplify_graph(graph);
                process = remove_branch(graph);
        } while (process > 0);

        seq = malloc(1);
        seq[0] = '\0';
        for (i = 0, n = 2*MEAN_LEN; i < graph->n; ++i){
                for (k = 0; k < graph->node[i].n_child; ++k){
                        if (graph->node[i].child[k].len > n){
                                free(seq);
                                seq = graph->node[i].child[k].seq;
                                n = graph->node[i].child[k].len;
                                seq[n] = '\0';
                        } else {
                                free(graph->node[i].child[k].seq);
                        }
                }
        }
       
        free(graph->node);
        free(graph);

        if (seq[0] != '\0')
                get_gene_pair(seq);
        free(seq);
	pthread_exit(NULL);
	return NULL;
}