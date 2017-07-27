#include "EM.h"

/******************************************************************************
 * ************************************************************************** *
 * *                               HDF5 WRITE                               * *
 * ************************************************************************** * 
 ******************************************************************************/

/******************************************************************************
 * Write array of integer to hdf5 file                                        *
 * @ len   : length of array                                                  *
 * @ vec   : array of integer                                                 *
 * @ folder: folder to store array in hdf5 file                               *
 * @ name  : name of array                                                    *
 *****************************************************************************/
void write_h5_int(unsigned int len, int *vec, hid_t folder, char *name)
{
	herr_t status;
	hsize_t dims[1] = {len};
	hid_t prop_id, dataspace, dataset;

	prop_id = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(prop_id, 1, dims);
	status = H5Pset_deflate(prop_id, COMPRESS);
	dataspace = H5Screate_simple(1, dims, NULL);
	dataset = H5Dcreate(folder, name, H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT, prop_id, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_INT,
			H5S_ALL, H5S_ALL, H5P_DEFAULT, vec);
	status = H5Pclose(prop_id);
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);
}

/******************************************************************************
 * Write array of double to hdf5 file                                         *
 * @ len   : length of array                                                  *
 * @ vec   : array of double                                                  *
 * @ folder: folder to store array in hdf5 file                               *
 * @ name  : name of array                                                    *
 *****************************************************************************/
void write_h5_double(unsigned int len, double *vec, hid_t folder, char *name)
{
	herr_t status;
	hsize_t dims[1] = {len};
	hid_t prop_id, dataspace, dataset;

	prop_id = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(prop_id, 1, dims);
	status = H5Pset_deflate(prop_id, COMPRESS);
	dataspace = H5Screate_simple(1, dims, NULL);
	dataset = H5Dcreate(folder, name, H5T_NATIVE_DOUBLE,
			dataspace, H5P_DEFAULT, prop_id, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE,
			H5S_ALL, H5S_ALL, H5P_DEFAULT, vec);
	status = H5Pclose(prop_id);
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);
}

/******************************************************************************
 * Write set of strings with equal length to hdf5 file                        *
 * @ len   : total length of strings set                                      *
 * @ max   : length of each string in set                                     *
 * @ string: concatenation of strings in set                                  *
 * @ folder: folder to store string set in hdf5                               *
 * @ name  : name of string set                                               *
 *****************************************************************************/
void write_h5_string(unsigned int len, unsigned int max,
                                        char *string, hid_t folder, char *name)
{
	herr_t status;
	hsize_t dims[1] = {len};
	hid_t prop_id, dataspace, dataset, type;

	prop_id = H5Pcreate(H5P_DATASET_CREATE);
	type = H5Tcopy(H5T_C_S1);
	status = H5Pset_chunk(prop_id, 1, dims);
	status = H5Tset_size(type, max);
	status = H5Pset_deflate(prop_id, COMPRESS);
	dataspace = H5Screate_simple(1, dims, NULL);
	dataset = H5Dcreate(folder, name, type,
                                  dataspace, H5P_DEFAULT, prop_id, H5P_DEFAULT);
	status = H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, string);

	status = H5Pclose(prop_id);
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);
}

/******************************************************************************
 * Create distribution table of reads from all transcripts                    *
 * The occurence of each transcript in table is proportional                  *
 * with its read density                                                      *
 *****************************************************************************/
void create_distribution()
{
	unsigned int tmp[REF_INF->n], i, k;
	DIST = calloc(1, sizeof(Discrete_dist));

	tmp[0] = REF_INF->count[0];
	for (i = 1; i < REF_INF->n; ++i)
		tmp[i] = REF_INF->count[i] + tmp[i - 1];

	DIST->n = tmp[i - 1];
	DIST->table = calloc(DIST->n, sizeof(unsigned int));

	for (i = 0, k = 0; i < REF_INF->n; ++i)
		while (k < tmp[i])
			DIST->table[k++] = i;
}

/******************************************************************************
 * 1. Divide total read into transcript according to its distribution         *
 * @ DIST->n    : total number of mapped reads                                * 
 * @ DIST->table: table of read distribution                                  *
 * 2. Estimate read count                                                     *
 * @ em_val->overlap: store estimate read count for each transcript           * 
 * 3. Write result to file                                                    *
 *****************************************************************************/
void *bootstrap_sample(void *data)
{
	Thread_data3 *bs_data = (Thread_data3*) data;
	EM_val *em_val;
	unsigned int count[REF_INF->n];
	unsigned long i_start = 0;
	unsigned int n = DIST->n;
	unsigned int *table = DIST->table;
	char name[10];
	
	memset(count, 0, REF_INF->n * sizeof(int));
	while (--n > 0)
		++count[table[(rand_r(&seed)) % n]];

	em_val = estimate_count(count, bs_data->order);

	sprintf(name, "bs%u%c", bs_data->order, '\0');

	pthread_mutex_lock(&LOCK);
	write_h5_double(REF_INF->n, em_val->overlap, *bs_data->bs, name);
	pthread_mutex_unlock(&LOCK);

	destroy_emVal(em_val);
	free(data);

	pthread_exit(NULL);
	return NULL;
}

void write_fusion()
{
	unsigned long i;
	unsigned int k, rc;
	char *seq;
	pthread_t thr[NTHREAD];
	memset(&thr, '\0', NTHREAD);
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	fprintf(OUT_FUSION, "#left_chr\tleft_break\tleft_gene\tright_chr\tright_break\tright_gene\tfusion_seq\n");
	for (i = 0; i < FUSION->n;){
		for (k = 0; k < NTHREAD && i < FUSION->n; ++i) {
			if (FUSION->detail[i].n < 2)
				continue;
			rc = pthread_create(&thr[k], &attr,
					assembly_fusion, (void*)i);
			++k;
		}

		for (; k > 0; --k)
			rc = pthread_join(thr[k - 1], NULL);
	}
	pthread_attr_destroy(&attr);
}

void write_result(char *out_dir, unsigned int n_bs, EM_val *em_val)
{
	unsigned int i, k, l, v, g;
	float sum[2] = {0.0, 0.0};
	FILE *out;
	hid_t file_id, root, aux, bs;
	herr_t status;

	l = strlen(out_dir);
        char file_path[l + 50];
        v = 1; 
	if (out_dir[l - 1] == '/')
		--l;
	memcpy(file_path, out_dir, l);

	// Abundance
	memcpy(file_path + l, "/abundance.tsv\0", 15);
	out = fopen(file_path, "w");

	fprintf(out, "#target_id\tunique_map\tlength\teff_length\test_counts\ttpm\n");
	for (i = 0; i < REF_INF->n; ++i)
		fprintf(out, "%s\t%u\t%u\t%f\t%f\t%f\n",
			REF_INF->ref_name + i*REF_INF->name_len,
				REF_INF->count[i], REF_INF->len[i], 
				   REF_INF->eff_len[i], em_val->overlap[i],
				   			 em_val->grama[i]);
	fclose(out);

	// Genome abundance TODO: Estimate directly from read count
	memcpy(file_path + l, "/abundance.gene.tsv\0", 20);
	out = fopen(file_path, "w");
	fprintf(out, "#gene_id\tgene_name\tlength\test_counts\ttpm\n");
	for (i = 0; i < REF_INF->n; ++i){
		g = GENE_MAP->gene[i];
		sum[0] = em_val->overlap[i];
		sum[1] = em_val->grama[i];
		while (i < REF_INF->n && g == GENE_MAP->gene[i + 1]){
			 ++i;
			 sum[0] += em_val->overlap[i];
			 sum[1] += em_val->grama[i];
		}
		fprintf(out, "%s\t%s\t%f\t%f\n",
			REF_INF->ref_name + i*REF_INF->name_len + 16,
				GENE_MAP->gene_name + g*GENE_MAP->l_gene + 16,
							      sum[0], sum[1]);
	}
	fclose(out);

	// Fusion
	if (FUSION->n > 0){
		memcpy(file_path + l, "/fusion.bedpe\0", 14);
		OUT_FUSION = fopen(file_path, "w");
		write_fusion();
	}

	// HDF5
	if (n_bs > 0) {
		memcpy(file_path + l, "/abundance.h5\0", 14);
		file_id = H5Fcreate(file_path, H5F_ACC_TRUNC,
							H5P_DEFAULT, H5P_DEFAULT);
		root = H5Gopen(file_id, "/", H5P_DEFAULT);
		aux = H5Gcreate(file_id, "/aux", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		if (MAPPED == 0)
			write_h5_string(1, 11, "single-end\0", aux, "description");
		else
			write_h5_string(1, 11, "paired-end\0", aux, "description");

		write_h5_string(1, 7, "0.43.0\0", aux, "kallisto_version");
		write_h5_string(1, 8, "Unknown\0", aux, "start_time");
		write_h5_int(1, &v, aux, "index_version");

		write_h5_double(REF_INF->n, em_val->overlap, root, "est_counts");
		write_h5_double(REF_INF->n, REF_INF->eff_len, aux, "eff_lengths");
		write_h5_double(1, &MEAN_FRAG, aux, "mean_fragment_length");
		write_h5_double(1, &MEAN_LEN, aux, "mean_read_length");
		write_h5_int(REF_INF->n, REF_INF->len, aux, "lengths");
		write_h5_int(1, &n_bs, aux, "num_bootstrap");
		write_h5_string(REF_INF->n, REF_INF->name_len, 
					REF_INF->ref_name, aux, "ids");
		destroy_emVal(em_val);
		unsigned int rc, n;
		pthread_t thr[NTHREAD];
		memset(&thr, '\0', NTHREAD);
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setstacksize(&attr, STACK_SIZE);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		create_distribution();
		bs = H5Gcreate(file_id, "/bootstrap", H5P_DEFAULT,
                                                  H5P_DEFAULT, H5P_DEFAULT);

		for (i = 0, n = 0; i < n_bs;) {
			for (k = 0; k < NTHREAD && i < n_bs; ++k, ++i) {
				Thread_data3 *data = calloc(1,
                                                      sizeof(Thread_data3));
				data->order = i;
				data->bs = &bs;
				rc = pthread_create(&thr[k], &attr,
                                                    bootstrap_sample, data);
			}

			for (; k > 0; --k, ++n)
				rc = pthread_join(thr[k - 1], NULL);
		}
		pthread_attr_destroy(&attr);

		status = H5Gclose(bs);
		status = H5Fclose(file_id);
		status = H5Gclose(root);
		status = H5Gclose(aux);
	}
}

/******************************************************************************
 * ************************************************************************** *
 * *                             ESTIMATE COUNT                             * *
 * ************************************************************************** * 
 ******************************************************************************/

double *M_step(EM_val *em_val, unsigned int *count, unsigned int *stop)
{
	unsigned int i, k;
	double total, max_err, err;
	double *new_grama = calloc(REF_INF->n, sizeof (double));
	for (i = 0; i < CLASS->n; ++i){
                for (k = 0, total = 0.0; k < CLASS->cls[i].n; ++k)
			total += em_val->grama[CLASS->cls[i].ref[k]];
                if (total == 0)
			continue;

		for (k = 0; k < CLASS->cls[i].n; ++k)
			em_val->overlap[CLASS->cls[i].ref[k]] +=
	                        em_val->grama[CLASS->cls[i].ref[k]] *
					 CLASS->cls[i].count / total;
        }

	for (i = 0, max_err = 0.0, *stop = 1; i < REF_INF->n; i++) {
		new_grama[i] = (double) (count[i] + em_val->overlap[i]) /
       	                                                 REF_INF->len[i];
		em_val->overlap[i] = 0.0;
		if (new_grama[i] > 1e-3){
			err = ABS(new_grama[i] - em_val->grama[i])/new_grama[i];
			max_err = err > max_err? err: max_err;
			if (max_err > 0.01)
				*stop = 0;
		}
	}

	em_val->grama = new_grama;
	return new_grama;
}

void E_step(unsigned int *count, EM_val *em_val)
{
	unsigned int i;

	em_val->grama = calloc(REF_INF->n, sizeof (double));
	em_val->overlap = calloc(REF_INF->n, sizeof (double));

	for (i = 0; i < REF_INF->n; i++)
		em_val->grama[i] = (double) 1/REF_INF->n;
}

EM_val *estimate_count(unsigned int *count, unsigned int order)
{
	DEBUG_PRINT("Running EM ...");
	fflush(stdout);

	unsigned int i, k, stop, minRound, maxRound;
	double *u, *u1, *u2, sum, error, alpha, r[REF_INF->n], v[REF_INF->n];
	EM_val *em_val = malloc(sizeof (EM_val));

	E_step(count, em_val);

	stop = 0;
	minRound = 50;
	maxRound = 1000;
	k = 0;
	while (k < minRound || (k < maxRound && stop == 0)) {
		u = em_val->grama;
		u1 = M_step(em_val, count, &stop);
		u2 = M_step(em_val, count, &stop);

		for (i = 0; i < REF_INF->n; ++i) {
			r[i] = u1[i] - u[i];
			v[i] = (u2[i] - u1[i]) - r[i];
		}

		alpha = -Norm(r, REF_INF->n) / Norm(v, REF_INF->n);

		if (alpha > -1) alpha = -1;
		free(u2);

		for (i = 0; i < REF_INF->n; ++i) {
			u1[i] = u[i] - 2 * alpha * r[i] + alpha * alpha * v[i];
			if (u1[i] < 0) u1[i] = 0;
		}

		em_val->grama = u1;
		M_step(em_val, count, &stop);

		DEBUG_PRINT("\rRunning EM ...%u round", k);
		fflush(stdout);
		free(u);
		free(u1);
		++k;
	}

	for (i = 0, sum = 0.0; i < REF_INF->n; ++i)
		sum += em_val->grama[i];

	for (i = 0; i < REF_INF->n; ++i) {
		em_val->overlap[i] = em_val->grama[i] * REF_INF->len[i];
		em_val->grama[i] = em_val->grama[i]*1000000 / sum;
	}

	if (order == 0)
		DEBUG_PRINT("\rFinish EM with %u rounds\n", k);
	else
		DEBUG_PRINT("\rFinish bootstrap number %u with %u rounds",
                                                                 order, k);

	return em_val;
}

destroy_emVal(EM_val *em_val)
{
	free(em_val->overlap);
	free(em_val->grama);
	free(em_val);
}
