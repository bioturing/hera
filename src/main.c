#include "EM.h"
#include "bam_write.h"

void print_usage()
{
	fprintf(stdout, "VERSION:\n\t%s\n", VERSION);
	fprintf(stdout, "QUANTILITY:\n\t ./hera quant -i path/to/index_directory [options] -1 <R1.fq> -2 <R2.fq>\n");
	fprintf(stdout, "\t\t-1:\t Input left read-files, separated by comma\n");
	fprintf(stdout, "\t\t-2:\t Input right read-files, separated by comma (optional)\n");
	fprintf(stdout, "\t\t-o:\t Output directory (default: ./)\n");
	fprintf(stdout, "\t\t-t:\t Number of threads (default: 1)\n");
	fprintf(stdout, "\t\t-z:\t Compress level (1 - 9) (default: -1)\n");
	fprintf(stdout, "\t\t-b:\t Number of bootstraps (default: 0)\n");
	fprintf(stdout, "\t\t-w:\t Output bam file 0:true, 1: false (default: 0)\n");
	fprintf(stdout, "\t\t-f:\t Genome fasta file (if not define, genome mapping will be ignore\n");
	fprintf(stdout, "\t\t-p:\t Output prefix (default: '')\n");
	exit(EXIT_FAILURE);
}

void directory_exist(const char* sFilePath)
{
	struct stat oBuffer;

	if (stat(sFilePath, &oBuffer) == 0)
		if (oBuffer.st_mode & S_IFDIR)
			return;

	mkdir(sFilePath, 0755);
}

void index_ref(char *transcriptome, char *genome, char *out_file, int full)
{
	// Initialize
	init_refInf();
	init_hash();

	//Read reference sequence + construct hash for transcriptome    
	printf("Index transcripts in: %s\n", transcriptome);
	get_ref_seq(transcriptome);

	// Write to file
	printf("Number of sequences\t: %u\n", REF_INF->n);
	write_to_file(out_file);

	// Index genome 
	if (full == 1){
		printf("Index genome in: %s\n", genome);
		indexGenome(genome, out_file);
	}
}

void get_file_name(char *file_list, char *file_name, unsigned int *idx)
{
	unsigned int i = *idx;
	while (file_list[*idx] != ',' && file_list[*idx] != '\0')
		++*idx;
	memcpy(file_name, file_list + i, *idx - i);
	file_name[*idx - i] = '\0';
}

void process_read_files(char *left_read, char *right_read)
{
	char left_file[1000], right_file[1000];
	unsigned int i, l_len, r_len;

	l_len = r_len = 0;

	while (1){
		get_file_name(left_read, left_file, &l_len);
		get_file_name(right_read, right_file, &r_len);
		get_alignment_pair(left_file, right_file);
		if (left_read[l_len] == '\0' || right_read[r_len] == '\0')
			break;
		++l_len;
		++r_len;
	}

	MEAN_LEN /= 2 * N_READ;
	MEAN_FRAG = MEAN_FRAG / PAIRED + MEAN_LEN;
	for (i = 0; i < REF_INF->n; ++i) {
		if (REF_INF->len[i] < MEAN_FRAG)
			REF_INF->eff_len[i] = 0;
		else
			REF_INF->eff_len[i] = REF_INF->len[i] - MEAN_FRAG + 1;
	}
	fprintf(stdout, "\rNumber of aligned pairs\t: %u/%u\n", MAPPED, N_READ);
	fprintf(stdout, "Mean read length\t: %f\n", MEAN_LEN);
	fprintf(stdout, "Mean fragment length\t: %f\n", MEAN_FRAG);

	fprintf(SUMMARY, "\nMAPPING RESULT:\n");
	fprintf(SUMMARY, "\tTotal number of read pairs\t\t\t: %u\n", N_READ);
	fprintf(SUMMARY, "\t\t- Number of proper mapped pairs\t: %u (%f)\n",														PAIRED, (float) PAIRED/N_READ);
	fprintf(SUMMARY, "\t\t- Number of unproper mapped pairs\t: %u (%f)\n",
	  		MAPPED -  PAIRED, (float)(MAPPED -  PAIRED)/N_READ);
	fprintf(SUMMARY, "\t\t- Number of unmapped pairs\t\t: %u (%f)\n",
			N_READ - MAPPED, (float) (N_READ - MAPPED)/N_READ);
	fprintf(SUMMARY, "\tMean read length\t\t\t: %f\n", MEAN_LEN);
	fprintf(SUMMARY, "\tMean fragment length\t\t\t: %f\n", MEAN_FRAG);
}

void process_read_file(char *left_read)
{
	char left_file[1000];
	unsigned int i, l_len = 0;
	while (1){
		get_file_name(left_read, left_file, &l_len);
		get_alignment(left_file);
		if (left_read[l_len] == '\0')
			break;

		++l_len;							
	}

	MEAN_LEN /= N_READ;
	MEAN_FRAG = MEAN_LEN;
	for (i = 0; i < REF_INF->n; ++i) {
		if (REF_INF->len[i] < MEAN_FRAG)
			REF_INF->eff_len[i] = 0;
		else
			REF_INF->eff_len[i] = REF_INF->len[i] - MEAN_FRAG + 1;
	}
	fprintf(stdout, "\rNumber of aligned reads\t: %u/%u\n", PAIRED, N_READ);
	fprintf(stdout, "Mean read lenght\t: %f\n", MEAN_LEN);

	fprintf(SUMMARY, "\nMAPPING RESULT:\n");
	fprintf(SUMMARY, "\tTotal number of reads\t\t\t: %u\n", N_READ);
	fprintf(SUMMARY, "\t\t- Number of mapped reads\t: %u (%f)\n",
			PAIRED, (float) PAIRED/N_READ);
	fprintf(SUMMARY, "\t\t- Number of unmapped reads\t: %u (%f)\n",
			N_READ - PAIRED, (float) (N_READ - PAIRED)/N_READ);
	fprintf(SUMMARY, "\tMean read length\t\t\t: %f\n", MEAN_LEN);
}

void quantility(int argc, char *argv[]) {
	char *idx_dir, *left_read, *right_read;
	char *out_dir = "./", *genome = NULL, *prefix = "\0";
	EM_val *em_val;
	int c, temp_nthread, temp_writebam, temp_nbs;
	unsigned int n_bs = 0;
	unsigned long i_start = 0;

	// Get parameter
	left_read = right_read = NULL;
	while ((c = getopt(argc - 1, argv + 1, "1:2:t:i:o:b:z:h:w:f:p:")) != -1) {
		switch (c) {
			case '1':
				left_read = optarg;
				break;
			case '2':
				right_read = optarg;
				break;
			case 't':
				temp_nthread = atoi(optarg);
				if (temp_nthread <= 0) {
					fprintf(stdout, "Invalid number of thread!!!\n");
					exit(EXIT_FAILURE);     
				}
				NTHREAD = temp_nthread;
				break;
			case 'b':
				temp_nbs = atoi(optarg);
				if (temp_nbs < 0) {
				       fprintf(stdout, "Invalid number of bootstraps!!!\n");
				       exit(EXIT_FAILURE);     
				}
				n_bs = temp_nbs;
				break;
			case 'i':
				idx_dir = optarg;
				break;
			case 'p':
				prefix = optarg;
				break;
			case 'f':
				genome = optarg;
				break;
			case 'o':
				out_dir = optarg;
				break;
			case 'z':
				COMPRESS_LEVEL = atoi(optarg);
				if (COMPRESS_LEVEL < -1 || COMPRESS_LEVEL > 9) {
					fprintf(stdout, "Invalid compress level!!!\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 'w':
				temp_writebam = atoi(optarg);
				if (temp_writebam < 0 || temp_writebam > 1) {
					fprintf(stdout, "Invalid output bam file parameter!!!\n");
					exit(EXIT_FAILURE);     
				}
				WRITE_BAM = temp_writebam;
				break;
			case 'h':
				print_usage();
				break;
			case '?':
				print_usage();
				break;
		}
	}
	directory_exist(out_dir);

	// Retrieve index
	init_refInf();
	init_hash();
	load_index(idx_dir, genome);
	init_output(idx_dir, out_dir, argc, argv, prefix);

	//Get alignment
	if (left_read != NULL && right_read != NULL) {
		fprintf(stdout, "Process paired-end reads in:\n\t 1. %s\n\t 2. %s\n",
					     		left_read, right_read);
		process_read_files(left_read, right_read);
	} else if (left_read != NULL) {
		fprintf(stdout, "Process single-end reads in: %s\n", left_read);
		process_read_file(left_read);
	} else {
		fprintf(stdout, "Input file is not defined\n");
		print_usage();
	}

	// Flush and end file
	if (WRITE_BAM == 0)
		bam_close(BAM);

	// Correct estimate count
	em_val = estimate_count(REF_INF->count, 0);

	// Finish
	write_result(out_dir, n_bs, em_val, prefix);
	fclose(SUMMARY);
}

int main(int argc, char *argv[]) {
	printf("Hera is a program developed by BioTuring for RNA-Seq analysis.\nPlease contact info@bioturing.com if you need further support\n");
	if (argc < 3)
		print_usage();
	else if (strcmp(argv[1], "index") == 0) {
		index_ref(argv[2], argv[3], argv[4], atoi(argv[5]));
	} else if (strcmp(argv[1], "quant") == 0 && 
		   strcmp(argv[2], "-i") == 0)
		quantility(argc, argv);
	else
		print_usage();

	fflush(stderr);
	fflush(stdout);
	pthread_exit(0);
	return 0;
}