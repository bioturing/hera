#include "EM.h"
#include "bam_write.h"

void print_usage()
{
	fprintf(stdout, "VERSION:\n\t%s\n", VERSION);
	fprintf(stdout, "QUANTILITY:\n\t ./hera quant -i path/to/index_directory [options] <R1.fq> <R2.fq>\n");
	fprintf(stdout, "\t\t-o:\t Output directory (default: ./)\n");
	fprintf(stdout, "\t\t-t:\t Number of threads (default: 1)\n");
	fprintf(stdout, "\t\t-z:\t Compress level (1 - 9) (default: -1)\n");
	fprintf(stdout, "\t\t-b:\t Number of bootstraps (default: 0)\n");
	fprintf(stdout, "\t\t-w:\t Output bam file 0:true, 1: false (default: 0)\n");
	fprintf(stdout, "\t\t-f:\t Genome fasta file (if not define, genome mapping will be ignore\n");
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
	DEBUG_PRINT("Index transcripts in: %s\n", transcriptome);
	get_ref_seq(transcriptome);

	// Write to file
	DEBUG_PRINT("Number of sequences\t: %u\n", REF_INF->n);
	write_to_file(out_file);

	// Index genome 
	if (full == 1){
		DEBUG_PRINT("Index genome in: %s\n", genome);
		indexGenome(genome, out_file);
	}
}

void quantility(int argc, char *argv[]) {
	char *idx_dir, *out_dir = "./", *genome = NULL;
	EM_val *em_val;
	int c, temp_nthread, temp_writebam, temp_nbs;
	unsigned int n_bs = 0;
	unsigned long i_start = 0;

	// Get parameter
	while ((c = getopt(argc - 1, argv + 1, "t:i:o:b:z:h:w:f:")) != -1) {
		switch (c) {
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
	init_output(idx_dir, out_dir, argc, argv);

	//Get alignment
	if (optind + 2 < argc) {
		DEBUG_PRINT("Process paired-end reads in:\n\t 1. %s\n\t 2. %s\n",
					     argv[optind + 1], argv[optind + 2]);
		get_alignment_pair(argv[optind + 1], argv[optind + 2]);
	} else {
		DEBUG_PRINT("Process single-end reads in: %s\n",
							       argv[optind + 1]);
		get_alignment(argv[optind + 1]);
	}

	// Flush and end file
	if (WRITE_BAM == 0)
		bam_close(BAM);

	// Correct estimate count
	em_val = estimate_count(REF_INF->count, 0);

	// Finish
	write_result(out_dir, n_bs, em_val);
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
