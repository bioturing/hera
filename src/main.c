#include "EM.h"
#include "bam_write.h"
#include "argument.h"
#include "log.h"

void directory_exist(const char* sFilePath)
{
	struct stat oBuffer;

	if (stat(sFilePath, &oBuffer) == 0)
		if (oBuffer.st_mode & S_IFDIR)
			return;

	if (mkdir(sFilePath, 0755) == -1)
		_error("Cannot create output directory: %s", sFilePath);
}

void process_read_pair(struct arg_quant_data *args)
{
	int32_t i;

	for (i = 0; i < args->nLeft; ++i) {
		if (i)
			_verbose("\n");
		get_alignment_pair(args->left_file[i], args->right_file[i],
				   args->threads);
	}

	MEAN_LEN /= 2 * N_READ;
	MEAN_FRAG = MEAN_FRAG / PAIRED + MEAN_LEN;

	for (i = 0; i < REF_INF->n; ++i) {
		if (REF_INF->len[i] < MEAN_FRAG)
			REF_INF->eff_len[i] = 0;
		else
			REF_INF->eff_len[i] = REF_INF->len[i] - MEAN_FRAG + 1;
	}
	
	_verbose("\n");
	_verbose("Mapping result:\n");
	_verbose("    Number of aligned read pairs   : %u/%u\n", MAPPED, N_READ);
	_verbose("    Mean read length               : %f\n", MEAN_LEN);
	_verbose("    Mean fragment length           : %f\n", MEAN_FRAG);

	_log("MAPPING RESULT:\n");
	_log("  Total number of read pairs            : %u\n", N_READ);
	_log("      - Number of proper mapped pairs   : %u (%f%%)\n",
	      PAIRED, (float) PAIRED/N_READ);
	_log("      - Number of unproper mapped pairs : %u (%f%%)\n",
	      MAPPED -  PAIRED, (float)(MAPPED -  PAIRED)/N_READ);
	_log("      - Number of unmapped pairs        : %u (%f%%)\n",
	      N_READ - MAPPED, (float) (N_READ - MAPPED)/N_READ);
	_log("  Mean read length                      : %f\n", MEAN_LEN);
	_log("  Mean fragment length                  : %f\n", MEAN_FRAG);
	_log("\n");
}

void process_read_single(struct arg_quant_data *args)
{
	int32_t i;

	for (i = 0; i < args->nLeft; ++i) {
		if (i)
			_verbose("\n");
		get_alignment_single(args->left_file[i], args->threads);
	}

	MEAN_LEN /= N_READ;
	MEAN_FRAG = MEAN_LEN;
	for (i = 0; i < REF_INF->n; ++i) {
		if (REF_INF->len[i] < MEAN_FRAG)
			REF_INF->eff_len[i] = 0;
		else
			REF_INF->eff_len[i] = REF_INF->len[i] - MEAN_FRAG + 1;
	}

	_verbose("\n");
	_verbose("Mapping result:\n");
	_verbose("    Number of aligned reads        : %u/%u\n", PAIRED, N_READ);
	_verbose("    Mean read lenght               : %f\n", MEAN_LEN);

	_log("MAPPING RESULT:\n");
	_log("  Total number of reads                 : %u\n", N_READ);
	_log("      - Number of mapped reads          : %u (%f%%)\n",
	      PAIRED, (float) PAIRED/N_READ);
	_log("      - Number of unmapped reads        : %u (%f%%)\n",
	      N_READ - PAIRED, (float) (N_READ - PAIRED)/N_READ);
	_log("  Mean read length                      : %f\n", MEAN_LEN);
	_log("\n");
}

void show_quant_args_info(struct arg_quant_data *args)
{
	_log("************* Used parameters ***************\n");
	_log("No. of threads             : %d\n", args->threads);
	_log("No. of bootstrap samples   : %d\n", args->bootstrap);
	_log("Output directory           : %s\n", args->out_dir);
	
	_log("Prefix                     : ", args->prefix);
	if (strlen(args->prefix))
		_log("%s\n", args->prefix);
	else
		_log("None\n");
	
	_log("Genome mapping             : ");
	if (args->genome)
		_log("%s\n", args->genome);
	else
		_log("None\n");
	
	_log("Output bam file            : ");
	if (args->bam) {
		_log("Yes\n");
		if (args->compress == -1)
			_log("Bam compress level         : default\n");
		else
			_log("Bam compress level         : %d\n", 
				  args->compress);	
	} else {
		_log("No\n");
	}
	
	if (args->nRight) {
		_log("Quantify mode              : paired-end\n");
		_log("No. of input pair file     : %d\n", args->nRight);
	} else {
		_log("Quantify mode              : single-end\n");
		_log("No. of input file          : %d\n", args->nLeft);
	}
	_log("\n");

	if (!args->verbose) return;

	_verbose("************* Used parameters ***************\n");
	_verbose("No. of threads             : %d\n", args->threads);
	_verbose("No. of bootstrap samples   : %d\n", args->bootstrap);
	_verbose("Output directory           : %s\n", args->out_dir);
	
	_verbose("Prefix                     : ", args->prefix);
	if (strlen(args->prefix))
		_verbose("%s\n", args->prefix);
	else
		_verbose("None\n");
	
	_verbose("Genome mapping             : ");
	if (args->genome)
		_verbose("%s\n", args->genome);
	else
		_verbose("None\n");
	
	_verbose("Output bam file            : ");
	if (args->bam) {
		_verbose("Yes\n");
		if (args->compress == -1)
			_verbose("Bam compress level         : default\n");
		else
			_verbose("Bam compress level         : %d\n", 
				  args->compress);	
	} else {
		_verbose("No\n");
	}
	
	if (args->nRight) {
		_verbose("Quantify mode              : paired-end\n");
		_verbose("No. of input pair files    : %d\n", args->nRight);
	} else {
		_verbose("Quantify mode              : single-end\n");
		_verbose("No. of input files         : %d\n", args->nLeft);
	}
	_verbose("\n");
}

void init_output(char *out_dir, char *prefix)
{
	int32_t len, p_len, pos;
	char *file_path;

	// Make file path
	len = strlen(out_dir);
	p_len = strlen(prefix);
	pos = 0;

	file_path = calloc(len + p_len + 17, sizeof(char));
	memcpy(file_path + pos, out_dir, len);
	pos += len;
	memcpy(file_path + pos, "/", 1);
	pos += 1;
	memcpy(file_path + pos, prefix, p_len);
	pos += p_len;

	if (p_len > 0) {
		memcpy(file_path + pos, "_", 1);
		pos += 1;
	}

	if (WRITE_BAM == 0) {
		memcpy(file_path + pos, "alignment.bam\0", 14);
		BAM = bam_open(file_path, "wb");
	}

	memcpy(file_path + pos, "hera.log\0", 9);
	init_log_file(file_path);

	free(file_path);
}

void quantility(int32_t pos, int32_t argc, char *argv[])
{
	int32_t i;

	// Get parameter
	struct arg_quant_data *args = arg_get_quant(pos, argc, argv);

	// FIXME: remove extern global
	COMPRESS_LEVEL = args->compress;
	if (args->bam)
		WRITE_BAM = 0;
	else
		WRITE_BAM = 1;

	directory_exist(args->out_dir);
	init_output(args->out_dir, args->prefix);

	_log("VERSION: %s\n\n", HERA_VERSION);
	_log("COMMAND:\n  ");
	for (i = 0; i < argc; ++i)
		_log("%s ", argv[i]);
	_log("\n\n");

	show_quant_args_info(args);
	
	// Retrieve index
	_log("************* Process Info ***************\n");
	_verbose("Loading index ...\n");
	_log("Loading index ...\n");
	init_refInf();
	init_hash();
	load_index(args->idx_dir, args->genome);	
	init_bam_header(args->idx_dir, argc, argv);

	//Get alignment
	if (args->nRight)
		process_read_pair(args);
	else
		process_read_single(args);

	// FIXME: why close in this function?
	if (args->bam) 
		bam_close(BAM);

	_verbose("Quantifying the abundances ...\n");
	_log("Quantifying the abundances ...\n");
	// Correct estimate count
	EM_val *em_val = estimate_count(REF_INF->count, 0);

	// Finish
	write_result(args->out_dir, args->bootstrap, em_val,
		     args->prefix, args->threads);

	close_log_file();
}

void index_reference(int32_t pos, int32_t argc, char *argv[])
{
	// Get parameter
	struct arg_index_data *args = arg_get_index(pos, argc, argv);

	// Initialize
	init_refInf();
	init_hash();

	//Read reference sequence + construct hash for transcriptome    
	_verbose("Index transcript in: %s\n", args->transcript);
	get_ref_seq(args->transcript);

	// Write to file
	_verbose("Number of sequences: %u\n", REF_INF->n);
	write_to_file(args->out_dir);

	// Index genome 
	if (args->genome) {
		_verbose("Index genome in: %s\n", args->genome);
		indexGenome(args->genome, args->out_dir);
	}

	_verbose("Finish!\n");
}

int32_t main(int32_t argc, char *argv[])
{
	// For better looking verbose
	_verbose("\n");

	if (argc < 2)
		quant_usage();
	else if (!strcmp(argv[1], "index"))
		index_reference(2, argc, argv);
	else if (!strcmp(argv[1], "quant"))
		quantility(2, argc, argv);
	else
		quant_usage();

	fflush(stderr);
	return 0;
}