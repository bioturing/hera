#include "argument.h"

static struct arg_quant_data *args_quant;
static struct arg_index_data *args_index;

#define __arg_error(fmt, args...) {						\
	fprintf(stderr, "ERROR: " fmt "\n"					\
		        "Using -h to view usage\n", ##args);			\
	exit(EXIT_FAILURE);							\
}

static void check_num(int32_t argc, int32_t pos, char *argv, char *s)
{
	if (argv[0] == '-' || pos == argc)
		__arg_error("Missing data for argument %s", s);

	int32_t i;
	for (i = 0; argv[i] != '\0'; ++i)
		if (argv[i] < '0' || argv[i] > '9')
			__arg_error("Invalid data for argument %s: %s",
				    s, argv);	
}

static void check_str(int32_t argc, int32_t pos, char *argv, char *s)
{
	if (argv[0] == '-' || pos == argc)
		__arg_error("Missing data for argument %s", s);
}

static int32_t get_left_file(int32_t pos, int32_t argc, char *argv[])
{
	int32_t i = pos;

	while (i < argc) {
		if (argv[i][0] != '-')
			++args_quant->nLeft;
		else
			break;
		++i;
	}
	
	if (!args_quant->nLeft)
		__arg_error("Missing data for argument -1");

	args_quant->left_file = malloc(args_quant->nLeft * sizeof(char*));
	for (i = 0; i < args_quant->nLeft; ++i)
		args_quant->left_file[i] = argv[i + pos];

	return args_quant->nLeft + pos;
}

static int32_t get_right_file(int32_t pos, int32_t argc, char *argv[])
{
	int32_t i = pos;

	while (i < argc) {
		if (argv[i][0] != '-')
			++args_quant->nRight;
		else
			break;
		++i;
	}

	if (!args_quant->nRight)
		__arg_error("Missing data for argument -2");

	args_quant->right_file = malloc(args_quant->nRight * sizeof(char*));
	for (i = 0; i < args_quant->nRight; ++i)
		args_quant->right_file[i] = argv[i + pos];

	return args_quant->nRight + pos;
}

static void init_quant_argument()
{
	args_quant->nLeft = args_quant->nRight = 0;
	args_quant->left_file = args_quant->right_file = NULL;
	args_quant->threads = 1;
	args_quant->bootstrap = 0;
	args_quant->idx_dir = NULL;
	args_quant->out_dir = "./";
	args_quant->genome = NULL;
	args_quant->bam = false;
	args_quant->compress = -1;
	args_quant->prefix = "\0";
	args_quant->verbose = false;
}

static void check_valid_quant_argument()
{
	if (args_quant->nRight && args_quant->nLeft != args_quant->nRight)
		__arg_error("Number of paired-end reads are not equal");

	if (args_quant->threads < 1)
		__arg_error("Invalid number of threads: %d", 
			  args_quant->threads);

	if (args_quant->bootstrap < 0)
		__arg_error("Invalid number of bootstrap: %d",
			  args_quant->bootstrap);

	if (args_quant->compress == 0 || args_quant->compress < -1 ||
	    args_quant->compress > 9)
		__arg_error("Invalid compress level: %d", 
			    args_quant->bootstrap);

	if (args_quant->idx_dir == NULL)
		__arg_error("Missing argument -i");

	if (args_quant->nLeft == 0)
		__arg_error("Input file is not defined");
}

struct arg_quant_data *arg_get_quant(int32_t pos, int32_t argc, char *argv[])
{
	if (pos == argc)
		quant_usage();

	int32_t i = pos;
	args_quant = malloc(sizeof(struct arg_quant_data));
	init_quant_argument();

	while (i < argc) {
		if (!strcmp(argv[i], "-1")) {
			i = get_left_file(i + 1, argc, argv);
		} else if (!strcmp(argv[i], "-2")) {
			i = get_right_file(i + 1, argc, argv);
		} else if (!strcmp(argv[i], "-i")) {
			++i;
			check_str(argc, i, argv[i], "-i");
			args_quant->idx_dir = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-o")) {
			++i;
			check_str(argc, i, argv[i], "-o");
			args_quant->out_dir = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-t")) {
			++i;
			check_num(argc, i, argv[i], "-t");
			args_quant->threads = atoi(argv[i]);
			++i;
		} else if (!strcmp(argv[i], "-b")) {
			++i;
			check_num(argc, i, argv[i], "-b");
			args_quant->bootstrap = atoi(argv[i]);
			++i;
		} else if (!strcmp(argv[i], "-f")) {
			++i;
			check_str(argc, i, argv[i], "-f");
			args_quant->genome = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-p")) {
			++i;
			check_str(argc, i, argv[i], "-p");
			args_quant->prefix = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-w")) {
			args_quant->bam = true;
			++i;
		} else if (!strcmp(argv[i], "-v")) {
			args_quant->verbose = true;
			++i;
		} else if (!strcmp(argv[i], "-z")) {
			++i;
			check_num(argc, i, argv[i], "-z");
			args_quant->compress = atoi(argv[i]);
			++i;
		} else if (!strcmp(argv[i], "-h")) {
			quant_usage();
		}  else {
			__arg_error("Invalid argument %s", argv[i]);
		}
	}

	check_valid_quant_argument();
	return args_quant;
}

static void init_index_argument()
{
	args_index->transcript = NULL;
	args_index->genome = NULL;
	args_index->out_dir = NULL;
}

static void check_valid_index_argument()
{
	if (!args_index->transcript)
		__arg_error("Missing argument -t");

	if (!args_index->out_dir)
		__arg_error("Missing argument -o");
}

struct arg_index_data *arg_get_index(int32_t pos, int32_t argc, char *argv[])
{
	if (pos == argc)
		index_usage();

	int32_t i = pos;
	args_index = malloc(sizeof(struct arg_index_data));
	init_index_argument();

	while (i < argc) {
		if (!strcmp(argv[i], "-t")) {
			++i;
			check_str(argc, i, argv[i], "-t");
			args_index->transcript = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-o")) {
			++i;
			check_str(argc, i, argv[i], "-o");
			args_index->out_dir = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-g")) {
			++i;
			check_str(argc, i, argv[i], "-g");
			args_index->genome = argv[i];
			++i;
		} else if (!strcmp(argv[i], "-h")) {
			index_usage();
		} else {
			__arg_error("Invalid argument %s", argv[i]);
		}
	}

	check_valid_index_argument();
	return args_index;
}