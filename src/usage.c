#include "usage.h"

void quant_usage()
{
	fprintf(stderr, "Hera is a program developed by BioTuring for RNA-Seq analysis\n");
	fprintf(stderr, "Please contact info@bioturing.com if you need further support\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Hera version: %s\n", HERA_VERSION);
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: ./hera quant [arguments]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Required arguments:\n");
	fprintf(stderr, "  -1 <read-files>    Input left read-files, separated by space\n");
	fprintf(stderr, "  -2 <read-files>    Input right read-files, separated by space\n");
	fprintf(stderr, "                     (using -1 only if quantify for single-end reads)\n");
	fprintf(stderr, "  -i <STRING>        path to hera index directory\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "  -o <STRING>        Output directory (default: ./)\n");
	fprintf(stderr, "  -t <INT>           Number of threads (default: 1)\n");
	fprintf(stderr, "  -b <INT>           Number of bootstrap samples (default: 0)\n");
	fprintf(stderr, "  -f <genome-file>   Genome mapping, need full index to use this option\n");
	fprintf(stderr, "                     (if not define, genome mapping will be ignore)\n");
	fprintf(stderr, "  -p <STRING>        Output prefix (default: '')\n");
	fprintf(stderr, "  -w                 Output bam file\n");
	fprintf(stderr, "  -z <INT>           Bam compress level (1 - 9) (default: -1)\n");
	fprintf(stderr, "  -v                 Verbose mode\n");
	fprintf(stderr, "  -h                 Print help\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Example:\n");
	fprintf(stderr, "  ./hera quant -i hera_index/ -w -t 32 -b 100 -o hera_output/\n");
	fprintf(stderr, "    -1 read_1.fq.gz -2 read_2.fq.gz\n");
	fprintf(stderr, "  (Output bam file, 32 threads, 100 bootstrap samples, paired-end mode)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  ./hera quant -i hera_index/ -t 32 -f GRCh37_75_homo_sapiens.fa\n");
	fprintf(stderr, "    -o hera_output/ -1 read_lane_1.fq.gz read_lane_2.fq.gz\n");
	fprintf(stderr, "  (No bam, 32 threads, genome mapping, single-end mode with multiple file)\n");
	
	exit(EXIT_FAILURE);
}

void index_usage()
{
	fprintf(stderr, "Hera is a program developed by BioTuring for RNA-Seq analysis\n");
	fprintf(stderr, "Please contact info@bioturing.com if you need further support\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Hera version: %s\n", HERA_VERSION);
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: ./hera index [arguments]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Required arguments:\n");
	fprintf(stderr, "  -t <STRING>        Transcriptome file to be index\n");
	fprintf(stderr, "  -o <STRING>        Output file name\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "  -g <STRING>        Genome file for full index\n");
	fprintf(stderr, "  -h                 Print help\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Example:\n");
	fprintf(stderr, "  ./hera index -t transcript.fa -o index -g Homo_sapiens.fa\n");
	
	exit(EXIT_FAILURE);
}