#include <stdio.h>
#include "align.h"
#include "em.h"
#include "bow_parser.h"

void print_help()
{
        fprintf(stderr, "Usage: Hera-EM alignment.bam n_thread\n");
}

const struct result run_em(const struct ref_info_t *ref, struct align_batch_t *a, int thread_cnt)
{
        struct model m;

	uint32_t total_cnt = 0;
	for (int i = 0; i < thread_cnt; i++)
		total_cnt += a[i].read_cnt;

	model_uniform_init(&m, ref, total_cnt, 1);

	squarem(&m, a, thread_cnt, LOOP_BEFORE_OFF);

	//calculate mean fragment length and deviation
	int i;
	double mean = 0, var = 0;
	for (i = m.f_dis.min_l; i <= m.f_dis.max_l; i++) {
		mean += i * m.f_dis.pdf[i];
	}

	for (i = m.f_dis.min_l; i <= m.f_dis.max_l; i++) {
		var += pow(i - mean, 2) * m.f_dis.pdf[i];
	}

	float temp = (float)mean, sd = sqrt(var);

	__VERBOSE("Mean fragment length %lf (+= %lf sd)\n", temp, sd);

	const struct result r = model_to_result(&m);

	model_release(&m);

        return r;
}

int main(int argc, char **argv)
{
        if (argc != 3) {
                print_help();
                return 0;
        }

        //Read transcript
        struct ref_info_t ref;

        int thread_cnt = atoi(argv[2]);
        struct align_batch_t *a = parse_bam(argv[1], thread_cnt, &ref);

	double start = realtime();

        const struct result r = run_em(&ref, a, thread_cnt);

	double end = realtime();

	fprintf(stderr, "EM running time (real time): %lf s\n", end - start);

        write_abundance("result", &r);
	result_release(&r);

        int i;
	free(a);
        release_ref_info(&ref); 

        return 0;
}