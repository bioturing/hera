#include "em.h"
#include "timing.h"
#include "utility.h"
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <pthread.h>

#define FULL_STATE 1
#define SIMP_STATE 2

#define align_cnt(a) ((a)->size)
#define noise_cnt(a) ((a)->noise)
#define get_align(a, id) ((a)->align[id])

struct global_thread_param
{
	void *m;
	void *master_m;

	struct model *last_full_m;

	int thread_running;
	int thread_done;

	int running_full;
	int convert_simp;
	int running_simp;
	int state;
	int filtering;
} m_param;

#define set_em_param(m1, m2, f) (m_param.m = (m1),                      \
				 m_param.master_m = (m2),               \
				 m_param.filtering = (f))

#define set_align_list_t(a) (m_param.align = (a))

#define running_full() (m_param.running_full =                          \
			m_param.running_simp = m_param.convert_simp = 1)
#define running_simp() (m_param.running_full = 0, m_param.running_simp = 1)

#define only_simp()    (m_param.running_full = m_param.convert_simp = 0,\
			m_param.running_simp = 1)
#define stop_em()      (m_param.running_full = m_param.convert_simp =   \
			m_param.running_simp = 0)

#define is_converge()                                                   \
	((m_param.state == FULL_STATE &&                                \
	  model_converge((struct model*)m_param.m,                      \
			 (struct model*)m_param.master_m)) ||           \
	 (m_param.state == SIMP_STATE &&                                \
	s_model_converge((struct s_model*)m_param.m,                    \
			 (struct s_model*)m_param.master_m)))

static pthread_mutex_t join_mutex       = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t wait_mutex       = PTHREAD_MUTEX_INITIALIZER;

static pthread_cond_t finish_cond       = PTHREAD_COND_INITIALIZER;
static pthread_cond_t thread_done_cond  = PTHREAD_COND_INITIALIZER;

static void unlock_synchronize() {
	m_param.thread_done = 0;
	pthread_cond_broadcast(&finish_cond);
}

static void thread_done() {
	pthread_mutex_lock(&wait_mutex);
	m_param.thread_done++;
	pthread_cond_signal(&thread_done_cond);
	pthread_cond_wait(&finish_cond, &wait_mutex);
	pthread_mutex_unlock(&wait_mutex);
}

#define e_step(name, type, init,                                        \
		j_prob, n_prob, update, n_update,                       \
		join, release, filtering, use_noise)                    \
static inline void e_step_## name (struct align_batch_t *a,		\
				   double *prob,			\
				   type *m, type *new_m)                \
{                                                                       \
	init(new_m, m);                                                 \
									\
	uint64_t shorten = 0;						\
									\
	/* update nosie read */                                         \
	n_update(new_m, NULL, noise_cnt(a));                            \
									\
	uint64_t from, to = 0;						\
	struct s_read_t *read;						\
									\
	uint64_t i;							\
	while (to < align_cnt(a)) {					\
		if (filtering) {					\
			get_align(a, shorten) = get_align(a, to);	\
			read = &get_align(a, shorten).read;		\
			++shorten;					\
		} else {						\
			read = &get_align(a, to).read;			\
		}							\
									\
		from = to + 1;						\
		to += get_align(a, to).read.cnt;			\
									\
		double total = read->prob;				\
									\
		for (i = from; i < to; i++)				\
			total += prob[i - from]				\
			       = j_prob(m, read, &get_align(a,i));	\
									\
		double noise = 0;                                       \
		if (use_noise) { 					\
			total += noise = n_prob(m, read);               \
			n_update(new_m, read, noise / total); 		\
		}							\
									\
		if (total < EPSILON) {					\
			noise_cnt(a) += filtering;			\
			continue;					\
		}							\
									\
		uint32_t cnt = 0;					\
		for (i = from; i < to; i++) {				\
			int mode = update(new_m, read, &get_align(a, i),\
					  prob[i - from] / total);	\
			if (filtering) {				\
				cnt += (mode == 0);			\
				if (mode == 0)				\
					get_align(a, shorten++)		\
						= get_align(a, i);	\
				else if (mode == 1)			\
					read->prob += prob[i - from];	\
			}						\
		}							\
									\
		if (filtering) {					\
			if (cnt == 0) {					\
				--shorten;				\
				noise_cnt(a) += noise / total;		\
			} else {					\
				read->cnt = cnt + 1;			\
			}						\
		}							\
	}								\
									\
	if (filtering) {						\
		align_cnt(a) = shorten;					\
		batch_resize(a);					\
	}								\
									\
	pthread_mutex_lock(&join_mutex);				\
	join(m_param.master_m, new_m);					\
	pthread_mutex_unlock(&join_mutex);				\
	release(new_m);							\
									\
	thread_done();							\
}

e_step(full_no_filter,  struct model, model_init,
			joint_prob, noise_prob,
			model_update, model_update_noise,
			model_join, model_release, 0, USE_NOISE);

e_step(full_filter,     struct model, model_init,
			joint_prob, noise_prob,
			model_update, model_update_noise,
			model_join, model_release, 1, USE_NOISE);

e_step(nf_no_filter,    struct model, model_init,
			joint_no_frag_prob, noise_prob,
			model_no_frag_update, model_update_noise,
			model_join, model_release, 0, USE_NOISE);

e_step(nf_filter,       struct model, model_init,
			joint_no_frag_prob, noise_prob,
			model_no_frag_update, model_update_noise,
			model_join, model_release, 1, USE_NOISE);

e_step(simp_non_filter, struct s_model, s_model_init,
			s_joint_prob, s_noise_prob,
			s_model_update, s_model_update_noise,
			s_model_join, s_model_release, 0, USE_NOISE);

e_step(simp_filter,     struct s_model, s_model_init,
			s_joint_prob, s_noise_prob,
			s_model_update, s_model_update_noise,
			s_model_join, s_model_release, 1, USE_NOISE);

#define to_s_algn(name, cond_prob)                                      \
static void to_s_align_## name (struct align_batch_t *a,		\
				   const struct model *m)		\
{									\
	uint64_t i;							\
	uint64_t from, to = 0;						\
									\
	while (to < align_cnt(a)) {					\
		struct s_read_t *read = &get_align(a, to).read;		\
									\
		from = to + 1;						\
		to += read->cnt;					\
									\
		for (i = from; i < to; i++)				\
			get_align(a, i).simp.c_prob			\
				= cond_prob(m, read, &get_align(a, i));	\
	}								\
}

to_s_algn(full, cond_prob);
to_s_algn(no_frag, cond_no_frag_prob);

static inline void run_full_step(struct align_batch_t *a, double *prob, struct model *m, struct model *new_m)
{
	if (m->use_frag) {
		if (m_param.filtering)
			e_step_full_filter(a, prob, m, new_m);
		else
			e_step_full_no_filter(a, prob, m, new_m);
	} else {
		if (m_param.filtering)
			e_step_nf_filter(a, prob, m, new_m);
		else
			e_step_nf_no_filter(a, prob, m, new_m);
	}
}

static void *thread_func(void *parameter)
{
	struct align_batch_t a = *(struct align_batch_t*)parameter;
	thread_done();

	double *prob = malloc(sizeof(double) * MAX_ALIGNMENT);
	void *new_m = malloc(sizeof(struct model));

	while (m_param.running_full)
		run_full_step(&a, prob, m_param.m, new_m);

	if (m_param.convert_simp) {
		if (m_param.last_full_m->use_frag)
			to_s_align_full(&a, m_param.last_full_m);
		else
			to_s_align_no_frag(&a, m_param.last_full_m);
	}

	new_m = realloc(new_m, sizeof(struct s_model));

	while (m_param.running_simp)
		if (m_param.filtering)
			e_step_simp_filter(&a, prob,
					(struct s_model*) m_param.m,
					(struct s_model*) new_m);
		else
			e_step_simp_non_filter(&a, prob,
					(struct s_model*) m_param.m,
					(struct s_model*) new_m);

	free(prob);
	free(new_m);

	align_batch_release(&a);

	return NULL;
}

static void m_finish(void *m)
{
	if (m_param.running_full) {
		model_finish((struct model*)m);
		m_param.state = FULL_STATE;
	} else if (m_param.running_simp) {
		s_model_finish((struct s_model*)m);
		m_param.state = SIMP_STATE;
	}/* else {
		assert(0);
	}*/
}

static void m_init()
{
	if (m_param.running_full)
		model_master_init((struct model*)m_param.master_m,
				  (struct model*)m_param.m);
	else if (m_param.running_simp)
		s_model_master_init((struct s_model*)m_param.master_m,
				    (struct s_model*)m_param.m);
}

static void m_release(void *m)
{
	if (m_param.state == FULL_STATE)
		model_release((struct model*)m);
	else if (m_param.state == SIMP_STATE)
		s_model_release((struct s_model*)m);
}

static void m_accelerate(void *m1, void *m2, void *m3)
{
	if (m_param.state == FULL_STATE)
		model_accelerate((struct model*)m1,
				 (struct model*)m2,
				 (struct model*)m3);
	else if (m_param.state == SIMP_STATE)
		s_model_accelerate((struct s_model*)m1,
				   (struct s_model*)m2,
				   (struct s_model*)m3);
}

static void wait_syncrhonize() {
	pthread_mutex_lock(&wait_mutex);
	while (m_param.thread_done < m_param.thread_running)
		pthread_cond_wait(&thread_done_cond, &wait_mutex);
	pthread_mutex_unlock(&wait_mutex);
}

static inline void em()
{
	/* Init */
	m_init();

	/* Start E step */
	unlock_synchronize();

	/* Wait E step done */
	wait_syncrhonize();

	/* M step */
	m_finish(m_param.master_m);
}

static inline void a_em(int *times)
{
	++*times;
	em();

	/* analysis */
#ifdef VERBOSE
#ifdef ANALYSIS
	m->diff(m->m, m->new_m);
#endif
#endif

#ifdef PRINT_ROUND
	fprintf(stderr, "[EM] done %d round\r", *times);
#endif

	/* save result */
#ifdef PRINT_PROGRESS
#ifdef SAVE_MODEL
	model_save(new, "model");
#endif
#endif
}

pthread_t *em_thread;

void em_thread_init(struct align_batch_t *es)
{
	em_thread = malloc(sizeof(pthread_t) * m_param.thread_running);

	int i;
	for (i = 0; i < m_param.thread_running; i++)
		pthread_create(&em_thread[i], NULL, &thread_func, &es[i]);

	wait_syncrhonize();
}

void em_thread_stop()
{
	stop_em();
	unlock_synchronize();

	int i;
	for (i = 0; i < m_param.thread_running; i++)
		pthread_join(em_thread[i], NULL);

	free(em_thread);
}

//return save model to m1
int squarem_loop(void *m1, void *m2, void *m3, int *times)
{
	set_em_param(m1, m2, 0);
	a_em(times);

	set_em_param(m2, m3, 1);
	a_em(times);

	#ifdef ACCELERATE
	#ifdef VERBOSE
	#endif
	m_accelerate(m1, m2, m3);
	#ifdef VERBOSE
	#endif
	#endif

	m_release(m1);
	m_release(m2);

	set_em_param(m3, m1, 0);
	a_em(times);

	int converged = is_converge();
	m_release(m3);

	return converged;
}

//return is converged
int f_squarem(struct model *m1, int *times, int m_loop)
{
	struct model m2, m3;
	running_full();
	while(m_loop && !squarem_loop(m1, &m2, &m3, times))
		--m_loop;

	return m_loop != 0;
};

int s_squarem(struct s_model *m1, int *times, int m_loop)
{
	struct s_model m2, m3;
	running_simp();

	while(m_loop && !squarem_loop(m1, &m2, &m3, times))
		--m_loop;
	return m_loop!= 0;
};

void squarem(struct model *m0, struct align_batch_t *a, int thread_cnt, int full_loop)
{
	start_timing();

	m_param.thread_running = thread_cnt;
	model_global_init();

	int times = 0;
	em_thread_init(a);

	fprintf(stderr, CYAN "[EM] start running\n" RESET);

	if (!f_squarem(m0, &times, full_loop)) {
		fprintf(stderr, GREEN "[EM] running simplified model\n" RESET);

		/*if not converged do simplified model em */
		struct s_model sm0;
		m_param.last_full_m = m0;
		char *locking = copy_to_s_model(&sm0, m0);

		s_squarem(&sm0, &times, -1);

		merge_to_f_model(m0, &sm0);
		m_release(&sm0);
		free(locking - 1);
	}
	em_thread_stop();
}
