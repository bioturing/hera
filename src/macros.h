#ifndef MARCOS_H_
#define MARCOS_H_

#define _verbose(fmt, args...) fprintf(stderr, fmt, ##args)

#define _error(fmt, args...) {							\
	fprintf(stderr, "ERROR: " fmt, ##args);					\
	exit(EXIT_FAILURE);							\
}

#define _warning(fmt, args...) fprintf(stderr, "WARNING: " fmt, ##args)

#define _abs(x) ((x) < 0 ? -(x) : (x))

#define _min(a, b) ((a) < (b) ? (a) : (b))

#define _min3(a, b, c) _min(_min(a, b), c)

#define _max(a, b) ((a) > (b) ? (a) : (b))

#define _max3(a, b, c) _max(_max(a, b), c)

#endif /* __MARCOS_H */