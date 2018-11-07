#ifndef _TIMING_
#define _TIMING_
#include <time.h>

#define start_timing() clock_t __start_time__ = clock();
#define restart_timing() __start_time__ = clock();


#define end_timing() ((double) (clock() - __start_time__) / CLOCKS_PER_SEC)

#endif