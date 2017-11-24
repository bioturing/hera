#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>
#include <stdarg.h>

void init_log_file(char *filename);

void close_log_file();

void _log(const char * format, ...);

#endif /* LOG_H_ */