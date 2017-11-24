#include "log.h"

static FILE *_logFile = NULL;

void init_log_file(char *filename)
{
	if (_logFile) {
		fprintf(stderr, "BUG: Log file already exist\n");
		return;
	}

	_logFile = fopen(filename, "w");
}

void close_log_file()
{
	if (!_logFile) {
		fprintf(stderr, "BUG: Log file does not exist\n");
		return;
	}

	fclose(_logFile);
}

void _log(const char * format, ...)
{
	if (!_logFile) {
		fprintf(stderr, "BUG: Log file does not exist\n");
		return;
	}

	char buffer[256];
	va_list args;
	va_start(args, format);
	vsprintf(buffer, format, args);
	fprintf(_logFile, "%s", buffer);
	va_end(args);
}