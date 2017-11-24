#include "read.h"

void free_read_inf(struct read_inf *read)
{
	free(read->seq);
	free(read->qual);
	free(read->name);
	free(read);
}

enum read_exit_code get_read_from_fq_buffer(struct read_inf *read,
					    char *buffer,
					    int32_t *pos)
{
	int32_t i = *pos;
	int32_t prev;

	// name part
	prev = i;
	if (buffer[i] != '@')
		return READ_FAIL;
	for (; buffer[i] != '\0' && buffer[i] != '\n'; ++i);
	if (buffer[i] == '\0')
		return READ_FAIL;
	read->name = buffer + prev + 1; // skip @ character
	buffer[i] = '\0';
	++i;

	// sequence part
	prev = i;
	for (; buffer[i] != '\0' && buffer[i] != '\n'; ++i);
	if (buffer[i] == '\0')
		return READ_FAIL;
	read->seq = buffer + prev;
	read->len = i - prev;
	buffer[i] = '\0';
	++i;

	if (read->len == 0)
		return READ_FAIL;

	// optionally part
	if (buffer[i] != '+')
		return READ_FAIL;
	for (; buffer[i] != '\0' && buffer[i] != '\n'; ++i);
	if (buffer[i] == '\0')
		return READ_FAIL;
	++i;

	// quality part
	prev = i;
	for (; buffer[i] != '\0' && buffer[i] != '\n'; ++i);
	if (i - prev != read->len)
		return READ_FAIL;
	read->qual = buffer + prev;
	if (buffer[i] == '\0')
		return READ_END;
	buffer[i] = '\0';
	if (buffer[i + 1] == '\0')
		return READ_END;
	*pos = ++i;
	return READ_SUCCESS;
}

enum read_exit_code get_read_from_fa_buffer(struct read_inf *read,
					    char *buffer,
					    int32_t *pos)
{
	int32_t i = *pos;
	int32_t prev, pos_cpy;
	read->qual = NULL;
	read->len = 0;

	// name part
	prev = i;
	if (buffer[i] != '>')
		return READ_FAIL;
	for (; buffer[i] != '\0' && buffer[i] != '\n'; ++i);
	if (buffer[i] == '\0')
		return READ_FAIL;
	read->name = buffer + prev + 1; // skip > character
	buffer[i] = '\0';
	++i;

	// sequence part
	pos_cpy = i;
	read->seq = buffer + i;

	while (1) {
		prev = i;
		if (buffer[i] == '>' || buffer[i] == '\0') {
			if (read->len == 0)
				return READ_FAIL;
			*pos = i;
			buffer[pos_cpy] = '\0';
			return READ_SUCCESS;
		}
		for (; buffer[i] != '\0' && buffer[i] != '\n'; ++i);
		memcpy(buffer + pos_cpy, buffer + prev, i - prev);
		pos_cpy += i - prev;
		read->len += i - prev;

		if (buffer[i] == '\0' || buffer[i + 1] == '\0') {
			buffer[pos_cpy] = '\0';
			return READ_END;
		}
		++i;
	}
}