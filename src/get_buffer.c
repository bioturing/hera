#include "get_buffer.h"

static int32_t buffer_size = 1048576;

enum gb_file_type get_format(char *file_name)
{
	char tmp_buffer[10];
	gzFile tmp_file;
	int32_t n;
	enum gb_file_type ret;
	
	tmp_buffer[0] = 0;
	tmp_file = gzopen(file_name, "r");
	if (!tmp_file)
		_error("Cannot open file: %s\n", file_name);
	n = gzread(tmp_file, tmp_buffer, 1);
	if (tmp_buffer[0] == '@')
		ret = TYPE_FASTQ;
	else if (tmp_buffer[0] == '>')
		ret = TYPE_FASTA;
	else
		_error("Unsupport format of file: %s\n", file_name);
	gzclose(tmp_file);

	return ret;
}

struct gb_pair_data *gb_init_pair(char *file_name1, char *file_name2)
{
	struct gb_pair_data *data = malloc(sizeof(struct gb_pair_data));
	if (!strcmp(file_name1, file_name2))
		_warning("Two identical files\n");
	data->format = get_format(file_name1);
	if (get_format(file_name2) != data->format)
		_error("Format in two files are not equal\n");
	data->left = malloc(sizeof(struct gb_file_inf));
	data->right = malloc(sizeof(struct gb_file_inf));
	data->left->file_name = file_name1;
	data->right->file_name = file_name2;
	data->left->file = gzopen(file_name1, "r");
	data->right->file = gzopen(file_name2, "r");
	data->left->buffer = malloc(buffer_size + 1);
	if (!data->left->buffer)
		_error("Not enough memory for allocate");
	data->right->buffer = malloc(buffer_size + 1);
	if (!data->right->buffer)
		_error("Not enough memory for allocate");
	data->left->end = data->right->end = 0;
	data->left->finish_read = data->right->finish_read = false;
	data->finish_flag = false;
	data->warning_flag = false;
	
	return data;
}

struct gb_single_data *gb_init_single(char *file_name)
{
	struct gb_single_data *data = malloc(sizeof(struct gb_single_data));
	data->format = get_format(file_name);
	data->single = malloc(sizeof(struct gb_file_inf));
	data->single->file_name = file_name;
	data->single->file = gzopen(file_name, "r");
	if (!data->single->file)
		_error("Cannot open file: %s\n", file_name);
	data->single->buffer = malloc(buffer_size + 1);
	data->single->end = 0;
	data->single->finish_read = false;
	data->finish_flag = false;

	return data;
}

void gb_destroy_pair_data(struct gb_pair_data *data)
{
	free(data->left);
	free(data->right);
	free(data);
}

void gb_destroy_single_data(struct gb_single_data *data)
{
	free(data->single);
	free(data);
}

/*
 * return 0 if read not found
 * otherwise return the next read's position
 */
static int32_t get_buffer_read_fq(struct gb_file_inf *data,
				  int32_t prev_i)
{
	int32_t id;
	int32_t i = prev_i;
	char *buffer = data->buffer;
	int32_t end = data->end;

	id = 0;
	while (1) {
		for (; buffer[i] != '\n' && i < end; ++i);
		if (i == end) {
			if (!data->finish_read) 
				break; 
		} 
		else { 
			++i; 
		}
		id = (id + 1) & 3;
		if (!id)
			return i;
		if (i == end)
			break;
	}

	if (end - prev_i >= buffer_size)
		_error("Read too long of file: %s\n", data->file_name);

	return 0;
}

/*
 * return 0 if read not found
 * otherwise return the next read's position
 */
static int32_t get_buffer_read_fa(struct gb_file_inf *data,
				  int32_t prev_i)
{
	int32_t i;
	char *buffer = data->buffer;
	int32_t end = data->end;

	if (prev_i == end)
		return 0;
	else
		i = prev_i + 1;
	
	for (; buffer[i] != '>' && i < end; ++i);
	if (i == end) {
		if (end - prev_i >= buffer_size)
			_error("Read too long of file: %s\n", data->file_name);
		if (!data->finish_read)
			return 0;
		else
			return end;
	}
	else { 
		return i;
	}
}

static void load_buffer(struct gb_file_inf *data)
{
	int32_t n_byte;

	n_byte = gzread(data->file, data->buffer + data->end,
			buffer_size);
	data->end += n_byte;
	if (n_byte < buffer_size)
		data->finish_read = true;
}

static void split_buffer(struct gb_file_inf *data, 
			 char *old_buffer,
			 int32_t prev_i)
{
	int32_t padding;

	padding = data->end - prev_i;
	data->buffer = malloc(padding + buffer_size + 1);
	if (!data->buffer)
		_error("Not enough memory for allocate");
	memcpy(data->buffer, old_buffer + prev_i, padding);
	data->end = padding;
	data->buffer[padding] = 0;
	old_buffer[prev_i] = 0;
}

void gb_destroy_pair_buf(struct gb_pair_buf *data)
{
	free(data->left);
	free(data->right);
	free(data);
}

void gb_destroy_single_buf(struct gb_single_buf *data)
{
	free(data->single);
	free(data);
}

struct gb_pair_buf *gb_get_pair(struct gb_pair_data *data)
{	
	// default ret->left and ret->right are NULL
	struct gb_pair_buf *ret = calloc(1, sizeof(struct gb_pair_buf));

	if (data->finish_flag)
		return ret;

	int32_t prev_i1, prev_i2, new_pos1, new_pos2;
	int32_t padding1, padding2;

	load_buffer(data->left);
	load_buffer(data->right);

	prev_i1 = prev_i2 = 0;
	while (1) {
		if (data->format == TYPE_FASTQ) {
			new_pos1 = get_buffer_read_fq(data->left, prev_i1);
			new_pos2 = get_buffer_read_fq(data->right, prev_i2);
		} else {
			new_pos1 = get_buffer_read_fa(data->left, prev_i1);
			new_pos2 = get_buffer_read_fa(data->right, prev_i2);
		}
		if (new_pos1)
			prev_i1 = new_pos1;
		if (new_pos2)
			prev_i2 = new_pos2;
		if (!new_pos1 || !new_pos2) 
			break;
	}
	if (!prev_i1 && !prev_i2) {
		data->finish_flag = true;
		return ret;
	}
	if (!prev_i1) {
		data->finish_flag = true;
		ret->left = data->left->buffer;
		return ret;
	}
	if (!prev_i2) {
		data->finish_flag = true;
		ret->right = data->right->buffer;
		return ret;
	}
	ret->left = data->left->buffer;
	ret->right = data->right->buffer;
	split_buffer(data->left, ret->left, prev_i1);
	split_buffer(data->right, ret->right, prev_i2);
	
	return ret;
}

struct gb_single_buf *gb_get_single(struct gb_single_data *data)
{	
	// default ret->single is NULL
	struct gb_single_buf *ret = calloc(1, sizeof(struct gb_single_buf));

	if (data->finish_flag)
		return ret;

	int32_t prev_i, new_pos;
	int32_t padding;

	load_buffer(data->single);

	prev_i = 0;
	while (1) {
		if (data->format == TYPE_FASTQ)
			new_pos = get_buffer_read_fq(data->single, prev_i);
		else
			new_pos = get_buffer_read_fa(data->single, prev_i);

		if (!new_pos) 
			break;

		prev_i = new_pos;
	}
	if (!prev_i) {
		data->finish_flag = true;
		return ret;
	}

	ret->single = data->single->buffer;
	split_buffer(data->single, ret->single, prev_i);
	
	return ret;
}