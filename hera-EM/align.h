#ifndef _ALIGN_
#define _ALIGN_

#include <stdint.h>

struct f_align_t {
        uint32_t trans_id;
        uint32_t pos;
        uint16_t f_len;
        uint16_t dis;
        uint8_t ori;
};

struct s_align_t {
        uint32_t trans_id;
        double c_prob;
};

struct s_read_t {
        uint16_t len;
        uint16_t cnt;
        double prob; //lock prob
};

union ualign_t {
        struct s_read_t read;
        struct f_align_t full;
        struct s_align_t simp;
};

struct align_batch_t {
        double noise;
        uint64_t size;
        uint64_t read_cnt;
        uint64_t cap;
        union ualign_t *align;
};
void align_batch_release(struct align_batch_t *a);
void batch_resize(struct align_batch_t *a);

#endif