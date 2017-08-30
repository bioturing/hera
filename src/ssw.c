/* The MIT License
   Copyright (c) 2012-1015 Boston College.
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Mengyao Zhao <zhangmp@bc.edu> */

/*
 *  ssw.c
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 1.2.3
 *	Last revision by Mengyao Zhao on 2017-06-26.
 *      Modified by Tuan Tran 2017-08-29.
 *
 */

#include "ssw.h"
#include "hash_align.h"

int8_t mat[25] = {
	 1, -2, -2, -2, 0,
	-2,  1, -2, -2, 0,
	-2, -2,  1, -2, 0,
    	-2, -2, -2,  1, 0,
    	 0,  0,  0,  0, 0
}; 	

int8_t ascii2int[128] = {
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

// Put the largest number of the 16 numbers in vm into m.
#define _mm_max16(m, vm) 					    		\
	(vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 8)); 			\
	(vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 4)); 			\
	(vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 2)); 			\
	(vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 1)); 			\
	(m) = _mm_extract_epi16((vm), 0)

// Find max of two number
#define _max2(x, y) (x > y ? x : y)

// Find max of three number
#define _max3(x, y, z) _max2(x, _max2(y, z))

static inline uint8_t get_value_pH(int32_t x_pos, int32_t y_pos,
    int32_t seg_len, __m128i *pvH)
{
	int32_t shift_pos = x_pos / seg_len;		
	x_pos = x_pos % seg_len;			
	__m128i vH = pvH[y_pos * seg_len + x_pos];			
	uint8_t *ret = (uint8_t*) &vH;
	return ret[shift_pos];
}

Cigar *ssw_align(char *ref_seq, int32_t ref_len, char *read_seq, int32_t read_len,
    int *score, int *skip)
{
	uint8_t weight_gapO = 4;
	uint8_t max = 0;
	uint8_t bias = 2;
	uint8_t read[read_len], ref[ref_len]; 
	int32_t end_read = read_len - 1;
	int32_t end_ref = -1;
	int32_t seg_len = (read_len + 15) >> 4; 
	int32_t edge, pos_max, cmp, i, j;
	uint8_t temp;
	uint8_t *t;

	for (i = 0; i < read_len; ++i)
		read[i] = ascii2int[read_seq[i]];

	for (i = 0; i < ref_len; ++i)
		ref[i] = ascii2int[ref_seq[i]];

	Cigar *cigar = init_cigar();
	__m128i vZero = _mm_set1_epi32(0);
	__m128i vProfile[MAT_N * seg_len];
	__m128i pvH[(ref_len + 1) * seg_len];
	__m128i vGapO = _mm_set1_epi8(weight_gapO);
	__m128i vBias = _mm_set1_epi8(bias);
	__m128i vMaxScore = vZero; 
	__m128i vMaxMark = vZero; 
	__m128i vTemp, vE, vF, vH, vMaxColumn;
	__m128i* vP;


	memset(vProfile, 0, MAT_N * seg_len * sizeof(__m128i));
	memset(pvH, 0, seg_len * sizeof(__m128i));
	for (i = 0; i < read_len; ++i){
		j = read[i] * seg_len + i % seg_len;
		t = (uint8_t*) &vProfile[j];
		t[i / seg_len] = 3;
	}

		
	for (i = 1; i <= ref_len; ++i) {
		vF = vZero;
		vMaxColumn = vZero; 
		vH = pvH[(i - 1) * seg_len + seg_len - 1];
		vH = _mm_slli_si128(vH, 1); 
		vP = vProfile + ref[i - 1] * seg_len;
	
		// Iner loop
		for (j = 0; j < seg_len; ++j) {
			vH = _mm_adds_epu8(vH, vP[j]);
			vH = _mm_subs_epu8(vH, vBias); 

			vE = pvH[(i - 1) * seg_len + j];

			// Get max from vH, vE and vF
			vH = _mm_max_epu8(vH, _mm_subs_epu8(vE, vGapO));
			vH = _mm_max_epu8(vH, vF);

			vMaxColumn = _mm_max_epu8(vMaxColumn, vH);
			pvH[i * seg_len + j] = vH;

			vF = _mm_subs_epu8(vH, vGapO);

			// Load the next vH
			vH = vE;
		}

		// Lazy_F loop
	        j = 0;
	        vH = pvH[i * seg_len + j];
	        vF = _mm_slli_si128(vF, 1);
	        vTemp = _mm_subs_epu8(vF, vH);
		vTemp = _mm_cmpeq_epi8(vTemp, vZero);
		cmp  = _mm_movemask_epi8(vTemp);

	        while (cmp != 0xffff) {
			vH = _mm_max_epu8(vH, vF);
			vMaxColumn = _mm_max_epu8(vMaxColumn, vH);
			pvH[i * seg_len + j] = vH;
			vF = _mm_subs_epu8(vH, vGapO);
			++j;

			if (j >= seg_len) {
				j = 0;
				vF = _mm_slli_si128(vF, 1);
			}

			vH = pvH[i * seg_len + j];
			vTemp = _mm_subs_epu8(vF, vH);
			vTemp = _mm_cmpeq_epi8(vTemp, vZero);
			cmp  = _mm_movemask_epi8(vTemp);
	        }

		vMaxScore = _mm_max_epu8(vMaxScore, vMaxColumn);
		vTemp = _mm_cmpeq_epi8(vMaxMark, vMaxScore);
		cmp = _mm_movemask_epi8(vTemp);

		if (cmp != 0xffff) {
			vMaxMark = vMaxScore;
			_mm_max16(temp, vMaxScore);

			if (temp > max) {
				max = temp;
				end_ref = i - 1;
				pos_max = i;
			}
		}
	}

	int32_t begin_ref = end_ref;
	*score = *skip = 0;
	
	if (max < read_len >> 2) {
		*score = read_len;
		add_cigar(cigar, SCLIP, read_len, NULL);
		goto end;
	}

	// Trace the alignment ending position on read.
	t = (uint8_t*) (pvH + pos_max * seg_len);
	int32_t column_len = seg_len * 16;
	int32_t tmp;
	
	for (i = 0; i < column_len; ++i, ++t) {
		if (*t == max) {
			tmp = (i >> 4) + (i & 15) * seg_len;
			if (tmp < end_read) end_read = tmp;
		}
	}
	
	// Trace the cigar
	int32_t x_pos = end_read;
	int32_t y_pos = end_ref + 1;
	int32_t up_value, left_value, diag_value, cur_value, max_value;
	uint8_t dir = 0;
	int32_t begin_read = end_read;
	int32_t right_sclip = read_len - end_read - 1;

	if (right_sclip > 0) {
		if (right_sclip < 5) {
			for (i = read_len - 1, j = end_ref + right_sclip; 
			    i > end_read; i--, j--) {
				if (j >= ref_len) {
					add_cigar(cigar, SCLIP, 1, NULL);
					++*score;
				} else {
					if (read[i] == ref[j]) {
						add_cigar(cigar, MATCH, 1, NULL);
					} else {
						add_cigar(cigar, DIFF, 1,
						    ref_seq + j);
						++*score;
					}
				}
			}
		} else {
			*score += right_sclip;
			add_cigar(cigar, SCLIP, right_sclip, NULL);
		}
	}

	cur_value = get_value_pH(x_pos, y_pos, seg_len, pvH);

	do {	
		switch (dir) {
		case 0:
			if (read[begin_read] == ref[begin_ref]) {
				add_cigar(cigar, MATCH, 1, NULL);
			} else {
				add_cigar(cigar, DIFF, 1,
				    ref_seq + begin_ref);
				++*score;
			}
			break;
		case 1:
			++*score;
			add_cigar(cigar, INSERT, 1, NULL);
			break;
		case 2:
			++*score;
			add_cigar(cigar, DELETE, 1, ref_seq + begin_ref);
			break;
		}
		
		if (x_pos == 0 || y_pos == 1) break;
		
		diag_value = get_value_pH(x_pos - 1, y_pos - 1, seg_len, pvH);
		up_value = get_value_pH(x_pos - 1, y_pos, seg_len, pvH);
		left_value = get_value_pH(x_pos, y_pos - 1, seg_len, pvH);

		if (cur_value == diag_value +
		    mat[read[begin_read] * MAT_N + ref[begin_ref]]) {
			if (!diag_value) break;
			cur_value = diag_value;
			--begin_read;
			--begin_ref;
			--x_pos;
			--y_pos;
			dir = 0;
		} else if (cur_value == left_value - weight_gapO) {
			if (!left_value) break;
			cur_value = left_value;
			--y_pos;
			--begin_ref;
			dir = 1;
		} else {
			if (!up_value) break;
			cur_value = up_value;
			--x_pos;
			--begin_read;
			dir = 2;
		}
	} while (1);

	if (begin_read < 5) {
		for (--begin_read, --begin_ref;
			begin_read >= 0; --begin_read, --begin_ref) {
			if (begin_ref < 0) {
				add_cigar(cigar, SCLIP, 1, NULL);
				++*score;
			} else {
				if (read[begin_read] == 
					ref[begin_ref]) {
					add_cigar(cigar, MATCH, 1, NULL);
				} else {
					add_cigar(cigar, DIFF, 1,
						ref_seq + begin_ref);
					++*score;
				}
			}
		}
	} else if (begin_read > 0) {
		*score += begin_read;
		add_cigar(cigar, SCLIP, begin_read, NULL);
	}
	
end:
	*skip = begin_ref;
	return cigar;
}	
