/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *      Version 1.2.3
 *      Modified by Tuan Tran 2017-08-29.
 *
 */

#ifndef SSW_H
#define SSW_H
#include <stdint.h>
#include <emmintrin.h>
#include "hash_align.h"

#define MAT_N 5

// If the query sequence is: ACGTATC, the sequence that read points to can be: 0123031
//
// Then if the penalty for match is 1 and for mismatch is -2, the substitution
// matrix of parameter mat will be:
//               A  C  G  T  N
//               1 -2 -2 -2  0  |A
//              -2  1 -2 -2  0  |C
//              -2 -2  1 -2  0  |G
//              -2 -2 -2  1  0  |T
//               0  0  0  0  0  |N
//
// mat is the pointer to the array
//               { 1, -2, -2, -2,  0,
//                -2,  1, -2, -2,  0,
//                -2, -2,  1, -2,  0,
//                -2, -2, -2,  1,  0,
//                 0,  0,  0,  0,  0 }
//
// mat has size (n * n)

// Do Striped Smith-Waterman alignment.
Cigar *ssw_align(char *ref, int32_t ref_len, char *read,
    int32_t read_len, int *score, int *skip);

#endif  // SSW_H