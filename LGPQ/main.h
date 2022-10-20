/*
 * Copyright (C) 2011-2021 Intel Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of Intel Corporation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef _MAIN_H_
#define _MAIN_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include <unordered_map> 
#include <unordered_set> 
using namespace std;

#include "DataRead.h"
#include "BF/BloomFilter.h"

#include "sgx_error.h"       /* sgx_status_t */
#include "sgx_eid.h"     /* sgx_enclave_id_t */

#ifndef TRUE
# define TRUE 1
#endif

#ifndef FALSE
# define FALSE 0
#endif

#if   defined(__GNUC__)
# define TOKEN_FILENAME   "enclave.token"
# define ENCLAVE_FILENAME "enclave.signed.so"
#endif

extern sgx_enclave_id_t global_eid;    /* global enclave id */

#if defined(__cplusplus)
extern "C" {
#endif

// For BF
typedef unordered_set<long> BFvalue;
typedef unordered_map<int, BFvalue> BF_value;
// For nodes of tree pattern
typedef unordered_map<VertexID, int> Degree;
typedef unordered_set<VertexLabel> PatternLabel;
// typedef unordered_map<VertexID, int> MapMatrix;			
typedef unordered_map<VertexID, VertexLabel> VLabels;

void ecall_libcxx_functions(void);

void SGX_init_Query(int*, int*, int, int);

double SGX_LGPQ_BF(BF_value *, ofstream &, ofstream &, int, int, int, int, VertexLabel, int*);

int constructBF(BF_value *, DIGRAPH<VertexLabel, EdgeLabel> *, VertexID, int, Degree *, Degree *, VertexLabel *, long *, int);

void BF(BF_value **, BloomFilter *, int, VertexID, int, DIGRAPH<VertexLabel, EdgeLabel>*, DIGRAPH<VertexLabel, EdgeLabel>*, int*);

void LongToChar(char *C_data, long L_data, int len, int base);
void MultiServers(int, int, unordered_map<int, double> &, int, int *, ofstream &);

#if defined(__cplusplus)
}
#endif

#endif /* !_MAIN_H_ */
