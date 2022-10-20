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
#include <stdio.h>
#include "main.h"
#include "Enclave_u.h"
#include <thread>

/* ecall_libcxx_functions:
 *   Invokes standard C++11 functions.
 */

// This function is part of mutex demo
void demo_counter_without_mutex()
{
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_mutex_demo_no_protection(global_eid);
	if (ret != SGX_SUCCESS)
		abort();
}

// This function is part of mutex demo
void demo_counter_mutex()
{
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_mutex_demo(global_eid);
	if (ret != SGX_SUCCESS)
		abort();
}

// This function is used by processing thread of condition variable demo
void demo_cond_var_run()
{
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_condition_variable_run(global_eid);
	if (ret != SGX_SUCCESS)
		abort();
}

// This function is used by the loader thread of condition variable demo
void demo_cond_var_load()
{
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_condition_variable_load(global_eid);
	if (ret != SGX_SUCCESS)
		abort();
}

// Examples for C++11 library and compiler features
void ecall_libcxx_functions(void)
{
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;

	// Example for lambda function feature:
	ret = ecall_lambdas_demo(global_eid);
	if (ret != SGX_SUCCESS)
		abort();
}

void SGX_init_Query(int *QueryMatrix, int *QueryLabel, int querysize, int labelsize)
{
	/*Initializing Query into SGX*/
	sgx_status_t ret;
	size_t QueryLen = querysize * querysize * sizeof(int), LabelLen = querysize * sizeof(int);
	ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_init_Query(global_eid, QueryMatrix, QueryLen, QueryLabel, LabelLen, querysize, labelsize);
	if (ret != SGX_SUCCESS)
		abort();
}

double SGX_LGPQ_BF(BF_value *BallBF, ofstream &BFRuntime, ofstream &BFOutFile, int pointer, int querysize, int sizeball, int BFtype, VertexLabel CenterLabel, int* result_BF)
{
	sgx_status_t ret;
	clock_t startTime, endTime;

	char buffer[ENTRY_VALUE_LEN + 1];
	size_t resultsize = sizeof(char);
	int bf_size = 0;
	char *BFresult = new char[1];
	double thisBallTime = 0;

	startTime = clock();
	for (auto it = BallBF->begin(); it != BallBF->end(); it++)
	{
		bf_size += (it->second.size() - 1);
	}

	// uint64_t vector_size = 315000000;//hold up 15 mil k,v // about 40MB
	// uint64_t vector_size = 460000000;//hold up 22 mil k,v // about 55MB
	uint64_t vector_size = 15000; // bits (bool) to char //
	// uint64_t vector_size = 830000000;//hold up 40 mil k,v // about 110MB
	uint8_t numHashs = 23;
	BloomFilter *myBloomFilter;
	size_t len = ENTRY_VALUE_LEN + 1;

	//construct BF
	if ((bf_size > 0) && (BFtype == 1))
	{
		myBloomFilter = new BloomFilter(vector_size, numHashs);

		for (auto it1 = BallBF->begin(); it1 != BallBF->end(); it1++)
		{
			for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
			{
				if (*it2 == -1) //-1 is used to initialize
					continue;
				// if(*it2==532632423)

				// sprintf(buffer, "%ld", *it2);
				LongToChar(buffer, *it2, ENTRY_VALUE_LEN, 10);
				myBloomFilter->add((uint8_t *)buffer, len);
			}
		}
	}
	else
	{
		myBloomFilter = new BloomFilter(vector_size, numHashs);
	}



	endTime = clock();
	thisBallTime += (double)(endTime - startTime) / CLOCKS_PER_SEC;


	result_BF[pointer] = 1;

	// BF Test in the enclave
	if ((bf_size > 0) && (BFtype == 1))
	{
		if (myBloomFilter->getbitsLen() % 8 != 0)
			cout << "The BloomFilter size is not the times of 8!" << endl;
		char *bits = new char[myBloomFilter->getbitsLen() / 8];
		if (bf_size > 10000)
			cout << "BF size: " << bf_size << endl;
		myBloomFilter->unloadBloomFilter(bits, myBloomFilter->getbitsLen());
		size_t bitsize = myBloomFilter->getbitsLen() / 8 * sizeof(char);
		BFresult[0] = '1';
		startTime = clock();

		ret = SGX_ERROR_UNEXPECTED;
		ret = ecall_BF(global_eid, bits, bitsize, BFresult, resultsize, numHashs, querysize, CenterLabel);
		if (ret != SGX_SUCCESS)
			abort();

		// if((BFresult[0]=='1')&&(result_BF[pointer]!=1))
		//	cout<<">> Ball " << pointer << ": The results of BF inside and outside SGX are different!"<<endl;
		if (BFresult[0] != '1')
		{
			result_BF[pointer] = 0;
		}
		endTime = clock();
		thisBallTime += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		BFRuntime << sizeball << ", " << thisBallTime * 1000 << endl;
		delete[] bits;
	}
	else
	{
		BFRuntime << sizeball << ", " << thisBallTime * 1000 << endl;
	}


	delete[] BFresult;
	delete myBloomFilter;
	return thisBallTime;

	// ret = SGX_ERROR_UNEXPECTED;

	// Example for lambda function feature:
	// Serialization
	// ret = ecall_lambdas_demo(global_eid);
	// if (ret != SGX_SUCCESS)
	// abort();
}

void LongToChar(char *C_data, long L_data, int len, int base)
{
	long tempdata = L_data;
	int templen = len - 1;
	while (tempdata > 0)
	{
		C_data[templen] = ('0' + (tempdata % base));
		tempdata = tempdata / base;
		templen--;
	}

	while (templen >= 0)
	{
		C_data[templen] = '\0';
		templen--;
	}
}