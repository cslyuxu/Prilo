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


//#include <iostream>
//#include "cgbe.h"
#include <stdio.h>
#include <string.h>
#include <string> 
#include <assert.h>
#include <algorithm>

# include <unistd.h>
# include <pwd.h>
# define MAX_PATH FILENAME_MAX

#include "sgx_urts.h"
#include "main.h"
#include "Enclave_u.h"


#include "DataRead.h"
#include "LGPQ.h"
#include <fstream>
#include "cgbe.h"
using namespace std;


/* Global EID shared by multiple threads */
sgx_enclave_id_t global_eid = 0;

typedef struct _sgx_errlist_t {
    sgx_status_t err;
    const char *msg;
    const char *sug; /* Suggestion */
} sgx_errlist_t;

/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] = {
    {
        SGX_ERROR_UNEXPECTED,
        "Unexpected error occurred.",
        NULL
    },
    {
        SGX_ERROR_INVALID_PARAMETER,
        "Invalid parameter.",
        NULL
    },
    {
        SGX_ERROR_OUT_OF_MEMORY,
        "Out of memory.",
        NULL
    },
    {
        SGX_ERROR_ENCLAVE_LOST,
        "Power transition occurred.",
        "Please refer to the sample \"PowerTransition\" for details."
    },
    {
        SGX_ERROR_INVALID_ENCLAVE,
        "Invalid enclave image.",
        NULL
    },
    {
        SGX_ERROR_INVALID_ENCLAVE_ID,
        "Invalid enclave identification.",
        NULL
    },
    {
        SGX_ERROR_INVALID_SIGNATURE,
        "Invalid enclave signature.",
        NULL
    },
    {
        SGX_ERROR_OUT_OF_EPC,
        "Out of EPC memory.",
        NULL
    },
    {
        SGX_ERROR_NO_DEVICE,
        "Invalid SGX device.",
        "Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards."
    },
    {
        SGX_ERROR_MEMORY_MAP_CONFLICT,
        "Memory map conflicted.",
        NULL
    },
    {
        SGX_ERROR_INVALID_METADATA,
        "Invalid enclave metadata.",
        NULL
    },
    {
        SGX_ERROR_DEVICE_BUSY,
        "SGX device was busy.",
        NULL
    },
    {
        SGX_ERROR_INVALID_VERSION,
        "Enclave version was invalid.",
        NULL
    },
    {
        SGX_ERROR_INVALID_ATTRIBUTE,
        "Enclave was not authorized.",
        NULL
    },
    {
        SGX_ERROR_ENCLAVE_FILE_ACCESS,
        "Can't open enclave file.",
        NULL
    },
    {
        SGX_ERROR_NDEBUG_ENCLAVE,
        "The enclave is signed as product enclave, and can not be created as debuggable enclave.",
        NULL
    },
    {
        SGX_ERROR_MEMORY_MAP_FAILURE,
        "Failed to reserve memory for the enclave.",
        NULL
    },
};

/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret)
{
    size_t idx = 0;
    size_t ttl = sizeof sgx_errlist/sizeof sgx_errlist[0];

    for (idx = 0; idx < ttl; idx++) {
        if(ret == sgx_errlist[idx].err) {
            if(NULL != sgx_errlist[idx].sug)
                printf("Info: %s\n", sgx_errlist[idx].sug);
            printf("Error: %s\n", sgx_errlist[idx].msg);
            break;
        }
    }
    
    if (idx == ttl)
        printf("Error: Unexpected error occurred.\n");
}

/* Initialize the enclave:
 *   Call sgx_create_enclave to initialize an enclave instance
 */
int initialize_enclave(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    
    /* Call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */
    ret = sgx_create_enclave(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, NULL, NULL, &global_eid, NULL);
    if (ret != SGX_SUCCESS) {
        print_error_message(ret);
        return -1;
    }

    return 0;
}

/* OCall functions */
void ocall_print_string(const char *str)
{
    /* Proxy/Bridge will check the length and null-terminate 
     * the input string to prevent buffer overflow. 
     */
    printf("%s", str);
}



int SGX_CDECL main(int argc, char *argv[])	{
	//used in the DIGraph.h for the random label distribution
	//srand((int)time(0));
	//use the same randomness
	srand(0);

	//string QueryName = argv[1];

    (void)(argc);
    (void)(argv);
    

    /* Initialize the enclave */
    if(initialize_enclave() < 0){
        printf("Enter a character before exit ...\n");
        getchar();
        return -1; 
    }
 
	//cout<<"Query Name:" << QueryName<<endl;
    
   

	
	cout << "\n\n*************************************" << endl;
	cout << "********   Loading Data    **********" << endl;
	cout << "*************************************" << endl;
	string QueryName = argv[1];
	//data_read_with_label.h
	//Read<VertexLabel, EdgeLabel> Query(QueryName);		
	DataRead<VertexLabel, EdgeLabel> Query(QueryName);				
	cout << "Query: "<<QueryName;
	cout << "\t|V|: " << Query.graph_ptr->getVcnt() << "\t\t|E|: " << Query.graph_ptr->getEcnt() <<endl;
	//cout<< QueryName << endl;
	string GraphName = argv[2];
	int pos = QueryName.find_last_of('/');
	string temp(QueryName.substr(pos+1));
	replace(temp.begin(), temp.end(), '/', '-');
	string OutFile = "LGPQ/results/" + temp;
	//cout<< OutFile<< endl;

	//data_read.h
	DataRead<VertexLabel, EdgeLabel> Graph(GraphName);
	cout << "Graph: "<<GraphName;
	cout << "\t|V|: " << Graph.graph_ptr->getVcnt() << "\t|E|: " << Graph.graph_ptr->getEcnt() <<endl;
	cout << endl;

	int labelsize = stoi(argv[3]);	
	int twigletlength = stoi(argv[4]);
	int mode = stoi(argv[5]);
    int threshold = stoi(argv[6]);
    int NLhop, pathlength;
    NLhop = twigletlength;
    pathlength = twigletlength;
	LGPQ<VertexLabel, EdgeLabel> ssim(&Query, &Graph, OutFile, NLhop, pathlength, labelsize, pathlength, mode, threshold);
    //For testing LGPQ
	//ssim.Match();
	


    
    /* Utilize trusted libraries */  
    //For testing SGX+BloomFilter
    //unordered_map<int, double>* BFTime = new unordered_map<int, double>;
    //BFTime->clear();
    //ecall_LGPQ_BF(&Query, &Graph, ssim.result_BF, labelsize, BFTime, mode, OutFile, threshold);

     ssim.Match();
    
    /* Destroy the enclave */
    sgx_destroy_enclave(global_eid);
    
    printf("Info: LGPQEnclave successfully returned.\n");
   
   

	cout << "Service done!" << endl;
	return 0;
	

}



