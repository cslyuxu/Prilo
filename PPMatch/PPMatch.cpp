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

 //This function is part of mutex demo
void demo_counter_without_mutex()
{
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_mutex_demo_no_protection(global_eid);
	if (ret != SGX_SUCCESS)
		abort();
}

//This function is part of mutex demo
void demo_counter_mutex()
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    ret = ecall_mutex_demo(global_eid);
    if (ret != SGX_SUCCESS)
        abort();
}

//This function is used by processing thread of condition variable demo
void demo_cond_var_run()
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    ret = ecall_condition_variable_run(global_eid);
    if (ret != SGX_SUCCESS)
        abort();
}

//This function is used by the loader thread of condition variable demo
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

void ecall_PPMatch_BF(DataRead<VertexLabel, EdgeLabel>& Query, DataRead<VertexLabel, EdgeLabel>& Graph, int* result_BF, int labelsize, unordered_map<int, double>* BFTime, int mode, string OutFileName, int threshold)
{
	sgx_status_t ret;
    int querysize = Query.graph_ptr->getVcnt();
    int graphsize = Graph.graph_ptr->getVcnt();
    //cout<<"query size: "<< querysize<<endl;
    //cout<<"graph size: "<< graphsize<<endl;
	int* BallCount = new int[graphsize]; //used for multiple servers
	for(int i = 0; i <graphsize;i++)
		BallCount[i] = 0;
    clock_t startTime, endTime;

	//output the runtimes
	string RuntimeName = OutFileName + "-BFruntime";	
	string OverallName = OutFileName + "-BF";
	ofstream Runtime(RuntimeName);	
	ofstream OutFile(OverallName);
	Runtime << "# |V_Q| = " << querysize << endl;
	Runtime << "# ball size, exact, 2-iter, NL, Path" << endl;

	char buffer[ENTRY_VALUE_LEN+1];
	char *BFresult = new char[1];
	size_t resultsize = sizeof(char);

     /*intialize matrix M_Q*/
    Query.graph_ptr->getDiameter();
    int dia = Query.graph_ptr->diameter;
	cout << "The diameter of the Query is " << dia << endl;

    int **Q= new int *[querysize];
	size_t QueryLen = querysize*querysize*sizeof(int), LabelLen = querysize*sizeof(int);
	int* QueryMatrix = new int[querysize*querysize];
	int* QueryLabel = new int[querysize];
	for (int i = 0; i < querysize; i++)
	{
		QueryLabel[i] = Query.graph_ptr->getVLabel(Query.graph_ptr->matrix[i]);
		Q[i] = new int[querysize];		
		for (int j = 0; j < querysize; j++)
		{
			Q[i][j] = 1;			
			if (Query.graph_ptr->isEdge(Query.graph_ptr->matrix[i], Query.graph_ptr->matrix[j]))
			{
				Q[i][j] = 0;				
			}
			QueryMatrix[i*querysize+j]=Q[i][j];
			//cout << Q[i][j] << " ";		
		}
		//cout << endl;
	}
    cout << endl;

	
	/*Initializing Query into SGX*/
	ret = SGX_ERROR_UNEXPECTED;
	ret = ecall_init_Query(global_eid, QueryMatrix, QueryLen, QueryLabel, LabelLen, querysize, labelsize);
    if (ret != SGX_SUCCESS)
        abort();
	
	delete[] QueryLabel;
	delete[] QueryMatrix;


    VertexLabel LabelofQ;
	int MaxLabelNum = 0;
	Graph.graph_ptr->CountLabelNum(Query.graph_ptr);
	for (auto it = Graph.graph_ptr->getLabelCount().begin(); it != Graph.graph_ptr->getLabelCount().end(); it++)
	{
		if (it->second > MaxLabelNum)
		{
			MaxLabelNum = it->second;
			LabelofQ = it->first;
		}
	}
	cout << "The Label " << LabelofQ << " has a maximum number of " << MaxLabelNum << ",  Total Label: "<<labelsize<<endl;
    cout << endl;


    // precompute the tree pattern for bloomfilter
	unordered_set<long> templong = {-1};
	long *weight = new long[6];
	weight[0] = 1;
	for (int i = 1; i < 6; i++)
	{
		weight[i] = weight[i - 1] * labelsize;
		//cout << weight[i] << endl;
	}
    BF_value** QueryBF;
    QueryBF = new BF_value *[querysize];
	for (int i = 0; i < querysize; i++)
	{
        QueryBF[i] = new BF_value;
		//QueryBF[i] = new BF_value;
		//output the labels of Q's vertices
		//cout << "The " << i << "th vertex of Q!" <<this->query->getVLabel(this->query->matrix[i])<< endl;
		for (int j = 1; j < 5; j++)
		{
			QueryBF[i]->emplace(j, templong);
		}
		constructBF(QueryBF[i], Query.graph_ptr, Query.graph_ptr->matrix[i], -1, nullptr, nullptr, nullptr, weight, threshold);
	}
    //cout<<">> The Tree Patterns for Query have been built!"<<endl;


    // pointer for short test
	int pointer = 0;
    int bf_size;

	// Num_Ball_Computed
	int Num_Ball_Computed = 0, Num_Ball_Ignored = 0;
	unordered_map<int, double> localBFTime;
	localBFTime.clear();

    double BF_Time = 0, thisBallTime;
	int BFtype;

	for (typename VLabels::iterator it1 = Graph.graph_ptr->getVLabel().begin(); it1 != Graph.graph_ptr->getVLabel().end(); it1++)
	{

        VertexID s = it1->first, t;

		/*if the center of the ball's label is not in the label set of query, then continue*/
		if (!(Query.graph_ptr->isLabel(Graph.graph_ptr->getVLabel(s))))
		{
			result_BF[pointer] = -1;
			pointer++;
			continue;
		}

        //For Graph Homomorphism
        if(mode == 1){
			if (Graph.graph_ptr->getVLabel(s) != LabelofQ)
			{
				result_BF[pointer] = -1;
				pointer++;
				continue;
			}
		}


        unordered_set<VertexID> s_neighbor;
		Graph.graph_ptr->getDNeighbor(s, dia, s_neighbor, Query.graph_ptr);
		int size_ball = s_neighbor.size();


        //mapping # to vertex
        VertexID *Matrix_ball = new VertexID[size_ball];
		unordered_map<VertexID, int> Matrix_Ball;
		int tempnum = 0;
		for (unordered_set<VertexID>::iterator it2 = s_neighbor.begin();
			 it2 != s_neighbor.end(); it2++)
		{
			t = *it2;
			Matrix_ball[tempnum] = t;
			Matrix_Ball[t] = tempnum;
			tempnum++;
			//	cout << "The vertex " << t << " is contained in the neighbor of vertex " << s << endl;
		}


        bool Prune_flag;	
		for (int i = 0; i < querysize; i++)
		{
			Prune_flag = true;
			/*Pruning*/
			for (int j = 0; j < size_ball; j++)
			{
				if (Query.graph_ptr->getVLabel(Query.graph_ptr->matrix[i]) == Graph.graph_ptr->getVLabel(Matrix_ball[j]))
				{
					Prune_flag = false;
				}
			}
			if (Prune_flag)
			{
				break;
			}
		}

        if (Prune_flag)
		{
			result_BF[pointer] = 0;
			delete[] Matrix_ball;
			Matrix_Ball.clear();
			s_neighbor.clear();
			pointer++;
			continue;
		}







        // build the bloomfilter
		BF_value *BallBF = new BF_value;
		for (int i = 1; i < 5; i++)
		{
			BallBF->emplace(i, templong);
		}

		// use for the ball structure in graph structure
		DIGRAPH<VertexLabel, EdgeLabel> *Ball = new DIGRAPH<VertexLabel, EdgeLabel>;
		for (int i = 0; i < size_ball; i++)
		{
			Ball->insertVertex(Matrix_ball[i], Graph.graph_ptr->getVLabel(Matrix_ball[i]));
		}

		for (auto it3 = Ball->getVLabel().begin(); it3 != Ball->getVLabel().end(); it3++)
		{
			for (auto it4 = Ball->getVLabel().begin(); it4 != Ball->getVLabel().end(); it4++)
			{
				if (it3->first == it4->first)
					continue;
				if ((Graph.graph_ptr->isEdge(it3->first, it4->first)) && (!Ball->isEdge(it3->first, it4->first)))
					Ball->insertEdge(it3->first, it4->first, 0);
				
			}
		}
		//cout << "Ball Vcnt:" << Ball->getVcnt() << endl;
		//cout << "Ball Ecnt:" << Ball->getEcnt() << endl;
		//cout << "Ball center degree:" << Ball->getOutDegree(s) + Ball->getInDegree(s) << endl;
		
        startTime = clock();
        BFtype = constructBF(BallBF, Ball, s, -1, nullptr, nullptr, nullptr, weight, threshold);
		

		bf_size = 0;
		localBFTime[pointer]=0;

		for (auto it = BallBF->begin(); it != BallBF->end(); it++)
		{
			// cout << it->second.size() - 1 << endl;
			bf_size += (it->second.size() - 1);
		}

		//uint64_t vector_size = 315000000;//hold up 15 mil k,v // about 40MB
    	//uint64_t vector_size = 460000000;//hold up 22 mil k,v // about 55MB
		uint64_t vector_size = 15000;// bits (bool) to char //
   		//uint64_t vector_size = 830000000;//hold up 40 mil k,v // about 110MB    
    	uint8_t numHashs = 23;
		BloomFilter *myBloomFilter;		
		size_t len = ENTRY_VALUE_LEN+1;

		if ((bf_size > 0)&&(BFtype==1))
		{	
			myBloomFilter = new BloomFilter(vector_size, numHashs);	
			
			for (auto it1 = BallBF->begin(); it1 != BallBF->end(); it1++)
			{
				for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
				{
					if (*it2 == -1) //-1 is used to initialize
						continue;
					//if(*it2==532632423)
					
					//sprintf(buffer, "%ld", *it2);
					LongToChar(buffer, *it2, ENTRY_VALUE_LEN, 10);
					myBloomFilter->add((uint8_t*)buffer, len);
				}
			}
		}
		else
		{
			myBloomFilter = new BloomFilter(vector_size, numHashs);	
		}

		endTime = clock();
		thisBallTime = 0;
		BF_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		thisBallTime += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		//cout<<"The number of Hashes: " << myBloomFilter->getnumHashes()<<endl;
		//cout<<"The length of bits: " << myBloomFilter->getbitsLen()<<endl;  

		
		BallCount[pointer] = 1;
		result_BF[pointer] = 1;

        // BloomFilter
		
		if ((bf_size > 0)&&(BFtype==1)){
			if(myBloomFilter->getbitsLen()%8 != 0)
				cout<<"The BloomFilter size is not the times of 8!"<<endl;
			char *bits=new char[myBloomFilter->getbitsLen()/8];
			if(bf_size>10000)
				cout<<"BF size: "<<bf_size<<endl;
			myBloomFilter->unloadBloomFilter(bits, myBloomFilter->getbitsLen());
			size_t bitsize = myBloomFilter->getbitsLen()/8*sizeof(char);
			BFresult[0]='1';
			startTime = clock();
			/*Initializing Query into SGX*/
			ret = SGX_ERROR_UNEXPECTED;
			ret = ecall_BF(global_eid, bits, bitsize, BFresult, resultsize, numHashs, querysize, Graph.graph_ptr->getVLabel(s));
    		if (ret != SGX_SUCCESS)
        		abort();
			
			///////BFTest(QueryBF, BallBF, query_size, s, pointer);
			//BF(QueryBF, myBloomFilter, querysize, s, pointer, Query.graph_ptr, Graph.graph_ptr, result_BF);	
			
			//if((BFresult[0]=='1')&&(result_BF[pointer]!=1))
			//	cout<<">> Ball " << pointer << ": The results of BF inside and outside SGX are different!"<<endl;
			if(BFresult[0]!='1'){
				result_BF[pointer]=0;
				BallCount[pointer] = -1;
			}
			endTime = clock();
			BF_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			thisBallTime += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			localBFTime[pointer] = thisBallTime;
			Runtime << size_ball << ", " << thisBallTime * 1000 << endl; 
			delete[] bits;
		}else{
			localBFTime[pointer] = thisBallTime;
			Runtime << size_ball << ", " << thisBallTime * 1000 << endl; 
		}


        /*free the space*/
		delete[] Matrix_ball;
		Matrix_Ball.clear();
		s_neighbor.clear();	
		delete BallBF;
		delete Ball;
		
		delete myBloomFilter;
		Num_Ball_Computed++;
		pointer++;
		
    }


    // cout number
    int BF_num=0;
	for (int i = 0; i < graphsize; i++)
	{
		if (result_BF[i] == 1)
			BF_num++;
    }

    cout <<"Total computed balls:"<<Num_Ball_Computed<<endl;
	cout << BF_num << " BF matches:" << BF_Time << " s" << endl;

	OutFile << "#ALL BF#	BFTime(s)	Player#	BFTime’(ms)" << endl;
	OutFile << Num_Ball_Computed << endl;
	OutFile << BF_num << endl;
	OutFile << BF_Time << endl;

    //cout << BF_num << " BF matches:" << BF_Time << " s" << endl;
	MultiServers(10, 4, localBFTime, graphsize, BallCount, OutFile);


    //ret = SGX_ERROR_UNEXPECTED;
    
    // Example for lambda function feature:
    //Serialization
    //ret = ecall_lambdas_demo(global_eid);
    //if (ret != SGX_SUCCESS)
    // abort();
    
   

    for (int i = 0; i < querysize; i++)
	{
		delete QueryBF[i];
	}
	delete QueryBF;
	delete weight;
	delete BFresult;
	delete[] BallCount;
	Runtime.close();
	OutFile.close();
	localBFTime.clear();

}

int constructBF(BF_value *BF, DIGRAPH<VertexLabel, EdgeLabel> *Graph, VertexID center, int type, Degree *left, Degree *right, VertexLabel *temptree, long *weight, int threshold)
{
	if (type < -1)
	{ // compute the weighted value
		long value = 0;
		for (int i = 0; i < 6; i++)
			value += temptree[i] * weight[i];

		for (auto it = BF->begin(); it != BF->end(); it++)
		{
			if (it->first != type + 7)
				continue;
			if (it->second.find(value) == it->second.end())
				it->second.insert(value);
			// cout << "center " << center << " with " << it->first << " pattern" << endl;
			// for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++)
			// cout << *it1 << endl;
		}

		return 1;
	}

	if (type == 6) // BF[1]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		VertexLabel *labels = new VertexLabel[6];
		// unordered_set<VertexID> *neighbor=new unordered_set<VertexID>;
		labels[4] = 0;
		labels[5] = 0;
		labels[3] = 0;
		patternlabel->insert(Graph->getVLabel(center));
		for (auto it1 = right->begin(); it1 != right->end(); it1++)
		{
			labels[1] = Graph->getVLabel(it1->first);
			patternlabel->insert(labels[1]);
			for (auto it2 = left->begin(); it2 != left->end(); it2++)
			{
				if (patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				labels[0] = Graph->getVLabel(it2->first);
				patternlabel->insert(labels[0]);
				leftchildlabel->clear();
				for (auto it3 = Graph->getOutEdge()[it2->first].begin(); it3 != Graph->getOutEdge()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}

				for (auto it3 = Graph->getInVertex()[it2->first].begin(); it3 != Graph->getInVertex()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}

				for (auto it3 = leftchildlabel->begin(); it3 != leftchildlabel->end(); it3++)
				{
					labels[2] = *it3;
					constructBF(BF, Graph, center, -6, nullptr, nullptr, labels, weight, threshold); // compute the weighted value
				}
				patternlabel->erase(Graph->getVLabel(it2->first));
			}
			patternlabel->erase(Graph->getVLabel(it1->first));
		}

		delete patternlabel;
		delete leftchildlabel;
		delete labels;
		// delete neighbor;
		return 1;
	}

	if (type == 7) // BF[2]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		VertexLabel *labels = new VertexLabel[6];
		Degree *degree0 = new Degree;
		Degree *degree1 = new Degree;		
		labels[4] = 0;
		labels[5] = 0;
		patternlabel->insert(Graph->getVLabel(center));
		for (auto it1 = right->begin(); it1 != right->end(); it1++)
		{
			labels[1] = Graph->getVLabel(it1->first);
			patternlabel->insert(labels[1]);
			for (auto it2 = left->begin(); it2 != left->end(); it2++)
			{
				if (patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				labels[0] = Graph->getVLabel(it2->first);
				patternlabel->insert(labels[0]);

				//all cases				
				degree1->clear();
				degree0->clear();
				degree1->emplace(it2->first, it2->second);
				degree0->emplace(it1->first, it1->second);
				constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight, threshold);

				leftchildlabel->clear();
				for (auto it3 = Graph->getOutEdge()[it2->first].begin(); it3 != Graph->getOutEdge()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}

				for (auto it3 = Graph->getInVertex()[it2->first].begin(); it3 != Graph->getInVertex()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}

				for (auto it3 = leftchildlabel->begin(); it3 != leftchildlabel->end(); it3++)
				{
					for (auto it4 = leftchildlabel->begin(); it4 != leftchildlabel->end(); it4++)
					{
						if (*it3 <= *it4)
							continue;
						labels[2] = *it3;
						labels[3] = *it4;
						constructBF(BF, Graph, center, -5, nullptr, nullptr, labels, weight, threshold); // compute the weighted value
					}
				}
				patternlabel->erase(Graph->getVLabel(it2->first));
			}
			patternlabel->erase(Graph->getVLabel(it1->first));
		}

		delete patternlabel;
		delete leftchildlabel;
		delete labels;
		delete degree0;
		delete degree1;		
		return 1;
	}

	if (type == 9) // BF[3]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		PatternLabel *rightchildlabel = new PatternLabel;
		Degree *degree0 = new Degree;
		Degree *degree1 = new Degree;
		VertexLabel *labels = new VertexLabel[6];
		labels[5] = 0;
		patternlabel->insert(Graph->getVLabel(center));
		for (auto it1 = right->begin(); it1 != right->end(); it1++)
		{
			labels[1] = Graph->getVLabel(it1->first);
			patternlabel->insert(labels[1]);
			rightchildlabel->clear();
			for (auto it2 = Graph->getOutEdge()[it1->first].begin(); it2 != Graph->getOutEdge()[it1->first].end(); it2++)
			{
				if (patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				rightchildlabel->insert(Graph->getVLabel(it2->first));
			}

			for (auto it2 = Graph->getInVertex()[it1->first].begin(); it2 != Graph->getInVertex()[it1->first].end(); it2++)
			{
				if (patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				rightchildlabel->insert(Graph->getVLabel(it2->first));
			}

			for(auto it2 = left->begin(); it2 != left->end(); it2++){
				if(patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				
                

				degree1->clear();
				degree0->clear();
				degree0->emplace(it1->first, it1->second);
				degree1->emplace(it2->first, it2->second);
				constructBF(BF, Graph, center, 6, degree0, degree1, nullptr, weight, threshold);
				constructBF(BF, Graph, center, 7, degree1, degree0, nullptr, weight, threshold);

				/*
				if((Graph->isEdge(it2->first, it1->first))||(Graph->isEdge(it1->first, it2->first))){
					degree1->clear();
					degree0->clear();
					degree0->emplace(it1->first, it1->second-1);
					degree1->emplace(it2->first, it2->second-1);
					//constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight);
					if(it2->second==2)
						constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight);
					else
						constructBF(BF, Graph, center, 7, degree1, degree0, nullptr, weight);
					//continue;
				}*/



				labels[0] = Graph->getVLabel(it2->first);
				patternlabel->insert(labels[0]);
				leftchildlabel->clear();

				for (auto it3 = Graph->getOutEdge()[it2->first].begin(); it3 != Graph->getOutEdge()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}

				for (auto it3 = Graph->getInVertex()[it2->first].begin(); it3 != Graph->getInVertex()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}


				for(auto it3 = rightchildlabel->begin(); it3 != rightchildlabel->end(); it3++){
					if(patternlabel->find(*it3)!=patternlabel->end())
						continue;
					labels[4]=*it3;
					patternlabel->insert(*it3);

					for (auto it4 = leftchildlabel->begin(); it4 != leftchildlabel->end(); it4++)
					{
						for (auto it5 = leftchildlabel->begin(); it5 != leftchildlabel->end(); it5++)
						{
							if (*it4 <= *it5)
								continue;
							labels[2] = *it4;
							labels[3] = *it5;
							constructBF(BF, Graph, center, -4, nullptr, nullptr, labels, weight, threshold); // compute the weighted value
						}
					}

					patternlabel->erase(labels[4]);
				}

				patternlabel->erase(labels[0]);

			}			
			patternlabel->erase(Graph->getVLabel(it1->first));
		}

		delete patternlabel;
		delete leftchildlabel;
		delete rightchildlabel;		
		delete degree0;
		delete degree1;
		delete labels;
		return 1;
	}

	if (type == 10) // BF[4]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		PatternLabel *rightchildlabel = new PatternLabel;
		Degree *degree0 = new Degree;
		Degree *degree1 = new Degree;
		Degree *degree2 = new Degree;
		VertexLabel *labels = new VertexLabel[6];
		patternlabel->insert(Graph->getVLabel(center));
		for (auto it1 = right->begin(); it1 != right->end(); it1++)
		{
			labels[1] = Graph->getVLabel(it1->first);
			patternlabel->insert(labels[1]);
			rightchildlabel->clear();
			for (auto it2 = Graph->getOutEdge()[it1->first].begin(); it2 != Graph->getOutEdge()[it1->first].end(); it2++)
			{
				if (patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				rightchildlabel->insert(Graph->getVLabel(it2->first));
			}

			for (auto it2 = Graph->getInVertex()[it1->first].begin(); it2 != Graph->getInVertex()[it1->first].end(); it2++)
			{
				if (patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;
				rightchildlabel->insert(Graph->getVLabel(it2->first));
			}



			for (auto it2 = left->begin(); it2 != left->end(); it2++)
			{
				if(patternlabel->find(Graph->getVLabel(it2->first)) != patternlabel->end())
					continue;


				degree1->clear();
				degree0->clear();
				degree0->emplace(it1->first, it1->second);
				degree1->emplace(it2->first, it2->second);
				//constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight);
				constructBF(BF, Graph, center, 9, degree1, degree0, nullptr, weight, threshold);
				/*
				if((Graph->isEdge(it2->first, it1->first))||(Graph->isEdge(it1->first, it2->first))){
					degree1->clear();
					degree0->clear();
					degree0->emplace(it1->first, it1->second);
					degree1->emplace(it2->first, it2->second);
					//constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight);
					if((it1->second>2)&&(it2->second==2))
					constructBF(BF, Graph, center, 9, degree1, degree0, nullptr, weight);
				}*/

				labels[0] = Graph->getVLabel(it2->first);
				patternlabel->insert(labels[0]);
				leftchildlabel->clear();
				
				
				for (auto it3 = Graph->getOutEdge()[it2->first].begin(); it3 != Graph->getOutEdge()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}

				for (auto it3 = Graph->getInVertex()[it2->first].begin(); it3 != Graph->getInVertex()[it2->first].end(); it3++)
				{
					if (patternlabel->find(Graph->getVLabel(it3->first)) != patternlabel->end())
						continue;
					leftchildlabel->insert(Graph->getVLabel(it3->first));
				}


				for(auto it3 = rightchildlabel->begin(); it3 != rightchildlabel->end(); it3++){
					if(patternlabel->find(*it3)!=patternlabel->end())
						continue;
					patternlabel->insert(*it3);
					for(auto it4 = rightchildlabel->begin(); it4 != rightchildlabel->end(); it4++){
						if(patternlabel->find(*it4)!=patternlabel->end())
							continue;
						if (*it3 <= *it4)
							continue;
						labels[4] = *it3;
						labels[5] = *it4;
						patternlabel->insert(*it4);

						for (auto it5 = leftchildlabel->begin(); it5 != leftchildlabel->end(); it5++)
						{
							for (auto it6 = leftchildlabel->begin(); it6 != leftchildlabel->end(); it6++)
							{
								if (*it5 <= *it6)
									continue;
								labels[2] = *it5;
								labels[3] = *it6;
								constructBF(BF, Graph, center, -3, nullptr, nullptr, labels, weight, threshold); // compute the weighted value
							}
						}

						patternlabel->erase(*it4);
					}
					patternlabel->erase(*it3);
				}
				patternlabel->erase(labels[0]);
			}
			patternlabel->erase(labels[1]);
		}

		delete patternlabel;
		delete leftchildlabel;
		delete rightchildlabel;
		delete degree0;
		delete degree1;
		delete degree2;
		delete labels;
		return 1;
	}

	if (type == -1) // initial, undirected tree pattern
	{
		Degree *degree0 = new Degree;
		Degree *degree1 = new Degree;
		Degree *degree2 = new Degree;
		// int degree = 0;

		int tempdegree;
		PatternLabel *patternlabel = new PatternLabel;
		patternlabel->insert(Graph->getVLabel(center));
		for (auto it = Graph->getOutEdge()[center].begin(); it != Graph->getOutEdge()[center].end(); it++)
		{ // compute the type of children vertices
			if (patternlabel->find(Graph->getVLabel(it->first)) == patternlabel->end())
			{
				// degree++;
				patternlabel->insert(Graph->getVLabel(it->first));
				for (auto it1 = Graph->getOutEdge()[it->first].begin(); it1 != Graph->getOutEdge()[it->first].end(); it1++)
				{
					if (patternlabel->find(Graph->getVLabel(it1->first)) == patternlabel->end())
					{
						patternlabel->insert(Graph->getVLabel(it1->first));
					}
					else
						continue;
				}

				for (auto it2 = Graph->getInVertex()[it->first].begin(); it2 != Graph->getInVertex()[it->first].end(); it2++)
				{
					if (patternlabel->find(Graph->getVLabel(it2->first)) == patternlabel->end())
					{
						patternlabel->insert(Graph->getVLabel(it2->first));
					}
					else
						continue;
				}

				// degree0
				if (patternlabel->size() == 2)
					degree0->emplace(it->first, 0);
				// degree0[it->first]=Graph->getVLabel(it->first);

				// degree1
				if (patternlabel->size() == 3)
					degree1->emplace(it->first, 1);
				// degree1[it->first]=Graph->getVLabel(it->first);

				// degree2
				if (patternlabel->size() > 3)
					degree2->emplace(it->first, patternlabel->size()-2);
				// degree2[it->first]=Graph->getVLabel(it->first);
			}
			else
				continue;

			patternlabel->clear();
			patternlabel->insert(Graph->getVLabel(center));
		}

		for (auto it = Graph->getInVertex()[center].begin(); it != Graph->getInVertex()[center].end(); it++)
		{ // for parents
			if (patternlabel->find(Graph->getVLabel(it->first)) == patternlabel->end())
			{
				// degree++;
				patternlabel->insert(Graph->getVLabel(it->first));
				for (auto it1 = Graph->getOutEdge()[it->first].begin(); it1 != Graph->getOutEdge()[it->first].end(); it1++)
				{
					if (patternlabel->find(Graph->getVLabel(it1->first)) == patternlabel->end())
					{
						patternlabel->insert(Graph->getVLabel(it1->first));
					}
					else
						continue;
				}

				for (auto it2 = Graph->getInVertex()[it->first].begin(); it2 != Graph->getInVertex()[it->first].end(); it2++)
				{
					if (patternlabel->find(Graph->getVLabel(it2->first)) == patternlabel->end())
					{
						patternlabel->insert(Graph->getVLabel(it2->first));
					}
					else
						continue;
				}

				// degree0
				if (patternlabel->size() == 2)
					degree0->emplace(it->first, 0);
				// degree0[it->first]=Graph->getVLabel(it->first);

				// degree1
				if (patternlabel->size() == 3)
					degree1->emplace(it->first, 1);
				// degree1[it->first]=Graph->getVLabel(it->first);

				// degree2
				if (patternlabel->size() > 3)
					degree2->emplace(it->first, patternlabel->size()-2);
				// degree2[it->first]=Graph->getVLabel(it->first);
			}
			else
				continue;

			patternlabel->clear();
			patternlabel->insert(Graph->getVLabel(center));
		}

		/* output the labels for debug
		for(auto it6 = degree0->begin(); it6!=degree0->end();it6++){
			cout<<"degree 0 labels:"<< it6->second<<endl;
		}
			for(auto it6 = degree1->begin(); it6!=degree1->end();it6++){
			cout<<"degree 1 labels:"<< it6->second<<endl;
		}
			for(auto it6 = degree2->begin(); it6!=degree2->end();it6++){
			cout<<"degree 2 labels:"<< it6->second<<endl;
		}*/
		if(degree2->size()>threshold){
			//cout<<"Ball with center "<<center<< " have >=20 neighbors that have more than two children!"<<endl;
			return 0;
		}

		if ((degree0->size() > 0) && (degree1->size() > 0))
			constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight, threshold);

		if ((degree0->size() > 0) && (degree2->size() > 0))
			constructBF(BF, Graph, center, 7, degree2, degree0, nullptr, weight, threshold);

		if ((degree1->size() > 0) && (degree2->size() > 0))
			constructBF(BF, Graph, center, 9, degree2, degree1, nullptr, weight, threshold);

		if (degree2->size() > 1)
			constructBF(BF, Graph, center, 10, degree2, degree2, nullptr, weight, threshold);

		delete degree0;
		delete degree1;
		delete degree2;
		delete patternlabel;
		return 1;
	}
}


void BF(BF_value **QueryBF, BloomFilter *ballBF, int querysize, VertexID center, int pointer, DIGRAPH<VertexLabel, EdgeLabel>* query, DIGRAPH<VertexLabel, EdgeLabel>* graph, int* result_BF)
{
	bool flag1 = false;
	bool flag2;
	
	char buffer[ENTRY_VALUE_LEN+1];
	size_t len = ENTRY_VALUE_LEN+1;
	for (int i = 0; i < querysize; i++)
	{
		if (query->getVLabel(query->matrix[i]) != graph->getVLabel(center))
			continue;
		flag2 = true;
		for (auto it1 = QueryBF[i]->begin(); it1 != QueryBF[i]->end(); it1++)
		{
			for (auto it = it1->second.begin(); it != it1->second.end(); it++)
			{
				if (*it == -1)
					continue;
				//sprintf(buffer, "%ld", *it);
				LongToChar(buffer, *it, ENTRY_VALUE_LEN, 10);
				if (!ballBF->possiblyContains((uint8_t*)buffer, len))
				{				
					flag2 = false;
					break;
				}
			}
			if (!flag2)
				break;
		}
		if (!flag2)
		{
			continue;
		}
		else
		{
			flag1 = true;
			break;
		}
	}

	if (!flag1)
	{
		result_BF[pointer] = 0;
	}	
}

void LongToChar(char *C_data, long L_data, int len, int base){
    long tempdata = L_data;
    int templen = len-1;
    while(tempdata>0){
        C_data[templen]=('0'+(tempdata%base));
        tempdata=tempdata/base;
        templen--;
    }

    while(templen>=0){
        C_data[templen] = '\0';
        templen--;
    }
}



void MultiServers(int times, int servernum, unordered_map<int, double> & BFTime, int graphsize, int *BallCount, ofstream & OutFile)
{
		int Num_Computed_Balls = 0, Num_positive=0, Num_negative=0;
	for(int i =0; i<graphsize;i++){
		if(BallCount[i]!=0)
			Num_Computed_Balls++;
		if(BallCount[i]==1)
			Num_positive++;
		if(BallCount[i]==-1)
			Num_negative++;
	}

	cout<< "time: " <<times <<" servernum: "<<servernum<<" Num: "<<Num_Computed_Balls<<" Positive: "<<Num_positive<<endl;

	int* ComputedBalls = new int[Num_Computed_Balls];
	int* Positive = new int[Num_positive];
	int* Negative = new int[Num_negative];
	int* Mark = new int[Num_Computed_Balls];
	int* Mark1 = new int[Num_positive];
	int* Mark0 = new int[Num_negative];

	int count = 0,count1 = 0, count2 = 0, position;
	for(int i =0; i<graphsize;i++){
		if(BallCount[i]==1){
			Positive[count1] = i;
			count1++;
		}		
		if(BallCount[i]==-1){
			Negative[count2] = i;
			count2++;
		}	
		if(BallCount[i]!=0){
			ComputedBalls[count] = i;
			count++;
		}			
	}


	double bfTime = 0, tempBF = 0, BF_Avg=0;


	int *** ballassign = new int**[times];
	int *** order = new int**[times];
	int *earlyposition = new int[servernum];
	int size = Num_Computed_Balls/servernum + 1;
	int mod1 = Num_Computed_Balls%servernum;
	int mod2 = Num_positive%servernum;
	for(int i=0;i< times;i++){
		ballassign[i] = new int*[servernum];
		order[i] = new int*[servernum];
		for(int j=0; j< servernum;j++){
			ballassign[i][j] = new int[size];
			for(int k = 0;k<size;k++)
				ballassign[i][j][k] = -1;
			order[i][j] = new int[2*size];
			for(int k = 0;k<2*size;k++)
				order[i][j][k] = -1;
		}
	}

	

	// conduct multiple times
	for(int i=0; i<times;i++){
		//baseline
		count = 0;
		for(int j=0;j<Num_Computed_Balls;j++){
			Mark[j] = 0;
		}

		//cout <<"hello2!"<< " mod1:" <<mod1<< " mod2:"<<mod2<<" size:"<<size<<endl;

		for(int j=0; j<servernum;j++){
			count = 0;
			while(1){				
				if((count==size)&&(j<mod1))
					break;
				if((count==(size-1))&&(j>=mod1))
					break;	
				position=rand()%Num_Computed_Balls;
				//cout<<position<<"	"<<Mark[position]<<"	"<<count<<endl;
				while(Mark[position]!=0){
					position = (position+1)%Num_Computed_Balls;
				}
				//cout<<position<<"	"<<Mark[position]<<"	"<<count<<endl;
				Mark[position] = 1;
				ballassign[i][j][count] = ComputedBalls[position];
				
				count++;
			}
		}

		

		//compute the baseline time
		bfTime =0;
		for(int j = 0; j< servernum;j++){
			tempBF = 0;
			for(int k = 0; k< size; k++){
				if(ballassign[i][j][k]!=-1){
					tempBF += BFTime[ballassign[i][j][k]];
				}
			}
			if(bfTime<tempBF)
				bfTime = tempBF;
		}


		//cout<<"The baselineTime: " << baselineTime <<"s !	The SOP Time: " << orderTime <<" s!"<<endl;
		//cout<<"The TwigTime: "<< twigTime<<endl;
		//cout<<"The TwigBallReadTime: "<< ReadTime<<endl;

		BF_Avg += bfTime;

	}

	cout << "******************************************************************" << endl;
	cout << "m = " << servernum << endl;
	cout << "The average BFTime: " << (BF_Avg/times)*1000 << " ms!"<<endl;


	OutFile << servernum << endl;
	OutFile << (BF_Avg/times)*1000 << endl;



	for(int i=0;i< times;i++){		
		for(int j=0; j< servernum;j++)
			delete[] ballassign[i][j];
		delete[] ballassign[i];
	}


	delete[] Mark;
	delete[] Mark0;
	delete[] Mark1;
	delete[] Positive;
	delete[] Negative;
	delete[] ComputedBalls;
	delete[] earlyposition;


}