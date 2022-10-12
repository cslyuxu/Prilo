#pragma once

#ifndef PPMATCH_H_
#define PPMATCH_H_

#include <iostream>
#include "DIGraph.h"
#include "DataRead.h"
#include <time.h>
#include <sstream>
#include "cgbe.h"
#include <stdio.h>
typedef mpz_t Ciphertext;

// BloomFilter
//#include <bf/all.hpp>
//using namespace bf;
#include "BF/BloomFilter.h"

using namespace std;

int cmp(const void *a, const void *b)
{

	return *(int *)a - *(int *)b; // ascending

	// return *(int *)b - *(int *)a; //descending
}

template <class VLabelType, class ELabelType>
class PPMATCH
{
public:
	int *result, *result_EncSSim, *OneIterNL, *result_Path, *EncSpecial, *result_BF, *result_Twig, *result_GH, *result_PPMatch;
	int query_size, graph_size, query_selected_label_size, GlobalCountForBuildingReplaceTable, hoplength, pathlength, LabelTreeLength, twiglength, labelsize, mode;
	double exact_solution_num, OneIter_num, TwoIter_num, NeighborLabel_num, Path_num, BF_num, portion, Twig_num;
	int* BallCount;
	string OutFileName;
	//string PathFileName;
	//string Index1, Index2;
	double StrongSimulation_Time, Precompute_Q_Time, Enc_OneIter_Time, Enc_TwoIter_Time, NeighborLabel_Time, Path_Time, Enc_Time, Ball_Time, Decrypt_OneIter_Time, Decrypt_TwoIter_Time, Decrypt_NeighborLabel, Decrypt_Path, Index_Time, Twig_Time, Decrypt_Twig, BF_Time, Temp_Time, Decrypt_GH;
	int NL_Improved, Path_Improved, Twig_Improved, BF_Improved;
	DIGRAPH<VLabelType, ELabelType> *query;
	DIGRAPH<VLabelType, ELabelType> *graph;

	//
	typedef unordered_map<int, Ciphertext> P_Column_Value;
	typedef unordered_map<int, P_Column_Value> P_Row;

	// For Path
	typedef unordered_map<int, Ciphertext> EncPath;
	// typedef unordered_map<int, EncPath> Path_Index;
	typedef unordered_map<VertexID, EncPath> PathIndex;
	typedef unordered_map<int, PathIndex> Path_Index;

	// typedef unordered_map<VLabelType, int> Path_label;
	// typedef unordered_map<int, Path_label> Path_Label;
	typedef unordered_map<int, VLabelType> PathLabel;
	typedef unordered_map<int, PathLabel> Path_Label;

	// typedef unordered_map<int, VLabelType> Path_num;
	// typedef unordered_map<int, Path_num> Path_Num;
	// typedef unordered_map<VertexID, Path_Num> PathNum;
	typedef unordered_map<VLabelType, int> PathNum;
	typedef unordered_map<int, PathNum> Path_Num;

	//For GH
	typedef unordered_map<int, int> Vertex_Map;
	typedef unordered_map<int, Vertex_Map> GH;

	// For NL
	typedef unordered_map<VLabelType, int> LABEL;
	typedef unordered_map<VertexID, LABEL> Index;
	typedef unordered_map<int, Index> B_Index; //<hop, Vertex, Label, 1=exists>

	// For Two_Iter
	typedef unordered_map<int, Ciphertext> Ciphertext_Match;
	typedef unordered_map<int, Ciphertext_Match> Replace_Table; // <# for the vertex, # for the case, Ciphertext>

	// For BF
	typedef unordered_set<long> BFvalue;
	typedef unordered_map<int, BFvalue> BF_value;

	// typedef unordered_map<VertexID, int> MapMatrix;				//record the map for the matrix
	typedef unordered_map<VertexID, VLabelType> VLabels;

	// For nodes of tree pattern
	typedef unordered_map<VertexID, int> Degree;
	typedef unordered_set<VLabelType> PatternLabel;

	/*For query with random labels*/
	// PPMATCH(Read<VLabelType, ELabelType> *, Read<VLabelType, ELabelType> *, string);

	/*For query with labels*/
	PPMATCH(DataRead<VLabelType, ELabelType> *, DataRead<VLabelType, ELabelType> *, string, int, int, int, int, int, int);
	~PPMATCH();
	void Match();
	void Exact(int **, int ***, int **, int, int, int);
	void Enc_OneIter(int **, Ciphertext **, P_Row &, int, int, CGBE *, Ciphertext &, Ciphertext &, int);
	void Enc_TwoIter(int **, Ciphertext **, Ciphertext **, P_Row &, int, int, CGBE *, Ciphertext &, Ciphertext &, Replace_Table &, Replace_Table &, Replace_Table &, Replace_Table &, int);
	void Enc_TwoIter_Random(int **, Ciphertext **, Ciphertext **, P_Row &, int, int, CGBE *, Ciphertext &, Ciphertext &, Replace_Table &, Replace_Table &, Replace_Table &, Replace_Table &, int, Ciphertext &, Ciphertext &, double);
	void Enc_GH(int **, Ciphertext **, int, GH &, Vertex_Map &, int, int, CGBE *, Ciphertext &, Ciphertext &, int, int, bool &, int, int, int &, clock_t &);
	void Enc_NeighborLabel(int **, Ciphertext ***, Ciphertext ***, P_Row &, P_Row &, int, int, CGBE *, B_Index &, B_Index &, VLabelType *, Ciphertext &, Ciphertext &, int, VertexID, VertexID *);
	void NLTest(int **, int ***, int ***, int **, int, int, CGBE *, B_Index &G, B_Index &, VLabelType *, int, VertexID, VertexID *);
	void Enc_Path(Path_Label &, Path_Num &, Path_Index &, Path_Index &, int, int **, VertexID *, int, CGBE *, double &, int, VertexID, int, int, Ciphertext &);
	void Enc_Twig(Path_Label &, Path_Num &, Path_Index &, Path_Index &, int, int **, VertexID *, int, CGBE *, double &, unordered_map<int, double>&, unordered_map<int, double>&, int, VertexID, int, int, Ciphertext &);
	void BF(BF_value **, BloomFilter *, int, VertexID, int);
	void BFTest(BF_value **, BF_value *, int, VertexID, int);
	void constructBF(BF_value *, DIGRAPH<VLabelType, ELabelType> *, VertexID, int, Degree *, Degree *, VLabelType *, long *, int&);
	// void Enc_Path_Center(Path_Label &, Path_Label &, Path_Index &, int **, int, int, VertexID *, int, CGBE *);
	void PathMatch(DIGRAPH<VLabelType, ELabelType> *, Path_Num &, unordered_map<VertexID, unordered_map<int, int>> &, unordered_map<VertexID, unordered_map<int, int>> &, VLabelType *, int, int, int, VertexID, VertexID, int);
	void TwigMatchAll(DIGRAPH<VLabelType, ELabelType> *, Path_Num &, unordered_map<VertexID, unordered_map<int, int>> &, unordered_map<VertexID, unordered_map<int, int>> &, VLabelType *, VLabelType *, int, int, int, VertexID, VertexID, int);
	// bool PathMatch(DIGRAPH<VLabelType, ELabelType>*, int *, int, int, int, VertexID);
	// bool PathMatch_Center(DIGRAPH<VLabelType, ELabelType>*, int *, int, int, int, int, VertexID);
	void Q_K_HOP(int ***, Ciphertext ***, Ciphertext ***, int ***, int ***, VLabelType *, CGBE *);
	void Q_Replace_Table_Child(int ***, Ciphertext **, Replace_Table &, Replace_Table &, CGBE *, Ciphertext &, Ciphertext &, int, int, int, Ciphertext &);
	void Q_Replace_Table_Parent(int ***, Ciphertext **, Replace_Table &, Replace_Table &, CGBE *, Ciphertext &, Ciphertext &, int, int, int, Ciphertext &);
	void BuildPathIndex(Path_Label &, Path_Num &, Path_Index &, Path_Index &, CGBE *, int, int, VLabelType *);
	void BuildTwigIndex(Path_Label &, Path_Num &, Path_Index &, Path_Index &, CGBE *, int, int, VLabelType *);
	void FindQueryPath(Path_Label &, Path_Num &, Path_Index &, Path_Index &, CGBE *, int, int, VLabelType *, VertexID, VertexID, int);
	void FindQueryTwig(Path_Label &, Path_Num &, Path_Index &, Path_Index &, CGBE *, int, int, VLabelType *, VLabelType *, VertexID, VertexID, int);
	int CantorExpansion(int, int *);
	void CantorExpansionDecode(int, int *, int);
	int ConstructBallMatrix(int **, VertexID *, int);
	int MaxDegree(int **, int);
	void ConstructBallIndex(B_Index &, B_Index &, int **, VertexID *, int);
	void ConstructNLIndex(B_Index &, B_Index &, DIGRAPH<VLabelType, ELabelType> *, VertexID, VertexID, int, unordered_map<VertexID, int> &);
	void LoadNLIndex(B_Index &, B_Index &, string, string);
	void MultiServers(int, int, unordered_map<int, double> &, unordered_map<int, double> &, unordered_map<int, double> &, ofstream &);
};

/*For query with random labels*/
/*
template <class VLabelType, class ELabelType>
PPMATCH<VLabelType, ELabelType>::PPMATCH(Read<VLabelType, ELabelType> *Q, Read<VLabelType, ELabelType> *G, string OutFile)
{
*/

/*For query with labels*/

template <class VLabelType, class ELabelType>
PPMATCH<VLabelType, ELabelType>::PPMATCH(DataRead<VLabelType, ELabelType> *Q, DataRead<VLabelType, ELabelType> *G, string OutFile, int khop, int pathlength, int label_size, int LabelTreeLength, int twiglength, int mode)
{
	this->portion = 0.25; // Random TwoIter
	this->query = Q->graph_ptr;
	this->graph = G->graph_ptr;
	this->graph_size = this->graph->getVcnt();
	this->query_size = this->query->getVcnt();
	//sec. 4.3	
	this->BallCount = new int[this->graph_size];
	for(int i = 0; i<this->graph_size; i++)
		this->BallCount[i] = 0;
	this->query_selected_label_size = this->query->VLabelSet.size();
	// if (this->query_selected_label_size > K_LABEL)
	//	this->query_selected_label_size = K_LABEL;
	this->exact_solution_num = 0;
	this->OneIter_num = 0;
	this->TwoIter_num = 0;
	this->NeighborLabel_num = 0;
	this->BF_num = 0;
	this->Path_num = 0;
	this->Twig_num = 0;
	this->result = new int[graph_size];
	this->result_EncSSim = new int[graph_size];
	this->result_PPMatch = new int [graph_size];
	this->OneIterNL = new int[graph_size];
	this->EncSpecial = new int[graph_size];
	this->result_Path = new int[graph_size];
	this->result_BF = new int[graph_size];
	this->result_Twig = new int[graph_size];
	this->result_GH = new int[graph_size];
	this->StrongSimulation_Time = 0;
	this->Precompute_Q_Time = 0;
	this->Enc_OneIter_Time = 0;
	this->Enc_TwoIter_Time = 0;
	this->NeighborLabel_Time = 0;
	this->Path_Time = 0;
	this->Twig_Time = 0;
	this->Temp_Time = 0;
	this->BF_Time = 0;
	this->Enc_Time = 0;
	this->Ball_Time = 0;
	this->GlobalCountForBuildingReplaceTable = 0;
	this->Decrypt_OneIter_Time = 0;
	this->Decrypt_TwoIter_Time = 0;
	this->Decrypt_NeighborLabel = 0;
	this->Decrypt_Path = 0;
	this->Decrypt_Twig = 0;
	this->Decrypt_GH = 0;
	this->Index_Time = 0;
	this->OutFileName = OutFile;
	this->NL_Improved = 0;
	this->Path_Improved = 0;
	this->Twig_Improved = 0;
	this->BF_Improved = 0;
	this->hoplength = khop;
	this->labelsize = label_size;
	this->mode = mode;
	this->pathlength = pathlength;
	this->LabelTreeLength = LabelTreeLength;
	this->twiglength = twiglength;
	cout << "*************************************" << endl;
	cout << "The Strong Simulation has been built!" << endl;
	cout << "*************************************" << endl;
	cout << endl;
}

template <class VLabelType, class ELabelType>
PPMATCH<VLabelType, ELabelType>::~PPMATCH()
{
	delete[] this->result;
	delete[] this->result_EncSSim;
	delete[] this->result_PPMatch;
	delete[] this->OneIterNL;
	delete[] this->EncSpecial;
	delete[] this->result_Path;
	delete[] this->result_BF;
	delete[] this->result_Twig;
	delete[] this->result_GH;
	delete[] this->BallCount;
	cout << "The ptr in Strong Simulation has been deleted!" << endl;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Match()
{
	clock_t startTime, endTime, startTime1, endTime1, startTimeTest, endTimeTest;
	clock_t startOne, endOne, startTwo, endTwo, startThree, endThree, startFour, endFour, startFive, endFive, startSix, endSix, startSeven, EndSeven;


	/*Output the data*/
	ofstream OutFile(OutFileName);
	string RuntimeName = OutFileName + "-runtime";
	//string PrecisionName = OutFileName + "-Accuracy";
	//string NumberName = OutFileName + "-Number";
	ofstream Runtime(RuntimeName);
	//ofstream Precision(PrecisionName);
	//ofstream Number(NumberName);
	Runtime << "# |V_Q| = " << query_size << endl;
	Runtime << "# ball size, exact, 2-iter, NL, Path" << endl;
	//Precision << "# name, dataset" << endl;
	//Number << "# name, dataset" << endl;

	char buffer[ENTRY_VALUE_LEN+1];

	///////////////////////////////////////
	/////////////   Client   //////////////
	///////////////////////////////////////
	// compute the diameter of Query
	this->query->getDiameter();
	int dia = this->query->diameter;
	cout << "The diameter of the Query is " << dia << endl;

	/*construct the matrix $\overline{M_Q}$ and encrypted one for the query*/
	startTime = clock();
	// Ciphertext ***EM_Q = new Ciphertext **[K_HOP];
	// used for one-hop neighbor
	int ***Q = new int **[hoplength];
	Q[0] = new int *[query_size];
	int ***Q_Child = new int **[hoplength];
	int ***Q_Parent = new int **[hoplength];
	Ciphertext **EM_Q = new Ciphertext *[query_size];
	Ciphertext **EM_Q_Two = new Ciphertext *[query_size];
	Ciphertext ***EM_Q_Child = new Ciphertext **[hoplength];
	Ciphertext ***EM_Q_Parent = new Ciphertext **[hoplength];
	CGBE *cgbe = new CGBE();

	Ciphertext cipher_one, cipher_zero; // for SP to compute in ciphertext
	mpz_init(cipher_one);
	mpz_init(cipher_zero);
	cgbe->setvalue(cipher_one, 1);
	cgbe->setvalue(cipher_zero, cgbe->encoding);
	cgbe->encrypt(cipher_one, cipher_one);
	cgbe->encrypt(cipher_zero, cipher_zero);

	Ciphertext cipher_one_Two;
	Ciphertext cipher_one_Four_Vq;
	mpz_init(cipher_one_Two);
	mpz_init(cipher_one_Four_Vq);
	cgbe->mul(cipher_one_Two, cipher_one, cipher_one);
	cgbe->mul(cipher_one_Two, cipher_one_Two, cipher_one_Two);
	cgbe->setvalue(cipher_one_Four_Vq, 1);
	for (int i = 0; i < query_size; i++)
	{
		cgbe->mul(cipher_one_Four_Vq, cipher_one_Four_Vq, cipher_one_Two);
	}
	cgbe->mul(cipher_one_Two, cipher_one, cipher_one);

	/*intialize matrix M_Q*/
	for (int i = 0; i < query_size; i++)
	{
		Q[0][i] = new int[query_size];
		EM_Q[i] = new Ciphertext[query_size];
		EM_Q_Two[i] = new Ciphertext[query_size];
		for (int j = 0; j < query_size; j++)
		{
			Q[0][i][j] = 1;
			mpz_init(EM_Q[i][j]);
			mpz_init(EM_Q_Two[i][j]);
			cgbe->setvalue(EM_Q[i][j], 1);
			if (this->query->isEdge(this->query->matrix[i], this->query->matrix[j]))
			{
				Q[0][i][j] = 0;
				cgbe->setvalue(EM_Q[i][j], cgbe->encoding);
			}
			cout << Q[0][i][j] << " ";
			cgbe->encrypt(EM_Q[i][j], EM_Q[i][j]);
			cgbe->setvalue(EM_Q_Two[i][j], EM_Q[i][j]);
			cgbe->mul(EM_Q_Two[i][j], EM_Q_Two[i][j], EM_Q_Two[i][j]);
		}
		cout << endl;
	}

	/*preprocessing for Q*/
	/*construct K_HOP indices for Q*/
	// NL

	VLabelType *Column = new VLabelType[this->query_selected_label_size]; // label for column	
	Q_K_HOP(Q, EM_Q_Child, EM_Q_Parent, Q_Child, Q_Parent, Column, cgbe);

	// choose label with maximum (optional) num
	VLabelType LabelofQ;
	int MaxLabelNum = 0;
	this->graph->CountLabelNum(this->query);
	for (auto it = this->graph->getLabelCount().begin(); it != this->graph->getLabelCount().end(); it++)
	{
		if (it->second > MaxLabelNum)
		{
			MaxLabelNum = it->second;
			LabelofQ = it->first;
		}
	}
	cout << "The Label " << LabelofQ << " has a maximum number of " << MaxLabelNum << endl;

	// For testing the Label choosing
	/*for(auto it = this->graph->getLabelCount().begin(); it != this->graph->getLabelCount().end(); it++){
		cout << it->first << "	" << it->second << endl;
	}*/

	endTime = clock();
	this->Enc_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

	startTime = clock();
	// load the Path_Index
	Path_Index PathIndex1, PathIndex2;
	Path_Label PathLabel;
	Path_Num PathNum;
	for (int i = 3; i <= this->pathlength; i++)
		BuildPathIndex(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, i, 0, nullptr);

	Path_Index TwigIndex1, TwigIndex2;
	Path_Label TwigLabel;
	Path_Num TwigNum;
	for (int i = 3; i <= this->twiglength; i++)
		BuildTwigIndex(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, i, 0, nullptr);

	endTime = clock();
	this->Enc_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	/*
	//check the index
	cout << "#######################" << endl;
	cout << "#######################" << endl;
	for(auto it = TwigLabel[3].begin();it != TwigLabel[3].end(); it++){
		cout<<it->first<<": "<<it->second<<endl;
	}

	for(auto it = TwigIndex2[3].begin();it != TwigIndex2[3].end(); it++){
		cout << "Vertex " << it->first <<":" << endl;
		for(auto it2 = TwigIndex2[3][it->first].begin(); it2 != TwigIndex2[3][it->first].end(); it2++){
			cgbe->decrypt(it2->second, it2->second);
			if(cgbe->isZero(it2->second))
				cout << "Twig " << it2->first << ": 0" <<endl;
		}
		for(auto it2 = PathIndex2[3][it->first].begin(); it2 != PathIndex2[3][it->first].end(); it2++){
			cgbe->decrypt(it2->second, it2->second);
			if(cgbe->isZero(it2->second))
				cout << "Path " << it2->first << ": 0" <<endl;
		}

	}*/

	int MaximumBFNum = 0,tempBFNum = 0;

	startTime = clock();
	clock_t starttemp, endtemp;
	starttemp = clock();
	// precompute the tree pattern for bloomfilter
	BF_value **QueryBF = new BF_value *[query_size];
	unordered_set<long> templong = {-1};
	long *weight = new long[6];
	weight[0] = 1;
	for (int i = 1; i < 6; i++)
	{
		weight[i] = weight[i - 1] * this->labelsize;
		cout << weight[i] << endl;
	}
	for (int i = 0; i < query_size; i++)
	{
		QueryBF[i] = new BF_value;
		//output the labels of Q's vertices
		//cout << "The " << i << "th vertex of Q!" <<this->query->getVLabel(this->query->matrix[i])<< endl;
		for (int j = 1; j < 5; j++)
		{
			QueryBF[i]->emplace(j, templong);
		}
		constructBF(QueryBF[i], this->query, this->query->matrix[i], -1, nullptr, nullptr, nullptr, weight, tempBFNum);
		if(tempBFNum > MaximumBFNum)
			MaximumBFNum = tempBFNum;
	}
	endtemp = clock();
	cout<<"The encodings are constructed in " << (double)(endtemp - starttemp) / CLOCKS_PER_SEC<<" seconds!"<<endl;
	
	
	/*quick test for the maximum num of encodings of query*/
	//cout<<"The maximum number of encodings is " << MaximumBFNum <<endl;
	//return;
	



	// Replace_Table
	Replace_Table Origin_Child, Replace_Child, Origin_Parent, Replace_Parent;
	Q_Replace_Table_Child(Q, EM_Q, Origin_Child, Replace_Child, cgbe, cipher_one, cipher_zero, 0, -1, 0, cipher_zero);
	Q_Replace_Table_Parent(Q, EM_Q, Origin_Parent, Replace_Parent, cgbe, cipher_one, cipher_zero, 0, -1, 0, cipher_zero);

	endTime = clock();
	this->Enc_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	cout << endl;
	cout << "The runtime for preprocessing Q on Client is " << Enc_Time << endl;
	cout << "\n*************************************" << endl;

	// OutFile << "The runtime for preprocessing Q on Client is " << Enc_Time << endl;
	// OutFile << "\n*************************************" << endl;

	/////////////////
	/////////////////
	////////////////
	// Here
	// OutFile << "Ball#"<< "	"<< "BallSize"<< "	"<< "runtime"<< "	"<< "Time-Ball" "	"<< "Time-Index" << "	"<< "Time-Neighbor"<< "	"<< "Time-TwoIter"<< "	"<< "Time-Path" << "	"<< "Time-Exact" << endl;
	/////////////////
	/////////////////
	/////////////////

	///////////////////////////////////////
	/////////////   Server   //////////////
	///////////////////////////////////////

	// load the graph index
	/**/
	/**/
	/**/

	// The new part

	/**/
	/**/
	/**/
	//startTimeTest = clock();
	//B_Index Graph_Index, Graph_Index_Reverse;
	//LoadNLIndex(Graph_Index, Graph_Index_Reverse, Index1, Index2);
	//endTimeTest = clock();
	//cout << "Indexes are loaded in " << (double)(endTimeTest - startTimeTest) / CLOCKS_PER_SEC << " seconds!" << endl;

	// pointer for short test
	int pointer = 0;

	// Num_Ball_Computed
	int Num_Ball_Computed = 0, Num_Ball_Ignored = 0;

	int bf_size;







	//sec. 4.3
	unordered_map<int, double> BallTime;
	unordered_map<int, double> TwigTime;
	unordered_map<int, double> TwigBallReadTime;
	BallTime.clear();
	TwigTime.clear();
	TwigBallReadTime.clear();


	//count the maximum degree;
	//int maxdegree = 0, tempdegree;




	double each_ball_time;
	// build each candidate
	for (typename VLabels::iterator it1 = this->graph->getVLabel().begin();
		 it1 != this->graph->getVLabel().end(); it1++)
	{
		//if(pointer>10000){
		//	pointer++;
		//	continue;
		//}

		startOne = clock();
		// No need to short test now
		//if (pointer != 195354){
		//	pointer++;
		//	continue;
		//}
		//	break;

		// if (Num_Ball_Computed == 1)
		//	break;

		VertexID s = it1->first, t;

		/*if the center of the ball's label is not in the label set of query, then continue*/
		if (!(this->query->isLabel(this->graph->getVLabel(s))))
		{
			this->result_EncSSim[pointer] = -1;
			this->OneIterNL[pointer] = -1;
			this->EncSpecial[pointer] = -1;
			this->result_Path[pointer] = -1;
			this->result[pointer] = -1;
			//this->result_BF[pointer] = -1;
			this->result_Twig[pointer] = -1;
			this->result_GH[pointer] = -1;

			pointer++;
			continue;
		}

		////////////////////////////////////////////////////////////////////
		/////		For GraphHomomorphism: Choose arbitrary label in Q /////
		////////////////////////////////////////////////////////////////////
		if(mode == 1){
			if (this->graph->getVLabel(s) != LabelofQ)
			{
				this->result_EncSSim[pointer] = -1;
				this->OneIterNL[pointer] = -1;
				this->EncSpecial[pointer] = -1;
				this->result_Path[pointer] = -1;
				this->result[pointer] = -1;
				//this->result_BF[pointer] = -1;
				this->result_Twig[pointer] = -1;
				this->result_GH[pointer] = -1;

				pointer++;
				continue;
			}
		}

		each_ball_time = 0;
		startTime = clock();
		startTwo = startTime;

		unordered_set<VertexID> s_neighbor;
		this->graph->getDNeighbor(s, dia, s_neighbor, query);
		int size_ball = s_neighbor.size();
		int edge_ball = 0;

		// If the ball is too large, continue
		/*if (size_ball > 3000)
		{
			OutFile << pointer << "	"<< size_ball << endl;
			pointer++;
			Num_Ball_Ignored++;
			continue;
		}*/

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

		/*construct the matrix M_B for the ball*/
		int **M_B = new int *[size_ball];
		edge_ball = ConstructBallMatrix(M_B, Matrix_ball, size_ball);


		/*
		//for counting the maximum degree
		tempdegree = MaxDegree(M_B, size_ball);
		if(tempdegree>maxdegree)
			maxdegree = tempdegree;		
		//free the space
		for (int i = 0; i < size_ball; i++)
		{
			delete[] M_B[i];
		}
		delete[] M_B;
		Matrix_Ball.clear();
		s_neighbor.clear();
		pointer++;
		continue;
		*/


		/////////////////////Pruning Techniques///////////////////////
		// If there is one row in P equals to zero vector, then prune//
		//////////////////////////////////////////////////////////////

		/*construct the matrix P for two methods*/
		int **answer = new int *[query_size];
		bool Prune_flag;
		int outposition;
		for (int i = 0; i < query_size; i++)
		{
			answer[i] = new int[size_ball];
			Prune_flag = true;
			/*intialize matrix P*/
			for (int j = 0; j < size_ball; j++)
			{
				answer[i][j] = 0;
				if (this->query->getVLabel(this->query->matrix[i]) == this->graph->getVLabel(Matrix_ball[j]))
				{
					answer[i][j] = 1;
					Prune_flag = false;
				}
			}
			if (Prune_flag)
			{
				outposition = i;
				break;
			}
		}

		if (Prune_flag)
		{
			this->result_EncSSim[pointer] = 0;
			this->OneIterNL[pointer] = 0;
			this->EncSpecial[pointer] = 0;
			this->result_Path[pointer] = 0;
			this->result[pointer] = 0;
			//this->result_BF[pointer] = 0;
			this->result_Twig[pointer] = 0;
			this->result_GH[pointer] = 0;

			/*free the space*/
			for (int i = 0; i < size_ball; i++)
			{
				delete[] M_B[i];
			}
			delete[] M_B;

			for (int i = 0; i < outposition; i++)
			{
				delete[] answer[i];
			}

			delete[] answer;
			delete[] Matrix_ball;
			Matrix_Ball.clear();
			s_neighbor.clear();
			pointer++;
			continue;
		}

		P_Row P_OneIter, P_NeighborLabel, P_TwoIter, P_Replace;
		GH P_GH;
		Vertex_Map GH_Temp;
		// NLtest
		int **M_P = new int *[query_size];

		for (int i = 0; i < query_size; i++)
		{
			M_P[i] = new int[size_ball];
			/*intialize matrix P*/
			for (int j = 0; j < size_ball; j++)
			{
				M_P[i][j] = 0;
				if (this->query->getVLabel(this->query->matrix[i]) == this->graph->getVLabel(Matrix_ball[j]))
				{
					answer[i][j] = 1;
					M_P[i][j] = 1;
					P_GH[i][j] = 1;
					mpz_init(P_OneIter[i][j]);
					mpz_init(P_NeighborLabel[i][j]);
					mpz_init(P_Replace[i][j]);
					mpz_init(P_TwoIter[i][j]);
					cgbe->setvalue(P_OneIter[i][j], 1);
					// cgbe->encrypt(P_OneIter[i][j], P_OneIter[i][j]);
					cgbe->setvalue(P_NeighborLabel[i][j], 1);
					// cgbe->encrypt(P_NeighborLabel[i][j], P_NeighborLabel[i][j]);
					cgbe->setvalue(P_TwoIter[i][j], 1);
				}
			}
		}

		// build the bloomfilter
		BF_value *BallBF = new BF_value;
		for (int i = 1; i < 5; i++)
		{
			BallBF->emplace(i, templong);
		}

		// use for the ball structure in graph structure
		DIGRAPH<VLabelType, ELabelType> *Ball = new DIGRAPH<VLabelType, ELabelType>;
		for (int i = 0; i < size_ball; i++)
		{
			Ball->insertVertex(Matrix_ball[i], this->graph->getVLabel(Matrix_ball[i]));
		}

		for (auto it3 = Ball->getVLabel().begin(); it3 != Ball->getVLabel().end(); it3++)
		{
			
			for(auto it4 = this->graph->getOutEdge()[it3->first].begin(); it4 != this->graph->getOutEdge()[it3->first].end();it4++){
				if(Ball->isVertex(it4->first))
					Ball->insertEdge(it3->first, it4->first, 0);
			}

			/*
			for (auto it4 = Ball->getVLabel().begin(); it4 != Ball->getVLabel().end(); it4++)
			{
				if (it3->first == it4->first)
					continue;
				if ((this->graph->isEdge(it3->first, it4->first)) && (!Ball->isEdge(it3->first, it4->first)))
					Ball->insertEdge(it3->first, it4->first, 0);
				// if ((this->graph->isEdge(it4->first, it3->first))&&(!Ball->isEdge(it4->first, it3->first)))
				// Ball->insertEdge(it4->first, it3->first, 0);
			}*/
		}

		// cout << "Ball Vcnt:" << Ball->getVcnt() << endl;
		// cout << "Ball Ecnt:" << Ball->getEcnt() << endl;
		// cout << "Ball center degree:" << Ball->getOutDegree(s) + Ball->getInDegree(s) << endl;

		/* BloomFilter is transmitted into SGX
		constructBF(BallBF, Ball, s, -1, nullptr, nullptr, nullptr, weight);

		bf_size = 0;
		for (auto it = BallBF->begin(); it != BallBF->end(); it++)
		{
			// cout << it->second.size() - 1 << endl;
			bf_size += (it->second.size() - 1);
		}


		//uint64_t vector_size = 315000000;//hold up 15 mil k,v // about 40MB
    	uint64_t vector_size = 460000000;//hold up 22 mil k,v // about 55MB
   		//uint64_t vector_size = 830000000;//hold up 40 mil k,v // about 110MB    
    	uint8_t numHashs = 23;
		BloomFilter *myBloomFilter;
		
		size_t len = ENTRY_VALUE_LEN+1;
		if (bf_size > 0)
		{	
			myBloomFilter = new BloomFilter(vector_size, numHashs);	
			
			for (auto it1 = BallBF->begin(); it1 != BallBF->end(); it1++)
			{
				for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
				{
					if (*it2 == -1) //-1 is used to initialize
						continue;
					//if(*it2==532632423)
					
					sprintf(buffer, "%ld", *it2);
					myBloomFilter->add((uint8_t*)buffer, len);
				}
			}
		}
		else
		{
			myBloomFilter = new BloomFilter(vector_size, numHashs);	
		}*/


		/* This first version of BloomFilter library in github
		bloom_filter *bf;
		if (bf_size > 0)
		{
			bf = new basic_bloom_filter(0.8, bf_size);
			for (auto it1 = BallBF->begin(); it1 != BallBF->end(); it1++)
			{
				for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
				{
					if (*it2 == -1) //-1 is used to initialize
						continue;
					if(*it2==532632423)
						cout<<"*it2: "<<*it2;
					bf->add(*it2);
				}
			}
		}
		else
		{
			bf = new basic_bloom_filter(0.8, 1);
		}*/
		endTime = clock();
		endTwo = endTime;
		each_ball_time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		this->Ball_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		//cout<<"BloomFilter for ball " << pointer <<" has been built!"<<endl;







		// Build the index for the NL technique

		/*Construct the indices for ball B */
		// B_Index Ball_Index, Ball_Index_Reverse;

		/* using the matrix multiplication */
		// ConstructBallIndex(Ball_Index, Ball_Index_Reverse, M_B, Matrix_ball, size_ball);

		/* using the graph traverse */
		// startTime1 = clock();
		// Ball = new DIGRAPH<VLabelType, ELabelType>;
		// Ball->ConstructInducedGraph(s_neighbor, this->graph);
		// ConstructNLIndex(Ball_Index, Ball_Index_Reverse, Ball, -1, -1, -1, Matrix_Ball);
		// endTime1 = clock();
		// this->Index_Time += (double)(endTime1 - startTime1) / CLOCKS_PER_SEC;

		this->result_EncSSim[pointer] = 1;
		this->OneIterNL[pointer] = 1;
		this->EncSpecial[pointer] = 1;
		this->result_Path[pointer] = 1;
		this->result[pointer] = 1;
		//this->result_BF[pointer] = 1;
		this->result_Twig[pointer] = 1;
		this->result_GH[pointer] = 1;

		/*Enc_One_Iter Algorithm*/
		// Enc_OneIter(M_B, EM_Q, P_OneIter, size_ball, pointer, cgbe, cipher_one, cipher_zero, Matrix_Ball[s]);

		/*Enc_Two_Iter Algorithm*/ // For strong simulation
		startFour = clock();
		
		endFour = clock();


		// BloomFilter
		startSix = clock();
		//if (bf_size > 0){
			///////BFTest(QueryBF, BallBF, query_size, s, pointer);
			//BF(QueryBF, myBloomFilter, query_size, s, pointer);
		//}
		endSix = clock();


		/*Enc_NieghborLabel Algorithm*/
		startThree = clock();
		// Enc_NeighborLabel(M_B, EM_Q_Child, EM_Q_Parent, P_NeighborLabel, P_Replace, size_ball, pointer, cgbe, Graph_Index, Graph_Index_Reverse, Column, cipher_one, cipher_zero, Matrix_Ball[s], s, Matrix_ball);
		endThree = clock();

		// NLTest(M_B, Q_Child, Q_Parent, M_P, size_ball, pointer, cgbe, Graph_Index, Graph_Index_Reverse, Column, Matrix_Ball[s], s, Matrix_ball);
		/*Enc_Path Algorithm*/
		startFive = clock();
		double PathTimeTemp = 0;
		for (int i = 3; i <= this->pathlength; i++)
			Enc_Path(PathLabel, PathNum, PathIndex1, PathIndex2, pointer, M_B, Matrix_ball, size_ball, cgbe, PathTimeTemp, Matrix_Ball[s], s, i, 50, cipher_one);
		endFive = clock();

		// Twig-Pruning
		double TwigTimeTemp = 0;
		startSeven = clock();
		TwigTime[pointer]=0;
		for (int i = 3; i <= this->pathlength; i++)
			Enc_Twig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, pointer, M_B, Matrix_ball, size_ball, cgbe, TwigTimeTemp, TwigTime, TwigBallReadTime, Matrix_Ball[s], s, i, 50, cipher_one);
		
		//Enc_Twig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, pointer, M_B, Matrix_ball, size_ball, cgbe, TwigTimeTemp, TwigTime, TwigBallReadTime, Matrix_Ball[s], s, twiglength, 50, cipher_one);
		EndSeven = clock();
		



		
		startTime = clock();

		/*Exact strong simulation Algorithm*/
		if(mode==2){
			//Exact(M_B, Q, answer, size_ball, pointer, Matrix_Ball[s]);
			if (this->portion == 1)
			{
					Enc_TwoIter(M_B, EM_Q, EM_Q_Two, P_TwoIter, size_ball, pointer, cgbe, cipher_one, cipher_zero, Origin_Child, Replace_Child, Origin_Parent, Replace_Parent, Matrix_Ball[s]);
				// cout << "hello portion 1 " << this->portion << endl;
			}
			else
					Enc_TwoIter_Random(M_B, EM_Q, EM_Q_Two, P_TwoIter, size_ball, pointer, cgbe, cipher_one, cipher_zero, Origin_Child, Replace_Child, Origin_Parent, Replace_Parent, Matrix_Ball[s], cipher_one_Two, cipher_one_Four_Vq, this->portion);
		

		}

		/*Enc_GraphHomomorphism Algorithm*/
		int count_GH = 0;
		bool flag_GH;		
		if(mode==1)
			Enc_GH(M_B, EM_Q, query_size, P_GH, GH_Temp, size_ball, pointer, cgbe, cipher_one, cipher_zero,  Matrix_Ball[s], -1, flag_GH, 0, 50, count_GH, startTime);
		endTime = clock();
		each_ball_time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		this->StrongSimulation_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		//For debug
		if((this->result_Twig[pointer]==1)&&(this->result_Path[pointer]==0)){
			cout<<"The " << pointer <<" ball error with center label:"<< this->graph->getVLabel(s)<<endl;
		}


		
		//sec. 4.3
		BallTime[pointer]=(double)(endTime - startTime) / CLOCKS_PER_SEC;
		this->BallCount[pointer] = -1;


		/*free the space*/
		for (int i = 0; i < size_ball; i++)
		{
			delete[] M_B[i];
		}
		delete[] M_B;

		for (int i = 0; i < query_size; i++)
		{
			delete[] answer[i];
			delete[] M_P[i];
		}
		delete[] answer;
		delete[] M_P;
		delete[] Matrix_ball;
		Matrix_Ball.clear();
		s_neighbor.clear();
		P_OneIter.clear();
		P_TwoIter.clear();
		P_NeighborLabel.clear();
		P_Replace.clear();
		P_GH.clear();
		GH_Temp.clear();
		//delete bf;
		//delete myBloomFilter;
		delete BallBF;
		delete Ball;

		// Ball_Index.clear();
		// Ball_Index_Reverse.clear();
		// Ball->~DIGRAPH();
		Num_Ball_Computed++;
		pointer++;
		if (each_ball_time > 1)
		{
			//cout << "The " << pointer << " ball is finished in " << each_ball_time << " seconds! The size is: " << size_ball <<  endl;
			cout << "The " << pointer << " ball is finished in " << each_ball_time << " seconds! The size is: " << size_ball << "The GH enumerates "<< count_GH <<" candidates!" << endl;
		}

		endOne = clock();

		Runtime << size_ball << ", " << (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC << ", " << (double)(endFour - startFour) * 1000 / CLOCKS_PER_SEC << ", " << (double)(endThree - startThree) * 1000 / CLOCKS_PER_SEC << ", " << PathTimeTemp << ", " << (double)(endSix - startSix) * 1000 / CLOCKS_PER_SEC<< ", " << TwigTimeTemp << ", " << edge_ball << endl; // total, twoIter, NL, Path, BF, Twig

		// OutFile << pointer << "	"<< size_ball<< "	"<< (double)(endOne - startOne) / CLOCKS_PER_SEC<< "	"<< (double)(endTwo - startTwo) / CLOCKS_PER_SEC << "	" << (double)(endTime1 - startTime1) / CLOCKS_PER_SEC << "	"<< (double)(endThree - startThree) / CLOCKS_PER_SEC << "	"<< (double)(endFour - startFour) / CLOCKS_PER_SEC << "	"<< (double)(endFive - startFive) / CLOCKS_PER_SEC << "	" << (double)(endTime - startTime) / CLOCKS_PER_SEC << endl;
	}

	if (this->query_selected_label_size <= 2)
		cout << "There are little labels to conduct the NeighborLabel Pruning!\n"
			 << endl;

	int ppmatch_num=0;

	// cout path number
	for (int i = 0; i < this->graph_size; i++)
	{
		if (this->result_Path[i] == 1)
			Path_num++;
		if (this->result_Twig[i] == 1)
			Twig_num++;
		if (this->result_BF[i] == 1)
			BF_num++;

		if((this->result_Twig[i]==1)&&(this->result_BF[i]==1))
			ppmatch_num++;
		//sec. 4.3 initial pruning number
		if ((this->result_Twig[i]==1)&&(this->result_BF[i]==1)){
			this->BallCount[i] = 1;
		}
	}





	/*print the results */
	cout << OutFileName << endl;
	cout << "******************************************************************" << endl;
	cout << "\nTotal computed Ball: " << Num_Ball_Computed << endl;
	cout << "The runtime for constructing the traversed balls and Matrices Ps is " << Ball_Time << endl;
	//cout << "Index time " << Index_Time << endl;
	cout << "******************************************************************" << endl;
	if(mode==2){
		cout << exact_solution_num << " exact strong simulation matches:" << StrongSimulation_Time << " s" << endl;
		cout << TwoIter_num << " TwoIter-" << this->portion << " matches:" << Enc_TwoIter_Time << " s" << endl;
		cout << "Decryption Time: " << Decrypt_TwoIter_Time << endl;
		cout << "******************************************************************" << endl;
	}
	if(mode==1){
		cout << exact_solution_num << " exact homo matches:" << StrongSimulation_Time << " s" << endl;
		cout << "******************************************************************" << endl;
	}
	
	cout << Path_num << " Path matches:" << Path_Time << " s" << endl;
	cout << "Decryption Time: " << Decrypt_Path << endl;
	//cout << "The improvement of Path is " << Path_Improved << endl;
	//cout << "******************************************************************" << endl;
	//cout << BF_num << " BF matches:" << BF_Time << " s" << endl;
	//cout << "The improvement of BF is " << BF_Improved << endl;
	cout << "******************************************************************" << endl;
	//Temp_Time
	cout << Twig_num << " Twig matches:" << Twig_Time << " s" << endl;
	cout << "Decryption Time: " << Decrypt_Twig << endl;
	cout << "The improvement of Twig is " << Twig_Improved << endl;
	cout << "******************************************************************" << endl;
	cout << "PPMatch matches:" << ppmatch_num << endl;

	int falseCount10 = 0;
	for (int i = 0; i < pointer; i++)
	{
		if ((result[i] == 1) && (result_BF[i] != 1))
		{
			// cout << "Some Strong Simulations are missing!" << endl;
			falseCount10++;
			cout << "Ball " << i << "is missed!" << endl;
		}
	}
	cout << falseCount10 << " matches are missed by BF!" << endl;


	int falseCountGH_BF = 0, falseCountGH_Twig=0;
	for (int i = 0; i < pointer; i++)
	{
		if ((result[i] == 1) && (result_BF[i] != 1))
		{
			// cout << "Some Strong Simulations are missing!" << endl;
			falseCountGH_BF++;
			cout << "Ball " << i << "is missed by BF!" << endl;
		}
		if ((result[i] == 1) && (result_Twig[i] != 1))
		{
			// cout << "Some Strong Simulations are missing!" << endl;
			falseCountGH_Twig++;
			cout << "Ball " << i << "is missed by Twig!" << endl;
		}
	}
	cout << falseCountGH_BF << " matches are missed by BF!" << endl;
	cout << falseCountGH_Twig << " matches are missed by Twig!" << endl;





	OutFile << "#ALL	LGPM#	LPGM(s) Path#	Path(s)	D(Path)(s) Twig#	Twig(s)	D(Twig)(s) PP#	PPCR	PPCRTree	PPCRTwig Player# AvgBase(ms)	AvgPP(ms) ImpTime AvgTwig(ms) AvgTwigRead(ms)" << endl;
	OutFile << "#The runtime for constructing the traversed balls and Matrices Ps is " << Ball_Time << endl;
	OutFile << Num_Ball_Computed << endl;
	if(mode==2){
		//OutFile << exact_solution_num << " exact strong simulation matches:" << StrongSimulation_Time << " s" << endl;
		OutFile << TwoIter_num << endl;
		OutFile << StrongSimulation_Time << endl;
		//OutFile << "Decryption Time: " << Decrypt_TwoIter_Time << endl;
		
	}
	if(mode==1){
		OutFile << exact_solution_num << endl;
		OutFile << StrongSimulation_Time << endl;
		
	}
	if(mode==0){
		OutFile << "0" << endl;
		OutFile << "0" << endl;
	}
	
	OutFile << Path_num << endl;
	OutFile << Path_Time << endl;
	OutFile << Decrypt_Path << endl;	
	OutFile << Twig_num << endl;
	OutFile << Twig_Time << endl;
	OutFile << Decrypt_Twig << endl;
	OutFile << ppmatch_num << endl;
	OutFile << ((double) (ppmatch_num))/ ((double) (Num_Ball_Computed)) << endl;
	OutFile << BF_num / Num_Ball_Computed << endl;
	OutFile << Twig_num / Num_Ball_Computed << endl;


	// if (exact_solution_num != 0)
	//	cout << "The false positive is " << OneIter_num - exact_solution_num << endl;

	int falseCount1 = 0;
	for (int i = 0; i < pointer; i++)
	{
		if ((result[i] == 1) && (result_EncSSim[i] != 1))
		{
			// cout << "Some Strong Simulations are missing!" << endl;
			falseCount1++;
		}
	}

	int trueCount1 = 0;
	for (int i = 0; i < pointer; i++)
	{
		if (result_EncSSim[i] == 1)
		{
			trueCount1++;
		}
	}

	int falseCount2 = 0;
	for (int i = 0; i < pointer; i++)
	{
		if ((result[i] == 1) && (OneIterNL[i] != 1))
		{
			// cout << "Some Strong Simulations are missing!" << endl;
			falseCount2++;
		}
	}

	int OneCount = 0;
	for (int i = 0; i < pointer; i++)
	{
		if (OneIterNL[i] == 1)
		{
			OneCount++;
		}
	}

	int EncCount = 0;
	for (int i = 0; i < pointer; i++)
	{
		if (EncSpecial[i] == 1)
		{
			EncCount++;
		}
	}

	//OutFile << "There are total " << OneCount << " obtained by One+EncSSim with missing " << falseCount2 << endl;

	/*
	Precision << "# Name, Dataset" << endl;
	Precision << "OneIter\t" << (1 - (OneIter_num - exact_solution_num) / Num_Ball_Computed) << endl;
	Precision << "TwoIter\t" << (1 - (TwoIter_num - exact_solution_num) / Num_Ball_Computed) << endl;
	Precision << "NL\t" << (1 - (NeighborLabel_num - exact_solution_num) / Num_Ball_Computed) << endl;
	Precision << "Path\t" << (1 - (Path_num - exact_solution_num) / Num_Ball_Computed) << endl;
	// Precision << "All,\t" << (1-(trueCount1-exact_solution_num)/Num_Ball_Computed) << endl;
	Precision << "Enc\t" << (1 - (OneCount - exact_solution_num) / Num_Ball_Computed) << endl;
	Precision << "Enc*\t" << (1 - (EncCount - exact_solution_num) / Num_Ball_Computed) << endl;

	Number << "# Name, Dataset" << endl;
	Number << "All\t" << Num_Ball_Computed << endl;
	Number << "Exact\t" << exact_solution_num << endl;
	Number << "OneIter\t" << OneIter_num << endl;
	Number << "TwoIter\t" << TwoIter_num << endl;
	Number << "NL\t" << NeighborLabel_num << endl;
	Number << "Path\t" << Path_num << endl;
	// Number << "All,\t" << trueCount1 << endl;
	Number << "Enc\t" << OneCount << endl;
	Number << "Enc*\t" << EncCount << endl;
	*/



	//sec 4.3
	int times = 10;
	int servers = 4;
	startTime = clock();
	//MultiServers(times, servers, BallTime, TwigTime, TwigBallReadTime, OutFile);
	MultiServers(times, 4, BallTime, TwigTime, TwigBallReadTime, OutFile);
	MultiServers(times, 8, BallTime, TwigTime, TwigBallReadTime, OutFile);
	MultiServers(times, 16, BallTime, TwigTime, TwigBallReadTime, OutFile);
	endTime = clock();
	cout<< "The Dealer need " << (double)(endTime - startTime) / CLOCKS_PER_SEC << " seconds to generate the ball suborders!"<< endl;

	OutFile << MaximumBFNum << endl;


	/*free the space for Queries */
	for (int i = 0; i < hoplength; i++)
	{
		for (int j = 0; j < query_size; j++)
		{
			delete[] Q[i][j];
			delete[] EM_Q_Child[i][j];
			delete[] EM_Q_Parent[i][j];
			delete[] Q_Child[i][j];
			delete[] Q_Parent[i][j];
		}
		delete[] Q[i];
		delete[] EM_Q_Child[i];
		delete[] EM_Q_Parent[i];
		delete[] Q_Child[i];
		delete[] Q_Parent[i];
	}

	for (int i = 0; i < query_size; i++)
	{
		delete[] EM_Q[i];
		delete[] EM_Q_Two[i];
	}

	delete[] Q;
	delete[] EM_Q;
	delete[] EM_Q_Two;
	delete[] EM_Q_Child;
	delete[] EM_Q_Parent;
	delete[] Q_Parent;
	delete[] Q_Child;
	delete[] Column;
	delete cgbe;
	for (int i = 0; i < query_size; i++)
	{
		delete QueryBF[i];
	}
	delete QueryBF;
	delete weight;
	mpz_clear(cipher_one);
	mpz_clear(cipher_zero);
	mpz_clear(cipher_one_Two);
	mpz_clear(cipher_one_Four_Vq);
	Origin_Child.clear();
	Origin_Parent.clear();
	Replace_Child.clear();
	Replace_Parent.clear();
	PathIndex1.clear();
	PathIndex2.clear();
	PathLabel.clear();
	PathNum.clear();
	TwigIndex1.clear();
	TwigIndex2.clear();
	TwigLabel.clear();
	TwigNum.clear();
	OutFile.close();
	Runtime.close();
	BallTime.clear();
	TwigTime.clear();
	TwigBallReadTime.clear();
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Exact(int **M_B, int ***Q, int **answer, int size_ball, int pointer, int center)
{
	bool flag;
	int temp_value;

	/* For test of OneIter in Plaintext */

	// int **answer_fake = new int *[query_size];
	// for (int i = 0; i < query_size; i++){
	// 	answer_fake[i] = new int[size_ball];
	// 	for (int j =0; j< size_ball; j++)
	// 		answer_fake[i][j] = 0;
	// }

	// for (int iteration = 0; iteration < 1; iteration++)
	// {
	// 	flag = true;
	// 	for (int i = 0; i < query_size; i++)
	// 	{
	// 		for (int j = 0; j < size_ball; j++)
	// 		{
	// 			if(answer[i][j]==0)
	// 				continue;
	// 			answer_fake[i][j] = answer[i][j];
	// 			temp_value = answer[i][j];
	// 			int follow = 1, parent = 1, sum1, sum2;
	// 			for (int k = 0; k < query_size; k++)
	// 			{
	// 				sum1 = 0;
	// 				sum2 = 0;
	// 				for (int l = 0; l < size_ball; l++)
	// 				{
	// 					sum1 += (M_B[j][l] * answer[k][l]);
	// 					sum2 += (M_B[l][j] * answer[k][l]);
	// 				}
	// 				sum1 += Q[0][i][k];
	// 				sum2 += Q[0][k][i];
	// 				if ((sum1 != 0) && (sum2 != 0))
	// 					answer_fake[i][j] = answer_fake[i][j] * 1;
	// 				else
	// 					answer_fake[i][j] = 0;
	// 				if (answer[i][j] != temp_value)
	// 					flag = false;
	// 			}
	// 		}
	// 	}
	// 	if (flag)
	// 		break;
	// }

	// for (int i = 0; i < query_size; i++)
	// 	for (int j= 0; j< size_ball; j++)
	// 		answer[i][j] = answer_fake[i][j];

	/* since each element is updated immediately, the pruning is much better */
	for (int iteration = 0; iteration < size_ball; iteration++)
	{
		flag = true;
		for (int i = 0; i < query_size; i++)
		{
			for (int j = 0; j < size_ball; j++)
			{
				if (answer[i][j] == 0)
					continue;
				temp_value = answer[i][j];
				int follow = 1, parent = 1, sum1, sum2;
				for (int k = 0; k < query_size; k++)
				{
					sum1 = 0;
					sum2 = 0;
					for (int l = 0; l < size_ball; l++)
					{
						sum1 += (M_B[j][l] * answer[k][l]);
						sum2 += (M_B[l][j] * answer[k][l]);
					}
					sum1 += Q[0][i][k];
					sum2 += Q[0][k][i];
					if ((sum1 != 0) && (sum2 != 0))
						answer[i][j] = answer[i][j] * 1;
					else
						answer[i][j] = 0;
					if (answer[i][j] != temp_value)
						flag = false;
				}
			}
		}
		if (flag)
			break;
	}

	/*compute the result for this ball*/
	int product = 1, sum = 0;
	int match_num;
	for (int i = 0; i < query_size; i++)
	{
		sum += answer[i][center];
		match_num = 0;
		for (int j = 0; j < size_ball; j++)
		{
			match_num += answer[i][j];
		}
		product *= match_num;
		if (product > 0)
			product = 1;
	}

	if (product > 0 && sum > 0)
	{
		exact_solution_num++;
		this->result[pointer] = 1;
		// cout << "The " << pointer << " ball is a strong simulation match!" << endl;
	}
	else
	{
		if (product < 0)
			cout << "The product in Exact Algorithm is wrong!" << endl;
		this->result[pointer] = 0;
		// cout << "The " << pointer << " ball is not a strong simulation match!" << endl;
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_OneIter(int **M_B, Ciphertext **EM_Q, P_Row &P_OneIter, int size_ball, int pointer, CGBE *cgbe, Ciphertext &cipher_one, Ciphertext &cipher_zero, int center)
{
	//"pointer" here in case of saving the result for Enc_OneIter
	clock_t startTime, endTime;
	P_Row P_temp;

	/*For Decryption with different combined private key */
	// OneIter
	cgbe->setCombinedPrivateKey_OneIter(2 * query_size);

	/*Violation Detector for one time*/
	mpz_t sum1, sum2, follow, parent;
	mpz_init(sum1);
	mpz_init(sum2);
	mpz_init(follow);
	mpz_init(parent);

	startTime = clock();
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_OneIter[i].find(j) == P_OneIter[i].end())
				continue;
			mpz_init(P_temp[i][j]);
			cgbe->setvalue(follow, 1);
			cgbe->setvalue(parent, 1);
			for (int k = 0; k < query_size; k++)
			{
				cgbe->setvalue(sum1, 0);
				cgbe->setvalue(sum2, 0);
				for (int l = 0; l < size_ball; l++)
				{
					if (M_B[j][l] == 1)
					{
						if (P_OneIter[k].find(l) != P_OneIter[k].end())
							cgbe->add(sum1, sum1, P_OneIter[k][l]);
					}
					if (M_B[l][j] == 1)
					{
						if (P_OneIter[k].find(l) != P_OneIter[k].end())
							cgbe->add(sum2, sum2, P_OneIter[k][l]);
					}
				}
				if (!(cgbe->isZero(sum1)))
					cgbe->setvalue(sum1, cipher_one);
				if (!(cgbe->isZero(sum2)))
					cgbe->setvalue(sum2, cipher_one);

				cgbe->add(sum1, sum1, EM_Q[i][k]);
				cgbe->add(sum2, sum2, EM_Q[k][i]);
				cgbe->mul(follow, follow, sum1);
				cgbe->mul(parent, parent, sum2);
			}
			cgbe->mul(P_temp[i][j], P_OneIter[i][j], follow);
			// cgbe->mul(P_temp[i][j], P_OneIter[i][j], parent);   This is wrong!!!! Below is correct!!!! only 2*|V_Q| for combinedPrivateKey
			cgbe->mul(P_temp[i][j], P_temp[i][j], parent);
		}
	}

	// Update the value of matrix P
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_OneIter[i].find(j) == P_OneIter[i].end())
				continue;
			cgbe->setvalue(P_OneIter[i][j], P_temp[i][j]);
		}
	}
	endTime = clock();
	this->Enc_OneIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

	////////////////////////////////////
	///////////    Result    ///////////
	////////////////////////////////////
	bool Result_flag, center_flag = false;
	for (int i = 0; i < query_size; i++)
	{
		startTime = clock();
		Result_flag = true;
		cgbe->setvalue(sum1, 0);
		for (int j = 0; j < size_ball; j++)
		{
			if (P_OneIter[i].find(j) == P_OneIter[i].end())
				continue;
			cgbe->add(sum1, sum1, P_OneIter[i][j]);
		}
		endTime = clock();
		this->Enc_OneIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		// decryption
		///////////////////////////////////////
		/////////////   Client   //////////////
		///////////////////////////////////////
		startTime = clock();
		cgbe->decryption_OneIter(sum1, sum1);
		cgbe->decryption_OneIter(P_OneIter[i][center], P_OneIter[i][center]);
		endTime = clock();
		this->Decrypt_OneIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		if (cgbe->isZero(sum1))
		{
			Result_flag = false;
			break;
		}
		if (!(cgbe->isZero(P_OneIter[i][center])))
		{
			center_flag = true;
		}
	}

	if (Result_flag && center_flag)
	{
		OneIter_num++;
		this->OneIterNL[pointer] = 1;
	}
	else
	{
		this->OneIterNL[pointer] = 0;
	}

	mpz_clear(sum1);
	mpz_clear(sum2);
	mpz_clear(parent);
	mpz_clear(follow);
	P_temp.clear();
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_TwoIter(int **M_B, Ciphertext **EM_Q, Ciphertext **EM_Q_Two, P_Row &P_TwoIter, int size_ball, int pointer, CGBE *cgbe, Ciphertext &cipher_one, Ciphertext &cipher_zero, Replace_Table &Origin_Child, Replace_Table &Replace_Child, Replace_Table &Origin_Parent, Replace_Table &Replace_Parent, int center)
{
	//"pointer" here in case of saving the result for Enc_TwoIter
	clock_t startTime, endTime;
	P_Row P_temp;

	/*For Decryption with different combined private key */
	// TwoIter .... to be determined
	cgbe->setCombinedPrivateKey_TwoIter(4 * query_size);

	/*Violation Detector for the first time*/
	mpz_t sum1, sum2, follow, parent;
	mpz_init(sum1);
	mpz_init(sum2);
	mpz_init(follow);
	mpz_init(parent);

	startTime = clock();
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			mpz_init(P_temp[i][j]);
			cgbe->setvalue(follow, 1);
			cgbe->setvalue(parent, 1);
			for (int k = 0; k < query_size; k++)
			{
				cgbe->setvalue(sum1, 0);
				cgbe->setvalue(sum2, 0);
				for (int l = 0; l < size_ball; l++)
				{
					if (M_B[j][l] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][l]);
					}
					if (M_B[l][j] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][l]);
					}
				}
				if (!(cgbe->isZero(sum1)))
					cgbe->setvalue(sum1, cipher_one);
				else
					cgbe->setvalue(sum1, cipher_zero);

				if (!(cgbe->isZero(sum2)))
					cgbe->setvalue(sum2, cipher_one);
				else
					cgbe->setvalue(sum2, cipher_zero);

				cgbe->add(sum1, sum1, EM_Q[i][k]);
				cgbe->add(sum2, sum2, EM_Q[k][i]);
				cgbe->mul(follow, follow, sum1);
				cgbe->mul(parent, parent, sum2);
			}

			// Replacement of Ciphertext
			for (auto it = Origin_Child[i].begin(); it != Origin_Child[i].end(); it++)
			{
				if (cgbe->isEqual(it->second, follow))
				{
					cgbe->setvalue(follow, Replace_Child[i][it->first]);
					// cout << "Child Replacement Suceed! " <<endl;
					break;
				}
			}

			for (auto it = Origin_Parent[i].begin(); it != Origin_Parent[i].end(); it++)
			{
				if (cgbe->isEqual(it->second, parent))
				{
					cgbe->setvalue(parent, Replace_Parent[i][it->first]);
					// cout << "Parent Replacement Suceed! " <<endl;
					break;
				}
			}

			cgbe->mul(P_temp[i][j], P_TwoIter[i][j], follow);
			cgbe->mul(P_temp[i][j], P_temp[i][j], parent);
		}
	}

	// Update the value of matrix P
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->setvalue(P_TwoIter[i][j], P_temp[i][j]);
		}
	}

	// mpz_t temp;
	// mpz_init(temp);
	/*Violation Detector for the second time*/
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->setvalue(follow, 1);
			cgbe->setvalue(parent, 1);
			for (int k = 0; k < query_size; k++)
			{
				cgbe->setvalue(sum1, 0);
				cgbe->setvalue(sum2, 0);
				for (int l = 0; l < size_ball; l++)
				{
					if (M_B[j][l] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][l]);
					}
					if (M_B[l][j] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][l]);
					}
				}
				// cgbe->setvalue(temp, EM_Q[i][k]);
				// cgbe->mul(temp, temp, temp);
				// cgbe->add(sum1, sum1, temp);
				cgbe->add(sum1, sum1, EM_Q_Two[i][k]);
				// cgbe->setvalue(temp, EM_Q[k][i]);
				// cgbe->mul(temp, temp, temp);
				// cgbe->add(sum2, sum2, temp);
				cgbe->add(sum2, sum2, EM_Q_Two[k][i]);
				cgbe->mul(follow, follow, sum1);
				cgbe->mul(parent, parent, sum2);
			}

			cgbe->setvalue(P_temp[i][j], follow);
			cgbe->mul(P_temp[i][j], P_temp[i][j], parent);
		}
	}

	// Update the value of matrix P
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->setvalue(P_TwoIter[i][j], P_temp[i][j]);
		}
	}

	endTime = clock();
	this->Enc_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

	////////////////////////////////////
	///////////    Result    ///////////
	////////////////////////////////////
	bool Result_flag, center_flag = false;
	for (int i = 0; i < query_size; i++)
	{
		startTime = clock();
		Result_flag = true;
		cgbe->setvalue(sum1, 0);
		for (int j = 0; j < size_ball; j++)
		{
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->add(sum1, sum1, P_TwoIter[i][j]);
		}
		endTime = clock();
		this->Enc_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		// decryption
		///////////////////////////////////////
		/////////////   Client   //////////////
		///////////////////////////////////////
		startTime = clock();
		cgbe->decryption_TwoIter(sum1, sum1);
		endTime = clock();
		this->Decrypt_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		if (cgbe->isZero(sum1))
		{
			Result_flag = false;
			break;
		}

		if (P_TwoIter[i].find(center) != P_TwoIter[i].end())
		{
			startTime = clock();
			cgbe->decryption_TwoIter(P_TwoIter[i][center], P_TwoIter[i][center]);
			endTime = clock();
			this->Decrypt_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			if (!(cgbe->isZero(P_TwoIter[i][center])))
			{
				center_flag = true;
			}
		}
	}

	if (Result_flag && center_flag)
	{
		TwoIter_num++;
		this->result_EncSSim[pointer] = 1;
		//	this->result_EncSSim[pointer] = 1;
	}
	else
	{
		this->result_EncSSim[pointer] = 0;
		this->OneIterNL[pointer] = 0;
		this->EncSpecial[pointer] = 0;
	}

	mpz_clear(sum1);
	mpz_clear(sum2);
	mpz_clear(parent);
	mpz_clear(follow);
	// mpz_clear(temp);
	P_temp.clear();
}


template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_GH(int **M_B, Ciphertext **EM_Q, int size_query, GH &P_GH, Vertex_Map &Temp_map, int size_ball, int pointer, CGBE *cgbe, Ciphertext &cipher_one, Ciphertext &cipher_zero, int center, int iter, bool &flag, int center_match, int decryptionNum, int &count, clock_t &Start)
{
	if(iter==size_query){
		count++;
		clock_t startTime, endTime;
		Ciphertext GH_Result;	
		bool temp_flag = true;
		mpz_init(GH_Result);
		cgbe->setvalue(GH_Result, 1);
		int decodeNum = 0;
		for(int i=0; i<size_query;i++){
			for(int j=0; j<size_query;j++){
				if(i==j)
					continue;
				if(M_B[Temp_map[i]][Temp_map[j]]==0){
					decodeNum++;
					cgbe->mul(GH_Result, GH_Result, EM_Q[i][j]);
					//decryption
					if (decodeNum == decryptionNum)
					{
						startTime = clock();
						cgbe->decryption_GH(GH_Result, GH_Result);
						if (cgbe->isZero(GH_Result)){
							temp_flag = false;
							cgbe->setvalue(GH_Result, 1);							
						}						
						decodeNum = 0;
						endTime = clock();
						this->Decrypt_GH += (double)(endTime - startTime) / CLOCKS_PER_SEC;
					}
				}
			}
		}

		//decryption
		while (decodeNum != decryptionNum)
		{		
			decodeNum++;
			cgbe->mul(GH_Result, GH_Result, cipher_one);			
		}

		startTime = clock();
		cgbe->decryption_GH(GH_Result, GH_Result);
		if (cgbe->isZero(GH_Result))
			temp_flag = false;			
		endTime = clock();
		this->Decrypt_GH += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		
		if(temp_flag){
			flag = true;
		}

		mpz_clear(GH_Result);
		return;
	}


	//Initialization
	if(iter == -1){	
	// first determine the ball center mapping
		clock_t startTime, currentTime;
		startTime=clock();
		cgbe->setCombinedPrivateKey_GH(decryptionNum);
		flag = false;
		for(int i=0; i<size_query; i++){
			if(P_GH[i].find(center)==P_GH[i].end())
				continue;
			Temp_map.clear();
			Temp_map[i]= center;
			Enc_GH(M_B, EM_Q, size_query, P_GH, Temp_map, size_ball, pointer, cgbe, cipher_one, cipher_zero, center, 0, flag, i, decryptionNum, count, startTime);
			currentTime = clock();
			if(((double)(currentTime - startTime) / CLOCKS_PER_SEC)>10)
				break;
		}

		//results
		if(flag){
			exact_solution_num++;
		}else{
			this->result_GH[pointer] = 0;
			this->result[pointer] = 0;
		}
	}else{
		clock_t currentTime;
		int position = iter;
		if(position == center_match){
			Enc_GH(M_B, EM_Q, size_query, P_GH, Temp_map, size_ball, pointer, cgbe, cipher_one, cipher_zero, center, iter+1, flag, center_match, decryptionNum, count, Start);
		}else{
			for(auto it = P_GH[position].begin(); it !=P_GH[position].end(); it++){
				Temp_map[position]=it->first;
				Enc_GH(M_B, EM_Q, size_query, P_GH, Temp_map, size_ball, pointer, cgbe, cipher_one, cipher_zero, center, iter+1, flag, center_match, decryptionNum, count, Start);
				currentTime = clock();
				if(((double)(currentTime - Start) / CLOCKS_PER_SEC)>10)
					break;
			}
		}
	}
}



template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_TwoIter_Random(int **M_B, Ciphertext **EM_Q, Ciphertext **EM_Q_Two, P_Row &P_TwoIter, int size_ball, int pointer, CGBE *cgbe, Ciphertext &cipher_one, Ciphertext &cipher_zero, Replace_Table &Origin_Child, Replace_Table &Replace_Child, Replace_Table &Origin_Parent, Replace_Table &Replace_Parent, int center, Ciphertext &cipher_one_Two, Ciphertext &cipher_one_Four_Vq, double portion)
{
	//"pointer" here in case of saving the result for Enc_TwoIter
	clock_t startTime, endTime;
	P_Row P_temp;

	// int* randomarray = new int[size_ball];
	// for(int i = 0; i < size_ball; i++)
	//	randomarray[i] = 0;
	// int count = (portion - 1)*size_ball/portion;
	int count;
	count = size_ball - (size_ball * portion);
	count--;
	if (count < 1)
	{
		this->result_EncSSim[pointer] = 1;
		P_temp.clear();
		return;
	}
	////////////////////another improvement
	/*
	int place;
	unordered_set<int> rarray;
		for(int i =0; i < count; i++){
		place = rand()%size_ball;
		while(rarray.find(place)!=rarray.end()){
			place = (place + 1)%size_ball;
		}
		rarray.insert(place);
	}*/

	/*
	int place;
	for(int i =0; i < count; i++){
		place = rand()%size_ball;
		while(randomarray[place]!=0){
			place = (place + 1)%size_ball;
		}
		randomarray[place] = 1;
	}*/

	int *randomarray = new int[count];
	int place;
	bool flag;
	for (int i = 0; i < count; i++)
	{
		place = rand() % size_ball;
		flag = true;
		while (flag)
		{
			flag = false;
			for (int j = 0; j <= i; j++)
			{
				if (randomarray[j] == place)
					flag = true;
			}
			if (place == center) // center cannot be always equal to 1
				flag = true;
			if (flag)
				place = (place + 1) % size_ball;
		}
		randomarray[i] = place;
	}

	int *indicate = new int[size_ball];
	for (int i = 0; i < size_ball; i++)
	{
		indicate[i] = 0;
	}
	for (int i = 0; i < count; i++)
	{
		indicate[randomarray[i]] = 1;
	}

	/*For Decryption with different combined private key */
	// TwoIter .... to be determined
	cgbe->setCombinedPrivateKey_TwoIter(4 * query_size);

	/*Violation Detector for the first time*/

	mpz_t sum1, sum2, follow, parent;
	mpz_init(sum1);
	mpz_init(sum2);
	mpz_init(follow);
	mpz_init(parent);

	startTime = clock();
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;

			mpz_init(P_temp[i][j]);

			/*
			if (randomarray[j] == 1){
				cgbe->setvalue(P_temp[i][j], cipher_one_Two);
				continue;
			}*/

			/*if (rarray.find(j)!=rarray.end()){
				cgbe->setvalue(P_temp[i][j], cipher_one_Two);
				continue;
			}*/

			if (indicate[j] == 1)
			{
				cgbe->setvalue(P_temp[i][j], cipher_one_Two);
				continue;
			}

			cgbe->setvalue(follow, 1);
			cgbe->setvalue(parent, 1);
			for (int k = 0; k < query_size; k++)
			{
				cgbe->setvalue(sum1, 0);
				cgbe->setvalue(sum2, 0);

				/*
				for (auto it = rarray.begin(); it != rarray.end(); it++)
				{
					if (M_B[j][*it] == 1)
					{
						if (P_TwoIter[k].find(*it) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][*it]);
					}
					if (M_B[*it][j] == 1)
					{
						if (P_TwoIter[k].find(*it) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][*it]);
					}
				}*/

				for (int l = 0; l < size_ball; l++)
				{
					if (M_B[j][l] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][l]);
					}
					if (M_B[l][j] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][l]);
					}
				}

				/*
				for (int l = 0; l < count; l++)
				{
					if (M_B[j][randomarray[l]] == 1)
					{
						if (P_TwoIter[k].find(randomarray[l]) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][randomarray[l]]);
					}
					if (M_B[randomarray[l]][j] == 1)
					{
						if (P_TwoIter[k].find(randomarray[l]) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][randomarray[l]]);
					}
				}*/

				if (!(cgbe->isZero(sum1)))
					cgbe->setvalue(sum1, cipher_one);
				else
					cgbe->setvalue(sum1, cipher_zero);

				if (!(cgbe->isZero(sum2)))
					cgbe->setvalue(sum2, cipher_one);
				else
					cgbe->setvalue(sum2, cipher_zero);

				cgbe->add(sum1, sum1, EM_Q[i][k]);
				cgbe->add(sum2, sum2, EM_Q[k][i]);
				cgbe->mul(follow, follow, sum1);
				cgbe->mul(parent, parent, sum2);
			}

			// Replacement of Ciphertext
			for (auto it = Origin_Child[i].begin(); it != Origin_Child[i].end(); it++)
			{
				if (cgbe->isEqual(it->second, follow))
				{
					cgbe->setvalue(follow, Replace_Child[i][it->first]);
					// cout << "Child Replacement Suceed! " <<endl;
					break;
				}
			}

			for (auto it = Origin_Parent[i].begin(); it != Origin_Parent[i].end(); it++)
			{
				if (cgbe->isEqual(it->second, parent))
				{
					cgbe->setvalue(parent, Replace_Parent[i][it->first]);
					// cout << "Parent Replacement Suceed! " <<endl;
					break;
				}
			}

			cgbe->setvalue(P_temp[i][j], follow);
			cgbe->mul(P_temp[i][j], P_temp[i][j], parent);
		}
	}

	// Update the value of matrix P
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->setvalue(P_TwoIter[i][j], P_temp[i][j]);
		}
	}

	// mpz_t temp;
	// mpz_init(temp);
	/*Violation Detector for the second time*/
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;

			/*
			if (randomarray[j] == 1){
				cgbe->setvalue(P_temp[i][j], cipher_one_Four_Vq);
				continue;
			}*/

			/*if (rarray.find(j)!=rarray.end()){
				cgbe->setvalue(P_temp[i][j], cipher_one_Four_Vq);
				continue;
			}*/

			if (indicate[j] == 1)
			{
				cgbe->setvalue(P_temp[i][j], cipher_one_Four_Vq);
				continue;
			}

			cgbe->setvalue(follow, 1);
			cgbe->setvalue(parent, 1);
			for (int k = 0; k < query_size; k++)
			{
				cgbe->setvalue(sum1, 0);
				cgbe->setvalue(sum2, 0);

				/*for (auto it = rarray.begin(); it != rarray.end(); it++)
				{
					if (M_B[j][*it] == 1)
					{
						if (P_TwoIter[k].find(*it) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][*it]);
					}
					if (M_B[*it][j] == 1)
					{
						if (P_TwoIter[k].find(*it) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][*it]);
					}
				}*/

				for (int l = 0; l < size_ball; l++)
				{
					if (M_B[j][l] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][l]);
					}
					if (M_B[l][j] == 1)
					{
						if (P_TwoIter[k].find(l) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][l]);
					}
				}
				/*
				for (int l = 0; l < count; l++)
				{
					if (M_B[j][randomarray[l]] == 1)
					{
						if (P_TwoIter[k].find(randomarray[l]) != P_TwoIter[k].end())
							cgbe->add(sum1, sum1, P_TwoIter[k][randomarray[l]]);
					}
					if (M_B[randomarray[l]][j] == 1)
					{
						if (P_TwoIter[k].find(randomarray[l]) != P_TwoIter[k].end())
							cgbe->add(sum2, sum2, P_TwoIter[k][randomarray[l]]);
					}
				}*/

				// cgbe->setvalue(temp, EM_Q[i][k]);
				// cgbe->mul(temp, temp, temp);
				// cgbe->add(sum1, sum1, temp);
				cgbe->add(sum1, sum1, EM_Q_Two[i][k]);
				// cgbe->setvalue(temp, EM_Q[k][i]);
				// cgbe->mul(temp, temp, temp);
				// cgbe->add(sum2, sum2, temp);
				cgbe->add(sum2, sum2, EM_Q_Two[k][i]);
				cgbe->mul(follow, follow, sum1);
				cgbe->mul(parent, parent, sum2);
			}

			cgbe->setvalue(P_temp[i][j], follow);
			cgbe->mul(P_temp[i][j], P_temp[i][j], parent);
		}
	}

	// Update the value of matrix P
	for (int i = 0; i < query_size; i++)
	{
		for (int j = 0; j < size_ball; j++)
		{
			// No need to compute 0 element
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->setvalue(P_TwoIter[i][j], P_temp[i][j]);
		}
	}

	endTime = clock();
	this->Enc_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

	////////////////////////////////////
	///////////    Result    ///////////
	////////////////////////////////////
	bool Result_flag, center_flag = false;
	for (int i = 0; i < query_size; i++)
	{
		startTime = clock();
		Result_flag = true;
		cgbe->setvalue(sum1, 0);
		for (int j = 0; j < size_ball; j++)
		{
			if (P_TwoIter[i].find(j) == P_TwoIter[i].end())
				continue;
			cgbe->add(sum1, sum1, P_TwoIter[i][j]);
		}
		endTime = clock();
		this->Enc_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		// decryption
		///////////////////////////////////////
		/////////////   Client   //////////////
		///////////////////////////////////////
		startTime = clock();
		cgbe->decryption_TwoIter(sum1, sum1);
		endTime = clock();
		this->Decrypt_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
		if (cgbe->isZero(sum1))
		{
			Result_flag = false;
			break;
		}

		if (P_TwoIter[i].find(center) != P_TwoIter[i].end())
		{
			startTime = clock();
			cgbe->decryption_TwoIter(P_TwoIter[i][center], P_TwoIter[i][center]);
			endTime = clock();
			this->Decrypt_TwoIter_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			if (!(cgbe->isZero(P_TwoIter[i][center])))
			{
				center_flag = true;
			}
		}
	}

	if (Result_flag && center_flag)
	{
		TwoIter_num++;
		this->result_EncSSim[pointer] = 1;
		//	this->result_EncSSim[pointer] = 1;
	}
	else
	{
		this->result_EncSSim[pointer] = 0;
		this->OneIterNL[pointer] = 0;
		this->EncSpecial[pointer] = 0;
	}

	mpz_clear(sum1);
	mpz_clear(sum2);
	mpz_clear(parent);
	mpz_clear(follow);
	// mpz_clear(temp);
	P_temp.clear();
	delete[] randomarray;
	delete[] indicate;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_NeighborLabel(int **M_B, Ciphertext ***EM_Q_Child, Ciphertext ***EM_Q_Parent, P_Row &P_NeighborLabel, P_Row &P_Replace, int size_ball, int pointer, CGBE *cgbe, B_Index &Graph_Index, B_Index &Graph_Index_Reverse, VLabelType *Column, Ciphertext &cipher_one, Ciphertext &cipher_zero, int center, VertexID centerNo, VertexID *Matrix_ball)
{
	if (hoplength <= 2)
		return;

	/*For Decryption with different combined private key */
	// NeighborLabel
	// cgbe->setCombinedPrivateKey_NeighborLabel(2 * (K_HOP - 1) * this->query_selected_label_size); //from 2 hop

	// cgbe->setCombinedPrivateKey_NeighborLabel(2 * (hoplength) * this->query_selected_label_size); //from 1 hop
	cgbe->setCombinedPrivateKey_NeighborLabel(2 * this->query_selected_label_size);

	clock_t startTime, endTime;

	mpz_t sum;
	mpz_init(sum);

	////////////////////////////////////
	/////////    Index Tech    /////////
	////////////////////////////////////
	int column_size = this->query_selected_label_size;

	bool Result_flag = false, result_temp;

	for (int i = 0; i < query_size; i++)
	{
		/////////////
		// for (int j = 0; j < size_ball; j++){

		// if (P_NeighborLabel[i].find(j) == P_NeighborLabel[i].end())
		// continue;
		// No need to compute 0 element

		/*
			for (int l = 0; l < K_HOP; l++)
			{
				for (int m = 0; m < column_size; m++)
				{
					if (Graph_Index[l][Matrix_ball[j]].find(Column[m]) == Graph_Index[l][Matrix_ball[j]].end())
					{
						cgbe->mul(P_NeighborLabel[i][j], P_NeighborLabel[i][j], EM_Q_Child[l][i][m]);
					}
					else
					{
						cgbe->mul(P_NeighborLabel[i][j], P_NeighborLabel[i][j], cipher_one);
					}

					if (Graph_Index_Reverse[l][Matrix_ball[j]].find(Column[m]) == Graph_Index_Reverse[l][Matrix_ball[j]].end())
					{
						cgbe->mul(P_NeighborLabel[i][j], P_NeighborLabel[i][j], EM_Q_Parent[l][i][m]);
					}
					else
					{
						cgbe->mul(P_NeighborLabel[i][j], P_NeighborLabel[i][j], cipher_one);
					}
				}
			}*/

		/*****************************/
		if (P_NeighborLabel[i].find(center) == P_NeighborLabel[i].end())
			continue;
		// l = 1 for two-hop neighbor, l = 2 for three-hop neighbor

		result_temp = true;
		for (int l = 0; l < hoplength; l++) // each hop
		{
			cgbe->setvalue(P_Replace[i][center], P_NeighborLabel[i][center]);

			startTime = clock();
			for (int m = 0; m < column_size; m++)
			{

				if (Graph_Index[l][centerNo].find(Column[m]) == Graph_Index[l][centerNo].end())
				{
					cgbe->mul(P_Replace[i][center], P_Replace[i][center], EM_Q_Child[l][i][m]);
				}
				else
				{
					cgbe->mul(P_Replace[i][center], P_Replace[i][center], cipher_one);
				}

				if (Graph_Index_Reverse[l][centerNo].find(Column[m]) == Graph_Index_Reverse[l][centerNo].end())
				{
					cgbe->mul(P_Replace[i][center], P_Replace[i][center], EM_Q_Parent[l][i][m]);
				}
				else
				{
					cgbe->mul(P_Replace[i][center], P_Replace[i][center], cipher_one);
				}
			}
			endTime = clock();
			this->NeighborLabel_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

			////////////////////////////////////
			///////////    Result    ///////////
			////////////////////////////////////

			// theorectically right without overflow, but here should combine every hop but not all the hops which may lead to overflow. Revise when have time.

			startTime = clock();
			cgbe->decryption_NeighborLabel(P_Replace[i][center], P_Replace[i][center]);
			endTime = clock();
			this->Decrypt_NeighborLabel += (double)(endTime - startTime) / CLOCKS_PER_SEC;

			if ((cgbe->isZero(P_Replace[i][center])))
			{
				result_temp = false;
			}
		}
		if (result_temp)
			Result_flag = true;
	}

	if (Result_flag)
	{
		NeighborLabel_num++;
		// this->result_EncSSim[pointer] = 1;
	}
	else
	{
		if (this->result_EncSSim[pointer] == 1)
			NL_Improved++;
		this->result_EncSSim[pointer] = 0;
		this->OneIterNL[pointer] = 0;
		this->EncSpecial[pointer] = 0;
	}
	mpz_clear(sum);
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::NLTest(int **M_B, int ***Q_Child, int ***Q_Parent, int **M_P, int size_ball, int pointer, CGBE *cgbe, B_Index &Graph_Index, B_Index &Graph_Index_Reverse, VLabelType *Column, int center, VertexID centerNo, VertexID *Matrix_ball)
{
	if (hoplength <= 2)
		return;

	/*For Decryption with different combined private key */
	// NeighborLabel
	// cgbe->setCombinedPrivateKey_NeighborLabel(2 * (K_HOP - 1) * this->query_selected_label_size); //from 2 hop
	// cgbe->setCombinedPrivateKey_NeighborLabel(2 * (hoplength -2) * this->query_selected_label_size); //from 1 hop

	clock_t startTime, endTime;

	startTime = clock();
	mpz_t sum;
	mpz_init(sum);

	////////////////////////////////////
	/////////    Index Tech    /////////
	////////////////////////////////////
	int column_size = this->query_selected_label_size;

	for (int i = 0; i < query_size; i++)
	{
		/*****************************/
		if (M_P[i][center] == 0)
			continue;
		// l = 1 for two-hop neighbor, l = 2 for three-hop neighbor

		for (int l = 0; l < hoplength; l++)
		{
			for (int m = 0; m < column_size; m++)
			{
				if (Graph_Index[l][centerNo].find(Column[m]) == Graph_Index[l][centerNo].end())
				{
					if (Q_Child[l][i][m] == 0)
						M_P[i][center] = 0;
				}

				if (Graph_Index_Reverse[l][centerNo].find(Column[m]) == Graph_Index_Reverse[l][centerNo].end())
				{
					if (Q_Parent[l][i][m] == 0)
						M_P[i][center] = 0;
				}
			}
		}

		/////////////
		//}
	}

	endTime = clock();
	this->NeighborLabel_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	////////////////////////////////////
	///////////    Result    ///////////
	////////////////////////////////////

	// bool Result_flag = true;
	bool Result_flag = false;

	for (int i = 0; i < query_size; i++)
	{
		/*
		startTime = clock();
		cgbe->setvalue(sum, 0);
		for (int j = 0; j < size_ball; j++)
		{
			if (P_NeighborLabel[i].find(j) == P_NeighborLabel[i].end())
				continue;
			cgbe->add(sum, sum, P_NeighborLabel[i][j]);
		}
		endTime = clock();
		this->NeighborLabel_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;



		startTime = clock();
		cgbe->decryption_NeighborLabel(sum, sum);
		endTime = clock();
		this->Decrypt_NeighborLabel += (double)(endTime - startTime) / CLOCKS_PER_SEC;


		if (cgbe->isZero(sum)){
			Result_flag = false;
			break;
		}*/

		// decryption
		///////////////////////////////////////
		/////////////   Client   //////////////
		///////////////////////////////////////

		if (M_P[i][center] != 0)
		{
			Result_flag = true;
			break;
		}
	}

	if (Result_flag)
	{
		NeighborLabel_num++;
		// this->result_EncSSim[pointer] = 1;
	}
	else
	{
		if (this->OneIterNL[pointer] == 1)
			NL_Improved++;
		this->result_EncSSim[pointer] = 0;
		this->OneIterNL[pointer] = 0;
	}
	mpz_clear(sum);
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_Path(Path_Label &PathLabel, Path_Num &PathNum, Path_Index &PathIndex1, Path_Index &PathIndex2, int pointer, int **M_B, VertexID *Matrix_ball, int size_ball, CGBE *cgbe, double &Time, int center, VertexID centerNo, int pathlength, int decryptionNum, Ciphertext &cipher_one)
{
	int RevisedLengthOfM_B, TotalCantorNum;
	int decodeNum = 0;
	// VLabelType

	clock_t startTime, endTime;
	// startTime = clock();

	Ciphertext PathPruning;
	Ciphertext Result;
	Ciphertext FinalResult;
	mpz_init(PathPruning);
	mpz_init(Result);
	mpz_init(FinalResult);
	cgbe->setvalue(PathPruning, 1);
	cgbe->setvalue(Result, 1);
	cgbe->setvalue(FinalResult, 0);

	DIGRAPH<VLabelType, ELabelType> *RBall = new DIGRAPH<VLabelType, ELabelType>;
	for (int i = 0; i < size_ball; i++)
	{
		RBall->insertVertex(Matrix_ball[i], 0);
	}

	for (auto it3 = RBall->getVLabel().begin(); it3 != RBall->getVLabel().end(); it3++)
	{
		for (auto it4 = RBall->getVLabel().begin(); it4 != RBall->getVLabel().end(); it4++)
		{
			if (it3->first == it4->first)
				continue;
			if (this->graph->isEdge(it3->first, it4->first))
				RBall->insertEdge(it3->first, it4->first, 0);
			if (this->graph->isEdge(it4->first, it3->first))
				RBall->insertEdge(it4->first, it3->first, 0);
		}
	}

	startTime = clock();
	// Save the paths of center
	unordered_map<VertexID, unordered_map<int, int>> Path_Center_Child, Path_Center_Parent;

	// Find the paths of center
	PathMatch(RBall, PathNum, Path_Center_Child, Path_Center_Parent, nullptr, RevisedLengthOfM_B, pathlength, 0, centerNo, centerNo, -1);
	endTime = clock();
	this->Path_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;

	/*
	cout << "#######################" << endl;
	cout << "#######################" << endl;

	for(auto it = Path_Center_Child.begin();it != Path_Center_Child.end(); it++){
		cout << "Vertex " << it->first <<":" << endl;
		for(auto it2 = Path_Center_Child[it->first].begin(); it2 != Path_Center_Child[it->first].end(); it2++){
			cout <<"Sum:" << it2->first<<endl;
		}

	}*/

	// Match
	cgbe->setCombinedPrivateKey_Path(decryptionNum);
	for (auto i = this->query->getVLabel().begin(); i != this->query->getVLabel().end(); i++)
	{
		if (this->query->getVLabel(i->first) != this->graph->getVLabel(centerNo))
			continue;

		for (auto it = PathIndex1[pathlength][i->first].begin(); it != PathIndex1[pathlength][i->first].end(); it++)
		{
			startTime = clock();
			if (Path_Center_Child[centerNo].find(it->first) == Path_Center_Child[centerNo].end())
			{
				decodeNum++;			
				cgbe->mul(PathPruning, PathPruning, PathIndex1[pathlength][i->first][it->first]);
			}
			endTime = clock();
			this->Path_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;

			if (decodeNum == decryptionNum)
			{
				startTime = clock();
				cgbe->decryption_Path(PathPruning, PathPruning);
				if (cgbe->isZero(PathPruning))
					cgbe->mul(Result, Result, cgbe->zero);
				cgbe->setvalue(PathPruning, 1);
				decodeNum = 0;
				endTime = clock();
				this->Decrypt_Path += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			}
		}

		for (auto it2 = PathIndex2[pathlength][i->first].begin(); it2 != PathIndex2[pathlength][i->first].end(); it2++)
		{

			startTime = clock();
			if (Path_Center_Parent[centerNo].find(it2->first) == Path_Center_Parent[centerNo].end())
			{
				decodeNum++;
				cgbe->mul(PathPruning, PathPruning, PathIndex2[pathlength][i->first][it2->first]);				
			}
			endTime = clock();
			this->Path_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;

			if (decodeNum == decryptionNum)
			{
				startTime = clock();
				cgbe->decryption_Path(PathPruning, PathPruning);
				if (cgbe->isZero(PathPruning))
					cgbe->mul(Result, Result, cgbe->zero);
				cgbe->setvalue(PathPruning, 1);
				decodeNum = 0;
				endTime = clock();
				this->Decrypt_Path += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			}
		}

		while (decodeNum != decryptionNum)
		{
			startTime = clock();
			decodeNum++;
			cgbe->mul(PathPruning, PathPruning, cipher_one);
			endTime = clock();
			this->Path_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;
		}

		startTime = clock();
		cgbe->decryption_Path(PathPruning, PathPruning);
		if (cgbe->isZero(PathPruning))
			cgbe->mul(Result, Result, cgbe->zero);
		endTime = clock();
		this->Decrypt_Path += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		cgbe->add(FinalResult, FinalResult, Result);
		cgbe->setvalue(PathPruning, 1);
		cgbe->setvalue(Result, 1);
		decodeNum = 0;
	}

	// endTime = clock();
	// this->Path_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

	/*		Client		*/

	if ((cgbe->isZero(FinalResult)))
	{
		//if (this->result_EncSSim[pointer] == 1)
		//	Path_Improved++;
		//this->result_EncSSim[pointer] = 0;
		//this->OneIterNL[pointer] = 0;
		this->result_Path[pointer] = 0;
		//this->EncSpecial[pointer] = 0;
	}
	mpz_clear(PathPruning);
	mpz_clear(Result);
	mpz_clear(FinalResult);
	Path_Center_Child.clear();
	Path_Center_Parent.clear();
	delete RBall;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::PathMatch(DIGRAPH<VLabelType, ELabelType> *RBall, Path_Num &PathNum, unordered_map<VertexID, unordered_map<int, int>> &Path_Center_Child, unordered_map<VertexID, unordered_map<int, int>> &Path_Center_Parent, VLabelType *Path, int RevisedBallSize, int PathLength, int iter, VertexID center, VertexID pointer, int point)
{
	if (iter == PathLength)
	{
		int sum = 0;
		for (int i = 1; i < PathLength; i++)
		{
			int temp = 1;
			for (int j = 1; j < i; j++)
			{
				temp = temp * PathNum[PathLength][-1];
			}
			sum += temp * PathNum[PathLength][Path[i]];
		}
		if (point == 0)
		{
			Path_Center_Child[center][sum] = 1;
		}
		if (point == 1)
		{
			Path_Center_Parent[center][sum] = 1;
		}
		return;
	}

	if (iter == 0)
	{
		// for center only 2020.8.28
		// if (this->graph->getVLabel(pointer) != Path[iter])
		//	return false;
		// if (PathMatch(RBall, Path, RevisedBallSize, PathLength, iter + 1, pointer))
		//	return true;

		// for center only 2020.8.28

		VLabelType *temppath1 = new VLabelType[PathLength];
		VLabelType *temppath2 = new VLabelType[PathLength];
		// for child
		temppath1[0] = this->graph->getVLabel(center);
		PathMatch(RBall, PathNum, Path_Center_Child, Path_Center_Parent, temppath1, RevisedBallSize, PathLength, iter + 1, center, center, 0);

		// for parent
		temppath2[0] = this->graph->getVLabel(center);
		PathMatch(RBall, PathNum, Path_Center_Child, Path_Center_Parent, temppath2, RevisedBallSize, PathLength, iter + 1, center, center, 1);

		delete[] temppath1;
		delete[] temppath2;
	}
	else if (iter < PathLength)
	{

		if (point == 0)
		{
			for (auto it3 = RBall->getOutEdge()[pointer].begin(); it3 != RBall->getOutEdge()[pointer].end(); it3++)
			{
				bool flag = true;
				for (int i = 0; i < iter; i++)
				{
					if ((this->graph->getVLabel(it3->first) == Path[i]))
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					Path[iter] = this->graph->getVLabel(it3->first);
					PathMatch(RBall, PathNum, Path_Center_Child, Path_Center_Parent, Path, RevisedBallSize, PathLength, iter + 1, center, it3->first, 0);
				}
			}
		}

		if (point == 1)
		{
			for (auto it4 = RBall->getInVertex()[pointer].begin(); it4 != RBall->getInVertex()[pointer].end(); it4++)
			{
				bool flag = true;
				for (int i = 0; i < iter; i++)
				{
					if ((this->graph->getVLabel(it4->first) == Path[i]))
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					Path[iter] = this->graph->getVLabel(it4->first);
					PathMatch(RBall, PathNum, Path_Center_Child, Path_Center_Parent, Path, RevisedBallSize, PathLength, iter + 1, center, it4->first, 1);
				}
			}
		}
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Q_K_HOP(int ***Q, Ciphertext ***EM_Q_Child, Ciphertext ***EM_Q_Parent, int ***Q_Child, int ***Q_Parent, VLabelType *Column, CGBE *cgbe)
{

	for (int i = 1; i < hoplength; i++)
	{
		Q[i] = new int *[query_size];
		for (int j = 0; j < query_size; j++)
		{
			Q[i][j] = new int[query_size];
			for (int k = 0; k < query_size; k++)
			{
				Q[i][j][k] = 1;
			}
		}
	}

	for (int i = 1; i < hoplength; i++)
	{
		for (int j = 0; j < query_size; j++)
		{
			for (int k = 0; k < query_size; k++)
			{
				// no cycles
				// if(j==k){
				//	Q[i][j][k] = 1;
				//	continue;
				// }

				int sum_value_temp = 1;
				/*compute the neighbor hop information*/
				for (int l = 0; l < query_size; l++)
				{
					sum_value_temp *= (Q[i - 1][j][l] + Q[0][l][k]);
					if (sum_value_temp == 0)
					{
						Q[i][j][k] = 0;
					}
				}

				// cout << Q[i][j][k] << " ";
			}
			// cout << endl;
		}
		// cout << endl;
	}

	int column = this->query_selected_label_size; // actual value for column

	/*Randomly choose |K_Label| labels to build the MQ_(K_HOP)*/
	int length = this->query->VLabelSet.size();
	int *Randombit = new int[length];
	int location;
	if (length < K_LABEL) // use all the labels in Q
	{
		for (int i = 0; i < length; i++)
			Randombit[i] = 1;
	}
	else // use K_LABEL labels in Q
	{
		for (int i = 0; i < length; i++)
			Randombit[i] = 0;
		for (int i = 0; i < K_LABEL; i++)
		{
			location = rand() % length;
			while (Randombit[location] == 1)
				location = ((location + 1) % length);
			Randombit[location] = 1;
		}
	}

	/*obtain the randomly chosen labels of Q*/
	location = 0;
	int pointer = 0;
	for (auto it1 = this->query->VLabelSet.begin(); it1 != this->query->VLabelSet.end(); it1++)
	{
		if (Randombit[location] == 1)
		{
			Column[pointer] = *it1;
			pointer++;
		}
		location++;
	}

	for (int i = 0; i < hoplength; i++)
	{
		Q_Child[i] = new int *[query_size];
		Q_Parent[i] = new int *[query_size];
		EM_Q_Child[i] = new Ciphertext *[query_size];
		EM_Q_Parent[i] = new Ciphertext *[query_size];
		for (int j = 0; j < query_size; j++)
		{
			Q_Child[i][j] = new int[query_size];
			Q_Parent[i][j] = new int[query_size];
			EM_Q_Child[i][j] = new Ciphertext[column];
			EM_Q_Parent[i][j] = new Ciphertext[column];
			for (int k = 0; k < column; k++)
			{
				mpz_init(EM_Q_Child[i][j][k]);
				mpz_init(EM_Q_Parent[i][j][k]);
				cgbe->setvalue(EM_Q_Child[i][j][k], 1);
				cgbe->setvalue(EM_Q_Parent[i][j][k], 1);
				Q_Child[i][j][k] = 1;
				Q_Parent[i][j][k] = 1;

				for (int l = 0; l < query_size; l++)
				{

					if ((Q[i][j][l] == 1) && (Q[i][l][j] == 1)) //   \overline{M_Q}
						continue;
					if (this->query->getVLabel(this->query->matrix[l]) == Column[k])
					{
						if (Q[i][j][l] == 0)
						{
							cgbe->setvalue(EM_Q_Child[i][j][k], cgbe->encoding);
							Q_Child[i][j][k] = 0;
						}

						if (Q[i][l][j] == 0)
						{
							cgbe->setvalue(EM_Q_Parent[i][j][k], cgbe->encoding);
							Q_Parent[i][j][k] = 0;
						}
					}
				}
				cgbe->encrypt(EM_Q_Child[i][j][k], EM_Q_Child[i][j][k]);
				cgbe->encrypt(EM_Q_Parent[i][j][k], EM_Q_Parent[i][j][k]);
				// cout<<endl;
			}
			// cout<<"The next Matrix \n"<<endl;
		}
	}

	delete[] Randombit;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Q_Replace_Table_Child(int ***Q, Ciphertext **EM_Q, Replace_Table &Origin, Replace_Table &Replace, CGBE *cgbe, Ciphertext &cipher_one, Ciphertext &cipher_zero, int recursion_length, int vertex, int value, Ciphertext &cipher)
{

	if (vertex == -1)
	{
		Ciphertext temp;
		mpz_init(temp);
		for (int i = 0; i < query_size; i++)
		{
			GlobalCountForBuildingReplaceTable = 0;
			cgbe->setvalue(temp, 1);
			Q_Replace_Table_Child(Q, EM_Q, Origin, Replace, cgbe, cipher_one, cipher_zero, 0, i, 1, temp);
		}
		mpz_clear(temp);
		return;
	}

	// if(recursion_length == (query_size - 1)){
	if (recursion_length == query_size)
	{
		cgbe->setvalue(Origin[vertex][GlobalCountForBuildingReplaceTable], cipher);
		if (value == 1)
			cgbe->setvalue(Replace[vertex][GlobalCountForBuildingReplaceTable], cipher_one);
		else
			cgbe->setvalue(Replace[vertex][GlobalCountForBuildingReplaceTable], cipher_zero);
		// <# for the vertex, # for the case, Ciphertext>
		GlobalCountForBuildingReplaceTable++;
		return;
	}
	else
	{
		int location;
		Ciphertext temp1, temp2;
		mpz_init(temp1);
		mpz_init(temp2);
		// if(recursion_length < vertex)
		//	location = recursion_length;
		// else
		//	location = recursion_length + 1;
		location = recursion_length;

		int temp_value;
		temp_value = value * Q[0][vertex][location];
		if (temp_value > 1)
			temp_value = 1;

		// case 0 for the location-th vertex
		cgbe->add(temp1, EM_Q[vertex][location], cipher_zero);
		// cgbe->setvalue(temp1, EM_Q[vertex][location]);
		cgbe->mul(temp1, temp1, cipher);
		Q_Replace_Table_Child(Q, EM_Q, Origin, Replace, cgbe, cipher_one, cipher_zero, (recursion_length + 1), vertex, temp_value, temp1);

		// case 1 for the location-th vertex
		cgbe->add(temp2, EM_Q[vertex][location], cipher_one);
		cgbe->mul(temp2, temp2, cipher);
		Q_Replace_Table_Child(Q, EM_Q, Origin, Replace, cgbe, cipher_one, cipher_zero, (recursion_length + 1), vertex, value * 1, temp2);

		mpz_clear(temp1);
		mpz_clear(temp2);
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Q_Replace_Table_Parent(int ***Q, Ciphertext **EM_Q, Replace_Table &Origin, Replace_Table &Replace, CGBE *cgbe, Ciphertext &cipher_one, Ciphertext &cipher_zero, int recursion_length, int vertex, int value, Ciphertext &cipher)
{

	if (vertex == -1)
	{
		Ciphertext temp;
		mpz_init(temp);
		for (int i = 0; i < query_size; i++)
		{
			GlobalCountForBuildingReplaceTable = 0;
			cgbe->setvalue(temp, 1);
			Q_Replace_Table_Parent(Q, EM_Q, Origin, Replace, cgbe, cipher_one, cipher_zero, 0, i, 1, temp);
		}
		mpz_clear(temp);
		return;
	}

	// if(recursion_length == (query_size - 1)){
	if (recursion_length == query_size)
	{
		cgbe->setvalue(Origin[vertex][GlobalCountForBuildingReplaceTable], cipher);
		if (value == 1)
			cgbe->setvalue(Replace[vertex][GlobalCountForBuildingReplaceTable], cipher_one);
		else
			cgbe->setvalue(Replace[vertex][GlobalCountForBuildingReplaceTable], cipher_zero);
		// <# for the vertex, # for the case, Ciphertext>
		GlobalCountForBuildingReplaceTable++;
		return;
	}
	else
	{
		int location;
		Ciphertext temp1, temp2;
		mpz_init(temp1);
		mpz_init(temp2);
		// if(recursion_length < vertex)
		//	location = recursion_length;
		// else
		//	location = recursion_length + 1;
		location = recursion_length;

		int temp_value;
		temp_value = value * Q[0][location][vertex];
		if (temp_value > 1)
			temp_value = 1;

		// case 0 for the location-th vertex
		cgbe->add(temp1, EM_Q[location][vertex], cipher_zero);
		// cgbe->setvalue(temp1, EM_Q[location][vertex]);
		cgbe->mul(temp1, temp1, cipher);
		Q_Replace_Table_Parent(Q, EM_Q, Origin, Replace, cgbe, cipher_one, cipher_zero, (recursion_length + 1), vertex, temp_value, temp1);

		// case 1 for the location-th vertex
		cgbe->add(temp2, EM_Q[location][vertex], cipher_one);
		cgbe->mul(temp2, temp2, cipher);
		Q_Replace_Table_Parent(Q, EM_Q, Origin, Replace, cgbe, cipher_one, cipher_zero, (recursion_length + 1), vertex, value * 1, temp2);

		mpz_clear(temp1);
		mpz_clear(temp2);
	}
}

template <class VLabelType, class ELabelType>
int PPMATCH<VLabelType, ELabelType>::ConstructBallMatrix(int **M_B, VertexID *Matrix_ball, int size)
{
	int edge_num = 0;
	for (int i = 0; i < size; i++)
	{
		M_B[i] = new int[size];
		/*intialize matrix M_B*/
		for (int j = 0; j < size; j++)
		{
			M_B[i][j] = 0;
			if (this->graph->isEdge(Matrix_ball[i], Matrix_ball[j])){
				M_B[i][j] = 1;
				edge_num++;
			}
		}
	}
	return edge_num;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::BuildPathIndex(Path_Label &PathLabel, Path_Num &PathNum, Path_Index &PathIndex1, Path_Index &PathIndex2, CGBE *cgbe, int pathlength, int templength, VLabelType *temppath)
{
	// intial each path
	if (templength == pathlength)
	{
		int sum = 0;
		for (int i = 1; i < pathlength; i++)
		{
			int temp = 1;
			for (int j = 1; j < i; j++)
			{
				temp = temp * PathNum[pathlength][-1];
			}
			sum += temp * PathNum[pathlength][temppath[i]];
			// cout<<"Path["<<i<<"]:"<<temppath[i]<<endl;
		}

		for (auto it = this->query->getVLabel().begin(); it != this->query->getVLabel().end(); it++)
		{
			bool flag = true;
			for (int k = 0; k < pathlength; k++)
			{
				if (it->second == temppath[k])
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				mpz_init(PathIndex1[pathlength][it->first][sum]);
				cgbe->setvalue(PathIndex1[pathlength][it->first][sum], 1);
				mpz_init(PathIndex2[pathlength][it->first][sum]);
				cgbe->setvalue(PathIndex2[pathlength][it->first][sum], 1);
			}
		}
		return;
	}

	// encryption
	if (templength == -1)
	{
		for (auto it1 = PathIndex1[pathlength].begin(); it1 != PathIndex1[pathlength].end(); it1++)
		{
			for (auto it2 = PathIndex1[pathlength][it1->first].begin(); it2 != PathIndex1[pathlength][it1->first].end(); it2++)
			{
				cgbe->encrypt(it2->second, it2->second);
			}
		}

		for (auto it3 = PathIndex2[pathlength].begin(); it3 != PathIndex2[pathlength].end(); it3++)
		{
			for (auto it4 = PathIndex2[pathlength][it3->first].begin(); it4 != PathIndex2[pathlength][it3->first].end(); it4++)
			{
				cgbe->encrypt(it4->second, it4->second);
			}
		}

		return;
	}

	// initial all the possible path
	if (templength == 0)
	{
		int labelnum = 1;
		VLabelType templabel;
		for (auto it = this->query->getVLabel().begin(); it != this->query->getVLabel().end(); it++)
		{
			templabel = it->second;
			bool flag = true;
			for (int i = 1; i < labelnum; i++)
			{
				if (templabel == PathLabel[pathlength][i])
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				PathLabel[pathlength][labelnum] = templabel;
				PathNum[pathlength][templabel] = labelnum;
				labelnum++;
			}
		}

		// PathLabel[-1] = labelnum;
		PathNum[pathlength][-1] = labelnum;

		VLabelType *path = new VLabelType[pathlength];
		for (int i = 0; i < pathlength; i++)
		{
			path[i] = -1;
		}
		for (auto it = PathLabel[pathlength].begin(); it != PathLabel[pathlength].end(); it++)
		{
			// if(it->first == -1)
			// continue;
			path[0] = it->second;
			BuildPathIndex(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, 1, path);
		}

		// for(auto it1 = PathNum.begin();it1 !=PathNum.end();it1++){
		// cout<<"Vertex " << it1->first << ": " << it1->second << endl;
		// }

		FindQueryPath(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, 0, nullptr, 0, 0, 0);

		BuildPathIndex(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, -1, path);
		delete[] path;
	}
	else
	{
		for (auto it2 = PathLabel[pathlength].begin(); it2 != PathLabel[pathlength].end(); it2++)
		{
			bool flag = true;
			for (int i = 0; i < templength; i++)
			{
				if (it2->second == temppath[i])
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				temppath[templength] = it2->second;
				BuildPathIndex(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, templength + 1, temppath);
			}
		}
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::FindQueryPath(Path_Label &PathLabel, Path_Num &PathNum, Path_Index &PathIndex1, Path_Index &PathIndex2, CGBE *cgbe, int pathlength, int templength, VLabelType *temppath, VertexID vertex, VertexID tempvertex, int point)
{
	if (templength == pathlength)
	{
		int sum = 0;
		for (int i = 1; i < pathlength; i++)
		{
			int temp = 1;
			for (int j = 1; j < i; j++)
			{
				temp = temp * PathNum[pathlength][-1];
			}
			sum += temp * PathNum[pathlength][temppath[i]];
		}
		if (point == 0)
			cgbe->setvalue(PathIndex1[pathlength][vertex][sum], cgbe->encoding);
		if (point == 1)
			cgbe->setvalue(PathIndex2[pathlength][vertex][sum], cgbe->encoding);
		return;
	}

	if (templength == 0)
	{
		VLabelType *path1 = new VLabelType[pathlength]; // ?????may need two paths?
		VLabelType *path2 = new VLabelType[pathlength];
		for (int i = 0; i < pathlength; i++)
		{
			path1[i] = -1;
			path2[i] = -1;
		}

		for (auto it = this->query->getVLabel().begin(); it != this->query->getVLabel().end(); it++)
		{
			path1[0] = this->query->getVLabel(it->first);
			FindQueryPath(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, 1, path1, it->first, it->first, 0);
			path2[0] = this->query->getVLabel(it->first);
			FindQueryPath(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, 1, path2, it->first, it->first, 1);
		}
		delete[] path1;
		delete[] path2;
	}
	else
	{
		// for child path
		if (point == 0)
		{
			for (auto it = this->query->getOutEdge()[tempvertex].begin(); it != this->query->getOutEdge()[tempvertex].end(); it++)
			{
				bool flag = true;
				for (int i = 0; i < templength; i++)
				{
					if (this->query->getVLabel(it->first) == temppath[i])
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					temppath[templength] = this->query->getVLabel(it->first);
					FindQueryPath(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, templength + 1, temppath, vertex, it->first, 0);
				}
			}
		}

		// for parent path
		if (point == 1)
		{
			for (auto it = this->query->getInVertex()[tempvertex].begin(); it != this->query->getInVertex()[tempvertex].end(); it++)
			{
				bool flag = true;
				for (int i = 0; i < templength; i++)
				{
					if (this->query->getVLabel(it->first) == temppath[i])
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					temppath[templength] = this->query->getVLabel(it->first);
					FindQueryPath(PathLabel, PathNum, PathIndex1, PathIndex2, cgbe, pathlength, templength + 1, temppath, vertex, it->first, 1);
				}
			}
		}
	}
}

template <class VLabelType, class ELabelType>
int PPMATCH<VLabelType, ELabelType>::CantorExpansion(int length, int *array)
{
	int fac[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320}; // n!
	int cnt, sum;
	sum = 0;
	for (int i = 0; i < length; i++)
	{
		cnt = 0;
		for (int j = i + 1; j < length; j++)
			if (array[j] < array[i])
				cnt++;
		sum += cnt * fac[length - i - 1];
	}
	return sum;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::CantorExpansionDecode(int length, int *array, int value)
{
	int i, j, t;
	const int fac[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320}; // n!

	bool vis[10];
	memset(vis, 0, sizeof(vis));
	value--;
	for (i = 0; i < length; ++i)
	{
		t = value / fac[length - i - 1];
		for (j = 1; j <= length; j++)
			if (!vis[j])
			{
				if (t == 0)
					break;
				t--;
			}
		array[i] = j, vis[j] = true;
		value %= fac[length - i - 1]; // remainder
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::ConstructBallIndex(B_Index &Ball_Index, B_Index &Ball_Index_Reverse, int **M_B, VertexID *Matrix_ball, int size)
{

	/*intialize matrix temp for saving the k_hop neighbors*/
	int ***temp = new int **[hoplength];
	for (int i = 1; i < hoplength; i++)
	{
		temp[i] = new int *[size];
		for (int j = 0; j < size; j++)
		{
			temp[i][j] = new int[size];
			for (int k = 0; k < size; k++)
				temp[i][j][k] = 0;
		}
	}

	for (int i = 0; i < hoplength; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				// no cycles
				// if(j == k){
				//	continue;
				// }
				if (0 < i && i < hoplength)
				{
					for (int l = 0; l < size; l++)
					{
						if (i == 1)
						{
							temp[i][j][k] += (M_B[j][l] * M_B[l][k]);
						}
						else
						{
							temp[i][j][k] += (temp[i - 1][j][l] * M_B[l][k]);
						}

						if (temp[i][j][k] > 0)
						{
							temp[i][j][k] = 1;
							if (Ball_Index[i][j].find(this->graph->getVLabel(Matrix_ball[k])) == Ball_Index[i][j].end())
							{
								Ball_Index[i][j][this->graph->getVLabel(Matrix_ball[k])] = 1;
							}

							if (Ball_Index_Reverse[i][k].find(this->graph->getVLabel(Matrix_ball[j])) == Ball_Index_Reverse[i][k].end())
							{
								Ball_Index_Reverse[i][k][this->graph->getVLabel(Matrix_ball[j])] = 1; // both parent and children
							}
						}
						else
						{
							if (temp[i][j][k] != 0)
								cout << "The ball index constructing is wrong!" << endl;
						}
					}
				}
				else
				{
					if (M_B[j][k] > 0) // The "i = 0" is used for the 1-hop
					{
						if (Ball_Index[i][j].find(this->graph->getVLabel(Matrix_ball[k])) == Ball_Index[i][j].end())
						{
							Ball_Index[i][j][this->graph->getVLabel(Matrix_ball[k])] = 1;
						}

						if (Ball_Index_Reverse[i][k].find(this->graph->getVLabel(Matrix_ball[j])) == Ball_Index_Reverse[i][k].end())
						{
							Ball_Index_Reverse[i][k][this->graph->getVLabel(Matrix_ball[j])] = 1; // both parent and children
						}
					}
					else
					{
						if (M_B[j][k] != 0)
							cout << "The ball index constructing is wrong!" << endl;
					}
				}
			}
		}
	}

	for (int i = 1; i < hoplength; i++)
	{
		for (int j = 0; j < size; j++)
			delete[] temp[i][j];
		delete[] temp[i];
	}
	delete[] temp;
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::ConstructNLIndex(B_Index &Ball_Index, B_Index &Ball_Index_Reverse, DIGRAPH<VLabelType, ELabelType> *Ball, VertexID origin, VertexID place, int hop, unordered_map<VertexID, int> &Matrix_Ball)
{
	if (hop == hoplength)
		return;

	if (hop == -1)
	{
		for (auto it1 = Ball->getVLabel().begin(); it1 != Ball->getVLabel().end(); it1++)
		{
			for (auto it2 = Ball->getOutEdge()[it1->first].begin(); it2 != Ball->getOutEdge()[it1->first].end(); it2++)
			{
				if (Ball_Index[0][Matrix_Ball[it1->first]].find(Ball->getVLabel(it2->first)) == Ball_Index[0][Matrix_Ball[it1->first]].end())
				{
					Ball_Index[0][Matrix_Ball[it1->first]][Ball->getVLabel(it2->first)] = 1;
				}

				if (Ball_Index_Reverse[0][Matrix_Ball[it2->first]].find(Ball->getVLabel(it1->first)) == Ball_Index_Reverse[0][Matrix_Ball[it2->first]].end())
				{
					Ball_Index_Reverse[0][Matrix_Ball[it2->first]][Ball->getVLabel(it1->first)] = 1; // both parent and children
				}
				ConstructNLIndex(Ball_Index, Ball_Index_Reverse, Ball, it1->first, it2->first, 1, Matrix_Ball);
			}
		}
	}
	else
	{
		for (auto it = Ball->getOutEdge()[place].begin(); it != Ball->getOutEdge()[place].end(); it++)
		{
			if (Ball_Index[hop][Matrix_Ball[origin]].find(Ball->getVLabel(it->first)) == Ball_Index[hop][origin].end())
			{
				Ball_Index[hop][Matrix_Ball[origin]][Ball->getVLabel(it->first)] = 1;
			}

			if (Ball_Index_Reverse[hop][Matrix_Ball[it->first]].find(Ball->getVLabel(origin)) == Ball_Index_Reverse[hop][Matrix_Ball[it->first]].end())
			{
				Ball_Index_Reverse[hop][Matrix_Ball[it->first]][Ball->getVLabel(origin)] = 1; // both parent and children
			}
			ConstructNLIndex(Ball_Index, Ball_Index_Reverse, Ball, origin, it->first, hop + 1, Matrix_Ball);
		}
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::LoadNLIndex(B_Index &Graph_Index, B_Index &Graph_Index_Reverse, string Index1, string Index2)
{
	int hop;
	VertexID vertex;
	VLabelType label;
	char str[100];
	string temp;

	ifstream OpenFile1(Index1);
	if (!(OpenFile1.getline(str, 100)))
	{
		cout << "The Index cannot be loaded!" << endl;
		return;
	}
	while (!OpenFile1.eof())
	{
		OpenFile1 >> temp;
		if (temp == "$")
		{
			// cout<<temp<<endl;
			OpenFile1 >> temp;
			// cout<<temp<<endl;
			hop = stol(temp);
		}
		else if (temp == "#")
		{
			OpenFile1 >> temp;
			vertex = stol(temp);
		}
		else
		{
			// cout<<temp<<endl;
			label = stol(temp);
			Graph_Index[hop][vertex][label] = 1;
		}
	}
	OpenFile1.close();

	ifstream OpenFile2(Index2);
	if (!(OpenFile2.getline(str, 100)))
	{
		cout << "The Index cannot be loaded!" << endl;
		return;
	}
	while (!OpenFile2.eof())
	{
		OpenFile2 >> temp;
		if (temp == "$")
		{
			OpenFile2 >> temp;
			hop = stol(temp);
		}
		else if (temp == "#")
		{
			OpenFile2 >> temp;
			vertex = stol(temp);
		}
		else
		{
			label = stol(temp);
			Graph_Index_Reverse[hop][vertex][label] = 1;
		}
	}
	OpenFile2.close();
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::BF(BF_value **QueryBF, BloomFilter *ballBF, int querysize, VertexID center, int pointer)
{
	clock_t startTime, endTime;
	bool flag1 = false;
	bool flag2;
	startTime = clock();
	char buffer[ENTRY_VALUE_LEN+1];
	size_t len = ENTRY_VALUE_LEN+1;
	for (int i = 0; i < querysize; i++)
	{
		if (this->query->getVLabel(this->query->matrix[i]) != this->graph->getVLabel(center))
			continue;
		flag2 = true;
		for (auto it1 = QueryBF[i]->begin(); it1 != QueryBF[i]->end(); it1++)
		{
			for (auto it = it1->second.begin(); it != it1->second.end(); it++)
			{
				if (*it == -1)
					continue;
				sprintf(buffer, "%ld", *it);
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
		if (this->result_EncSSim[pointer] == 1)
			BF_Improved++;
		this->result_EncSSim[pointer] = 0;
		// this->OneIterNL[pointer] = 0;
		// this->EncSpecial[pointer] = 0;
		this->result_BF[pointer] = 0;
	}
	endTime = clock();
	this->BF_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
}



template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::BFTest(BF_value **QueryBF, BF_value *BallBF, int querysize, VertexID center, int pointer)
{
	clock_t startTime, endTime;
	bool flag1 = false;
	bool flag2;
	startTime = clock();
	for (int i = 0; i < querysize; i++)
	{
		if (this->query->getVLabel(this->query->matrix[i]) != this->graph->getVLabel(center))
			continue;
		flag2 = true;
		for (auto it1 = QueryBF[i]->begin(); it1 != QueryBF[i]->end(); it1++)
		{
			for (auto it = it1->second.begin(); it != it1->second.end(); it++)
			{
				if (*it == -1)
					continue;
				for(auto it2 = BallBF->begin();it2!=BallBF->end();it2++){
					if(it1->first!=it2->first)
						continue;
					if(it2->second.find(*it)==it2->second.end()){
						//cout<<"center: "<<center<<"	Missed value: "<<*it;
						flag2 = false;
						break;
					}
				}
				if(!flag2)
				{
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
		if (this->result_EncSSim[pointer] == 1)
			BF_Improved++;
		this->result_EncSSim[pointer] = 0;
		// this->OneIterNL[pointer] = 0;
		// this->EncSpecial[pointer] = 0;
		this->result_BF[pointer] = 0;
	}
	endTime = clock();
	this->BF_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
}




template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::constructBF(BF_value *BF, DIGRAPH<VLabelType, ELabelType> *Graph, VertexID center, int type, Degree *left, Degree *right, VLabelType *temptree, long *weight, int& num)  //Outside SGX for test
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
			if (it->second.find(value) == it->second.end()){
				it->second.insert(value);
				num++;
			}

			// cout << "center " << center << " with " << it->first << " pattern" << endl;
			// for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++)
			// cout << *it1 << endl;
		}

		return;
	}

	if (type == 6) // BF[1]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		VLabelType *labels = new VLabelType[6];
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
					constructBF(BF, Graph, center, -6, nullptr, nullptr, labels, weight, num); // compute the weighted value
				}
				patternlabel->erase(Graph->getVLabel(it2->first));
			}
			patternlabel->erase(Graph->getVLabel(it1->first));
		}

		delete patternlabel;
		delete leftchildlabel;
		delete labels;
		// delete neighbor;
		return;
	}

	if (type == 7) // BF[2]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		VLabelType *labels = new VLabelType[6];
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
				constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight, num);

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
						constructBF(BF, Graph, center, -5, nullptr, nullptr, labels, weight, num); // compute the weighted value
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
		return;
	}

	if (type == 9) // BF[3]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		PatternLabel *rightchildlabel = new PatternLabel;
		Degree *degree0 = new Degree;
		Degree *degree1 = new Degree;
		VLabelType *labels = new VLabelType[6];
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
				constructBF(BF, Graph, center, 6, degree0, degree1, nullptr, weight, num);
				constructBF(BF, Graph, center, 7, degree1, degree0, nullptr, weight, num);

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
							constructBF(BF, Graph, center, -4, nullptr, nullptr, labels, weight, num); // compute the weighted value
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
		return;
	}

	if (type == 10) // BF[4]
	{
		PatternLabel *patternlabel = new PatternLabel;
		PatternLabel *leftchildlabel = new PatternLabel;
		PatternLabel *rightchildlabel = new PatternLabel;
		Degree *degree0 = new Degree;
		Degree *degree1 = new Degree;
		Degree *degree2 = new Degree;
		VLabelType *labels = new VLabelType[6];
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
				constructBF(BF, Graph, center, 9, degree1, degree0, nullptr, weight, num);
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
								constructBF(BF, Graph, center, -3, nullptr, nullptr, labels, weight, num); // compute the weighted value
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

		return;
	}

	if (type == -1) // initial, undirected tree pattern
	{
		num = 0;
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

		if ((degree0->size() > 0) && (degree1->size() > 0))
			constructBF(BF, Graph, center, 6, degree1, degree0, nullptr, weight, num);

		if ((degree0->size() > 0) && (degree2->size() > 0))
			constructBF(BF, Graph, center, 7, degree2, degree0, nullptr, weight, num);

		if ((degree1->size() > 0) && (degree2->size() > 0))
			constructBF(BF, Graph, center, 9, degree2, degree1, nullptr, weight, num);

		if (degree2->size() > 1)
			constructBF(BF, Graph, center, 10, degree2, degree2, nullptr, weight, num);

		delete degree0;
		delete degree1;
		delete degree2;
		delete patternlabel;
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::Enc_Twig(Path_Label &TwigLabel, Path_Num &TwigNum, Path_Index &TwigIndex1, Path_Index &TwigIndex2, int pointer, int **M_B, VertexID *Matrix_ball, int size_ball, CGBE *cgbe, double &Time, unordered_map<int, double>& TwigTime, unordered_map<int, double>& TwigBallReadTime, int center, VertexID centerNo, int twiglength, int decryptionNum, Ciphertext &cipher_one)
{
	int RevisedLengthOfM_B, TotalCantorNum;
	int decodeNum = 0;
	// VLabelType

	clock_t startTime, endTime;
	// startTime = clock();

	Ciphertext TwigPruning;
	Ciphertext Result;
	Ciphertext FinalResult;
	mpz_init(TwigPruning);
	mpz_init(Result);
	mpz_init(FinalResult);
	cgbe->setvalue(TwigPruning, 1);
	cgbe->setvalue(Result, 1);
	cgbe->setvalue(FinalResult, 0);

	////////
	// Time for reading the balls  (todo)
	////////
	startTime = clock();
	DIGRAPH<VLabelType, ELabelType> *RBall = new DIGRAPH<VLabelType, ELabelType>;
	for (int i = 0; i < size_ball; i++)
	{
		RBall->insertVertex(Matrix_ball[i], 0);
	}

	for (auto it3 = RBall->getVLabel().begin(); it3 != RBall->getVLabel().end(); it3++)
	{
		for(auto it4 = this->graph->getOutEdge()[it3->first].begin(); it4!= this->graph->getOutEdge()[it3->first].end(); it4++){
			if(RBall->isVertex(it4->first))
				RBall->insertEdge(it3->first, it4->first, 0);
		}

		for(auto it4 = this->graph->getInVertex()[it3->first].begin(); it4!= this->graph->getInVertex()[it3->first].end(); it4++){
			if(RBall->isVertex(it4->first))
				RBall->insertEdge(it4->first, it3->first, 0);
		}

		/*
		for (auto it4 = RBall->getVLabel().begin(); it4 != RBall->getVLabel().end(); it4++)
		{
			if (it3->first == it4->first)
				continue;
			if (this->graph->isEdge(it3->first, it4->first))
				RBall->insertEdge(it3->first, it4->first, 0);
			if (this->graph->isEdge(it4->first, it3->first))
				RBall->insertEdge(it4->first, it3->first, 0);
		}*/
	}
	endTime = clock();
	TwigBallReadTime[pointer] = (double)(endTime - startTime) / CLOCKS_PER_SEC;

	startTime = clock();
	// Save the twigs of center
	unordered_map<VertexID, unordered_map<int, int>> Twig_Center_Child, Twig_Center_Parent;

	// Find the twigs of center
	TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, nullptr, nullptr, RevisedLengthOfM_B, twiglength, 0, centerNo, centerNo, -1);
	endTime = clock();
	this->Twig_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	//this->Temp_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;
	TwigTime[pointer] += (double)(endTime - startTime) / CLOCKS_PER_SEC;
	/*
	cout << "#######################" << endl;
	cout << "#######################" << endl;

	for(auto it = Twig_Center_Child.begin();it != Twig_Center_Child.end(); it++){
		cout << "Vertex " << it->first <<":" << endl;
		for(auto it2 = Twig_Center_Child[it->first].begin(); it2 != Twig_Center_Child[it->first].end(); it2++){
			cout <<"Sum:" << it2->first<<endl;
		}

	}*/

	// Match
	cgbe->setCombinedPrivateKey_Twig(decryptionNum);

	for (auto i = this->query->getVLabel().begin(); i != this->query->getVLabel().end(); i++)
	{

		if (this->query->getVLabel(i->first) != this->graph->getVLabel(centerNo))
			continue;

		for (auto it = TwigIndex1[twiglength][i->first].begin(); it != TwigIndex1[twiglength][i->first].end(); it++)
		{
			startTime = clock();
			if (Twig_Center_Child[centerNo].find(it->first) == Twig_Center_Child[centerNo].end())
			{
				decodeNum++;
				cgbe->mul(TwigPruning, TwigPruning, TwigIndex1[twiglength][i->first][it->first]);				
			}
			endTime = clock();
			this->Twig_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;
			TwigTime[pointer] += (double)(endTime - startTime) / CLOCKS_PER_SEC;

			if (decodeNum == decryptionNum)
			{
				startTime = clock();
				cgbe->decryption_Twig(TwigPruning, TwigPruning);
				if (cgbe->isZero(TwigPruning))
					cgbe->mul(Result, Result, cgbe->zero);
				cgbe->setvalue(TwigPruning, 1);
				decodeNum = 0;
				endTime = clock();
				this->Decrypt_Twig += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			}
		}

		for (auto it2 = TwigIndex2[twiglength][i->first].begin(); it2 != TwigIndex2[twiglength][i->first].end(); it2++)
		{

			startTime = clock();
			if (Twig_Center_Parent[centerNo].find(it2->first) == Twig_Center_Parent[centerNo].end())
			{
				decodeNum++;
				cgbe->mul(TwigPruning, TwigPruning, TwigIndex2[twiglength][i->first][it2->first]);				
			}
			endTime = clock();
			this->Twig_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;
			TwigTime[pointer] += (double)(endTime - startTime) / CLOCKS_PER_SEC;

			if (decodeNum == decryptionNum)
			{
				startTime = clock();
				cgbe->decryption_Twig(TwigPruning, TwigPruning);
				if (cgbe->isZero(TwigPruning))
					cgbe->mul(Result, Result, cgbe->zero);
				cgbe->setvalue(TwigPruning, 1);
				decodeNum = 0;
				endTime = clock();
				this->Decrypt_Twig += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			}
		}

		while (decodeNum != decryptionNum)
		{
			//precomputing the decryption key for different powers


			//startTime = clock();
			decodeNum++;
			cgbe->mul(TwigPruning, TwigPruning, cipher_one);
			//endTime = clock();
			//this->Twig_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;
			//Time += (double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC;
		}

		startTime = clock();
		cgbe->decryption_Twig(TwigPruning, TwigPruning);
		if (cgbe->isZero(TwigPruning))
			cgbe->mul(Result, Result, cgbe->zero);
		endTime = clock();
		this->Decrypt_Twig += (double)(endTime - startTime) / CLOCKS_PER_SEC;

		cgbe->add(FinalResult, FinalResult, Result);
		cgbe->setvalue(TwigPruning, 1);
		cgbe->setvalue(Result, 1);
		decodeNum = 0;
	}

	// endTime = clock();
	// this->Twig_Time += (double)(endTime - startTime) / CLOCKS_PER_SEC;

	/*		Client		*/

	if ((cgbe->isZero(FinalResult)))
	{
		if (this->result_EncSSim[pointer] == 1)
			Twig_Improved++;
		this->result_EncSSim[pointer] = 0;
		// this->OneIterNL[pointer] = 0;
		// this->EncSpecial[pointer] = 0;
		this->result_Twig[pointer] = 0;
	}

	mpz_clear(TwigPruning);
	mpz_clear(Result);
	mpz_clear(FinalResult);
	Twig_Center_Child.clear();
	Twig_Center_Parent.clear();
	delete RBall;
}











template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::TwigMatchAll(DIGRAPH<VLabelType, ELabelType> *RBall, Path_Num &TwigNum, unordered_map<VertexID, unordered_map<int, int>> &Twig_Center_Child, unordered_map<VertexID, unordered_map<int, int>> &Twig_Center_Parent, VLabelType *Twig, VLabelType *templabel, int RevisedBallSize, int TwigLength, int iter, VertexID center, VertexID pointer, int point)
{
	//cout<<"iter:"<<iter<<endl;
	if (iter == (TwigLength + 1))
	{
		int sum1 = 0; // twig
		int sum2 = 0; // only left child
		int sum3 = 0; // only right child
		// int sum4 = 0;	//no children
		bool flag = true;

		for (int i = 1; i < (TwigLength + 1); i++)
		{
			//int temp = 1;
			if (TwigNum[TwigLength][Twig[i]] == 0 && (i != TwigLength))
				cout << "label 0 happen when generating the index!" << endl;

			//for (int j = 1; j < i; j++)
			//{
			//	temp = temp * TwigNum[TwigLength][-1];
			//}

			// if(i<TwigLength-1){
			//	sum4 += temp*TwigNum[TwigLength][Twig[i]];
			// }

			if (i < (TwigLength - 1))
			{
				sum2 += TwigNum[twiglength][i*(-1)] * TwigNum[TwigLength][Twig[i]];
				sum3 += TwigNum[twiglength][i*(-1)] * TwigNum[TwigLength][Twig[i]];
			}

			if (i == (TwigLength - 1))
			{
				if (Twig[i + 1] == 0)
					flag = false;
				sum2 += TwigNum[twiglength][i*(-1)] * TwigNum[TwigLength][Twig[i]];
				sum3 += TwigNum[twiglength][i*(-1)] * TwigNum[TwigLength][Twig[i + 1]]; // the last label may be 0 due to repeated label
			}

			sum1 += TwigNum[twiglength][i*(-1)] * TwigNum[TwigLength][Twig[i]];
			//cout<<"Twig["<<i<<"]:"<<Twig[i]<<endl;
		}

		// cout << sum1 << "	" << sum2 << "	" << sum3 << endl;

		if (point == 0)
		{
			if (flag)
			{
				Twig_Center_Child[center][sum2] = 1;
				Twig_Center_Child[center][sum3] = 1;
			}
			Twig_Center_Child[center][sum1] = 1;
		}
		if (point == 1)
		{
			if (flag)
			{
				Twig_Center_Parent[center][sum2] = 1;
				Twig_Center_Parent[center][sum3] = 1;
			}
			Twig_Center_Parent[center][sum1] = 1;
		}
		return;
	}

	bool flag, flag1, flag2;

	if (iter == 0)
	{
		VLabelType *temptwig = new VLabelType[TwigLength + 1];
		VLabelType *tempLabel = new VLabelType[this->labelsize];
		//unordered_map<VLabelType, int> LabelFlag;
		//for(auto it = ) later to quick check label existence

		// for child
		//cout<<"center label:" << this->graph->getVLabel(center);
		temptwig[0] = this->graph->getVLabel(center);
		TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, temptwig, tempLabel, RevisedBallSize, TwigLength, iter + 1, center, center, 0);

		// for parent
		temptwig[0] = this->graph->getVLabel(center);
		TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, temptwig, tempLabel, RevisedBallSize, TwigLength, iter + 1, center, center, 1);

		delete[] temptwig;
		delete[] tempLabel;
	}
	else
	{
		if (point == 0)
		{

			if (iter < (TwigLength - 1))
			{
				for (auto it = RBall->getOutEdge()[pointer].begin(); it != RBall->getOutEdge()[pointer].end(); it++)
				{
					flag = true;
					for (int i = 0; i < iter; i++)
					{
						if (this->graph->getVLabel(it->first) == Twig[i]) //here
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						Twig[iter] = this->graph->getVLabel(it->first);
						TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 1, center, it->first, 0);
					}
				}
			}
			else
			{ // the twig node

				if (RBall->getOutEdge()[pointer].size() == 0)
					return;

				if (RBall->getOutEdge()[pointer].size() == 1)
				{ // only one child
					for (auto it = RBall->getOutEdge()[pointer].begin(); it != RBall->getOutEdge()[pointer].end(); it++)
					{
						flag = true;
						for (int i = 0; i < iter; i++)
						{
							if (this->graph->getVLabel(it->first) == Twig[i])
							{
								flag = false;
								break;
							}
						}
						if (flag)
						{

							Twig[iter] = this->graph->getVLabel(it->first);
							Twig[iter + 1] = 0;
							TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 2, center, it->first, 0);
						}
					}
				}
				else
				{ // more than two child
					int labellength=0;
					for(int i = 0; i<iter;i++){
						templabel[labellength] = Twig[i];
						labellength++;
					}
					for (auto it1 = RBall->getOutEdge()[pointer].begin(); it1 != RBall->getOutEdge()[pointer].end(); it1++)
					{
						flag1 = true;
						for(int i = 0;i<labellength;i++){
							if(this->graph->getVLabel(it1->first)==templabel[i]){
								flag1 = false;
								break;
							}
						}

						if(flag1){
							templabel[labellength] = this->graph->getVLabel(it1->first);
							labellength++;
						}
					}


					if(iter == labellength)
						return;

					if(iter == (labellength - 1)){
						Twig[iter] = templabel[iter];
						Twig[iter + 1] = 0;
						TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 2, center, 0, 0);						
					}else{
						for(int i = iter;i<labellength;i++)
						{
							for (int j = iter; j<labellength;j++)
							{
								if (TwigNum[twiglength][templabel[i]] <= TwigNum[twiglength][templabel[j]])
									continue;

								Twig[iter] = templabel[i];
								Twig[iter + 1] = templabel[j];
								TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 2, center, 0, 0);						
							}
						}
					}
				}
			}
		}

		if (point == 1)
		{

			if (iter < (TwigLength - 1))
			{
				for (auto it = RBall->getInVertex()[pointer].begin(); it != RBall->getInVertex()[pointer].end(); it++)
				{
					flag = true;
					for (int i = 0; i < iter; i++)
					{
						if (this->graph->getVLabel(it->first) == Twig[i])
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						Twig[iter] = this->graph->getVLabel(it->first);
						TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 1, center, it->first, 1);
					}
				}
			}
			else
			{ // the twig node
				if (RBall->getInVertex()[pointer].size() == 0)
					return;
				if (RBall->getInVertex()[pointer].size() == 1)
				{ // only one child
					for (auto it = RBall->getInVertex()[pointer].begin(); it != RBall->getInVertex()[pointer].end(); it++)
					{
						flag = true;
						for (int i = 0; i < iter; i++)
						{
							if (this->graph->getVLabel(it->first) == Twig[i])
							{
								flag = false;
								break;
							}
						}
						if (flag)
						{
							Twig[iter] = this->graph->getVLabel(it->first);
							Twig[iter + 1] = 0;
							TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 2, center, it->first, 1);
						}
					}
				}
				else
				{ // more than two child
					int labellength=0;					
					for(int i = 0; i<iter;i++){
						templabel[labellength] = Twig[i];
						labellength++;
					}
					for (auto it1 = RBall->getInVertex()[pointer].begin(); it1 != RBall->getInVertex()[pointer].end(); it1++)
					{
						flag1 = true;
						for(int i = 0;i<labellength;i++){
							if(this->graph->getVLabel(it1->first)==templabel[i]){
								flag1 = false;
								break;
							}
						}

						if(flag1){
							templabel[labellength] = this->graph->getVLabel(it1->first);
							labellength++;
						}
					}


					if(iter == labellength)
						return;
					
					if(iter == (labellength - 1)){
						Twig[iter] = templabel[iter];
						Twig[iter + 1] = 0;
						TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 2, center, 0, 1);						
					}else{
						for(int i = iter;i<labellength;i++)
						{
							for (int j = iter; j<labellength;j++)
							{
								if (TwigNum[twiglength][templabel[i]] <= TwigNum[twiglength][templabel[j]])
									continue;							

								Twig[iter] = templabel[i];
								Twig[iter + 1] = templabel[j];
								TwigMatchAll(RBall, TwigNum, Twig_Center_Child, Twig_Center_Parent, Twig, templabel, RevisedBallSize, TwigLength, iter + 2, center, 0, 1);						
							}
						}
					}
				}
			}
		}
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::BuildTwigIndex(Path_Label &TwigLabel, Path_Num &TwigNum, Path_Index &TwigIndex1, Path_Index &TwigIndex2, CGBE *cgbe, int twiglength, int templength, VLabelType *temptwig)
{
	// intial each twig
	if (templength == (twiglength + 1))
	{
		int sum1 = 0; // twig
		int sum2 = 0; // only left child
		int sum3 = 0; // only right child
		// int sum4 = 0;	//no children

		for (int i = 1; i < (twiglength + 1); i++)
		{
			//int temp = 1;
			if (TwigNum[twiglength][temptwig[i]] == 0||temptwig[i]==0)
				cout << "label 0 of Q's index happen!" << endl;

			//for (int j = 1; j < i; j++)
			//{
			//	temp = temp * TwigNum[twiglength][-1];
			//}

			// if(i<twiglength-1){
			//	sum4 += temp*TwigNum[twiglength][temptwig[i]];
			// }

			if (i < (twiglength - 1))
			{
				sum2 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
				sum3 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
			}

			if (i == (twiglength - 1))
			{
				sum2 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
				sum3 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i + 1]];
			}

			sum1 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
			// cout<<"Twig["<<i<<"]:"<<temptwig[i]<<endl;
		}

		// cout << sum1 << "	" << sum2 << "	" << sum3 << endl;

		for (auto it = this->query->getVLabel().begin(); it != this->query->getVLabel().end(); it++)
		{
			bool flag = true;
			for (int k = 0; k < (twiglength + 1); k++)
			{
				if (it->second == temptwig[k])
				{
					flag = false; // check repeated label
					break;
				}
			}
			if (flag)
			{ // sum is the value without the first label, so the values are set to 1 (encoding corresponding to 0)
				mpz_init(TwigIndex1[twiglength][it->first][sum1]);
				cgbe->setvalue(TwigIndex1[twiglength][it->first][sum1], 1);
				mpz_init(TwigIndex1[twiglength][it->first][sum2]);
				cgbe->setvalue(TwigIndex1[twiglength][it->first][sum2], 1);
				mpz_init(TwigIndex1[twiglength][it->first][sum3]);
				cgbe->setvalue(TwigIndex1[twiglength][it->first][sum3], 1);
				// mpz_init(TwigIndex1[twiglength][it->first][sum4]);
				// cgbe->setvalue(TwigIndex1[twiglength][it->first][sum4], 1);

				mpz_init(TwigIndex2[twiglength][it->first][sum1]);
				cgbe->setvalue(TwigIndex2[twiglength][it->first][sum1], 1);
				mpz_init(TwigIndex2[twiglength][it->first][sum2]);
				cgbe->setvalue(TwigIndex2[twiglength][it->first][sum2], 1);
				mpz_init(TwigIndex2[twiglength][it->first][sum3]);
				cgbe->setvalue(TwigIndex2[twiglength][it->first][sum3], 1);
				// mpz_init(TwigIndex2[twiglength][it->first][sum4]);
				// cgbe->setvalue(TwigIndex2[twiglength][it->first][sum4], 1);
			}
		}
		return;
	}

	// encryption
	if (templength == -1)
	{
		for (auto it1 = TwigIndex1[twiglength].begin(); it1 != TwigIndex1[twiglength].end(); it1++)
		{
			for (auto it2 = TwigIndex1[twiglength][it1->first].begin(); it2 != TwigIndex1[twiglength][it1->first].end(); it2++)
			{
				cgbe->encrypt(it2->second, it2->second);
			}
		}

		for (auto it3 = TwigIndex2[twiglength].begin(); it3 != TwigIndex2[twiglength].end(); it3++)
		{
			for (auto it4 = TwigIndex2[twiglength][it3->first].begin(); it4 != TwigIndex2[twiglength][it3->first].end(); it4++)
			{
				cgbe->encrypt(it4->second, it4->second);
			}
		}

		return;
	}

	// initial all the possible twig
	if (templength == 0)
	{
		int labelnum = 1;
		VLabelType templabel;
		for (auto it = this->query->getVLabel().begin(); it != this->query->getVLabel().end(); it++)
		{
			templabel = it->second;
			bool flag = true;
			for (int i = 1; i < labelnum; i++)
			{
				if (templabel == TwigLabel[twiglength][i])
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				TwigLabel[twiglength][labelnum] = templabel;
				TwigNum[twiglength][templabel] = labelnum;
				labelnum++;
			}
		}

		// TwigLabel[-1] = labelnum;
		TwigNum[twiglength][-1] = labelnum;
		for(int i=2; i<(twiglength + 1); i++){
			TwigNum[twiglength][i*(-1)]= TwigNum[twiglength][(i-1)*(-1)]*(labelnum);
		}
		//cout << "size: " << labelnum << endl;

		VLabelType *twig = new VLabelType[twiglength + 1]; // twig has length + 1
		for (int i = 0; i < (twiglength + 1); i++)
		{
			twig[i] = -1;
		}
		for (auto it = TwigLabel[twiglength].begin(); it != TwigLabel[twiglength].end(); it++)
		{
			// if(it->first == -1)
			// continue;
			twig[0] = it->second;
			BuildTwigIndex(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, 1, twig);
		}

		// for(auto it1 = TwigNum.begin();it1 !=TwigNum.end();it1++){
		// cout<<"Vertex " << it1->first << ": " << it1->second << endl;
		// }

		// find the twigs in Query and set the value with encoding
		FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, 0, nullptr, nullptr, 0, 0, 0);

		BuildTwigIndex(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, -1, twig); // Encryption
		delete[] twig;
	}
	else
	{
		for (auto it2 = TwigLabel[twiglength].begin(); it2 != TwigLabel[twiglength].end(); it2++)
		{
			bool flag = true;
			for (int i = 0; i < templength; i++)
			{
				if (it2->second == temptwig[i])
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				temptwig[templength] = it2->second;
				if ((templength == twiglength) && (temptwig[templength] > temptwig[templength - 1]))
				{ // twig: consider only one case ---- left child with larger label
					continue;
				}
				BuildTwigIndex(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 1, temptwig);
			}
		}
	}
}

template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::FindQueryTwig(Path_Label &TwigLabel, Path_Num &TwigNum, Path_Index &TwigIndex1, Path_Index &TwigIndex2, CGBE *cgbe, int twiglength, int templength, VLabelType *temptwig, VLabelType *templabel, VertexID vertex, VertexID tempvertex, int point)
{
	if (templength == (twiglength + 1))
	{
		int sum1 = 0; // twig
		int sum2 = 0; // only left child
		int sum3 = 0; // only right child
		// int sum4 = 0;	//no children
		bool flag = true;

		for (int i = 1; i < (twiglength + 1); i++)
		{
			//int temp = 1;
			if (TwigNum[twiglength][temptwig[i]] == 0 && (i != twiglength))
				cout << "label 0 happen when generating the index!" << endl;

			//for (int j = 1; j < i; j++)
			//{
			//	temp = temp * TwigNum[twiglength][-1];
			//}

			// if(i<twiglength-1){
			//	sum4 += temp*TwigNum[twiglength][temptwig[i]];
			// }

			if (i < (twiglength - 1))
			{
				sum2 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
				sum3 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
			}

			if (i == (twiglength - 1))
			{
				if (temptwig[i + 1] == 0)
					flag = false;
				sum2 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
				sum3 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i + 1]]; // the last label may be 0 due to repeated label
			}

			sum1 += TwigNum[twiglength][i*(-1)] * TwigNum[twiglength][temptwig[i]];
			// cout<<"Twig["<<i<<"]:"<<temptwig[i]<<endl;
		}


		//cout << sum1 << "	" << sum2 << "	" << sum3 << endl;

		if (point == 0)
		{
			if (flag)
			{
				cgbe->setvalue(TwigIndex1[twiglength][vertex][sum2], cgbe->encoding); // sum2 is covered by sum1 in this case
				cgbe->setvalue(TwigIndex1[twiglength][vertex][sum3], cgbe->encoding);
			}
			cgbe->setvalue(TwigIndex1[twiglength][vertex][sum1], cgbe->encoding);
		}
		if (point == 1)
		{
			if (flag)
			{
				cgbe->setvalue(TwigIndex2[twiglength][vertex][sum2], cgbe->encoding);
				cgbe->setvalue(TwigIndex2[twiglength][vertex][sum3], cgbe->encoding);
			}
			cgbe->setvalue(TwigIndex2[twiglength][vertex][sum1], cgbe->encoding);
		}
		return;
	}

	if (templength == 0)
	{
		VLabelType *twig1 = new VLabelType[twiglength + 1]; // ?????may need two twigs?
		VLabelType *twig2 = new VLabelType[twiglength + 1];
		VLabelType *tempLabel = new VLabelType[this->labelsize];
		for (int i = 0; i < (twiglength + 1); i++)
		{
			twig1[i] = -1;
			twig2[i] = -1;
		}

		for (auto it = this->query->getVLabel().begin(); it != this->query->getVLabel().end(); it++)
		{
			twig1[0] = this->query->getVLabel(it->first);
			FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, 1, twig1, tempLabel, it->first, it->first, 0);
			twig2[0] = this->query->getVLabel(it->first);
			FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, 1, twig2, tempLabel, it->first, it->first, 1);
		}
		delete[] twig1;
		delete[] twig2;
		delete[] tempLabel;
	}
	else
	{
		bool flag, flag1, flag2;
		// for child path
		if (point == 0)
		{
			if (templength < (twiglength - 1))
			{
				for (auto it = this->query->getOutEdge()[tempvertex].begin(); it != this->query->getOutEdge()[tempvertex].end(); it++)
				{
					flag = true;
					for (int i = 0; i < templength; i++)
					{
						if (this->query->getVLabel(it->first) == temptwig[i])
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						temptwig[templength] = this->query->getVLabel(it->first);
						FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 1, temptwig, templabel, vertex, it->first, 0);
					}
				}
			}
			else
			{ // the twig node
				if (this->query->getOutEdge()[tempvertex].size() == 0)
					return;
				if (this->query->getOutEdge()[tempvertex].size() == 1)
				{ // only one child
					for (auto it = this->query->getOutEdge()[tempvertex].begin(); it != this->query->getOutEdge()[tempvertex].end(); it++)
					{
						flag = true;
						for (int i = 0; i < templength; i++)
						{
							if (this->query->getVLabel(it->first) == temptwig[i])
							{
								flag = false;
								break;
							}
						}
						if (flag)
						{
							temptwig[templength] = this->query->getVLabel(it->first);
							temptwig[templength + 1] = 0;
							FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 2, temptwig, templabel, vertex, 0, 0);
						}
					}
				}
				else
				{ // more than two child
					int labellength=0;
					for(int i = 0; i<templength;i++){
						templabel[labellength] = temptwig[i];
						labellength++;
					}
					for (auto it1 = this->query->getOutEdge()[tempvertex].begin(); it1 != this->query->getOutEdge()[tempvertex].end(); it1++)
					{
						flag1 = true;
						for(int i = 0;i<labellength;i++){
							if(this->query->getVLabel(it1->first)==templabel[i]){
								flag1 = false;
								break;
							}
						}

						if(flag1){
							templabel[labellength] = this->query->getVLabel(it1->first);
							labellength++;
						}
					}


					if(templength == labellength)
						return;

					if(templength == (labellength - 1)){
						temptwig[templength] = templabel[templength];
						temptwig[templength + 1] = 0;
						FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 2, temptwig, templabel, vertex, 0, 0);
										
					}else{
						for(int i = templength;i<labellength;i++)
						{
							for (int j = templength; j<labellength;j++)
							{
								if (TwigNum[twiglength][templabel[i]] <= TwigNum[twiglength][templabel[j]])
									continue;

								temptwig[templength] = templabel[i];
								temptwig[templength + 1] = templabel[j];
								FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 2, temptwig, templabel, vertex, 0, 0);				
							}
						}
					}
				}
			}
		}

		// for parent path
		if (point == 1)
		{
			if (templength < (twiglength - 1))
			{
				for (auto it = this->query->getInVertex()[tempvertex].begin(); it != this->query->getInVertex()[tempvertex].end(); it++)
				{
					flag = true;
					for (int i = 0; i < templength; i++)
					{
						if (this->query->getVLabel(it->first) == temptwig[i])
						{
							flag = false;
							break;
						}
					}
					if (flag)
					{
						temptwig[templength] = this->query->getVLabel(it->first);
						FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 1, temptwig, templabel, vertex, it->first, 1);
					}
				}
			}
			else
			{ // the twig node
				if (this->query->getInVertex()[tempvertex].size() == 0)
					return;
				if (this->query->getInVertex()[tempvertex].size() == 1)
				{ // only one parent
					for (auto it = this->query->getInVertex()[tempvertex].begin(); it != this->query->getInVertex()[tempvertex].end(); it++)
					{
						flag = true;
						for (int i = 0; i < templength; i++)
						{
							if (this->query->getVLabel(it->first) == temptwig[i])
							{
								flag = false;
								break;
							}
						}
						if (flag)
						{
							temptwig[templength] = this->query->getVLabel(it->first);
							temptwig[templength + 1] = 0;
							FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 2, temptwig, templabel, vertex, 0, 1);
						}
					}
				}
				else
				{
					int labellength=0;
					for(int i = 0; i<templength;i++){
						templabel[labellength] = temptwig[i];
						labellength++;
					}
					for (auto it1 = this->query->getInVertex()[tempvertex].begin(); it1 != this->query->getInVertex()[tempvertex].end(); it1++)
					{
						flag1 = true;
						for(int i = 0;i<labellength;i++){
							if(this->query->getVLabel(it1->first)==templabel[i]){
								flag1 = false;
								break;
							}
						}

						if(flag1){
							templabel[labellength] = this->query->getVLabel(it1->first);
							labellength++;
						}
					}


					if(templength == labellength)
						return;

					if(templength == (labellength - 1)){
						temptwig[templength] = templabel[templength];
						temptwig[templength + 1] = 0;
						FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 2, temptwig, templabel, vertex, 0, 1);
										
					}else{
						for(int i = templength;i<labellength;i++)
						{
							for (int j = templength; j<labellength;j++)
							{
								if (TwigNum[twiglength][templabel[i]] <= TwigNum[twiglength][templabel[j]])
									continue;

								temptwig[templength] = templabel[i];
								temptwig[templength + 1] = templabel[j];
								FindQueryTwig(TwigLabel, TwigNum, TwigIndex1, TwigIndex2, cgbe, twiglength, templength + 2, temptwig, templabel, vertex, 0, 1);				
							}
						}
					}					
				}
			}
		}
	}
}


template <class VLabelType, class ELabelType>
void PPMATCH<VLabelType, ELabelType>::MultiServers(int times, int servernum, unordered_map<int, double>& BallTime, unordered_map<int, double> & TwigTime, unordered_map<int, double> & TwigBallReadTime, ofstream & OutFile)
{
	

	int Num_Computed_Balls = 0, Num_positive=0, Num_negative=0;
	for(int i =0; i<this->graph_size;i++){
		if(this->BallCount[i]!=0)
			Num_Computed_Balls++;
		if(this->BallCount[i]==1)
			Num_positive++;
		if(this->BallCount[i]==-1)
			Num_negative++;
	}

	cout<< "time: " <<times <<" servernum: "<<servernum<<" Num: "<<Num_Computed_Balls<<" Positive: "<<Num_positive<<endl;

	int* ComputedBalls = new int[Num_Computed_Balls];
	int* Positive = new int[Num_positive];
	int* Negative = new int[Num_negative];
	int* Mark = new int[Num_Computed_Balls];
	int* Mark1 = new int[Num_positive];
	int* Mark0 = new int[Num_negative];

	int count = 0,count1 = 0, count2 = 0, pointer, position;
	for(int i =0; i<this->graph_size;i++){
		if(this->BallCount[i]==1){
			Positive[count1] = i;
			count1++;
		}		
		if(this->BallCount[i]==-1){
			Negative[count2] = i;
			count2++;
		}	
		if(this->BallCount[i]!=0){
			ComputedBalls[count] = i;
			count++;
		}			
	}


	double baselineTime=0, orderTime=0, twigTime = 0, tempTime = 0, tempTwig = 0, ReadTime = 0, tempRead = 0, base_Avg = 0, order_Avg = 0, twig_Avg=0, twig_read_Avg=0;


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

		
		//order
		if((2*Num_positive)<Num_Computed_Balls){
			cout<<"Case 1!"<<endl;			
			for(int j=0;j<Num_positive;j++){
				Mark1[j] = 0;
			}
			for(int j=0;j<Num_negative;j++){
				Mark0[j] = 0;
			}

			for(int j=0;j<servernum;j++){
				earlyposition[j] = 0;
			}

			int k = 0;
			//insert positive values
			for(int j = 0;j<servernum;j++){				
				while(k<Num_positive){
					if((earlyposition[j]==(Num_positive/servernum+1))&&(j<mod2)){						
						break;
					}
					if((earlyposition[j]==(Num_positive/servernum))&&(j>=mod2)){						
						break;
					}

					position = rand()%Num_positive;
					while(Mark1[position]!=0){
						position = (position+1)%Num_positive;
					}
					Mark1[position] = 1;
					order[i][j][earlyposition[j]]= Positive[position];
					earlyposition[j]++;
					k++;
				}
			}

			
			//insert negative values
			for(int j = 0;j<servernum;j++){		
				count = 0;
				while(1){						
					if(count == earlyposition[j]){						
						break;
					}
					position = rand()%Num_negative;
					while(Mark0[position]!=0){
						position = (position+1)%Num_negative;
					}
					Mark0[position] = 1;
					order[i][j][earlyposition[j]+count] = Negative[position];
					count++;
				}				
			}

			//fill the rest balls (to do)
			


		}else{
			cout<<"Case 0!"<<endl;
			count= 0;
			for(int j=0;j<Num_Computed_Balls;j++){
				Mark[j] = 0;
			}
			

			for(int j=0; j<servernum;j++){
				count = 0;
				while(1){		
					if((count==size)&&(j<mod1)){
						earlyposition[j] = size;
						break;
					}
					if((count==size-1)&&(j>=mod1)){
						earlyposition[j] = size-1;
						break;	
					}
					position=rand()%Num_Computed_Balls;
					while(Mark[position]!=0){
						position = (position+1)%Num_Computed_Balls;
					}
					Mark[position] = 1;
					order[i][j][count] = ComputedBalls[position];
					count++;
				}
			}			
		}

		//permutation (to do)




		//compute the baseline time
		baselineTime = 0;
		twigTime =0;
		for(int j = 0; j< servernum;j++){
			tempTime = 0;
			tempTwig = 0;
			tempRead = 0;
			for(int k = 0; k< size; k++){
				if(ballassign[i][j][k]!=-1){
					tempTime += BallTime[ballassign[i][j][k]];
					tempTwig += TwigTime[ballassign[i][j][k]];
					tempRead += TwigBallReadTime[ballassign[i][j][k]];
				}
			}
			if(baselineTime < tempTime)
				baselineTime = tempTime;
			if(twigTime<tempTwig)
				twigTime = tempTwig;
			if(ReadTime<tempRead)
				ReadTime = tempRead;
		}

		//compute the order time
		orderTime = 0;
		
		for(int j = 0; j< servernum;j++){
			tempTime = 0;
			
			for(int k = 0; k< earlyposition[j]; k++){
				if(ballassign[i][j][k]!=-1)
					tempTime += BallTime[order[i][j][k]];
					
			}
			if(orderTime < tempTime)
				orderTime = tempTime;
			
		}

		//cout<<"The baselineTime: " << baselineTime <<"s !	The SOP Time: " << orderTime <<" s!"<<endl;
		//cout<<"The TwigTime: "<< twigTime<<endl;
		//cout<<"The TwigBallReadTime: "<< ReadTime<<endl;

		base_Avg += baselineTime;
		order_Avg += orderTime;
		twig_Avg += twigTime;
		twig_read_Avg += ReadTime;
	}

	cout << "******************************************************************" << endl;
	cout << "m = " << servernum << endl;
	cout << "The average baselineTime: " << base_Avg/times << "s!"<<endl;
	cout << "The average SOPTime: " << order_Avg/times << " s!"<<endl;
	cout << "The average TwigTime: " << (twig_Avg/times)*1000 << " ms!"<<endl;
	cout << "The average TwigReadTime: " << (twig_read_Avg/times)*1000 << " ms!"<<endl;


	OutFile << servernum << endl;
	OutFile << base_Avg/times*1000 << endl;
	OutFile << order_Avg/times*1000 << endl;
	OutFile << (base_Avg/times)/(order_Avg/times) << endl;
	OutFile << (twig_Avg/times)*1000 << endl;
	OutFile << (twig_read_Avg/times)*1000 << endl;


	for(int i=0;i< times;i++){		
		for(int j=0; j< servernum;j++){
			delete[] ballassign[i][j];
			delete[] order[i][j];
		}
		delete[] ballassign[i];
		delete[] order[i];
	}


	delete[] Mark;
	delete[] Mark0;
	delete[] Mark1;
	delete[] Positive;
	delete[] Negative;
	delete[] ComputedBalls;
	delete[] earlyposition;
	delete[] ballassign;
	delete[] order;
}

template <class VLabelType, class ELabelType>
int PPMATCH<VLabelType, ELabelType>::MaxDegree(int **M_B, int size)
{
	int *degree = new int[size];
	int maxdegree = 0;
	for(int i = 0;i < size; i++){
		degree[i]=0;
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if(M_B[i][j]==1){
				degree[i]++;
				degree[j]++;
			}
		}
	}
	for(int i = 0; i < size; i++){
		if(degree[i]>maxdegree)
			maxdegree = degree[i];
	}

	delete[] degree;
	return maxdegree;
}


#endif
