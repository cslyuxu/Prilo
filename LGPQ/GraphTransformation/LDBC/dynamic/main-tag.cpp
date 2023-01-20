#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <vector>
#include<stdlib.h>
#include <unordered_map>
#include<math.h>
using namespace std;

//#define	GRAPH_LABEL_NUMBER	40
typedef unordered_map<int, int> destination;
typedef unordered_map<int, destination> Data;

long countedge(Data&);


int main(int argc, char *argv[]) {
	string scale = argv[1];
	//int interval = stoi(argv[2]);
	//int num = stoi(argv[3]);
	//int maximum = stoi(argv[4]);
	//double unit = stod(argv[5]);



    Data data;

	
	string OutFile = "LDBC_sf"+scale;
	string OutMessage = "Message_Label";
	ofstream Graph(OutFile);
	ofstream OutMessageLabel(OutMessage);

	unordered_map<long, long> graph, person, message, forum, tagclass, post, country;
	unordered_map<long, int> label, messageflag, label_index1, label_index2, label_index3, label_index4;
	unordered_map<int, int> person_index, message_index, forum_index, tagclass_index, post_index, country_index;
	//intermediate link
	unordered_map<long, long> city_country, tag_tagclass;
	
	char str[30000];
	string str1, str2;
	int counter, label_count;
	long temp, index, temp1;

	index = 0;	

	////////////// Vertex

	//Read Person		
	index++;
	person_index[0] = index;
	ifstream Person("person_0_0.csv");
	counter = -1;
	while (!Person.eof()){
		if((counter < 0) && (Person.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (Person.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp = stol(str2);
			//cout<<str2<<endl;		
			counter++;
			graph[index]=temp;
			person[temp] = index;
			index++;
		}	
	}	
	person_index[1] = index-1;
	cout<<"person: "<< person_index[0] <<"-"<<person_index[1]<<": "<<person_index[1]- person_index[0]+1<<endl;
	Person.close();
	

	//Read Forum
	forum_index[0] = index;
	ifstream Forum("forum_0_0.csv");
	counter = -1;
	while (!Forum.eof()){
		if((counter < 0) && (Forum.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (Forum.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp = stol(str2);
			//cout<<str2<<endl;		
			counter++;
			graph[index]=temp;
			forum[temp] = index;
			index++;
		}	
	}
	forum_index[1] = index-1;
	cout<<"forum: "<< forum_index[0] <<"-"<<forum_index[1]<<": "<<forum_index[1]- forum_index[0]+1<<endl;
	Forum.close();


	//Read Message
	message_index[0] = index;
	ifstream Message("comment_0_0.csv");
	counter = -1;
	while (!Message.eof()){
		if((counter < 0) && (Message.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (Message.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp = stol(str2);
			//cout<<str2<<endl;		
			counter++;
			graph[index]=temp;
			message[temp] = index;
			index++;
		}	
	}
	message_index[1] = index-1;
	cout<<"comment: "<< message_index[0] <<"-"<<message_index[1]<<": "<<message_index[1]- message_index[0]+1<<endl;
	Message.close();


	/*
	//Tagclass 
	tagclass_index[0] = index;
	ifstream Tagclass("../static/tagclass_0_0.csv");
	counter = -1;
	while (!Tagclass.eof()){
		if((counter < 0) && (Tagclass.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (Tagclass.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp = stol(str2);
			//cout<<str2<<endl;		
			counter++;
			graph[index]=temp;
			tagclass[temp] = index;
			index++;
		}	
	}
	tagclass_index[1] = index-1;
	cout<<"tagclass: "<< tagclass_index[0] <<"-"<<tagclass_index[1]<<": "<<tagclass_index[1]- tagclass_index[0]+1<<endl;
	Tagclass.close();
	*/
	
	
	//Post 
	post_index[0] = index;
	ifstream Post("post_0_0.csv");
	counter = -1;
	while (!Post.eof()){
		if((counter < 0) && (Post.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (Post.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp = stol(str2);
			//cout<<str2<<endl;		
			counter++;
			graph[index]=temp;
			post[temp] = index;
			index++;
		}	
	}
	post_index[1] = index-1;
	cout<<"post: "<< post_index[0] <<"-"<<post_index[1]<<": "<<post_index[1]- post_index[0]+1<<endl;
	Post.close();
	
	
	
	
	//country
	country_index[0] = index;
	ifstream Country("../static/place_isPartOf_place_0_0.csv");
	counter = -1;
	while (!Country.eof()){
		if((counter < 0) && (Country.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (Country.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');	
			temp1 = stol(str2);
			getline(input, str2, '\n');	
			temp = stol(str2);
			city_country[temp1] = temp;
			if(country.find(temp)==country.end()){	
				//cout<<str2<<endl;		
				counter++;
				graph[index]=temp;
				country[temp] = index;
				index++;
			}
		}	
	}
	country_index[1] = index-1;
	cout<<"country: "<< country_index[0] <<"-"<<country_index[1]<<": "<<country_index[1]- country_index[0]+1<<endl;
	Country.close();







	////////////// Edge


	//person-country
	ifstream PC("person_isLocatedIn_place_0_0.csv");
	counter = -1;
	while (!PC.eof()){
		if((counter < 0) && (PC.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PC.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');
			temp1 = stol(str2);	
			getline(input, str2, '\n');	
			temp = stol(str2);
			data[person[temp1]][country[city_country[temp]]] = 1;
			//cout<<"Person: " << person[temp1] <<"		Country: "<< city_country[temp] <<endl;
		}	
	}
	PC.close();
	cout<<"After person->country, the number of current edges : "<<countedge(data)<<endl;

	//forum-post
	ifstream FP("forum_containerOf_post_0_0.csv");
	counter = -1;
	while (!FP.eof()){
		if((counter < 0) && (FP.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (FP.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			data[forum[temp1]][post[temp]] = 1;
		}	
	}	
	FP.close();
	cout<<"After forum->post, the number of current edges : "<<countedge(data)<<endl;


	//forum-person
	ifstream FP1("forum_hasModerator_person_0_0.csv");
	counter = -1;
	while (!FP1.eof()){
		if((counter < 0) && (FP1.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (FP1.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			data[forum[temp1]][person[temp]] = 1;
		}	
	}	
	FP1.close();
	cout<<"After forum->person-HasModerator, the number of current edges : "<<countedge(data)<<endl;

	ifstream FP2("forum_hasMember_person_0_0.csv");
	counter = -1;
	while (!FP2.eof()){
		if((counter < 0) && (FP2.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (FP2.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '|');			
			temp = stol(str2);
			data[forum[temp1]][person[temp]] = 1;
		}	
	}	
	FP2.close();
	cout<<"After forum->person-HasMember, the number of current edges : "<<countedge(data)<<endl;



	//post-person
	ifstream PP("post_hasCreator_person_0_0.csv");
	counter = -1;
	while (!PP.eof()){
		if((counter < 0) && (PP.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PP.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			data[post[temp1]][person[temp]] = 1;
		}	
	}	
	PP.close();
	cout<<"After post->person, the number of current edges : "<<countedge(data)<<endl;



	//person-post
	ifstream PP1("person_likes_post_0_0.csv");
	counter = -1;
	while (!PP1.eof()){
		if((counter < 0) && (PP1.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PP1.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '|');			
			temp = stol(str2);
			data[person[temp1]][post[temp]] = 1;
		}	
	}	
	PP1.close();
	cout<<"After person->post, the number of current edges : "<<countedge(data)<<endl;




	//person-comment
	ifstream PC1("person_likes_comment_0_0.csv");
	counter = -1;
	while (!PC1.eof()){
		if((counter < 0) && (PC1.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PC1.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '|');			
			temp = stol(str2);
			data[person[temp1]][message[temp]] = 1;
		}	
	}	
	PC1.close();
	cout<<"After person->comment (message), the number of current edges : "<<countedge(data)<<endl;




	//person-person
	ifstream PP2("person_knows_person_0_0.csv");
	counter = -1;
	while (!PP2.eof()){
		if((counter < 0) && (PP2.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PP2.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '|');			
			temp = stol(str2);
			data[person[temp1]][person[temp]] = 1;
			data[person[temp]][person[temp1]] = 1;
		}	
	}	
	PP2.close();
	cout<<"After person<->person, the number of current edges : "<<countedge(data)<<endl;


	//comment-post
	ifstream CP("comment_replyOf_post_0_0.csv");
	counter = -1;
	while (!CP.eof()){
		if((counter < 0) && (CP.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (CP.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			data[message[temp1]][post[temp]] = 1;
		}	
	}	
	CP.close();
	cout<<"After comment (message)->post, the number of current edges : "<<countedge(data)<<endl;


	//comment-comment
	ifstream CC("comment_replyOf_comment_0_0.csv");
	counter = -1;
	while (!CC.eof()){
		if((counter < 0) && (CC.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (CC.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			data[message[temp1]][message[temp]] = 1;
			messageflag[message[temp]] = 1; //for identifying the message label out of comments
		}	
	}	
	CC.close();
	cout<<"After comment (message)->comment (message), the number of current edges : "<<countedge(data)<<endl;




	//comment-person
	ifstream CP1("comment_hasCreator_person_0_0.csv");
	counter = -1;
	while (!CP1.eof()){
		if((counter < 0) && (CP1.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (CP1.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			data[message[temp1]][person[temp]] = 1;
		}	
	}	
	CP1.close();
	cout<<"After comment (message)->person, the number of current edges : "<<countedge(data)<<endl;





	int messagecount = 0;
	label_count = 1;
	////////////// Label (Tag)

	for(long i = 1; i <= country_index[1]; i++){
		label[i] = 0;
	}

	//tag-tagclass
	ifstream TT("../static/tag_hasType_tagclass_0_0.csv");
	counter = -1;
	while (!TT.eof()){
		if((counter < 0) && (TT.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (TT.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			tag_tagclass[temp1] = temp;
		}	
	}	
	TT.close();	

	
    //person-tagclass	
	ifstream PT2("person_hasInterest_tag_0_0.csv");
	counter = -1;
	while (!PT2.eof()){
		if((counter < 0) && (PT2.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PT2.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			if(label_index1.find(tag_tagclass[temp])==label_index1.end()){
				label_index1[tag_tagclass[temp]] = label_count;
				label_count++;
			}
			label[person[temp1]] = label_index1[tag_tagclass[temp]];			
		}	
	}	
	PT2.close();
	cout<<"After person-tag, the number of current labels : "<<label_count-1<<endl;
	cout<<"After person-tag, the number of current edges : "<<countedge(data)<<endl;



	//forum-tagclass
	ifstream FT("forum_hasTag_tag_0_0.csv");
	counter = -1;
	while (!FT.eof()){
		if((counter < 0) && (FT.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (FT.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			if(label_index2.find(tag_tagclass[temp])==label_index2.end()){
				label_index2[tag_tagclass[temp]] = label_count;
				label_count++;
			}
			label[forum[temp1]] = label_index2[tag_tagclass[temp]];	
		}	
	}	
	FT.close();
	cout<<"After forum-tag, the number of current labels : "<<label_count-1<<endl;
	cout<<"After forum-tag, the number of current edges : "<<countedge(data)<<endl;


	//post-tagclass
	ifstream PT3("post_hasTag_tag_0_0.csv");
	counter = -1;
	while (!PT3.eof()){
		if((counter < 0) && (PT3.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (PT3.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			if(label_index3.find(tag_tagclass[temp])==label_index3.end()){
				label_index3[tag_tagclass[temp]] = label_count;
				label_count++;
			}
			label[post[temp1]] = label_index3[tag_tagclass[temp]];	
		}	
	}	
	PT3.close();
	cout<<"After post-tag, the number of current labels : "<<label_count-1<<endl;
	cout<<"After post-tag, the number of current edges : "<<countedge(data)<<endl;

	//message-tagclass	
	ifstream MT("comment_hasTag_tag_0_0.csv");
	counter = -1;
	while (!MT.eof()){
		if((counter < 0) && (MT.getline(str, 30000))) {			
			counter++;							
			continue;
		}

		if((counter >=0) && (MT.getline(str, 30000))) {	
			str1 = str;
			stringstream input(str1);
			getline(input, str2, '|');			
			temp1 = stol(str2);
			getline(input, str2, '\n');			
			temp = stol(str2);
			if(label_index4.find(tag_tagclass[temp])==label_index4.end()){
				label_index4[tag_tagclass[temp]] = label_count;
				label_count++;
			}
			label[message[temp1]] = label_index4[tag_tagclass[temp]];	
		}	
	}	
	MT.close();
	cout<<"After message-tag, the number of current labels : "<<label_count-1<<endl;
	cout<<"After message-tag, the number of current edges : "<<countedge(data)<<endl;
	




	for(long i = 1; i <= country_index[1]; i++){
		if(i<=person_index[1]){			
			continue;
		}
		if(i<=forum_index[1]){			
			continue;
		}

		//Message is a message having reply comments
		if(i<=message_index[1]){			
			if(messageflag.find(i)!=messageflag.end()){				
				messagecount++;
			}
			continue;
		}

		if(i<=post_index[1]){			
			continue;
		}

		if(i<=country_index[1]){
			label[i] = label_count;
			continue;
		}
	}

	cout<< "The number of message is : "<<messagecount <<" out of " << message_index[1]-message_index[0]+1<< " comments!"<<endl;
	cout<<"After coutry-tag, the number of current labels : "<<label_count<<endl;




	//output dataset
	Graph<<"# Trasformed LDBC-graph: person_forum_comment_message_post_country"<<endl;
	Graph<<"# LDBC_sf"<<scale<<endl;
	Graph<<"# Nodes: "<<country_index[1]<<"	 Edges: "<<countedge(data)<<endl;
	Graph<<"# FromNodeId	Label	ToNodeId";

	OutMessageLabel<<"#	NodeId	Label";

	for(long i = 1; i <=country_index[1];i++){
		if(messageflag.find(i)!=messageflag.end())
			OutMessageLabel<<"\n"<<i<<"\t"<<label[i];
		Graph << "\n#\t" << i << "\t" << label[i];
		for(auto it = data[i].begin();it != data[i].end();it++){
			Graph <<"\t"<<it->first;
		}
	}

	OutMessageLabel.close();
	Graph.close();
	return 0;
}



long countedge(Data& edge){
	long count = 0;
	for(auto it = edge.begin();it !=edge.end();it++){
		for(auto it1 = it->second.begin();it1!= it->second.end();it1++)
			count++;
	}
	return count;
}



