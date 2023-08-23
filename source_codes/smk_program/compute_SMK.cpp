#include<cstdlib>
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<ctime>
#include<algorithm>
#include<map>
#include<utility>
#include<queue>
#include<cmath>
#include<pthread.h>
#include<unistd.h>
#include<queue>
#include "Node.h"
#include "functions.h"
using namespace std;

int main(int argc,char* argv[])
{
	clock_t time1 = clock();
	//clock_t time2 = 0; // count the buliding time of Compress ETree

	// Read argument string
	string IIN_file = argv[1];
	string Influence_file = argv[2];
	//string META_file = argv[3];
	int h = str2int(argv[3]);
	//int prune_IIN = str2int(argv[4]);
	string SMK_directory = argv[4];

	// read meta-structure file
	/*cout << "Read Meta-Structure File ... " << endl;
	MetaGraph* metaGraph = new MetaGraph(); // the meta-structure
	read_META(META_file,metaGraph); // in function.cpp
	vector<int> exist_type;*/

	// show meta structure
	/*for(int i=1;i<=metaGraph->ds;i++)
	{
		cout << "Layer " << i << ": " << endl;
		for(size_t j=0;j<metaGraph->Layers[i].size();j++)
		{
			cout << "   ID: " << metaGraph->Layers[i][j]->mid << " Type: " << metaGraph->Layers[i][j]->type;
			cout << " Layer: " << metaGraph->Layers[i][j]->layer << endl;
			bool exist_type_flag = false;
			for(size_t k=0;k<exist_type.size();k++)
			{
				if(metaGraph->Layers[i][j]->type == exist_type[k])
				{
					exist_type_flag = true;
					break;
				}
			}
			if(!exist_type_flag)
				exist_type.push_back(metaGraph->Layers[i][j]->type);
		}
	}
	cout << "Total Type: ";
	for(size_t i=0;i<exist_type.size();i++)
		cout << exist_type[i] << " ";
	cout << endl;*/

	// Read IIN File
	// num_obj, num_link, num_type, num_rel
	int stat[4] = {-1,-1,-1,-1};
	map<int,int> num_of_each_type;
	map<int,int> num_of_each_rel;
	map<int,Node*> all_Node; // all nodes in IIN (ID->pointer)
	cout << "Reading IIN File ... " << endl;
	//if(prune_IIN == 1)
	//	read_IIN_prune(IIN_file,stat,num_of_each_type,num_of_each_rel,all_Node,exist_type); // in function.cpp
	//else
	read_IIN(IIN_file,stat,num_of_each_type,num_of_each_rel,all_Node); // in function.cpp
	cout << "Done! Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;

	/*map<int,Node*>::iterator iter;
	for(iter = all_Node.begin();iter != all_Node.end();iter++)
		iter->second->printInfo();*/

	map<int,Node*>::iterator iter;

	// Filter out the individual.
	map<int,Node*> all_user; // Collect of persons
	for(iter = all_Node.begin();iter != all_Node.end();iter++)
	{
		Node* n = iter->second;
		// if type is 1 and >= one friend, collect it.
		if(n->type == 1)
		{
			map<Node*,int>::iterator iter_uu;
			int countfriend = 0;
			for(iter_uu = n->fanin.begin();iter_uu != n->fanin.end();iter_uu++)
			{
				if(iter_uu->first->type == 1)
				countfriend ++;
			}
			for(iter_uu = n->fanout.begin();iter_uu != n->fanout.end();iter_uu++)
			{
				if(iter_uu->first->type == 1)
				countfriend ++;
			}
			if(countfriend > 0)
				all_user[n->num] = n;
		}
	}
	cout << "Total Nodes: " << all_Node.size() << " ,Total Users: " << all_user.size() << endl;

	// Read influence file
	cout << "Reading Influence File ... " << endl;
	map<pair<int,int>,double> prob;
	read_Inf(Influence_file,prob,all_user); // in function.cpp
	cout << "Done! Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;


	cout << "Computing SMK Matrix" << endl;
	int count = 1;
	for(iter = all_user.begin();iter != all_user.end();iter++)
	{
		cout << count << "/" << all_user.size() << endl;
		count ++;
		int source_node = iter->first;
		map<int,Node*> VS;
		//check_in_short_distance(VS,all_Node,0,metaGraph->Layers.size(),iter->second); // in function.cpp
		//check_in_short_distance(VS,all_Node,0,h,iter->second); // in function.cpp
		map<int,double> omega;
		//compute_omega(omega,prob,VS,source_node,metaGraph->ds,all_user);
		compute_omega(omega,prob,VS,source_node,h,all_user);

		//clear flag (check_in_short_distance)
		map<int,Node*>::iterator it2;
		//for(it2 = all_Node.begin();it2 != all_Node.end();it2++)
		//	it2->second->flag = false;

		// the user larger than distance h, set omega to 0.
		queue<Node*> Q;
		Q.push(iter->second);
		iter->second->dist = 0;
		while(Q.size() != 0)
		{
			Node* n = Q.front();
			if(n->dist < h)
			{
				map<Node*,int>::iterator it3;
				for(it3 = n->fanout.begin();it3!=n->fanout.end();it3++)
				{
					if(it3->first->type == 1 && it3->first->dist > n->dist + 1)
					{
						it3->first->dist = n->dist + 1;
						Q.push(it3->first);
					}
				}
			}
			Q.pop();
		}


		map<int,double>::iterator it;
		string SMK_file = SMK_directory+"/SMK_n"+int2str(source_node)+".txt";
		fstream ff;
		ff.open(SMK_file.c_str(),ios::out);
		for(it = omega.begin();it != omega.end();it++)
		{
			if(all_user[it->first]->dist > h)
				ff << it->first << ",0" << endl;
			else
				ff << it->first << "," << it->second << endl;
		}
		ff.close();
		for(it2 = all_Node.begin();it2 != all_Node.end();it2++)
			it2->second->dist = 2147483647;

	}
	cout << "Time Usage: "<< double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds.(Contain the time of building compressed-ETree)" << endl;
}

