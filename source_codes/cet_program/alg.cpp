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
#include<iomanip>
#include "Node.h"
#include "functions.h"
using namespace std;

int main(int argc,char* argv[])
{
	clock_t time1 = clock();

	// Read argument string
	char buffer[10000];
	string IIN_file = argv[1];
	string META_file_List_File = argv[2];
	string Node_List_File = "";
	string location = argv[3];
	string Node_Location_File = argv[4];
	if (argv[5])
		Node_List_File = argv[5];

	// Read IIN File
	// num_obj, num_link, num_type, num_rel
	int stat[4] = {-1,-1,-1,-1};
	map<int,int> num_of_each_type;
	map<int,int> num_of_each_rel;
	map<int,Node*> all_Node; // all nodes in IIN (ID->pointer)
	map<int,Node*>::iterator iter;
	map<int,Node*> all_user; // Collect of individual

	cout << "Reading IIN File ... " << endl;
	read_IIN(IIN_file,stat,num_of_each_type,num_of_each_rel,all_Node); // in function.cpp
	cout << "Done! Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;
	read_Node_Location(Node_Location_File, all_Node);

	// Collect the the individual.
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

	// Read Node List
	vector<Node*> Node_List; // store Node pointer
	if(argv[5])
	{
		fstream f;
		f.open(argv[5],ios::in);
		f.getline(buffer,sizeof(buffer));
		string s = "";
		if(f.eof())
		{
			for(int i=0;buffer[i]!='\0';i++)
				s += buffer[i];
			int ss = str2int(s);
			if(all_Node.find(ss) == all_Node.end())
				cout << "Node " << ss << " is not in IIN Graph." << endl;
			if(all_Node[ss]->type != 1)
				cout << "Node " << ss << " is not an individual." << endl;
			Node_List.push_back(all_Node[ss]);
		}
		while(!f.eof())
		{	
			for(int i=0;buffer[i]!='\0';i++)
				s += buffer[i];
			int ss = str2int(s);
			if(all_Node.find(ss) == all_Node.end())
				cout << "Node " << ss << " is not in IIN Graph." << endl;
			if(all_Node[ss]->type != 1)
				cout << "Node " << ss << " is not an individual." << endl;
			Node_List.push_back(all_Node[ss]);
			s = "";
			f.getline(buffer,sizeof(buffer));
		}
		f.close();
	}
	else
	{
		//for(iter = all_user.begin();iter != all_user.end();iter++)
		//	Node_List.push_back(iter->second);
		for(iter = all_Node.begin();iter != all_Node.end();iter++)
			Node_List.push_back(iter->second);
	}


	// Read meta-structure file
	//cout << "Read Meta-Structure File ... " << endl;
	
	fstream f;
	vector<string> META_List;
	f.open(META_file_List_File.c_str(),ios::in);
	f.getline(buffer,sizeof(buffer));
	while(!f.eof())
	{
		string s = "";
		for(int i=0;buffer[i]!='\0';i++)
			s += buffer[i];
		META_List.push_back(s);
		s = "";
		f.getline(buffer,sizeof(buffer));
	}
	f.close();
	
	int meta_counter = 1;
	for(size_t met = 0;met < META_List.size();met++)
	{
		clock_t time_m = clock();
		MetaGraph* metaGraph = new MetaGraph(); // create a new meta-structure
		read_META(META_List[met],metaGraph);
		for(int i=1;i<=metaGraph->ds;i++)
		{
			cout << "Layer " << i << ": " << endl;
			for(size_t j=0;j<metaGraph->Layers[i].size();j++)
			{
				cout << "   ID: " << metaGraph->Layers[i][j]->mid << " Type: " << metaGraph->Layers[i][j]->type;
				cout << " Layer: " << metaGraph->Layers[i][j]->layer << endl;
			}
		}
		string META_name = "";
		for(int fn=META_List[met].size()-1;fn>0;fn--)
		{
			if(META_List[met][fn] == '/')
				break;
			META_name = META_List[met][fn] + META_name;
		}
		
		// for each node in Node_List, find CET.
		for(size_t i=0;i<Node_List.size();i++)
		{
			cout << "Node: " << i+1 << "/" << Node_List.size() << ", Meta-structure: " << met+1 << "/" << META_List.size() << endl;
			//cout << "Node: " << Node_List[i]->num << endl;
			if(Node_List[i]->type != metaGraph->Layers[1][0]->type)
				continue;
			CETNode* CET = new CETNode(); // root of Compress-ETree
			CET->addNode(Node_List[i]);
			CET->layer = 1;
			CET->represent_X_location = CET->content[0]->X_location;
			CET->represent_Y_location = CET->content[0]->Y_location;
			CET->root_X_location = CET->content[0]->X_location;
			CET->root_Y_location = CET->content[0]->Y_location;
			CET->distance_from_root = 0;
			int total_ins = 0;
			Traversal(metaGraph,1,CET,total_ins); // in function.cpp
			//if(total_ins == 0) // if there is no instance, skip SIMER computation
			//	continue;
			//cout << "Instance: " << total_ins << endl;

			// prune non-instance CET Node
			check_Instance_CET(CET,metaGraph->ds);
			
			// output CET
			fstream fc;
			// output d-table
			fstream dt;
			//string filename= location + "/CET_m"+int2str(meta_counter)+"_n"+int2str(Node_List[i]->num)+".txt";
			string filename = location + "/CET_"+META_name+"_n"+int2str(Node_List[i]->num)+".txt";
			fc.open(filename.c_str(),ios::out);

			string filename2 = location + "/CET_" + META_name + "_n" + int2str(Node_List[i]->num) + "_d_table" + ".txt";
			dt.open(filename2.c_str(), ios::out);

			if(CET->child.size() == 0)
			{
				//cout << "Empty CET: " << Node_List[i]->num << endl;
				fc.close();
				continue;
			}

			// Doing BFS to List CET
			int CETn_count = 1;
			int layer = 1;
			int rank = 0;
			int parent_layer = 0;
			int parent_rank = 0;
			int temp = layer;

			vector<D_TABLE*> d_table_list;
			//pair<int, int> location(layer,rank);
			//queue<pair<int,CETNode*> > Q;
			//queue<pair<pair<int, int>, CETNode*> > Q;
			queue<CET_BFS> Q;
			//pair<int,CETNode*> p(0,CET);
			//pair<pair<int, int>, CETNode*> p(location, CET);
			CET_BFS p(layer, rank, parent_layer, parent_rank, CET);
			Q.push(p);
			while(Q.size() != 0)
			{
				CETNode* cn = Q.front().CET;
				for(size_t j=0;j<cn->child.size();j++)
				{
					layer = cn->child[j]->layer;
					if (layer != temp) {
						rank = 0;
						temp = layer;
					}
					//pair<int, int> p2(layer, rank);
					//pair<pair<int, int>,CETNode*> p1(p2,cn->child[j]);
					CET_BFS p1(layer, rank, Q.front().layer, Q.front().rank, cn->child[j]);
					Q.push(p1);
					rank += 1;
				}
				//fc << CETn_count << "," << Q.front().first << "," << Q.front().second->layer;
				if (Q.front().layer != metaGraph->ds) 
				{
					fc << Q.front().layer << "-" << Q.front().rank << "," << Q.front().parent_layer << "-" << Q.front().parent_rank;
					for (size_t j = 0; j < cn->content.size(); j++)
						fc << "," << cn->content[j]->num;
					fc << endl;
				}
				else
				{	
					D_TABLE* d_table = new D_TABLE(cn->distance_from_root, int2str_0(Q.front().layer) + "-" + int2str_0(Q.front().rank) + "," + int2str_0(Q.front().parent_layer) + "-" + int2str_0(Q.front().parent_rank));
					for (size_t j = 0; j < cn->content.size(); j++)
						d_table->leaf_node += "," + int2str_0(cn->content[j]->num);
					d_table_list.push_back(d_table);
					/*dt << setw(10) << cn->distance_from_root <<setw(20) << int2str(Q.front().layer) + "-" + int2str(Q.front().rank) + "," + int2str(Q.front().parent_layer) + "-" + int2str(Q.front().parent_rank);
					for (size_t j = 0; j < cn->content.size(); j++)
						dt << "," << cn->content[j]->num;
					dt << setw(20) << "test" << setw(20) << "test";
					dt << endl;*/
				}
				Q.pop();
				CETn_count++;
			}
			d_table_compare compare;
			sort(d_table_list.begin(), d_table_list.end(), compare);
			for (size_t i = 0; i<d_table_list.size(); i++) {
				dt << d_table_list[i]->distance_from_root << " " << d_table_list[i]->leaf_node << " " << d_table_list[i]->entities << " " << d_table_list[i]->traversed << endl;
			}

			fc.close();
			dt.close();
		}
		meta_counter++;
		double time_usage_m = double(clock()-time_m)/double(CLOCKS_PER_SEC);
		string filename_t = location + "/CET_"+META_name+"_time.txt";
		fstream ft;
		ft.open(filename_t.c_str(),ios::out);
		ft << "Time Usage of meta-structure (" << META_List[met] << "): "<< time_usage_m << " seconds." << endl;
		ft.close();
	}

	cout << "Time Usage: "<< double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds." << endl;
}

