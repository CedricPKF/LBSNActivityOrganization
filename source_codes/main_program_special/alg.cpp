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
#include "Node.h"
#include "functions.h"
using namespace std;
//////////////////////////////
//int smth;
//////////////////////////////
//#define NUM_THREADS 8

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void count_remaining_time(double x)
{
	int X = int(x);
	cout << "Estimate Remaining Time: ";
	if(X/86400 > 0)
		cout << X/86400 << " day ";
	X = X - 86400*(X/86400);
	if(X/3600 > 0)
		cout << X/3600 << " hour ";
	X = X - 3600*(X/3600);
	if(X/60 > 0)
		cout << X/60 << " minute ";
	X = X - 60*(X/60);
	cout << X << " second " << endl;
}
struct ARG{
	int tid;
	//int tn; // target_node
	vector<int> tn;// target_node vector
	int un; // user_node
	int d; // delta
	MetaGraph* mt; // metaGraph
	vector<Node*>* user_vec; // division set of user nodes
	map<int,Node*>* an; // all_Node
	map<int,Node*>* au; // all_user
	//map<pair<int,int>,double>* pb; // prob
	string smk;
	//vector<int>* sv; // source_vec
	//vector<double>* simerv; // SIMER_vec
	vector<pair<int, double> >* Simerv;
	string cet;
	string meta_name;
	int vt_flag;
	double Spatial_Distance;
	string sim;
	string iin_name;
};

void* Mission(void *t)
{
	struct ARG  *my_data;
	my_data = (struct ARG *)t;
	int tid = my_data->tid;
	//int target_node = my_data->tn;
	vector<int> target_nodes = my_data->tn;
	int user_node = my_data->un;
	int delta = my_data->d;
	MetaGraph* metaGraph = my_data->mt;
	vector<Node*>* user_vec = my_data->user_vec;
	map<int,Node*>* all_Node = my_data->an;
	map<int,Node*>* all_user = my_data->au;
	//map<pair<int,int>,double>* prob = my_data->pb;
	string smk = my_data->smk;
	//vector<int>* source_vec = my_data->sv;
	//vector<double>* SIMER_vec = my_data->simerv;
	vector<pair<int, double> >* Simer_vec = my_data->Simerv;
	string cet = my_data->cet;
	string meta_name = my_data->meta_name;
	int vt_flag = my_data->vt_flag;
	double Spatial_Distance = my_data->Spatial_Distance;
	string sim = my_data->sim;
	string iin_name = my_data->iin_name;
	//for(iter = all_user.begin();iter != all_user.end();iter++)
	int count = 0;
	int real_count = 0;
	//clock_t t_p = clock();
	for(size_t kk=0;kk<user_vec->size();kk++)
	{
		count++;
		real_count++;
		// choose V_s in all_user
		int source_node = (*user_vec)[kk]->num;
		//int source_node = iter->first;
		if(source_node == user_node) // if V_s == the specific user, skip
			continue;
		if((*all_user)[source_node]->dist == 2147483647) // if V_s does not have path to specific user, skip
			continue;

		// Traversal Algorithm and get all instance
		/*CETNode* CET = new CETNode(); // root of Compress-ETree
		CET->addNode((*user_vec)[kk]);
		CET->layer = 1;
		int total_ins = 0;
		Traversal(metaGraph,1,CET,total_ins); // in function.cpp
		if(total_ins == 0) // if there is no instance, skip SIMER computation
			continue;
		cout << "Instance: " << total_ins << endl;*/


		// load CET
		string CET_file = "";
		string D_TABLE_file = "";
		string meta_name2 = "";
		for(int fn=meta_name.size();fn>=0;fn--)
		{
			if(meta_name[fn] =='\0')
				continue;
			if(meta_name[fn] == '/')
				break;
			meta_name2 = meta_name[fn] + meta_name2;
		}
		CET_file = cet + "/CET_"+ meta_name2 +"_n"+int2str(source_node)+".txt";
		D_TABLE_file = cet + "/CET_" + meta_name2 + "_n" + int2str(source_node) + "_d_table.txt";

		//cout << "Reading CET File: " << tid << endl;

   // cout << "Reading CET File: " << tid << endl;

		//cout << "Reading CET File: " << cet + "/CET_"+meta_name2+"_n"+int2str(source_node)+".txt" << endl;

		//reading CET file:

		fstream fcet;
		fcet.open(CET_file.c_str(), ios::in);
		char buffer[1000];
		fcet.getline(buffer, sizeof(buffer));
		//map<int,CETNode*> all_CETNode;
		map<pair<int, int>, CET_FILE_information*> cet_file_list;
		while (!fcet.eof())
		{
			string s = "";
			int count = 0;
			//int nodeid = -1;
			int node_layer = -1;
			int node_rank = -1;
			//int parent = -1;
			int parent_layer = -1;
			int parent_rank = -1;
			vector<int> contentid;
			bool flag_notinNode = false;
			for (int i = 0; buffer[i] != '\0'; i++)
			{
				if (buffer[i] == '-' && count == 0)
				{
					//nodeid = str2int(s);
					node_layer = str2int(s);
					s = "";
					count++;
				}
				else if (buffer[i] == ',' && count == 1)
				{
					//parent = str2int(s);
					node_rank = str2int(s);
					s = "";
					count++;
				}
				else if (buffer[i] == '-' && count == 2)
				{
					parent_layer = str2int(s);
					s = "";
					count++;
				}
				else if (buffer[i] == ',' && count == 3)
				{
					//parent = str2int(s);
					parent_rank = str2int(s);
					s = "";
					count++;
				}
				else if (buffer[i] == ',')
				{
					int temp = str2int(s);
					if (all_Node->find(temp) == all_Node->end())
					{
						flag_notinNode = true;
						break;
					}
					contentid.push_back(temp);
					s = "";
				}
				else
					s += buffer[i];
			}
			if (flag_notinNode)
			{
				fcet.getline(buffer, sizeof(buffer));
				continue;
			}
			int temp = str2int(s);
			if (all_Node->find(temp) == all_Node->end())
			{
				fcet.getline(buffer, sizeof(buffer));
				continue;
			}
			contentid.push_back(temp);	
			pair<int, int> node_pair, parent_pair;
			node_pair = make_pair(node_layer, node_rank);
			parent_pair = make_pair(parent_layer, parent_rank);
			CET_FILE_information* CFI = new CET_FILE_information(node_pair, parent_pair, contentid);
			cet_file_list[node_pair] = CFI;
			fcet.getline(buffer, sizeof(buffer));
		}
		fcet.close();

		//cout << "Reading D_TABLE File: " << cet + "/CET_" + meta_name2 + "_n" + int2str(source_node) + "_d_table.txt" << endl;

		fstream fdtb;
		fdtb.open(D_TABLE_file.c_str(), ios::in);
		char buffer_1[10000];
		fdtb.getline(buffer_1, sizeof(buffer_1));
		//map<pair<int, int>, D_TABLE_information*> d_table_list;
		vector<D_TABLE_information*> d_table_list;
		while (!fdtb.eof())
		{
			string s = "";
			int count = 0;
			double distance_from_root = 0;
			vector<int> contentid;
			vector<int> entity;
			vector<pair<int, int> > traversed;
			int cet_leafnode_layer = -1, cet_leafnode_rank = -1, cet_leafnode_parent_layer = -1, cet_leafnode_parent_rank = -1 ,traversed_layer = -1, traversed_rank = -1;
			//bool over_distance_limit = false;
			bool flag_notinNode = false;
			bool entity_is_NULL = false;
			bool traversed_is_NULL = false;
			//if (over_distance_limit)//d_table has be sorted for distance.
				//break;
			for (size_t i = 0; buffer_1[i] != '\0'; i++)
			{
				if (buffer_1[i] == ' ' && count == 0)
				{
					distance_from_root = str2double(s);
					s = "";
					count++;
					/*if (distance_from_root > Spatial_Distance)
					{
						over_distance_limit = true;
						break;
					}*/
				}
				else if (buffer_1[i] == '-' && count == 1)
				{
					cet_leafnode_layer = str2int(s);
					s = "";
					count++;
				}
				else if (buffer_1[i] == ',' && count == 2)
				{
					cet_leafnode_rank = str2int(s);
					s = "";
					count++;
				}
				else if (buffer_1[i] == '-' && count == 3)
				{
					cet_leafnode_parent_layer = str2int(s);
					s = "";
					count++;
				}
				else if (buffer_1[i] == ',' && count == 4)
				{
					cet_leafnode_parent_rank = str2int(s);
					s = "";
					count++;
				}
				else if(buffer_1[i] == ',' && count == 5)
				{
					int temp = str2int(s);
					if (all_Node->find(temp) == all_Node->end())
					{
						flag_notinNode = true;
						break;
					}
					contentid.push_back(temp);
					s = "";
				}
				else if (buffer_1[i] == ' ' && count == 5)
				{
					int temp = str2int(s);
					if (all_Node->find(temp) == all_Node->end())
					{
						flag_notinNode = true;
						break;
					}
					contentid.push_back(temp);
					s = "";
					count++;
				}
				else if (buffer_1[i] == ',' && count == 6)
				{
					int temp = str2int(s);
					entity.push_back(temp);
					s = "";
				}
				else if (buffer_1[i] == ' ' && count == 6)
				{
					if (s == "NULL") 
					{
						count++;
						s = "";
						entity_is_NULL = true;
					}
					else if (s == "empty") {
						count++;
						s = "";
					}
					else
					{
						int temp = str2int(s);
						entity.push_back(temp);
						count++;
						s = "";
					}
				}
				else if (buffer_1[i] == '-')
				{
					traversed_layer = str2int(s);
					s = "";
				}
				else if (buffer_1[i] == ',')
				{
					traversed_rank = str2int(s);
					traversed.push_back(make_pair(traversed_layer, traversed_rank));
					s = "";
				}
				else
					s += buffer_1[i];
			}
			if (flag_notinNode)
			{
				fdtb.getline(buffer_1, sizeof(buffer_1));
				continue;
			}
			if (s == "NULL") 
			{
				traversed_is_NULL = true;
			}
			else if (s == "empty") 
			{

			}
			else 
			{
				traversed_rank = str2int(s);
				traversed.push_back(make_pair(traversed_layer, traversed_rank));
			}

			//pair<int, int> leaf_node_pair(cet_leafnode_layer, cet_leafnode_rank);
			D_TABLE_information* D_TABLEn = new D_TABLE_information(distance_from_root,cet_leafnode_layer, cet_leafnode_rank, cet_leafnode_parent_layer, cet_leafnode_parent_rank, contentid, entity, traversed,entity_is_NULL,traversed_is_NULL);
			//d_table_list[leaf_node_pair] = D_TABLEn;
			d_table_list.push_back(D_TABLEn);
			fdtb.getline(buffer_1, sizeof(buffer_1));
		}
		fdtb.close();

		//update D_TABLE
		map<pair<int, int>, CET_FILE_information*>::iterator iter_cet_file_list;
		map<pair<int, int>, bool> untraversed;
		for (iter_cet_file_list = cet_file_list.begin(); iter_cet_file_list != cet_file_list.end(); iter_cet_file_list++) {
			untraversed[iter_cet_file_list->second->Node_pair] = true;
		}
		//check whether the parent root had be traversed
		for (size_t i = 0; i<d_table_list.size();i++) {
			if (d_table_list[i]->distance_from_root > Spatial_Distance) {
				break;
			}
			for (size_t j = 0; j < d_table_list[i]->traversed.size(); j++) {
				untraversed[d_table_list[i]->traversed[j]] = false;
			}
		}
		for (size_t i = 0; i<d_table_list.size(); i++) {
			if (d_table_list[i]->distance_from_root > Spatial_Distance) {
				break;
			}
			//run this when parent node is untraversed.
			if(untraversed[d_table_list[i]->parent_pair]){
				pair<int, int> temp = d_table_list[i]->parent_pair;
				while (untraversed[temp] && temp.first != 0) 
				{
					d_table_list[i]->traversed.push_back(temp);
					for (size_t j = 0; j < cet_file_list[temp]->Content_id.size(); j++) 
					{
						d_table_list[i]->entity.push_back(cet_file_list[temp]->Content_id[j]);
					}
					untraversed[temp] = false;
					temp = cet_file_list[temp]->Parent_pair;
				}
			}
			d_table_list[i]->entity_null = false;
			d_table_list[i]->traversed_null = false;
		}
		//cout << "Completed update D_TABLE: "<< cet + "/CET_" + meta_name2 + "_n" + int2str(source_node) + "_d_table.txt" << endl;
		D_TABLE_compare compare;
		sort(d_table_list.begin(), d_table_list.end(), compare);
		//update d_table_file
		fstream Cdtbb;
		Cdtbb.open(D_TABLE_file.c_str(), ios::out);
		for (size_t i = 0; i < d_table_list.size(); i++) {
			Cdtbb << d_table_list[i]->distance_from_root << " " << d_table_list[i]->node_pair.first << "-" << d_table_list[i]->node_pair.second << ","
				<< d_table_list[i]->parent_pair.first << "-" << d_table_list[i]->parent_pair.second;
			for (size_t j = 0; j < d_table_list[i]->content_id.size(); j++) {
				Cdtbb << "," << d_table_list[i]->content_id[j];
			}
			Cdtbb << " ";
			if (d_table_list[i]->entity_null) {
				Cdtbb << "NULL";
			}
			else if (d_table_list[i]->entity.size() > 0) {
				for (size_t k = 0; k < d_table_list[i]->entity.size(); k++) {
					if (k > 0)
						Cdtbb << ",";
					Cdtbb << d_table_list[i]->entity[k];
				}
			}
			else {
				Cdtbb << "empty";
			}
			Cdtbb << " ";
			if (d_table_list[i]->traversed_null) {
				Cdtbb << "NULL";
			}
			else if (d_table_list[i]->traversed.size() > 0) {
				for (size_t k = 0; k < d_table_list[i]->traversed.size(); k++) {
					if (k > 0)
						Cdtbb << ",";
					Cdtbb << d_table_list[i]->traversed[k].first << "-" << d_table_list[i]->traversed[k].second;
				}
			}
			else {
				Cdtbb << "empty";
			}
			Cdtbb << endl;
		}

		/* //old_code
		//cout << "Reading CET File: " << cet + "/CET_" + meta_name2 + "_n" + int2str(source_node) + ".txt" << endl;
		fstream fcet;
		fcet.open(CET_file.c_str(),ios::in);
		char buffer[1000];
		fcet.getline(buffer,sizeof(buffer));
		CETNode* CET = new CETNode();
		//map<int,CETNode*> all_CETNode;
		map<pair<int,int>, CETNode*> all_CETNode;
		while(!fcet.eof())
		{
			string s = "";
			int count = 0;
			//int nodeid = -1;
			int node_layer = -1;
			int node_rank = -1;
			//int parent = -1;
			int parent_layer = -1;
			int parent_rank = -1;
			vector<int> contentid;
			bool flag_notinNode = false;
			for(int i=0;buffer[i]!='\0';i++)
			{
				if(buffer[i] == '-' && count == 0)
				{
					//nodeid = str2int(s);
					node_layer = str2int(s);
					s = "";
					count++;
				}
				else if(buffer[i] == ',' && count == 1)
				{
					//parent = str2int(s);
					node_rank = str2int(s);
					s = "";
					count++;
				}
				else if (buffer[i] == '-' && count == 2)
				{
					parent_layer = str2int(s);
					s = "";
					count++;
				}
				else if (buffer[i] == ',' && count == 3)
				{
					//parent = str2int(s);
					parent_rank = str2int(s);
					s = "";
					count++;
				}
				else if(buffer[i] == ',')
				{
					int temp = str2int(s);
					if(all_Node->find(temp) == all_Node->end())
					{
						flag_notinNode = true;
						break;
					}
					contentid.push_back(temp);
					s = "";
				}
				else
					s += buffer[i];
			}
			if(flag_notinNode)
			{
				fcet.getline(buffer,sizeof(buffer));
				continue;
			}
			int temp = str2int(s);
			if(all_Node->find(temp) == all_Node->end())
			{
				fcet.getline(buffer,sizeof(buffer));
				continue;
			}
			contentid.push_back(temp);


			CETNode* CETn = new CETNode();
			pair<int, int> node_pair, parent_pair;
			node_pair = make_pair(node_layer, node_rank);
			parent_pair = make_pair(parent_layer, parent_rank);
			all_CETNode[node_pair] = CETn;
			if(node_layer == 1)//only a root
			{
				CET = CETn;
				CETn->layer = 1;
			}
			if(parent_layer != 0)
			{
				//all_CETNode[parent]->child.push_back(all_CETNode[nodeid]);
				all_CETNode[parent_pair]->child.push_back(all_CETNode[node_pair]);
				//all_CETNode[nodeid]->layer = all_CETNode[parent]->layer + 1;
				all_CETNode[node_pair]->layer = all_CETNode[parent_pair]->layer + 1;
			}
			for(size_t item=0;item<contentid.size();item++)
				all_CETNode[node_pair]->addNode((*all_Node)[contentid[item]]);


			fcet.getline(buffer,sizeof(buffer));
		}
		//cout << "Finish Reading CET" << endl;

  // cout << "Finish Reading single CET" << tid << endl;

		//map<int,CETNode*>::iterator itc;
		map<pair<int, int>, CETNode*>::iterator itc;
		for(itc = all_CETNode.begin();itc != all_CETNode.end();itc++)
			all_CETNode.erase(itc->first);

		fcet.close();
		*/

		// Collect the Node with distance <= ds from V_s
		//map<int,Node*> VS;
		//check_in_short_distance(VS,*all_Node,0,metaGraph->Layers.size(),(*all_Node)[source_node]); // in function.cpp
		// Computing Omega
		map<int,double> omega;
		//cout << "VS = " << VS.size() << " all_user = " << all_user->size() << endl;
		//compute_omega(omega,*prob,VS,source_node,metaGraph->Layers.size(),*all_user); // in function.cpp

		//read smk file with source_node user
		//cout << "Read SMK File" << endl;

    //cout << "Read SMK File" << tid << endl;

		string smk_file = smk + "/SMK_n"+int2str(source_node)+".txt";
		read_smk(smk_file,omega);
		map<int,double>::iterator it;

		//Read SIM File
		map<int, map<int, double> > sim_file;
		for (size_t i = 0; i < target_nodes.size(); i++) 
		{
			string sim_location = sim + "/SIM_" + iin_name + "_n" + int2str(target_nodes[i]) + ".txt";
			map<int, double> sim_content;
			read_smk(sim_location, sim_content);
			sim_file[target_nodes[i]] = sim_content;
		}

		// Computing SIMER
		double d = 0.0;
		//for(size_t j=0;j<target_nodes.size();j++)
			//d += SIMER(source_node,target_nodes[j],user_node,*all_user,*all_Node,metaGraph,omega,CET,delta,vt_flag); // in function.cpp
		for (size_t j = 0; j < target_nodes.size(); j++)
		{
			map<int, double> sim_list = sim_file[target_nodes[j]];
			d += SIMER(source_node, target_nodes[j], user_node, *all_user, *all_Node, metaGraph, omega, d_table_list, cet_file_list, delta, vt_flag, Spatial_Distance, sim_list);
		}
		//double d = SIMER(source_node,target_node,user_node,*all_user,*all_Node,metaGraph,omega,CET,delta); // in function.cpp
		//cout << "Finish SIMER" << endl;

		pthread_mutex_lock(&mutex);
		//source_vec->push_back(source_node);
		//SIMER_vec->push_back(d);
		Simer_vec->push_back(make_pair(source_node, d));
		pthread_mutex_unlock(&mutex);


		//cout << "SIMER of " << source_node << " = " << d << " " << endl;
		//cout << "SIMER size() = " << SIMER_vec->size() << endl;
		/*if(count > 100)
		{
			count = 0;
			double aaa = double(clock()-t_p)/double(CLOCKS_PER_SEC);
			double bbb = double(user_vec->size()-real_count)/double(real_count);
			cout << "TID: " << tid << " ";
			count_remaining_time(aaa*bbb);
		}*/

	}
	cout << "TID: " << tid << " Finish." << endl;
	return NULL;
	//pthread_exit(NULL);
}

int main(int argc,char* argv[])
{
	clock_t time1 = clock();

	// Read argument string
	string IIN_file = argv[1];
	string Influence_file = argv[2];
	string META_file = argv[3];
	int user_node = str2int(argv[4]);
	//int target_node = str2int(argv[5]);
	string target_node_file = argv[5];
	//int topK = str2int(argv[6]);
	string source_node_file = argv[6];
	double delta = str2double(argv[7]);
	//string topK_file = argv[8];
	string detail_file = argv[8];
	int NUM_THREADS = str2int(argv[9]);
	int prune_IIN = str2int(argv[10]);
	string SMK_directory = argv[11];
	string CET_location = argv[12]; // new add
	int Distance = str2int(argv[13]); // new add
	int vt_flag = str2int(argv[14]); // new add
	double Spatial_Distance = str2double(argv[15]);
	string SIM_directory = argv[16];

//reading target nodes
	fstream tgn;
	vector<int> target_nodes;
	cout << "target node file " << target_node_file << endl;
	tgn.open(target_node_file.c_str());
	char buffer[50];
	tgn.getline(buffer,sizeof(buffer));
	if(tgn.eof())
	{
		string s = "";
		for(int i=0;buffer[i] != '\0' && int(buffer[i])>31;i++)
			s += buffer[i];
		target_nodes.push_back(str2int(s));
	}
	while(!tgn.eof())
	{
		if(buffer[0] == '\0')
			tgn.getline(buffer,sizeof(buffer));
		string s = "";
		for(int i=0;buffer[i] != '\0' && int(buffer[i])>31;i++)
			s += buffer[i];
		target_nodes.push_back(str2int(s));
		tgn.getline(buffer,sizeof(buffer));
	}
	tgn.close();

//reading source node file
	fstream sn;
	vector<int> source_nodes;
	cout << "source node file " << source_node_file << endl;
	sn.open(source_node_file.c_str(), ios::in);
	sn.getline(buffer, sizeof(buffer));
	if (sn.eof())
	{
		string s = "";
		for (int i = 0; buffer[i] != '\0' && int(buffer[i])>31; i++)
			s += buffer[i];
		source_nodes.push_back(str2int(s));
	}
	while (!sn.eof())
	{
		if (buffer[0] == '\0')
			sn.getline(buffer, sizeof(buffer));
		string s = "";
		for (int i = 0; buffer[i] != '\0' && int(buffer[i])>31; i++)
			s += buffer[i];
		source_nodes.push_back(str2int(s));
		sn.getline(buffer, sizeof(buffer));
	}
	sn.close();


	cout << "u = " << user_node << ", vt = ";
	for(size_t i=0;i<target_nodes.size();i++)
		cout << target_nodes[i] << ",";
	cout << endl << " ,delta = " << delta << endl;
	cout << "NUM_THREADS = " << NUM_THREADS << endl;

///////////////////////////////////////////////////
  //cout << "Type in the 1 key" << endl;
	//cin >> smth;
//////////////////////////////////////////////////

	// read meta-structure file
	cout << "Read Meta-Structure File ... " << endl;
	MetaGraph* metaGraph = new MetaGraph(); // the meta-structure
	read_META(META_file,metaGraph); // in function.cpp
	vector<int> exist_type;
	for(int i=1;i<=metaGraph->ds;i++)
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
	cout << endl;

///////////////////////////////////////////////////
	//cout << "Type in the 1 key" << endl;
	//cin >> smth;
//////////////////////////////////////////////////

	// Read IIN File
	// num_obj, num_link, num_type, num_rel
	int stat[4] = {-1,-1,-1,-1};
	map<int,int> num_of_each_type;
	map<int,int> num_of_each_rel;
	map<int,Node*> all_Node; // all nodes in IIN (ID->pointer)
	cout << "Reading IIN File ... " << endl;
	if(prune_IIN == 1)
		read_IIN_prune(IIN_file,stat,num_of_each_type,num_of_each_rel,all_Node,exist_type); // in function.cpp
	else
		read_IIN(IIN_file,stat,num_of_each_type,num_of_each_rel,all_Node); // in function.cpp
	//cout << "Done! Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;

	string IIN_name = "";
	for (int fn = IIN_file.size() - 1; fn>0; fn--)
	{
		if (IIN_file[fn] == '/')
			break;
		IIN_name = IIN_file[fn] + IIN_name;
	}

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

///////////////////////////////////////////////////
	//cout << "Type in the 1 key" << endl;
  //cin >> smth;
//////////////////////////////////////////////////

	// Read influence file
	//cout << "Reading Influence File ... " << endl;
	//map<pair<int,int>,double> prob;
	//read_Inf(Influence_file,prob,all_user);
	//cout << "Done! Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;

	// Doing BFS, only consider the node connect with specific user.
	cout << "Doing BFS To Prune Nodes ..." << endl;
	BFS(all_user[user_node],0);
	vector<int> remove_user;
	for(iter = all_user.begin();iter != all_user.end();iter++)
	{
		if(iter->second->dist == 2147483647)
			remove_user.push_back(iter->first);
	}
	for(size_t i=0;i<remove_user.size();i++)
		all_user.erase(remove_user[i]);
	cout << "Done. Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;
	cout << "Total Nodes: " << all_Node.size() << " ,Total Users: " << all_user.size() << endl;

///////////////////////////////////////////////////
	//cout << "Type in the 1 key" << endl;
	//cin >> smth;
//////////////////////////////////////////////////

	// Main Program
	cout << "Start Compute SIMER.  " << endl;
	//vector<int> source_vec; // record ID of V_s
	//vector<double> SIMER_vec; // record SIMER value of V_s
	vector<pair<int, double> > SIMER_VEC; //record ID,SIMER value of V_s
	vector<vector<Node*> > user_vec;
	for(int i=0;i<NUM_THREADS;i++)
	{
		vector<Node*> vec;
		user_vec.push_back(vec);
	}
	int count_user = 0;
	int mission_part = 0;
	/*for(iter = all_user.begin();iter!=all_user.end();iter++)
	{
		count_user++;
		if(iter->first == user_node)
			continue;
		if(all_user[iter->first]->dist == 2147483647)
			continue;
		if(all_user[iter->first]->dist > Distance)
			continue;
		user_vec[mission_part].push_back(iter->second);
		mission_part++;
		if(mission_part == NUM_THREADS)
			mission_part = 0;
		//if(user_vec[mission_part].size() > (all_user.size()/NUM_THREADS + 1))
			//mission_part++;
	}*/
	for (size_t i =0;i<source_nodes.size();i++)
	{
		if (source_nodes[i] == user_node)
			continue;
		if (all_user[source_nodes[i]]->dist == 2147483647)
			continue;
		if (all_user[source_nodes[i]]->dist > Distance)
			continue;
		user_vec[mission_part].push_back(all_user[source_nodes[i]]);
		mission_part++;
		if (mission_part == NUM_THREADS)
			mission_part = 0;
	}
	//for(size_t i=0;i<user_vec.size();i++)
	//	cout << user_vec[i].size() << endl;

	cout << "Multithread Start" << endl;
	int rc;
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;
	void* status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
	struct ARG arg_set[NUM_THREADS];
	clock_t time_thread = clock();
	for(int i = 0;i < NUM_THREADS;i++)
	{
		arg_set[i].tid = i;
		arg_set[i].tn = target_nodes;
		arg_set[i].un = user_node;
		arg_set[i].d = delta;
		arg_set[i].mt = metaGraph;
		arg_set[i].user_vec = &user_vec[i];
		arg_set[i].an = &all_Node;
		arg_set[i].au = &all_user;
		//arg_set[i].pb = &prob;
		arg_set[i].smk = SMK_directory;
		//arg_set[i].sv = &source_vec;
		//arg_set[i].simerv = &SIMER_vec;
		arg_set[i].Simerv = &SIMER_VEC;
		arg_set[i].cet = CET_location;
		arg_set[i].meta_name = META_file;
		arg_set[i].vt_flag = vt_flag;
		arg_set[i].Spatial_Distance = Spatial_Distance;
		arg_set[i].sim = SIM_directory;
		arg_set[i].iin_name = IIN_name;

//		cout << "main() : creating thread, " << i << endl;
		rc = pthread_create(&threads[i], &attr, Mission, (void *)&arg_set[i] );
		if(rc)
		{
			cout << "Error: Unable to create thread, " << rc << endl;
			exit(-1);
		}
	}
	pthread_attr_destroy(&attr);
	for(int i = 0;i < NUM_THREADS;i++)
	{
		rc = pthread_join(threads[i],&status);
		if(rc)
		{
			cout << "Error: Unable to join, " << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread id: " << i;
		cout << " exiting with status: " << status << endl;
	}
	//cout << "Main: program exiting." << endl;
	cout << "Multithread End" << endl;
	double thread_time_all = double(clock()-time_thread)/double(CLOCKS_PER_SEC);
	double thread_time_avg = double(clock()-time_thread)/double(CLOCKS_PER_SEC)/NUM_THREADS;
	//cout << "Thread. Time Acc.: " << thread_time_all << " sec" << endl;
	//cout << "Thread. Time Acc.: " << thread_time_avg << " sec" << endl;

///////////////////////////////////////////////////
  	//cout << "Type in the 1 key" << endl;
		//cin >> smth;
//////////////////////////////////////////////////

	// sort candidate V_s with SIMER
	/*cout << "Sorting V_s with SIMER ... " << endl;
	size_t N = source_vec.size();
	for(size_t i=0;i+1<N;i++)
	{
		for(size_t j=0;j+1+i<N;j++)
		{
			if(SIMER_vec[j] < SIMER_vec[j+1])
			{
				int tmp1 = source_vec[j];
				source_vec[j] = source_vec[j+1];
				source_vec[j+1] = tmp1;
				double tmp2 = SIMER_vec[j];
				SIMER_vec[j] = SIMER_vec[j+1];
				SIMER_vec[j+1] = tmp2;
			}
		}
	}
	cout << "Done. Time Acc.: " << double(clock()-time1)/double(CLOCKS_PER_SEC) << " sec" << endl;*/
	Topk_compare topk_compare;
	sort(SIMER_VEC.begin(), SIMER_VEC.end(), topk_compare);

	// Write Files
	/*fstream f;
	f.open(topK_file.c_str(),ios::out);
	for(size_t i=0;i<size_t(topK) && i != source_vec.size();i++)
		f << source_vec[i] << endl;
	f.close();
	f.open(detail_file.c_str(),ios::out);
	for(size_t i=0;i<source_vec.size();i++)
		f << source_vec[i] << "," << SIMER_vec[i] << "," << all_user[source_vec[i]]->dist << endl;
	f << "Time Usage: "<< double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds. (take summation of multi-thread time)" << endl;
	f << "Time Usage: "<< double(clock()-time1)/double(CLOCKS_PER_SEC)-thread_time_all+thread_time_avg;
	f << " seconds. (take average of multi-thread time)" << endl;*/
	fstream f;
	f.open(detail_file.c_str(), ios::out);
	for (size_t i = 0; i < SIMER_VEC.size(); i++)
		f << SIMER_VEC[i].first << "," << SIMER_VEC[i].second <<"," << all_user[SIMER_VEC[i].first]->dist << endl;
	f << "Time Usage: " << double(clock() - time1) / double(CLOCKS_PER_SEC) << " seconds. (take summation of multi-thread time)" << endl;
	f << "Time Usage: " << double(clock() - time1) / double(CLOCKS_PER_SEC) - thread_time_all + thread_time_avg;
	f << " seconds. (take average of multi-thread time)" << endl;
	f.close();


	cout << "Time Usage: "<< double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds. (take summation of multi-thread time)" << endl;
	cout << "Time Usage: "<< double(clock()-time1)/double(CLOCKS_PER_SEC)-thread_time_all+thread_time_avg;
	cout << " seconds. (take average of multi-thread time)" << endl;
}
