#include<cstdlib>
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<ctime>
#include<algorithm>
#include<map>
#include<utility>
#include<cmath>
#include<pthread.h>
#include<unistd.h>
#include<iomanip>
#include "Node.h"
#include "functions.h"
using namespace std;

int main(int argc, char* argv[])
{
	clock_t time1 = clock();
	string IIN_file = argv[1];
	int Node_Num = str2int(argv[2]);
	string location = argv[3];

	// Read IIN File
	// num_obj, num_link, num_type, num_rel
	int stat[4] = { -1,-1,-1,-1 };
	map<int, int> num_of_each_type;
	map<int, int> num_of_each_rel;
	map<int, Node*> all_Node; // all nodes in IIN (ID->pointer)
	map<int, Node*>::iterator iter;
	map<int, Node*>::iterator iter_two;
	map<int, Node*> Node_need_to_calculate;

	cout << "Reading IIN File ... " << endl;
	read_IIN(IIN_file, stat, num_of_each_type, num_of_each_rel, all_Node); // in function.cpp
	cout << "Done! Time Acc.: " << double(clock() - time1) / double(CLOCKS_PER_SEC) << " sec" << endl;

	if (all_Node.find(Node_Num) == all_Node.end())
	{
		cout << "Target Node is not in INN map.";
		return 1;
	}

	Node* target_Node = all_Node[Node_Num];

	//Collect Node need to be calculated.
	cout << "Collecting Nodes which need to be cal sim." << endl;
	for (iter = all_Node.begin(); iter != all_Node.end(); iter++) 
	{
		Node* n = iter->second;
		if (n->type == target_Node->type) {
			Node_need_to_calculate[n->num] = n;
		}
	}
	cout << "Done!" << endl;

	//Calculating
	map<int,double> sim_tn_node;//target node's similarity with all node with same type in iin map.
	cout << "Calculating.." << endl;
	for (iter = Node_need_to_calculate.begin(); iter != Node_need_to_calculate.end(); iter++) {
		sim_tn_node[iter->first] = sim(target_Node, iter->second);
	}
	cout << "Done!" << endl;

	string IIN_name = "";
	for (int fn = IIN_file.size() - 1; fn>0; fn--)
	{
		if (IIN_file[fn] == '/')
			break;
		IIN_name = IIN_file[fn] + IIN_name;
	}
	cout << "Outputing!" << endl;

	string filename = location + "/SIM_" + IIN_name + "_n" + int2str(Node_Num) + ".txt";
	fstream fs;
	fs.open(filename.c_str(), ios::out);
	
	for (iter = Node_need_to_calculate.begin(); iter != Node_need_to_calculate.end(); iter++) {
		fs << iter->first << "," << sim_tn_node[iter->first] << endl;
	}

	fs.close();

	cout << "Time Usage: " << double(clock() - time1) / double(CLOCKS_PER_SEC) << " seconds." << endl;
}