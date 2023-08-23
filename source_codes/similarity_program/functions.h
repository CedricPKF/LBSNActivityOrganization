#ifndef FUNCTIONS_H
#define FUNCTIONS_H

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
#include "Node.h"

using namespace std;

int str2int(string l);

double str2double(string l);

string int2str(int a);

void get_num(char buffer[],map<int,int>& m);

void setup_metaGraph(MetaGraph* g,int l,int uid,int utype,int e,int vid,int vtype);

void Traversal(MetaGraph* graph,int layer,CETNode* CETn,int& total_ins);

void Cartesian(vector<vector<Node*> >& sigma,vector<vector<Node*> >& Ins,vector<Node*> temp, int index);

void read_IIN(string FileName,int stat[],map<int,int>& num_of_each_type,map<int,int>& num_of_each_rel,map<int,Node*>& all_Node);
void read_IIN_prune(string FileName,int stat[],map<int,int>& num_of_each_type,map<int,int>& num_of_each_rel,map<int,Node*>& all_Node,vector<int>& exist_type);

void read_META(string FileName,MetaGraph* metaGraph);

double Distinctiveness(int vs,CETNode* CETn,vector<vector<Node*> >& instance,int AS,int& Total_ins,vector<Node*>& All_N_Set,vector<Node*>& All_T_Set,int ds,int vt,int vt_flag);

void DFS(int vs,CETNode* CETn,vector<vector<Node*> >& instance,int& Total_ins,vector<Node*>& All_N_Set,vector<Node*>& All_T_Set,int ds,int vt,int vt_flag);

//double IMDR(int vs,int vt,MetaGraph* metaGraph,map<int,double>& omega,CETNode* CET,map<int,Node*>& all_Node,double delta,int vt_flag);
double IMDR(int vs, int vt, MetaGraph* metaGraph, map<int, double>& omega, vector<D_TABLE_information*>& d_table_list, map<pair<int, int>, CET_FILE_information*>& cet_file_list, map<int, Node*>& all_Node, double delta, int vt_flag, double Spatial_distance);

double sim(Node* a,Node* b);

//double SIMER(int vs,int vt,int vu,map<int,Node*>& all_user,map<int,Node*>& all_Node,MetaGraph* metaGraph,map<int,double>& omega,CETNode* CET,double delta,int vt_flag);
double SIMER(int vs, int vt, int vu, map<int, Node*>& all_user, map<int, Node*>& all_Node, MetaGraph* metaGraph, map<int, double>& omega, vector<D_TABLE_information*>& d_table_list, map<pair<int, int>, CET_FILE_information*>& cet_file_list, double delta, int vt_flag, double Spatial_distance);

void BFS(Node* n,int dist);

void read_Inf(string FileName,map<pair<int,int>,double>& prob,map<int,Node*>& all_user);
void check_in_short_distance(map<int,Node*>& VS,map<int,Node*>& all_Node,int cur_dist,int ds,Node* n);
void compute_omega(map<int,double>& omega,map<pair<int,int>,double>& prob,map<int,Node*>& VS,int source_node,int ds,map<int,Node*>& all_user);
void read_smk(string smk,map<int,double>& omega);

#endif
