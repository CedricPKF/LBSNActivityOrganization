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
#include "functions.h"
using namespace std;

int str2int(string l)
{
	int ans = 0;
	int digit = l.length();
	for(int i=0;i<digit;i++)
	{
		ans *= 10;
		ans += (l[i]-48);
	}
	return ans;
}

double str2double(string l)
{
	double ans = 0;
	bool digitflag = false;
	int digitcount = 0;
	int len = l.length();
	for(int i=0;i<len;i++)
	{
		if(digitflag)
			digitcount++;
		if(l[i] != '.')
		{
			ans *= 10.0;
			ans += double(l[i]-48);
		}
		else
			digitflag = true;
	}
	for(int i=0;i<digitcount;i++)
		ans /= 10.0;
	return ans;
}

string int2str(int a)
{
	if(a == 0)
	{return "";}
	else
	{
		int b = a%10;
		int c = a/10;
		char d = b+48;
		string s = int2str(c)+d;
		return s;
	}
}

void get_num(char buffer[],map<int,int>& m)
{
	string s = "";
	string ss = "";
	int key;
	int num_key;
	for(int i=0;buffer[i]!='\0';i++)
	{
		if(buffer[i] == ' ')
		{
			for(int j=0;s[j]!='\0';j++)
			{
				if(s[j] == ':')
				{
					key = str2int(ss);
					ss = "";
				}
				else
					ss += s[j];
			}
			num_key = str2int(ss);
			m[key] = num_key;
			s = "";
			ss = "";
		}
		else
			s += buffer[i];
	}
	for(int j=0;s[j]!='\0';j++)
	{
		if(s[j] == ':')
		{
			key = str2int(ss);
			ss = "";
		}
		else
			ss += s[j];
	}
	num_key = str2int(ss);
	m[key] = num_key;
}

void setup_metaGraph(MetaGraph* g,int l,int uid,int utype,int e,int vid,int vtype)
{
	// original format (g,l,u,e,v)
	// new format (g,l,uid,utype,e,vid,vtype)
	map<int,vector<MetaNode*> >::iterator iter;
	iter = g->Layers.find(l);
	if(iter == g->Layers.end())
	{
		vector<MetaNode*> vec;
		g->Layers[l] = vec;
	}
	iter = g->Layers.find(l+1);
	if(iter == g->Layers.end())
	{
		vector<MetaNode*> vec;
		g->Layers[l+1] = vec;
	}
	MetaNode* node_u = NULL;
	MetaNode* node_v = NULL;
	
	for(size_t i=0;i<g->Layers[l].size();i++)
	{
		//if(g->Layers[l][i]->type == u)
		if(g->Layers[l][i]->mid == uid)
		{
			node_u = g->Layers[l][i];
			break;
		}
	}
	for(size_t i=0;i<g->Layers[l+1].size();i++)
	{
		//if(g->Layers[l1][i]->type == v)
		if(g->Layers[l+1][i]->mid == vid)
		{
			node_v = g->Layers[l+1][i];
			break;
		}
	}
	if(node_u == NULL)
	{
		//node_u = new MetaNode(u,l);
		node_u = new MetaNode(utype,l,uid);
		g->Layers[l].push_back(node_u);
	}
	if(node_v == NULL)
	{
		//node_v = new MetaNode(v,l+1);
		node_v = new MetaNode(vtype,l+1,vid);
		g->Layers[l+1].push_back(node_v);
	}
	node_u->fanout[node_v] = e;
	node_v->fanin[node_u] = e;
}

void Traversal(MetaGraph* graph,int layer,CETNode* CETn,int& total_ins)
{
	// if layer == ds: return
	//cout << "Layer: " << layer << endl;
	if(layer == graph->ds)
	{total_ins++;return;}

	// Initialize Ins[]
	vector<vector<Node*> > Ins;
	// for n in S[layer+1] do
	map<int,vector<MetaNode*> >::iterator iter;
	iter = graph->Layers.find(layer+1);
	vector<MetaNode*> N = iter->second;
	for(size_t i=0;i<N.size();i++)
	{
		MetaNode* n = N[i];
		// C <- {}
		vector<vector<Node*> > C;
		vector<Node*> F;
		vector<Node*> E;

		map<MetaNode*,int>::iterator iter_m;
		map<int,vector<Node*> >::iterator iter_i;
		// for (n',n) in M do
		for(iter_m = n->fanin.begin();iter_m != n->fanin.end();iter_m++)
		{
			MetaNode* nf = iter_m->first;
			int edgetype = iter_m->second;
			vector<Node*> node_layer = CETn->content; // modify here
			F.clear();
			// F <- {v|phi(g[n'],v) = (n',n)}
			// C <- C union {F}
			for(size_t j=0;j<node_layer.size();j++)
			{
				if(node_layer[j]->type == nf->type)
				{
					map<Node*,int>::iterator iter_k;
					for(iter_k = node_layer[j]->fanout.begin();iter_k != node_layer[j]->fanout.end();iter_k++)
					{
						if(iter_k->second == edgetype)
						{
							F.push_back(iter_k->first);
						}
					}
				}
			}
			C.push_back(F);
		}

		// E <- intersection C
		if(C.size() == 1)
		{
			E = C[0];
		}
		else if(C.size() > 1)
		{
			for(size_t j=0;j<C[0].size();j++)
			{
				Node* a = C[0][j];
				bool no_same = true;
				for(size_t k=1;k<C.size();k++)
				{
					for(size_t kk=0;kk<C[k].size();kk++)
					{
						if(C[k][kk] == a)
							no_same = false;
					}
					if(no_same)
						break;
				}
				if(!no_same)
					E.push_back(a);
			}
		}
		// Ins[n] <- E
		Ins.push_back(E);
	}


	// sigma <- Cartesian Product of all Ins[n]
	vector<vector<Node*> > sigma;
	vector<Node*> t;
	Cartesian(sigma,Ins,t,0);

	// rtn <- {}
	// for combination in sigma do
	//   g' <- g union combination
	//   I <- Traversal(G,S,g',w',layer+1)
	//   rtn <- rtn union I
	// return rtn
	for(size_t i=0;i<sigma.size();i++)
	{
		// check duplicate item in sigma
		bool dup = false;
		for(size_t j=0;j+1<sigma[i].size();j++)
		{
			for(size_t k=j+1;k<sigma[i].size();k++)
			{
				if(sigma[i][j] == sigma[i][k])
				{
					dup = true;
					break;
				}
			}
			if(dup)
				break;
		}
		if(dup)
			continue;
		CETNode* CETChild = new CETNode();
		CETChild->content = sigma[i];
		CETChild->layer = CETn->layer + 1;
		CETn->child.push_back(CETChild);
		Traversal(graph,layer+1,CETChild,total_ins);
	}
	return;
}

void Cartesian(vector<vector<Node*> >& sigma,vector<vector<Node*> >& Ins,vector<Node*> temp, int index)
{
	if(index == int(Ins.size()))
	{
		sigma.push_back(temp);
		return;
	}
	for(size_t i=0;i<Ins[index].size();i++)
	{
		vector<Node*> temp2 = temp;
		temp2.push_back(Ins[index][i]);
		Cartesian(sigma,Ins,temp2,index+1);
	}
}

double Distinctiveness(int vs,CETNode* CETn,vector<vector<Node*> >& instance,int AS,int& Total_ins,vector<Node*>& All_N_Set,vector<Node*>& All_T_Set,int ds)
{
	double dist = 0.0;
	DFS(vs,CETn,instance,Total_ins,All_N_Set,All_T_Set,ds); // search C-ETree
	dist = double(All_N_Set.size())/(double(Total_ins)*double(AS)); // Definition of Distinctiveness
	return dist;
}
void DFS(int vs,CETNode* CETn,vector<vector<Node*> >& instance,int& Total_ins,vector<Node*>& All_N_Set,vector<Node*>& All_T_Set,int ds)
{
	// reach the C-ETree End
	if(CETn->child.size() == 0 && CETn->layer == ds)
	{
		Total_ins++;
		All_T_Set.push_back(instance[instance.size()-1][0]);
		/*if(instance[0][0]->num == vs && instance[instance.size()-1][0]->num == vt)
			Target_ins++;*/
		// All_U_Set
		for(size_t i=0;i<instance.size();i++)
		{
			for(size_t j=0;j<instance[i].size();j++)
			{
				bool isExist = false;
				for(size_t k=0;k<All_N_Set.size();k++)
				{
				if(instance[i][j] == All_N_Set[k])
					{
						isExist = true;
						break;
					}
				}
				if(!isExist)
					All_N_Set.push_back(instance[i][j]);
			}
		}
		// All_U_Set
	}
	// Keep search C-Etree to find instance
	for(size_t i=0;i<CETn->child.size();i++)
	{
		instance.push_back(CETn->child[i]->content);
		DFS(vs,CETn->child[i],instance,Total_ins,All_N_Set,All_T_Set,ds);
		instance.pop_back();
	}
	return;
}

double IMDR(int vs,int vt,MetaGraph* metaGraph,map<int,double>& omega,CETNode* CET,map<int,Node*>& all_Node,double delta)
{
	vector<vector<Node*> > instance; // Current Instance
	vector<Node*> temp;
	temp.push_back(all_Node[vs]);
	instance.push_back(temp);
	int AS = metaGraph->get_size(); // AS in the definition of meta-structre
	int Total_ins = 0;
	vector<Node*> All_N_Set; // Union i for (V_S)i in definition of Distinctiveness
	vector<Node*> All_T_Set; // total Node for (V_S)i but can repeat (in definition of IMDR last term).

	double d = Distinctiveness(vs,CET,instance,AS,Total_ins,All_N_Set,All_T_Set,metaGraph->ds);
	double imdr = 0.0;

	// compute sum of omega
	double sum_omega = 0.0;
	for(size_t i=0;i<All_N_Set.size();i++)
	{
		if(omega[All_N_Set[i]->num] >= delta)
			sum_omega += omega[All_N_Set[i]->num];
	}

	// compute sum of Sim()
	double sum_sim = 0.0;//
	for(size_t i=0;i<All_T_Set.size();i++)
		sum_sim += sim(All_T_Set[i],all_Node[vt]);

	imdr = sum_omega * sum_sim * d / double(All_N_Set.size()); // Definition of IMDR
	return imdr;
}

double sim(Node* a,Node* b)
{
	double similarity = 0.0;
	vector<Node*> a_nb;
	vector<Node*> b_nb;
	map<Node*,int>::iterator iter;
	for(iter = a->fanin.begin();iter != a->fanin.end();iter++)
		a_nb.push_back(iter->first);
	for(iter = b->fanin.begin();iter != b->fanin.end();iter++)
		b_nb.push_back(iter->first);

	// get the set of intersection and union
	vector<Node*> intersection;
	vector<Node*> unionset;
	for(size_t i=0;i<a_nb.size();i++)
	{
		for(size_t j=0;j<b_nb.size();j++)
		{
			if(a_nb[i] == b_nb[j])
			{
				bool isExist = false;
				for(size_t k=0;k<intersection.size();k++)
				{
					if(a_nb[i] == intersection[k])
					{
						isExist = true;
						break;
					}
				}
				if(!isExist)
					intersection.push_back(a_nb[i]);
			}
		}
	}
	for(size_t i=0;i<a_nb.size();i++)
	{
		bool isExist = false;
		for(size_t k=0;k<unionset.size();k++)
		{
			if(a_nb[i] == unionset[k])
			{
				isExist = true;
				break;
			}
		}
		if(!isExist)
			unionset.push_back(a_nb[i]);
	}
	for(size_t i=0;i<b_nb.size();i++)
	{
		bool isExist = false;
		for(size_t k=0;k<unionset.size();k++)
		{
			if(b_nb[i] == unionset[k])
			{
				isExist = true;
				break;
			}
		}
		if(!isExist)
			unionset.push_back(b_nb[i]);
	}

	similarity = double(intersection.size())/double(unionset.size());
	return similarity;
}

double SIMER(int vs,int vt,int vu,map<int,Node*>& all_user,map<int,Node*>& all_Node,MetaGraph* metaGraph,map<int,double>& omega,CETNode* CET,double delta)
{
	double simer = 0.0;
	simer = IMDR(vs,vt,metaGraph,omega,CET,all_Node,delta)/(all_user[vs]->dist);
	return simer;
}

void read_IIN(string FileName,int stat[],map<int,int>& num_of_each_type,map<int,int>& num_of_each_rel,map<int,Node*>& all_Node)
{
	fstream f;
	f.open(FileName.c_str(),ios::in);
	char buffer[200];
	int line = 0;
	f.getline(buffer,sizeof(buffer));
	while(!f.eof())
	{
		line += 1;
		if(line == 1)
		{
			// read numbers of object, link, type, relationship
			int index = 0;
			string s = "";
			for(int i=0;buffer[i]!='\0';i++)
			{
				if(buffer[i] == ' ')
				{
					stat[index] = str2int(s);
					index += 1;
					s = "";
				}
				else
					s += buffer[i];
			}
			stat[index] = str2int(s);
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		else if(line == 2)
		{
			// read the number of each type
			get_num(buffer,num_of_each_type);
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		else if(line == 3)
		{
			// read the number of each relationship
			get_num(buffer,num_of_each_rel);
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		
		// read the edge
		int temp[5] = {-1,-1,-1,-1,-1};
		for(int i=0;buffer[i]!='\0';i++)
		{
			int index = 0;
			string s = "";
			for(int i=0;buffer[i]!='\0';i++)
			{
				if(buffer[i] == ' ')
				{
					temp[index] = str2int(s);
					index += 1;
					if(index == 5)
						break;
					s = "";
				}
				else
					s += buffer[i];
			}
			temp[index] = str2int(s);
		}
		// Check whether in all_Node
		if(all_Node.find(temp[0]) == all_Node.end())
		{
			Node* n = new Node(temp[0],temp[1]);
			all_Node[temp[0]] = n;
		}
		if(all_Node.find(temp[3]) == all_Node.end())
		{
			Node* n = new Node(temp[3],temp[4]);
			all_Node[temp[3]] = n;
		}
		map<int,Node*>::iterator iter1 = all_Node.find(temp[0]);
		map<int,Node*>::iterator iter2 = all_Node.find(temp[3]);
		iter1->second->fanout[iter2->second] = temp[2];
		iter2->second->fanin[iter1->second] = temp[2];
			
		f.getline(buffer,sizeof(buffer));
	}
	f.close();
}
void read_IIN_prune(string FileName,int stat[],map<int,int>& num_of_each_type,map<int,int>& num_of_each_rel,map<int,Node*>& all_Node,vector<int>& exist_type)
{
	fstream f;
	f.open(FileName.c_str(),ios::in);
	char buffer[200];
	int line = 0;
	f.getline(buffer,sizeof(buffer));
	while(!f.eof())
	{
		line += 1;
		if(line == 1)
		{
			// read numbers of object, link, type, relationship
			int index = 0;
			string s = "";
			for(int i=0;buffer[i]!='\0';i++)
			{
				if(buffer[i] == ' ')
				{
					stat[index] = str2int(s);
					index += 1;
					s = "";
				}
				else
					s += buffer[i];
			}
			stat[index] = str2int(s);
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		else if(line == 2)
		{
			// read the number of each type
			get_num(buffer,num_of_each_type);
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		else if(line == 3)
		{
			// read the number of each relationship
			get_num(buffer,num_of_each_rel);
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		
		// read the edge
		int temp[5] = {-1,-1,-1,-1,-1};
		for(int i=0;buffer[i]!='\0';i++)
		{
			int index = 0;
			string s = "";
			for(int i=0;buffer[i]!='\0';i++)
			{
				if(buffer[i] == ' ')
				{
					temp[index] = str2int(s);
					index += 1;
					if(index == 5)
						break;
					s = "";
				}
				else
					s += buffer[i];
			}
			temp[index] = str2int(s);
		}
		bool exist_flag1 = false;
		bool exist_flag2 = false;
		for(size_t k=0;k<exist_type.size();k++)
		{
			if(temp[1] == exist_type[k])
			{
				exist_flag1 = true;
				break;
			}
		}
		if(!exist_flag1)
		{
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		for(size_t k=0;k<exist_type.size();k++)
		{
			if(temp[4] == exist_type[k])
			{
				exist_flag2 = true;
				break;
			}
		}
		if(!exist_flag2)
		{
			f.getline(buffer,sizeof(buffer));
			continue;
		}
			
		// Check whether in all_Node
		if(all_Node.find(temp[0]) == all_Node.end())
		{
			Node* n = new Node(temp[0],temp[1]);
			all_Node[temp[0]] = n;
		}
		if(all_Node.find(temp[3]) == all_Node.end())
		{
			Node* n = new Node(temp[3],temp[4]);
			all_Node[temp[3]] = n;
		}
		map<int,Node*>::iterator iter1 = all_Node.find(temp[0]);
		map<int,Node*>::iterator iter2 = all_Node.find(temp[3]);
		iter1->second->fanout[iter2->second] = temp[2];
		iter2->second->fanin[iter1->second] = temp[2];
			
		f.getline(buffer,sizeof(buffer));
	}
	f.close();
}

void read_META(string FileName,MetaGraph* metaGraph)
{
	fstream f;
	char buffer[200];
	int line = 0;
	f.open(FileName.c_str(),ios::in);
	f.getline(buffer,sizeof(buffer));
	while(!f.eof())
	{
		line += 1;
		if(line <= 3 )
		{
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		//int temp[3];
		int temp[5]; // with meta node id
		int index = 0;
		string s = "";
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ' ')
			{
				temp[index] = str2int(s);
				index += 1;
				s = "";
				//if(index == 3)
				if(index == 5)
				{
					//setup_metaGraph(metaGraph,line-3,temp[0],temp[1],temp[2]);
					setup_metaGraph(metaGraph,line-3,temp[0],temp[1],temp[2],temp[3],temp[4]); // with meta node id
					index = 0;
				}
			}
			else
				s += buffer[i];
		}
		temp[index] = str2int(s);
		//setup_metaGraph(metaGraph,line-3,temp[0],temp[1],temp[2]);
		setup_metaGraph(metaGraph,line-3,temp[0],temp[1],temp[2],temp[3],temp[4]); // with meta node id
		f.getline(buffer,sizeof(buffer));
	}
	f.close();
	metaGraph->setSource(metaGraph->Layers[1][0]);
	metaGraph->setSink(metaGraph->Layers[line-2][0]);
	metaGraph->setLayer();
}

void BFS(Node* n,int dist)
{
	if(n->dist > dist)
		n->dist = dist;
	else
		return;
	map<Node*,int>::iterator iter;
	for(iter = n->fanout.begin();iter != n->fanout.end();iter++)
	{
		if(iter->first->type == 1)
			BFS(iter->first,dist+1);
	}
}

void read_Inf(string FileName,map<pair<int,int>,double>& prob,map<int,Node*>& all_user)
{
	fstream f;
	f.open(FileName.c_str(),ios::in);
	char buffer[200];
	int line = 0;
	f.getline(buffer,sizeof(buffer));
	while(!f.eof())
	{
		line += 1;
		string s = "";
		int id1 = 0;
		int id2 = 0;
		double value = 0.0;
		int index = 0;
		for(int i=0;buffer[i]!='\0';i++)
		{
			if(buffer[i] == ',' && index == 0)
			{
				id1 = str2int(s);
				index++;
				s = "";
			}
			else if(buffer[i] == ',' && index == 1)
			{
				id2 = str2int(s);
				index++;
				s = "";
			}
			else
				s += buffer[i];
		}
		if(all_user.find(id1) == all_user.end() || all_user.find(id2) == all_user.end())
		{
			f.getline(buffer,sizeof(buffer));
			continue;
		}
		value = str2double(s);
		pair<int,int> p(id1,id2);
		prob[p] = value;
		f.getline(buffer,sizeof(buffer));
	}
	
}
void check_in_short_distance(map<int,Node*>& VS,map<int,Node*>& all_Node,int cur_dist,int ds,Node* n)
{
	// Doing DFS
	if(cur_dist > ds)
		return;
	if(n->flag)
		return;
	//if(VS.find(n->num) == VS.end() && n->type != 1)
	if(n->type != 1)
	{
		n->flag = true;
		//cout << "VS size() = " << VS.size() << endl;
		VS[n->num] = n;
	}
	map<Node*,int>::iterator iter;
	for(iter = n->fanout.begin();iter != n->fanout.end();iter++)
	{
		if(iter->first->flag)
			continue;
		if(iter->first->type == 1)
			continue;
		check_in_short_distance(VS,all_Node,cur_dist+1,ds,iter->first);
	}
}

void compute_omega(map<int,double>& omega,map<pair<int,int>,double>& prob,map<int,Node*>& VS,int source_node,int ds,map<int,Node*>& all_user)
{
	// Initializing omega and some influence between non-individual
	omega[source_node] = 1.0;
	map<int,Node*>::iterator iter;
	// for v in V(S) or U
	for(iter = VS.begin();iter != VS.end();iter++)
	{
		Node* v = iter->second;
		omega[v->num] = 0.0;
		// if v belongs to all_user
		if(all_user.find(v->num) != all_user.end())
		{
			// for u in VS intersect N(v)
			map<Node*,int>::iterator iter2;
			for(iter2 = v->fanout.begin();iter2 != v->fanout.end();iter2++)
			{
				Node* u = iter2->first;
				// intersection
				if(VS.find(u->num) != VS.end())
				{
					pair<int,int> p(u->num,v->num);
					prob[p] = 1.0;
				}
			}
		}
	}
	for(iter = all_user.begin();iter != all_user.end();iter++)
	{
		Node* v = iter->second;
		omega[v->num] = 0;
		if(all_user.find(v->num) != all_user.end())
		{
			map<Node*,int>::iterator iter2;
			for(iter2 = v->fanout.begin();iter2 != v->fanout.end();iter2++)
			{
				Node* u = iter2->first;
				if(VS.find(u->num) != VS.end())
				{
					pair<int,int> p(u->num,v->num);
					prob[p] = 1.0;
				}
			}
		}
	}
	omega[source_node] = 1.0;
	// I <- U
	map<double,vector<Node*> > I;
	vector<Node*> y;
	y.push_back(all_user[source_node]);
	I[omega[source_node]] = y;
	map<double,vector<Node*> >::reverse_iterator rit;
	while(I.size() != 0)
	{
		// find user with max omega
		double omega_tmp = -1.0;
		Node* mu;
		rit = I.rbegin();
		vector<Node*> x = rit->second;
		mu = x[x.size()-1];

		omega_tmp = rit->first;
		rit->second.pop_back();
		if(rit->second.size() == 0)
			I.erase(omega_tmp);
		map<Node*,int>::iterator iter2;
		for(iter2 = mu->fanout.begin();iter2 != mu->fanout.end();iter2++)
		{
			// itersect N(mu) with union(U and VS)
			Node* v = iter2->first;
			if(VS.find(v->num) != VS.end() || all_user.find(v->num) != all_user.end())
			{
				pair<int,int> p(v->num,mu->num);
				if(omega[mu->num]*prob[p] > omega[v->num])
				{
					omega[v->num] = omega[mu->num]*prob[p];
					if(I.find(omega[v->num]) != I.end())
						I[omega[v->num]].push_back(v);
					else
					{
						vector<Node*> xx;
						xx.push_back(v);
						I[omega[v->num]] = xx;
					}
				}
			}
		}
	}
}



