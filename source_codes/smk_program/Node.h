#ifndef NODE_H
#define NODE_H

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
using namespace std;

class Node
{
	public:
		Node()
		{num = -1;type = -1;dist = 2147483647;}
		Node(int n,int t)
		{
			num = n;
			type = t;
			dist = 2147483647;
			flag = false;
		}
		int num;
		int type;
		int dist;
		bool flag;

		// store <fanin/fanout, edgetype>
		map<Node*,int> fanin;
		map<Node*,int> fanout;
		void printInfo()
		{
			cout << "==============================" << endl;
			cout << "# object number = " << num << endl;
			cout << "# object type   = " << type << endl;
			map<Node*,int>::iterator iter;
			if(fanin.size() > 0)
				cout << "<<< Fanin <<<" << endl;
			for(iter = fanin.begin();iter != fanin.end();iter++)
			{
				cout << "fanin num " << iter->first->num << " with edge type " << iter->second << endl;
			}
			if(fanout.size() > 0)
				cout << ">>> Fanout >>>" << endl;
			for(iter = fanout.begin();iter != fanout.end();iter++)
			{
				cout << "fanout num " << iter->first->num << " with edge type " << iter->second << endl;
			}
		}
};

class MetaNode
{
	public:
		MetaNode()
		{type = -1;layer = -1;}
		MetaNode(int t,int l,int i)
		{
			type = t;
			layer = l;
			mid = i;
		}
		int type;
		int layer;
		int mid;
		map<MetaNode*,int> fanin;
		map<MetaNode*,int> fanout;
};

class MetaGraph
{
	public:
		MetaGraph(){}
		void setSource(MetaNode* s) {source = s;}
		void setSink(MetaNode* s) {sink = s;}
		void setLayer() {ds = sink->layer;}
		int get_size()
		{
			int a = 0;
			map<int,vector<MetaNode*> >::iterator iter;
			for(iter = Layers.begin();iter != Layers.end();iter++)
				a += int(iter->second.size());
			return a;
		}
		int get_type_num()
		{
			vector<int> a;
			map<int,vector<MetaNode*> >::iterator iter;
			for(iter = Layers.begin();iter != Layers.end();iter++)
			{
				for(size_t i=0;i<iter->second.size();i++)
				{
					bool isExist = false;
					for(size_t j=0;j<a.size();j++)
					{
						if(iter->second[i]->type == a[j])
						{
							isExist = true;
							break;
						}
					}
					if(!isExist)
						a.push_back(iter->second[i]->type);
				}
			}
			return int(a.size());
		}
		
		MetaNode* source;
		MetaNode* sink;
		int ds;
		// put the node with same layer together
		map<int,vector<MetaNode*> > Layers;
		
};

class Instance
{
	public:
		Instance()
		{
			w = 1.0;
			/*SC = 1.0;
			SCSE = 1.0;
			BSCSE = 1.0;*/
		}
		Instance(Node* n)
		{
			root = n;
			w = 1.0;
			/*SC = 1.0;
			SCSE = 1.0;
			BSCSE = 1.0;*/
		}
		
		void copy_instance(Instance* i)
		{
			root = i->root;
			Layers = i->Layers;
			w = i->w;
			//SC = i->SC;
			//SCSE = i->SCSE;
			//BSCSE = i->BSCSE;
		}

		void Show_Instance(fstream& f)
		{
			for(size_t j=1;j<=Layers.size();j++)
			{
				f << "Layer " << j << ": ";
				for(size_t k=0;k<Layers[j].size();k++)
					f << Layers[j][k]->num << " ";
				f << endl;
			}
		}
		map<int,vector<Node*> > Layers;
		Node* root;

		double w;
		//double SC;
		//double SCSE;
		//double BSCSE;

};

class CETNode
{
	public:
		CETNode(){}
		void addNode(Node* n)
		{
			content.push_back(n);
		}
		vector<CETNode*> child;
		vector<Node*> content;
		int layer;
};

#endif
