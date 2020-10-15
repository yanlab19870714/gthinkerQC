//########################################################################
//## Copyright 2018 Da Yan http://www.cs.uab.edu/yanda
//##
//## Licensed under the Apache License, Version 2.0 (the "License");
//## you may not use this file except in compliance with the License.
//## You may obtain a copy of the License at
//##
//## //http://www.apache.org/licenses/LICENSE-2.0
//##
//## Unless required by applicable law or agreed to in writing, software
//## distributed under the License is distributed on an "AS IS" BASIS,
//## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//## See the License for the specific language governing permissions and
//## limitations under the License.
//########################################################################

//test whether the candidate order matter.

#include "kcore.h"
#define hash_set __gnu_cxx::hash_set

using namespace std;
typedef int VertexID;

using get_time = chrono::high_resolution_clock;
using ms = chrono::microseconds;

bool COVER_VERTEX_PRUNE = true;
bool LOOKAHEAD_PRUNE = true;
bool UPPER_BOUND_PRUNE = true;
bool LOWER_BOUND_PRUNE = true;
bool CRITICAL_VERTEX_PRUNE = true;

int COMPACT_SIZE = 10;

using get_time = chrono::high_resolution_clock;
using ms = chrono::microseconds;

double _gamma;
int min_size;
int expand_count = 0;

void load(char* fname, map<VertexID, kc_value>& g_map) {
	ifstream in(fname);
	VertexID vid;
	int temp, num;
	while(in >> vid)
	{
		g_map[vid].del = false;
		in >> num;
		QCValue nbs;
		for(int i=0; i<num; i++)
		{
			in >> temp;
			g_map[vid].nbs.insert(temp);
		}
	}
    in.close();
}

//put the vertex whose degree is max at the beginning
void max_deg_ahead(vector<QCVertex*> & cand_exts)
{
	//first find the max degree index
	//using std swap
	int max_index = 0;
	int max_value = 0;
	for(int i = 0; i<cand_exts.size();i++)
	{
		int v_deg = cand_exts[i]->value.size();
		if(v_deg > max_value)
		{
			max_index = i;
			max_value = v_deg;
		}
	}
	iter_swap(cand_exts.begin(), cand_exts.begin()+max_index);
}

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        cout<<"arg1 = input file name, arg2 = degree ratio, arg3 = min_size"<<endl;
        return -1;
    }
    string fn = argv[1];
    map<VertexID, kc_value> g_map;
	load(argv[1], g_map);
    cout<<"loaded"<<endl;

	string report_path = fn + "_non_kcore";
	ofstream fout;
	fout.open(report_path);

	vector<pair<int, VertexID> > deg_sorter; //<size, id>
	int max_deg = 0;
	VertexID root_id;
	for(auto it = g_map.begin(); it != g_map.end(); ++it)
	{
		deg_sorter.push_back(make_pair(it->second.nbs.size(), it->first));
//		deg_sorter[it->second.nbs.size()] = it->first;
		if(it->second.nbs.size() > max_deg)
		{
			root_id = it->first;
			max_deg = it->second.nbs.size();
		}
	}


	map<VertexID, int> id2index;
	int g_size = deg_sorter.size();
	vector<VertexID> order(g_size); //index -> vid

	sort(deg_sorter.begin(), deg_sorter.end());
//	VertexID root_id = deg_sorter[g_size-1].second;
	order[0] = root_id;
	id2index[root_id] = 0;

	set<VertexID>& root_nbs = g_map.find(root_id)->second.nbs;
	int i = 1;
	int j = g_size - root_nbs.size();
	for(auto ds: deg_sorter)
	{
		VertexID vid = ds.second;
		if(vid == root_id)
			continue;
		if(root_nbs.find(vid) == root_nbs.end())
		{
			order[i] = vid;
			id2index[vid] = i;
			i++;
		}
		else
		{
			order[j] = vid;
			id2index[vid] = j;
			j++;
		}
	}

	vector<VertexID> cover_vec;

	//print recoded id
	for(auto v : order)
	{
		auto it = g_map.find(v);
		set<VertexID>& nbs = it->second.nbs;
		fout<<id2index[v]<<" "<<nbs.size();
		bool round_1st = true;
		for(auto nb : nbs)
		{
			if(round_1st)
			{
				fout<<"\t"<<id2index[nb];
				round_1st = false;
			}
			else
				fout<<" "<<id2index[nb];
		}
		fout<<endl;
	}

    fout.close();
    return 0;
}
