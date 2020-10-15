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

#include "../deprecated_version/quasiclique.h"

struct kc_value
{
	bool del;
	set<VertexID> nbs;
	set<VertexID> nbs2hop;
};

obinstream & operator>>(obinstream & m, kc_value & v)
{
    m >> v.del;
    m >> v.nbs;
    m >> v.nbs2hop;
    return m;
}

ibinstream & operator<<(ibinstream & m, const kc_value & v)
{
	m << v.del;
    m << v.nbs;
    m << v.nbs2hop;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, kc_value & v)
{
    m >> v.del;
    m >> v.nbs;
    m >> v.nbs2hop;
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const kc_value & v)
{
	m << v.del;
    m << v.nbs;
    m << v.nbs2hop;
    return m;
}

//2hop prune
void prune_2hop(VertexID id, kc_value & kc_val, map<VertexID, kc_value> & kc_g, int k, int min_size, vector<VertexID> & to_del){
	kc_val.del = true;
	to_del.push_back(id);
	set<VertexID> & nbs2hop = kc_val.nbs2hop;
	// 2hop is the superset of 1hop
	for(auto it = nbs2hop.begin(); it != nbs2hop.end(); ++it){
		int nb_id = *it;
		auto itr = kc_g.find(nb_id);
		if(itr != kc_g.end() && itr->second.del == false){
			kc_value & nb_val = itr->second;
			set<VertexID> & nb_nbs = nb_val.nbs;
			set<VertexID> & nb_nbs2hop = nb_val.nbs2hop;
			//delete v from its nb's 2hop set
			nb_nbs2hop.erase(id);
			//delete v from its nb's 1hop set
			auto nb_it = nb_nbs.find(id);
			if(nb_it != nb_nbs.end())//maybe not 1hop neighbor
				nb_nbs.erase(nb_it);
			//prune recursively
			if(nb_nbs.size() < k || nb_nbs2hop.size() < min_size - 1)
				if(!nb_val.del)
					prune_2hop(nb_id, nb_val, kc_g, k, min_size, to_del);
		}
	}
}

//1hop prune
void prune_1hop(VertexID id, kc_value & kc_val, map<VertexID, kc_value> & kc_g, int k, vector<VertexID> & to_del){
	kc_val.del = true;
	to_del.push_back(id);
	set<VertexID> & nbs = kc_val.nbs;
	for(auto it = nbs.begin(); it != nbs.end(); ++it){
		int nb_id = *it;
		if(kc_g.find(nb_id) != kc_g.end() //ignore nb_id that has not been pulled
				&& kc_g[nb_id].del == false){
			kc_value & nb_val = kc_g[nb_id];
			set<VertexID> & nb_nbs = nb_val.nbs;
			nb_nbs.erase(id);
			if(nb_nbs.size() < k)
				if(!nb_val.del) prune_1hop(nb_id, nb_val, kc_g, k, to_del);
		}
	}
}

void k_core(QCSubgraph & g, int min_deg, map<int, int> & result_map){
	int k = min_deg;
	vector<QCVertex>& vertices = g.vertexes;
	map<VertexID, kc_value> kc_g;
	VertexID root_id = vertices[0].id;
	for(int i = 0; i < vertices.size(); i++){
		QCVertex & v = vertices[i];
		QCValue & nbs = v.value;
		kc_value kc_val;
		kc_val.del = false;
		set<VertexID> & kc_nbs = kc_val.nbs;
		kc_nbs.insert(nbs.begin(), nbs.end());
		kc_g[v.id] = kc_val;
	}
	while(true){
		//prune vertices with degree < k (set field "del")
		vector<VertexID> to_del;
		for(auto it = kc_g.begin(); it != kc_g.end(); ++it)
		{
			kc_value & val = it->second;
			if(val.nbs.size() < k)
				if(!val.del) prune_1hop(it->first, val, kc_g, k, to_del);
		}
		//see whether root_id is still not pruned
		if(kc_g[root_id].del){
			/*cout<<"kcore is "<<k-1<<endl;
			cout<<"root id is deleted";*/
			return;
		}
		//do the actual deletes
		for(auto it = to_del.begin(); it != to_del.end(); ++it) kc_g.erase(*it);
		/*//-----------debug--------
		set<VertexID> & check_nbs = kc_g[root_id].nbs;
		cout<<k<<" : "<<kc_g.size()<<endl;
		auto it = check_nbs.begin();
		while (it != check_nbs.end()){
			cout<<*it<<" ";
			it++;
		}
		cout<<endl<<endl;
		for(auto itr = kc_g.begin(); itr != kc_g.end(); ++itr){
			cout<<itr->first<<" ";
		}
		cout<<endl<<"--------------------"<<endl;

		//-----------debug--------*/
		//save the results
		result_map[k] = kc_g.size();
		k++;
	}
}
