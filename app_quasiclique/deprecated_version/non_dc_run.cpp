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

//additional pruning rules:
//1. in spawn, v only pulls 1-hop nbs > v; in step 1, only pulls 2-hop nbs > rootID
//2. in spawn, if degree too low (w.r.t. gamma & min_size), prune (ToDo)
//3. for pulled vertices, prune it if degree too low

#include <chrono>

#include "../deprecated_version/kcore.h"

struct ContextValue
{
	int round; // the iteration number in compute(.)
	map<VertexID, kc_value> g_map;  // for k-core pruning: g_map[v] = (neighbors, v_deleted?)
};
obinstream & operator>>(obinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.g_map;
    return m;
}

ibinstream & operator<<(ibinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.g_map;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.g_map;
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.g_map;
    return m;
}

//-------------------------------------------

typedef Task<QCVertex, ContextValue> QCliqueTask; //context = step
double min_deg;

class CliqueComper:public Comper<QCliqueTask>
{
public:

	virtual bool task_spawn(VertexT * v)
	{
		bool result = false;
		//get the most time-consuming part
		if(v->value.size() >= min_deg){
			QCliqueTask * t = new QCliqueTask;
			t->context.round = 1;
			//#1. add a new v and add it at the end
			VertexID vid = v->id;
			for(int i=0; i<v->value.size(); i++)
			{
				VertexID nb = v->value[i];
				//only pull the vetex larger than rootID
				if(nb > vid)
					t->pull(nb);
			}
			QCVertex root_v;
			root_v.id = vid;
			t->subG.addVertex(root_v);//edges are not added until now
			result = t->is_bigtask();
			add_task(t);
		}
		return result;
	}

    virtual bool compute(SubgraphT & g, ContextT & context, vector<VertexT *> & frontier)
    {
    	VertexID rootID = g.vertexes[0].id;
    	map<VertexID, kc_value>& g_map = context.g_map;
    	if(context.round == 1)
    	{
    		//use k-core like algorithm for iterative degree pruning (count adj-items that belong to 2nd-hop in degree)
    		set<VertexID> first_hop; //vertices in 1-ego network
    		g_map[rootID].del = false; //adj is empty
    		first_hop.insert(rootID);
    		set<VertexID>  to_pull, del_v;
    		//#2. add degree-satisfied frontier vertices (empty adj-list) to graph
    		vector<bool> bitmap(frontier.size(), false);
    		for(int i = 0; i < frontier.size(); i++) {
    			first_hop.insert(frontier[i]->id); //build 1-ego network
    			//--------------------------------
				//add degree-satisfied pulled vertex to g
    			if(frontier[i]->value.size() >= min_deg){ //add to g_map
    				bitmap[i] = true;
    				g_map[rootID].nbs.insert(frontier[i]->id);
    				g_map[frontier[i]->id].del = false; //adj is empty
    			}
    			else //add to del_v
    				del_v.insert(frontier[i]->id);
			}
    		//#3. set adj-list for every vertex in g  &  prepare to_pull set
    		for(int i = 0; i < frontier.size(); i++) {
    			if(bitmap[i]){
    				QCValue & front_adj = frontier[i]->value;
    				for(int j = 0; j < front_adj.size(); j++){
    					VertexID front_id = frontier[i]->id;
    					VertexID nb_id = front_adj[j];
    					//only add the edge not delete (vertex in map and 2-hop)
    					if(nb_id >= rootID && del_v.find(nb_id) == del_v.end()) //del_v used here to prune adj-items
    						g_map[front_id].nbs.insert(nb_id);
    				}
    			}
    		}
    		//#4. prune vertices with degree < min_deg (set field "del")
    		//run this because some vertices' adj-lists are filtered using del_v
    		vector<VertexID> to_del;
    		for(auto it = g_map.begin(); it != g_map.end(); ++it)
			{
				kc_value & val = it->second;
				if(val.nbs.size() < min_deg)
					if(!val.del) prune(it->first, val, g_map, min_deg, to_del);
			}
    		//now, g_map is degree-pruned
			//#5_1. see whether root_id is still not pruned
			if(g_map[rootID].del) return false; //the task terminates
			//#5_2. do the actual deletes
			for(auto it = to_del.begin(); it != to_del.end(); ++it) g_map.erase(*it);
			//#6. move gmap's nbs, which is not in 1-hop, to to_pull
			for(auto it = g_map.begin(); it != g_map.end(); ++it){
				//!!!! the root vertex was already added, only need to add neighbor for root vertex
				set<VertexID> & nbs = it->second.nbs;
				set<VertexID> new_nbs;
				for(auto itr = nbs.begin(); itr != nbs.end(); ++itr){
					if(first_hop.find(*itr) != first_hop.end())//nb is in 1-ego net
						new_nbs.insert(*itr);
					else
						to_pull.insert(*itr); //need pulling
				}
				nbs.swap(new_nbs);
			}
			for(auto it = to_pull.begin(); it != to_pull.end(); it++) pull(*it);
    		context.round++;
    		//vector<QCVertex> & g_vec = g.vertexes;
    		return true;
    	}
    	else //context = 2
    	{
    		//#7. construct 2_hop vector
    		set<VertexID> two_hop;
    		for(auto it = g_map.begin(); it != g_map.end(); ++it)
    			two_hop.insert(it->first);
    		for(int i = 0; i < frontier.size(); i++)
    			two_hop.insert(frontier[i]->id);

    		//#8. add the frontier vertices to the map
    		for(int i = 0; i < frontier.size(); i++) {
				QCValue & nbs = frontier[i]->value;
				if(nbs.size() >= min_deg){
					VertexID vid = frontier[i]->id;
					g_map[vid].del = false;
					for(int j = 0; j < nbs.size(); j++){
						VertexID nb_id = nbs[j];
						if(nb_id >= rootID && two_hop.find(nb_id) != two_hop.end()){
							g_map[vid].nbs.insert(nb_id);
							g_map[nb_id].nbs.insert(vid);
						}
					}
				}
			}

    		//#9. kcore prune 2_hop g_map
    		vector<VertexID> to_del;
			for(auto it = g_map.begin(); it != g_map.end(); ++it)
			{
				kc_value & val = it->second;
				if(val.nbs.size() < min_deg)
					if(!val.del) prune(it->first, val, g_map, min_deg, to_del);
			}
			//now, g_map is degree-pruned
			//see whether root_id is still not pruned
			if(g_map[rootID].del) return false; //the task terminates
			//do the actual deletes
			for(auto it = to_del.begin(); it != to_del.end(); ++it) g_map.erase(*it);

			//#10. convert the g_map to graph
			g.vertexes.clear();
			g.vmap.clear();
			for(auto it = g_map.begin(); it != g_map.end(); ++it){
				QCVertex v;
				v.id = it->first;
				set<VertexID> & nbs = it->second.nbs;
				QCValue & adj = v.value;
				for(auto itr = nbs.begin(); itr != nbs.end(); ++itr)
					adj.push_back(*itr);
				g.addVertex(v);
			}
			g_map.clear();
			vector<QCVertex>& vertices = g.vertexes;

    		//the 2-hop graph is complete until now.
			//#11. prepare graph S and candidate for QC, call QCQ()
    		QCSubgraph gs;
    		QCVertex v;
    		v.id = rootID;
    		gs.addVertex(v);
    		vector<QCVertex*> cand_exts;
    		//exclude the root vertex
    		for(int i = 1; i < vertices.size(); i++){
    			assert(vertices[i].value.size() >= min_deg); //additional degree prune on 2-hop graph
    			cand_exts.push_back(&vertices[i]);
    		}

    		/*//for quasiclique performance analysis
    		int edge_num = 0; //g's edge number
			int degree; //vertex's degree
			for(int i = 0; i < vertices.size(); i++){
				degree = vertices[i].value.size();
				edge_num += degree;
			}
			edge_num = edge_num/2;

			task_fout<<rootID<<" "<<vertices.size()<<" "<<edge_num<<" "<<g.getVertex(rootID)->value.size()<<" ";
			//get k-core map
			map<int, int> kc_map;
			k_core(g, min_deg, kc_map);
			//print map
			map<int,int>::reverse_iterator rit = kc_map.rbegin();
			task_fout<<rit->first<<endl;*/

			QCQ(gs, g, cand_exts, fout);
    		return false;
    	}
    }
};

class CliqueWorker:public Worker<CliqueComper>
{
public:
    CliqueWorker(int num_compers) : Worker<CliqueComper>(num_compers){
    	min_deg = ceil(_gamma * (min_size - 1));
    }

    virtual VertexT* toVertex(char* line)
    {
        VertexT* v = new VertexT;
        char * pch;
        pch=strtok(line, " \t");
        v->id=atoi(pch);
        strtok(NULL," \t");
        QCValue & val = v->value;
        while((pch=strtok(NULL, " ")) != NULL)
        {
            val.push_back(atoi(pch));
        }
        sort(val.begin(), val.end());
        return v;
    }

    virtual void task_spawn(VertexT * v, vector<TaskT> & tcollector)
	{
    	if(v->value.size() >= min_deg){
    		TaskT t;
			tcollector.push_back(t);
			TaskT & task = tcollector.back();
			task.context.round = 1;
			VertexID vid = v->id;
			for(int i=0; i<v->value.size(); i++)
			{
				VertexID nb = v->value[i];
				if(nb > vid)
					task.pull(nb);
			}
			QCVertex root_v;
			root_v.id = vid;
			task.subG.addVertex(root_v);
    	}
	}
};

int main(int argc, char* argv[])
{
    init_worker(&argc, &argv);
    WorkerParams param;
    if(argc != 5){
    	cout<<"arg1 = input path in HDFS, arg2 = number of threads"
    			<<", arg3 = degree ratio, arg4 = min_size"<<endl;
    	return -1;
    }
    param.input_path = argv[1];  //input path in HDFS
    int thread_num = atoi(argv[2]);  //number of threads per process
    _gamma = atof(argv[3]);
    min_size = atoi(argv[4]);
    //min_deg = ceil(_gamma * (min_size - 1));

    param.force_write=true;
    param.native_dispatcher=false;
    //------
    CliqueWorker worker(thread_num);
    worker.run(param);
    //report ------
    //job finished, no need to lock "agg_rwlock"
    if(_my_rank == MASTER_RANK)
    {
    	cout<<"\n Finish."<<endl;
    }
    //-------------
    worker_finalize();
    return 0;
}
