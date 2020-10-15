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

//#include "kcore.h"
#include <algorithm>

#include "../basic_version/graph.h"
int COMPACT_SIZE = 10;

struct ContextValue
{
	int round; // the iteration number in compute(.)
	vector<VertexID> cand_id_vec; // Xcand (to be added to X)
	QCSubgraph gs; // subgraph induced by X
	map<VertexID, kc_value> g_map; // for k-core pruning: g_map[v] = (neighbors, v_deleted?)
	QCHopSubgraph g_2hop; //for 2hop-prune
};



obinstream & operator>>(obinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.gs;
    m >> c.cand_id_vec;
    m >> c.g_map;
    m >> c.g_2hop;
    return m;
}

ibinstream & operator<<(ibinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.gs;
    m << c.cand_id_vec;
    m << c.g_map;
    m << c.g_2hop;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.gs;
    m >> c.cand_id_vec;
    m >> c.g_map;
    m >> c.g_2hop;
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.gs;
    m << c.cand_id_vec;
    m << c.g_map;
    m << c.g_2hop;
    return m;
}

//-------------------------------------------
typedef Task<QCVertex, ContextValue> QCliqueTask;

double min_deg; // = ceil(_gamma * (min_size - 1))
float TIME_THRESHOLD; // if running time >= TIME_THRESHOLD, split into subtasks

class CliqueComper:public Comper<QCliqueTask>
{
public:

	//check whether task is bigtask
	virtual bool is_bigtask(QCliqueTask * task){
		if(task->context.cand_id_vec.size() > BIGTASK_THRESHOLD
				||task->to_pull.size() > BIGTASK_THRESHOLD)
			return true;
		else
			return false;
	}

	virtual bool task_spawn(VertexT * v)
	{
		bool result = false;
		if(v->value.size() >= gnmin_deg
				&& find(v->value.begin(), v->value.end(), max_v) == v->value.end()){
			QCliqueTask * t = new QCliqueTask;
			t->context.round = 1;
			//#1. add a new v and add it at the end; pull 1-hop neighbor
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
			result = is_bigtask(t);
			add_task(t);
		}
		return result;
	}

    virtual bool compute(SubgraphT & g, ContextT & context, vector<VertexT *> & frontier)
    {
    	VertexID rootID = g.vertexes[0].id;
    	map<VertexID, kc_value>& g_map = context.g_map; //used for k-core pruning
    	if(context.round == 1)
    	{
    		//use k-core like algorithm for iterative degree pruning (count adj-items that belong to 2nd-hop in degree)
    		set<VertexID> first_hop; //track vertices in 1-ego network, to avoid pulling them again
    		g_map[rootID].del = false; //add root_v to g_map, v's adj-list is now empty
    		first_hop.insert(rootID);
    		set<VertexID>  to_pull, del_v;
    		//#2. add degree-satisfied frontier vertices (empty adj-list) to graph
    		vector<bool> bitmap(frontier.size(), false); // to filter degree-unsatisfied vertices, not to add to g_map
    		//--- add qualified one-hop neighbors to g_map (just set bitmap) and g_map[root].adjlist
    		for(int i = 0; i < frontier.size(); i++) {
    			first_hop.insert(frontier[i]->id); //build 1-ego network
    			//--------------------------------
				//add degree-satisfied pulled vertex to g
    			if(frontier[i]->value.size() >= gnmin_deg){ //add to g_map
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
    				VertexID front_id = frontier[i]->id;
    				QCValue & front_adj = frontier[i]->value;
    				for(int j = 0; j < front_adj.size(); j++){
    					VertexID nb_id = front_adj[j];
    					//only add the edge not delete (vertex in map and 2-hop)
    					if(nb_id >= rootID && del_v.find(nb_id) == del_v.end()) //del_v used here to prune adj-items
    						g_map[front_id].nbs.insert(nb_id);
    				}
    			}
    		}
    		//#4. prune vertices with degree < min_deg (set field "del")
    		//run this because some vertices' adj-lists are filtered using del_v
    		vector<VertexID> to_del; //output of k-core pruning
    		for(auto it = g_map.begin(); it != g_map.end(); ++it)
			{
				kc_value & val = it->second;
				if(val.nbs.size() < gnmin_deg)
					if(!val.del) prune_1hop(it->first, val, g_map, gnmin_deg, to_del);
			}
    		//now, g_map is degree-pruned
			//#5_1. see whether root_id is still not pruned
			if(g_map[rootID].del) return false; //the task terminates
			//#5_2. do the actual deletes
			for(auto it = to_del.begin(); it != to_del.end(); ++it) g_map.erase(*it);
			//#6. check all adj-items in g_map
			//-- if item is not in 1-hop, move to to_pull
			//-- we remove the item from g_map, since it may be pruned in round 2
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
    		return true;
    	}
    	else if(context.round == 2) //context == 2
    	{
    		//#7. construct 2_hop vector, to filter 3-hop neighbors
			set<VertexID> two_hop;
			// two_hop.add(one_hop vertices in g_map)
			for(auto it = g_map.begin(); it != g_map.end(); ++it)
				two_hop.insert(it->first);
			// two_hop.add(two_hop vertices pulled)
			for(int i = 0; i < frontier.size(); i++)
				two_hop.insert(frontier[i]->id);

			//#8. add qualified frontier vertices to g_map
			for(int i = 0; i < frontier.size(); i++) {
				QCValue & nbs = frontier[i]->value;
				if(nbs.size() >= gnmin_deg){
					VertexID vid = frontier[i]->id;
					g_map[vid].del = false;
					for(int j = 0; j < nbs.size(); j++){
						VertexID nb_id = nbs[j];
						if(nb_id >= rootID &&
								two_hop.find(nb_id) != two_hop.end()){
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
				if(val.nbs.size() < gnmin_deg)
					if(!val.del) prune_1hop(it->first, val, g_map, gnmin_deg, to_del);
			}
			//now, g_map is degree-pruned
			//see whether root_id is still not pruned
			if(g_map[rootID].del) return false; //the task terminates
			//do the actual deletes
			for(auto it = to_del.begin(); it != to_del.end(); ++it) g_map.erase(*it);
			if(g_map.size() < min_size) return false;
			//#10. convert the g_map to graph

			Graph gograph;
			int num_of_vertices = gograph.Cliques(g_map, gfpout);
			context.round++;
			//no time delay version only need 2 round
			return false;
			//### since there is no vertex requested, step 3 will run directly
			//### so we only need round 3 to call QCQ to simplify code
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

    virtual void task_spawn(VertexT * v, vector<TaskT*> & tcollector)
	{
		if(v->value.size() >= gnmin_deg
				&& find(v->value.begin(), v->value.end(), max_v) == v->value.end()){
			TaskT* task = new TaskT;

			task->context.round = 1;
			VertexID vid = v->id;
			for(int i=0; i<v->value.size(); i++)
			{
				VertexID nb = v->value[i];
				if(nb > vid)
					task->pull(nb);
			}
			QCVertex root_v;
			root_v.id = vid;
			task->subG.addVertex(root_v);
			tcollector.push_back(task);
		}
	}
};

int main(int argc, char* argv[])
{
    init_worker(&argc, &argv);
    WorkerParams param;
    if(argc != 6){
    	cout<<"arg1 = input path in HDFS, arg2 = number of threads"
    			<<", arg3 = degree ratio, arg4 = min_size, arg5 = BIGTASK_THRESHOLD"<<endl;
    	return -1;
    }
    param.input_path = argv[1];  //input path in HDFS
    int thread_num = atoi(argv[2]);  //number of threads per process
    gdmin_deg_ratio = atof(argv[3]);
    gnmin_size = atoi(argv[4]);
    gnmax_size = INT_MAX;
    gnmin_deg = ceil(gdmin_deg_ratio * (gnmin_size - 1));
//    _gamma = atof(argv[3]);
//    min_size = atoi(argv[4]);
    BIGTASK_THRESHOLD = atoi(argv[5]);
    max_v = 0;

    //for old gthinker
    _gamma = gdmin_deg_ratio;
    min_size = gnmin_size;
    min_deg = gnmin_deg;

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
