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

#include <algorithm>

#include "../deprecated_version/kcore.h"
int COMPACT_SIZE = 10;

struct ContextValue
{
	int round; // the iteration number in compute(.)
	vector<VertexID> cand_id_vec; // Xcand (to be added to X)
	QCSubgraph gs; // subgraph induced by X
	map<VertexID, kc_value> g_map; // for k-core pruning: g_map[v] = (neighbors, v_deleted?)
};



obinstream & operator>>(obinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.gs;
    m >> c.cand_id_vec;
    m >> c.g_map;
    return m;
}

ibinstream & operator<<(ibinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.gs;
    m << c.cand_id_vec;
    m << c.g_map;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.gs;
    m >> c.cand_id_vec;
    m >> c.g_map;
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.gs;
    m << c.cand_id_vec;
    m << c.g_map;
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
		if(v->value.size() >= min_deg){
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

	// After pruning, only those vertices in X and Xcand remain in the pruned g
	void graph_shrink(QCSubgraph & pre_gs, vector<QCVertex*> & cand, QCSubgraph & g,
			QCSubgraph & shrinked_G, vector<VertexID>& cand_id_vec){
		//Input:
		//  pre_gs is the subgraph induced by X
		//  cand is Xcand
		//Output:
		//  shrinked_G is the output to set, initially empty
		//  cand_id_vec is the output Xcand
		//  task_gs is the adjlist-shrinked subgraph induced by X
		set<VertexID> remain_set;
		// we have updated X and candX by iterative_bounding, remain_set = X U candX
		vector<QCVertex>& vertices = pre_gs.vertexes;
		int new_gs_size = vertices.size();
		int cand_size = cand.size();

		// remain_set = remain_set U X
		for(int i = 0; i < new_gs_size; i++)
			remain_set.insert(vertices[i].id);
		// remain_set = remain_set U Xcand
		for(int i = 0; i < cand_size; i++)
			remain_set.insert(cand[i]->id);

		//For each vertex in X, delete the adj-items not in remain_set
		for(int i = 0; i < new_gs_size; i++){
			QCVertex* it = g.getVertex(vertices[i].id);
			QCVertex v;
			v.id = it->id;
			QCValue temp_nbs;
			for(VertexID adj : it->value)
				if(remain_set.find(adj) != remain_set.end())
					temp_nbs.push_back(adj);

			v.value.swap(temp_nbs);
			shrinked_G.addVertex(v);
		}

		//For each vertex in Xcand, delete the adj-items not in remain_set
		for(int i = 0; i < cand_size; i++){
			QCVertex v;
			v.id = cand[i]->id;

			QCValue temp_nbs;
			for(VertexID adj: cand[i]->value)
				if(remain_set.find(adj) != remain_set.end())
					temp_nbs.push_back(adj);

			cand_id_vec.push_back(v.id);
			v.value.swap(temp_nbs);
			shrinked_G.addVertex(v);
		}
	}


	//revise 5: shrink g
	void compact_graph(QCSubgraph & gs, vector<QCVertex*> & cand, QCSubgraph & g,
			QCSubgraph & new_g)
	{
		//get the set of gs and cand ==> 2-hop-set  //which will ingore checking v before gs
		//create new_g based on 2-hop-set using g
		//the cand is QCVertex*, so it will change in new_g... be careful.
		set<VertexID> s_cand_set; //Union of gs and cand
		for(auto v : gs.vertexes)
			s_cand_set.insert(v.id);
		for(auto v : cand)
			s_cand_set.insert(v->id);


		for(auto sv : gs.vertexes)
		{
			QCVertex* it = g.getVertex(sv.id);
			QCVertex v;
			v.id = it->id;
			QCValue temp_nbs;
			for(VertexID adj : it->value)
				if(s_cand_set.find(adj) != s_cand_set.end())
					temp_nbs.push_back(adj);
			v.value.swap(temp_nbs);
			new_g.addVertex(v);
		}

		for(auto cv : cand)
		{
			QCVertex v;
			v.id  = cv->id;

			QCValue temp_nbs;
			for(VertexID adj: cv->value)
				if(s_cand_set.find(adj) != s_cand_set.end())
					temp_nbs.push_back(adj);

			v.value.swap(temp_nbs);
			new_g.addVertex(v);
		}

		vector<QCVertex*> new_cand;
		for(auto v : cand)
			new_cand.push_back(new_g.getVertex(v->id));
		cand.swap(new_cand);
	}

	struct deg_sorter{
		QCVertex* v;
		int indeg = 0;
		int exdeg = 0;
	};
	static bool comp(const deg_sorter &x, const deg_sorter &y)
	{
		if(x.indeg == y.indeg)
		{
			if(x.exdeg == y.exdeg)
				return x.v->id < y.v->id;
			else
				return x.exdeg < y.exdeg;
		}
		else
			return x.indeg < y.indeg;
	}

	void get_cand_deg(vector<QCVertex*>& cand, QCSubgraph& X_g, vector<deg_sorter>& cand_deg)
	{
		set<VertexID> cand_set;
		for(auto v:cand)
			cand_set.insert(v->id);

		int n_cand = cand.size();
		hash_map<VertexID, int> & Xg_map = X_g.vmap;
//		cand_deg[0].v = cand[0];
	//	for (int j = 1; j < n_cand; j++)
		for (int j = 0; j < n_cand; j++)
		{
			cand_deg[j].v = cand[j];
			QCValue & cand_j_adj = cand[j]->value;
			int adj_size = cand_j_adj.size();
			for (int k = 0; k < adj_size; k++)
			{
				VertexID nb = cand_j_adj[k];
				if(Xg_map.find(nb) != Xg_map.end())
					cand_deg[j].indeg++;
				if(cand_set.find(nb) != cand_set.end())
					cand_deg[j].exdeg++;
			}
		}
	}

	void get_cand_deg_1st(vector<QCVertex*>& cand, vector<deg_sorter>& cand_deg)
	{
		int n_cand = cand.size();
//		cand_deg[0].v = cand[0];
//		for (int j = 1; j < n_cand; j++)
		for (int j = 0; j < n_cand; j++)
		{
			cand_deg[j].v = cand[j];
			cand_deg[j].exdeg = cand[j]->value.size();
		}
	}

	void sort_deg(vector<QCVertex*> & cand_exts, QCSubgraph & gs){
		vector<deg_sorter> cand_deg(cand_exts.size());
		if(!gs.vertexes.empty())
			get_cand_deg(cand_exts, gs, cand_deg);
		else
			get_cand_deg_1st(cand_exts, cand_deg);
		sort(cand_deg.begin(), cand_deg.end(), comp);
//		sort(cand_deg.begin()+1, cand_deg.end(),comp);
		cand_exts.clear();
		for(auto c: cand_deg)
			cand_exts.push_back(c.v);
	}

	//set the 2hop set in g_map
	void set_2hop(map<VertexID, kc_value>& g_map)
	{
		for(auto it = g_map.begin(); it != g_map.end(); ++it)
		{
			VertexID vid = it->first;
			kc_value & val = it->second;
			set<VertexID>& nbs = val.nbs;
			set<VertexID>& nbs_2hop = val.nbs2hop;
			for(auto nb : nbs)
			{
				nbs_2hop.insert(nb);
				auto nb_it = g_map.find(nb);
				assert(nb_it != g_map.end());
				set<VertexID>& nb_nbs = nb_it->second.nbs;
				if(!nb_nbs.empty())
					nbs_2hop.insert(nb_nbs.begin(), nb_nbs.end());
			}
		}
	}

	bool tddq_QCQ(QCSubgraph & gs, QCSubgraph & g, vector<QCVertex*> & cand_exts,
			ofstream & fout, chrono::steady_clock::time_point & init_time){
		sort_deg(cand_exts, gs);
		vector<QCVertex>& vertices = g.vertexes;
		int cand_size = cand_exts.size();
		//----
		int gs_size = gs.vertexes.size();
		bool bhas_qclq = false;
		int cover_size = COVER_VERTEX_PRUNE ? cover_prune(gs, g, cand_exts) : 0;
		for(int i = 0; i < cand_size - cover_size; i++){
			if ((gs_size + cand_size - i) < min_size)
				return false;

			//whether the union of subgraph with candidate is QCQ
			if(LOOKAHEAD_PRUNE){
				QCSubgraph union_q = gs;
				for (int k = i; k < cand_size; k++)
					add_vertex(union_q, cand_exts[k]);

				if (is_QC(union_q)) {
					output_qcq(union_q, fout);
					return false;
				}
			}

			QCVertex * v = cand_exts[i];
			QCSubgraph new_gs = gs;
			add_vertex(new_gs, v);

			//construct new candidate

			vector<QCVertex*> new_cand;
			new_cand.insert(new_cand.begin(), cand_exts.begin() + i + 1,
					cand_exts.end());
			filter_by_2hop(new_cand, v, g);
			int new_gs_size = new_gs.vertexes.size();
			int new_cand_size = new_cand.size();
			//-----------
			// 1. leaf node in search space
			if (new_cand_size == 0) {
				if (new_gs_size >= min_size && is_QC(new_gs)) {
					bhas_qclq = true;
					output_qcq(new_gs, fout);
				}
			} else {
				// 2. not leaf node
				bool ext_prune = iterative_bounding(new_cand, new_gs, fout, g);
				new_gs_size = new_gs.vertexes.size();
				new_cand_size = new_cand.size();
				//-----------------------------------
				auto end = chrono::steady_clock::now();
				float exe_time = (float)chrono::duration_cast<chrono::milliseconds>(end - init_time).count()/1000.0;
				if(exe_time > TIME_THRESHOLD){
					//split task if it already run too long
					if(!ext_prune && new_cand_size + new_gs_size >= min_size){
						QCliqueTask * t = new QCliqueTask;
						//New task's subgraph only include vertex in X or candidate
						t->context.gs = new_gs;
						graph_shrink(new_gs, new_cand, g, t->subG,
								t->context.cand_id_vec);
						t->context.round = 3;
						add_task(t);
					}

					//check whether current X is_QC
					if (new_gs_size >= min_size && is_QC(new_gs)){
						bhas_qclq = true;
						output_qcq(new_gs, fout);
					}

				} else {
					//just run QCQ recursively
					if(!ext_prune && new_cand_size + new_gs_size >= min_size){
						bool bhas_super_qclq;
						if((new_cand_size+new_gs_size)*2 < g.vertexes.size()
								&& g.vertexes.size() > COMPACT_SIZE)
						{
							QCSubgraph new_g;
							compact_graph(new_gs, new_cand, g, new_g);
							bhas_super_qclq = tddq_QCQ(new_gs, new_g, new_cand, fout, init_time);
						} else {
							bhas_super_qclq = tddq_QCQ(new_gs, g, new_cand, fout, init_time);
						}

						bhas_qclq = bhas_qclq || bhas_super_qclq;
						if (!bhas_super_qclq && new_gs_size >= min_size && is_QC(new_gs)) {
							bhas_qclq = true;
							output_qcq(new_gs, fout);
						}
					}
				}
			}
		}
		return bhas_qclq;
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
				if(val.nbs.size() < min_deg)
					if(!val.del) prune_1hop(it->first, val, g_map, min_deg, to_del);
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
				if(nbs.size() >= min_deg){
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

			//revise: add 2hop prune
			QCHopSubgraph g_2hop;
			set_2hop(g_map);

			//#9. kcore prune 2_hop g_map
    		vector<VertexID> to_del;
			for(auto it = g_map.begin(); it != g_map.end(); ++it)
			{
				kc_value & val = it->second;
				if(val.nbs.size() < min_deg || val.nbs2hop.size() < min_size - 1)
					if(!val.del) prune_2hop(it->first, val, g_map, min_deg, min_size, to_del);
			}
			//now, g_map is degree-pruned
			//see whether root_id is still not pruned
			if(g_map[rootID].del) return false; //the task terminates
			//do the actual deletes
			for(auto it = to_del.begin(); it != to_del.end(); ++it) g_map.erase(*it);
			if(g_map.size() < min_size) return false;
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

    		//the 2-hop graph is complete now.
			//#11. prepare graph induced by X, and Xcand
			//-- g_map[0] => g(X)
			//-- g_map [1, ..] => Xcand
			//-- add them to the context
    		QCSubgraph& gs = context.gs;
    		QCVertex v;
    		v.id = rootID;
    		gs.addVertex(v);
    		vector<VertexID> & cand_id_vec = context.cand_id_vec;
    		//exclude the root vertex
    		for(int i = 1; i < vertices.size(); i++){
    			assert(vertices[i].value.size() >= min_deg);
				cand_id_vec.push_back(vertices[i].id);
    		}

			context.round++;
			return true;
			//### since there is no vertex requested, step 3 will run directly
			//### so we only need round 3 to call QCQ to simplify code
    	}else{
    		//context.round == 3
			vector<QCVertex*> cand_exts;
			vector<VertexID>& cand_id_vec = context.cand_id_vec;
			int cand_vec_size = cand_id_vec.size();
			for(int i = 0; i < cand_vec_size; i++){
				cand_exts.push_back(g.getVertex(cand_id_vec[i]));
			}

			sort_deg(cand_exts, context.gs);
			auto init_time = chrono::steady_clock::now();

			QCSubgraph new_g;
			compact_graph(context.gs, cand_exts, g, new_g);

    		tddq_QCQ(context.gs, new_g, cand_exts, fout, init_time);
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
    if(argc != 7){
    	cout<<"arg1 = input path in HDFS, arg2 = number of threads"
    			<<", arg3 = degree ratio, arg4 = min_size, arg5 = time delay threshold, arg6 = BIGTASK_THRESHOLD"<<endl;
    	return -1;
    }
    param.input_path = argv[1];  //input path in HDFS
    int thread_num = atoi(argv[2]);  //number of threads per process
    _gamma = atof(argv[3]);
    min_size = atoi(argv[4]);
    TIME_THRESHOLD = atof(argv[5]);
    BIGTASK_THRESHOLD = atoi(argv[6]);
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
