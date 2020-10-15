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

//#include "maxclique.h"
#include "subg-dev.h"

typedef vector<VertexID> CliqueValue;
typedef Vertex<VertexID, CliqueValue> CliqueVertex;
typedef Subgraph<CliqueVertex> CliqueSubgraph;
typedef hash_set<VertexID> VSet;

struct ContextValue
{
	int round; // the iteration number in compute(.)
	CliqueValue context;
	vector<CliqueVertex> vvec;
	vector<int> N;
	int stop_loc; // stop location in expand() split
	bool to_split; // weather need to do the real split in 1st round
	// split task from vertices in [start, end) in 1st round.
	int start; // split start index
	int end; // split end index + 1
};



obinstream & operator>>(obinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.context;
    m >> c.vvec;
    m >> c.N;
    m >> c.stop_loc;
    m >> c.to_split;
    m >> c.start;
    m >> c.end;
    return m;
}

ibinstream & operator<<(ibinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.context;
    m << c.vvec;
    m << c.N;
    m << c.stop_loc;
    m << c.to_split;
    m << c.start;
    m << c.end;
    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
	m >> c.round;
	m >> c.context;
	m >> c.vvec;
	m >> c.N;
	m >> c.stop_loc;
	m >> c.to_split;
	m >> c.start;
	m >> c.end;
	return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
	m << c.round;
	m << c.context;
	m << c.vvec;
	m << c.N;
	m << c.stop_loc;
	m << c.to_split;
	m << c.start;
	m << c.end;
    return m;
}

int VNUM_ALLOWED_BEFORE_SPLIT; //### note for split ###: edge number is counted once for (u, v) and (v, u)
int SPLIT_SIZE;
float TIME_THRESHOLD; // if running time >= TIME_THRESHOLD, split into subtasks

typedef Task<CliqueVertex, ContextValue> CliqueTask;
//### note for split ###: context is the set of vertices already included in Q

void print_vec(vector<VertexID> & vec){ //@@@@@@@ for debug report
    for(int i = 0; i < vec.size(); i++)
    cout << vec[i] << "  ";
    cout << endl;
}

class MCQTrimmer:public Trimmer<CliqueVertex>
{
    virtual void trim(CliqueVertex & v) {
        CliqueValue & val = v.value;
        CliqueValue newval;
        for (int i = 0; i < val.size(); i++) {
            if (v.id < val[i])
            	newval.push_back(val[i]);
        }
        val.swap(newval);
    }
};

class CliqueAgg:public Aggregator<VSet, VSet, VSet>  //Context = int
{
private:
	VSet Q_max;

public:

    virtual void init(){}

    virtual void init_udf(VSet& prev)
    {
    	if(prev.size() > Q_max.size())
    		Q_max = prev; //do not swap, so that prev can be used for reporting any time
    		//it is called infrequently anyway, i.e. after agg_sync
    }

    virtual void aggregate_udf(VSet & Qmax)
    {
    	if(Qmax.size() > Q_max.size())
    		Q_max.swap(Qmax);
    }

    virtual void stepFinal_udf(VSet & part)
    {
    	if(part.size() > Q_max.size())
    		Q_max.swap(part);
    }

    virtual void finishPartial_udf(VSet & collector)
    {
    	collector = Q_max;
    }

    virtual void finishFinal_udf(VSet & collector)
    {
    	collector = Q_max;
    }
};

class CliqueComper:public Comper<CliqueTask, CliqueAgg>
{
public:

	//define comparator for sorting
	struct deg_comparator{
	    bool operator()(const CliqueVertex * u, const CliqueVertex * v) const
	    {
	        return u->value.size() > v->value.size();
	    }
	}; //for sorting by degree (descending order)

	void degree_sort(vector<CliqueVertex> & vertexes, vector<CliqueVertex> & vout, vector<int> & N)
	{
	    //in: vertexes
	    //out: N
	    //side effect: "vertexes" is sorted by degree
	    vector<CliqueVertex *> vvec(vertexes.size());
	    for(int i=0; i<vvec.size(); i++) vvec[i] = &vertexes[i];
	    //--------------
	    //sort by degree
		deg_comparator comper;
	    sort(vvec.begin(), vvec.end(), comper);
	    //get max-degree (we require "vertexes" not empty)
	    int max_deg = vvec[0]->value.size();
	    int size = vertexes.size();
	    //set vout
	    for(int i=0; i<size; i++)
	    {
	        vout.push_back(*vvec[i]);
	    }
	    //set N-array
	    N.resize(size);
	    if(size > max_deg)
	    {
	        int i = 0;
	        for(; i<max_deg; i++) N[i] = i + 1;
	        for(; i<size; i++) N[i] = max_deg + 1;
	    }
	    else
	        for(int i=0; i<size; i++) N[i] = max_deg + 1;
	}

	//reference paper: an efficient branch-and-bound algorithm for finding a maximum clique
	//NUMBER-SORT function
	void color_sort(vector<CliqueVertex> & vertexes, vector<int> & N)
	{
	    //{NUMBER}
	    int maxno = 1;
	    vector<hash_map<VertexID, CliqueVertex *> > C;
	    hash_map<VertexID, int> Nmap;
	    C.resize(3);
	    for(int cur = 0; cur < vertexes.size(); cur++)
	    {
	        CliqueVertex & p = vertexes[cur];
	        int k = 1;
	        CliqueValue & nbs = p.value;
	        //coloring
	        bool stop = false;
	        while(!stop)
	        {
	            int pos = 0; //position in neighbor list
	            for(; pos < nbs.size(); pos++)
	            {
	                VertexID nb = nbs[pos];
	                if(C[k].find(nb) != C[k].end()) break; //find a neighbor in Ck
	            }
	            if(pos == nbs.size()) stop = true;
	            else k++;
	        }
	        //check k
	        if(k > maxno)
	        {
	            maxno = k;
	            C.resize(k + 2);
	        }
	        Nmap[p.id] = k;
	        C[k][p.id] = &p;
	    }
	    //{SORT}
	    vector<CliqueVertex> temp;
	    for(int i=1; i<=maxno; i++)
	    {
	        hash_map<VertexID, CliqueVertex *> & set = C[i];
	        for(hash_map<VertexID, CliqueVertex *>::iterator it = set.begin(); it!=set.end(); it++)
	        {
	            temp.push_back(*(it->second));
	            N.push_back(Nmap[it->first]);
	        }
	    }
	    temp.swap(vertexes);
	}

	void nbs_prune(vector<CliqueVertex *> & Rp_pts, vector<CliqueVertex> & Rp)
	{
	    VSet Rp_set;
	    for(int i = 0; i<Rp_pts.size(); i++)  Rp_set.insert(Rp_pts[i]->id);
	    for(int i = 0; i<Rp_pts.size(); i++)
	    {
	        Rp.resize(Rp.size() + 1);
	        CliqueVertex & u = Rp.back();
	        u.id = Rp_pts[i]->id;
	        vector<VertexID> & nb_list = Rp_pts[i]->value;
	        for(int j = 0; j < nb_list.size(); j++)
	        {
	            if(Rp_set.find(nb_list[j]) != Rp_set.end())
	                u.value.push_back(nb_list[j]);
	        }
	    }
	}

	//EXPAND
	void expand(vector<CliqueVertex> & vertexes, vector<int> & N, VSet & Q, VSet & Qmax,
			CliqueValue & context, int stop_loc, chrono::steady_clock::time_point & init_time) //vertexes = R in the paper
	{
		while((int)vertexes.size() - stop_loc > SPLIT_SIZE)
		{
			CliqueTask * t = new CliqueTask;
			t->context.round = 2;
			t->context.vvec = vertexes;
			t->context.stop_loc = vertexes.size() - SPLIT_SIZE;
			t->context.N = N;
			t->context.context = context;
			t->context.context.insert(t->context.context.end(), Q.begin(), Q.end());
			add_task(t);
			int rest = vertexes.size() - SPLIT_SIZE;
			vertexes.resize(rest);
			N.resize(rest);
		}

		while(vertexes.size() > stop_loc)
		{
			if(Q.size() + N.back() > Qmax.size())
			{
				CliqueVertex & p = vertexes.back(); //***
				Q.insert(p.id);
				vector<CliqueVertex *> Rp_pts; //new search space
				//compute Rp = vertexes intersects nbs(p)
				CliqueValue & nbs = p.value;
				set<VertexID> nbs_set; //somehow, using VSet would be very slow for a sparse graph, due to hash_set setup time
				for(int pos = 0; pos < nbs.size(); pos++) nbs_set.insert(nbs[pos]);
				for(int i = 0; i < vertexes.size(); i++)
				{
					CliqueVertex & v = vertexes[i];
					if(nbs_set.find(v.id) != nbs_set.end()) Rp_pts.push_back(&v);
				}
				if(!Rp_pts.empty())
				{
					vector<CliqueVertex> Rp;
					nbs_prune(Rp_pts, Rp);
					vector<int> colorN;
					color_sort(Rp, colorN);
					if(Q.size() + colorN.back() > Qmax.size()) // check before expand
					{
						auto end = chrono::steady_clock::now();
						float exe_time = (float)chrono::duration_cast<chrono::milliseconds>(end - init_time).count()/1000.0;
						if(exe_time > TIME_THRESHOLD)
						{ // split tasks
							CliqueTask * t = new CliqueTask;
							t->context.round = 2;
							t->context.N.swap(colorN);
							t->context.context = context;
							t->context.context.insert(t->context.context.end(), Q.begin(), Q.end());
							t->context.vvec.swap(Rp);
							t->context.stop_loc = 0;
							add_task(t);
						} else {
							expand(Rp, colorN, Q, Qmax, context, 0, init_time);
						}
					}
				}
				else if(Q.size() > Qmax.size()) Qmax = Q;
				Q.erase(p.id);
				vertexes.pop_back(); //***
				N.pop_back();
			} else return;
		}
	}

	//check whether task is bigtask
	virtual bool is_bigtask(CliqueTask * task){
		if(task->context.vvec.size() > BIGTASK_THRESHOLD
				||task->to_pull.size() > BIGTASK_THRESHOLD)
			return true;
		else
			return false;
	}

	virtual bool task_spawn(VertexT * v)
	{
		CliqueAgg* agg = get_aggregator();
		VSet Qmax;
		agg->finishPartial(Qmax);//cannot directly use agg->Qmax without rdlock it first
		if(Qmax.size() >= 1 + v->value.size()) return false; //==========> pruning with Qmax right at spawning
		//cout<<v->id<<": in task_spawn"<<endl;//@@@@@@@@@@@@@
		CliqueTask * t = new CliqueTask;
		t->context.round = 1;
		t->context.stop_loc = 0;
		t->context.context.push_back(v->id); //====> this is Q = {v}
		for(int i=0; i<v->value.size(); i++)
		{
			VertexID nb = v->value[i];
			t->pull(nb);
		}
		bool result = is_bigtask(t);
		t->context.start = 0;
		t->context.end = t->to_pull.size();
		t->context.to_split = false;
		add_task(t);
		return result;
	}

    virtual bool compute(SubgraphT & g, ContextT & context, vector<VertexT *> & frontier)
    {
    	CliqueAgg* agg = get_aggregator();
		VSet Qmax;
		agg->finishPartial(Qmax);

    	if(context.round == 1)
		{
			//cout<<"context: ";//@@@@@@@@@@@@@
			//print_vec(context);//@@@@@@@@@@@@@
			for(int i = 0; i < frontier.size(); i++) {
				CliqueVertex v;
				v.id = frontier[i]->id;
				g.addVertex(v);
			}
			//-----
			//int E = 0;
			for(int i = 0; i < frontier.size(); i++) {
				CliqueVertex *v = g.getVertex(frontier[i]->id);
				CliqueValue &fval = frontier[i]->value;
				for (int j = 0; j < fval.size(); j++) {
					CliqueVertex *v1 = g.getVertex(fval[j]);
					if (v1 != NULL) {
						v->value.push_back(fval[j]);
						v1->value.push_back(v->id);
						//E++;
					}
				}
			}
			//====> now "g" contains nbs(Q), i.e. intersection of all nbs(v), v is in Q
			/*
			//@@@@@@ report graph g @@@@@@
			cout<<"********** g"<<context<<endl;
			for(int i=0; i<g.vertexes.size(); i++)
			{
				VertexT & v = g.vertexes[i];
				cout<<"v"<<v.id<<": ";
				for(int j=0; j<v.value.size(); j++) cout<<v.value[j]<<" ";
				cout<<endl;
			}
			//@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			*/
			vector<CliqueVertex> & vertices = g.vertexes;
			//------
			/*This method is not faster than sequential program.
			 * It will also easy to use up the memory when expand parameter go low
			 */
			if (context.end - context.start > VNUM_ALLOWED_BEFORE_SPLIT)
			{
				context.to_split = true;
				while(context.end - context.start > SPLIT_SIZE)
				{
					CliqueTask * t = new CliqueTask;
					t->context.context = context.context;
					int new_end = context.start + SPLIT_SIZE;
					t->context.start = context.start;
					t->context.end = new_end;
					t->context.round = 1;
					t->context.to_split = true;
					t->subG = g;
					add_task(t);
					context.start = new_end;
				}

			}

			if(context.to_split)
			{//split
				cout<<_my_rank<<">> split: V = "<<vertices.size()<<endl;//@@@@@@@@@@@@@
				int i = context.start;
				for(; i<context.end; i++)
				{
					CliqueTask * t = new CliqueTask;
					VertexID newRoot = vertices[i].id;
					t->context.context = context.context; //inherits old context
					t->context.context.push_back(newRoot); //### note for split ###: context expands for one more level
					for(int j=0; j<vertices[i].value.size(); j++)// add newRoot's nbs into subG
					{
						if (vertices[i].value[j] < newRoot) {
							continue;
						}

						CliqueVertex v;
						v.id = vertices[i].value[j];
						t->subG.addVertex(v);
					}
					for(int j=0; j<vertices[i].value.size(); j++)
					{
						VertexID nb = vertices[i].value[j];

						if (nb  < newRoot) {
							continue;
						}

						CliqueVertex* v_nb = t->subG.getVertex(nb);
						CliqueValue &fval = g.getVertex(nb)->value;
						for (int k = 0; k < fval.size(); k++) { // add newRoot's nbs's edge into subG
							CliqueVertex *v1 = t->subG.getVertex(fval[k]);
							if (v1 != NULL) {
								v_nb->value.push_back(fval[k]);
								v1->value.push_back(nb);
							}
						}
					}
					if(t->context.context.size() + t->subG.vertexes.size() > Qmax.size())
					{
						t->context.round = 1; // continue split with "vertices.size() > VNUM_ALLOWED_BEFORE_SPLIT"
						t->context.stop_loc = 0;
						t->context.start = 0;
						t->context.end = t->subG.vertexes.size();
						t->context.to_split = false;
						add_task(t); //### note for split ###: add only if bigger than Qmax
					}
					else delete t; //if you do not give the system to manage (and delete), you need to delete yourself to avoid memory leak
				}
				return false;
			}
			else if(!vertices.empty())
			{
				degree_sort(vertices, context.vvec, context.N);
			}
		}
    	//round == 2
		vector<CliqueVertex> & vertices = context.vvec;
		//run single-threaded mining code
		//--- init Qmax as current max for best pruning: max {prev_max, new_max_to_sync}
		if(vertices.size() + context.context.size() <= Qmax.size()) return false;//===========================> add a pruning
		//====> below: to construct an arbitrary set of size |Qmax| - 1, for pruning purpose when calling MCQ(.):
		int remain = Qmax.size() - context.context.size();
		VSet Qmax_dummy;
		if(remain > 0) //just to make a set with "remain" vertices
		{
			for(auto it = Qmax.begin(); it != Qmax.end(); it++)
			{
				Qmax_dummy.insert(*it);
				if(Qmax_dummy.size() == remain) break;
			}
		}
		size_t old_size = Qmax_dummy.size();
		if(vertices.size() > 0) //====> note that "g" does not contain "root" (stored in context), can be empty
		{
			auto init_time = chrono::steady_clock::now();
			VSet Q;
			expand(context.vvec, context.N, Q, Qmax_dummy, context.context, context.stop_loc, init_time);
		}
		if(Qmax_dummy.size() > old_size || context.context.size() > Qmax.size())
		{
			Qmax_dummy.insert(context.context.begin(), context.context.end()); //root unions maxclique in subgraph g
			CliqueAgg* agg = get_aggregator();
			agg->aggregate(Qmax_dummy);
		}
		return false;
    }
};

class CliqueWorker:public Worker<CliqueComper>
{
public:
    CliqueWorker(int num_compers) : Worker<CliqueComper>(num_compers){}

    virtual VertexT* toVertex(char* line)
    {
        VertexT* v = new VertexT;
        char * pch;
        pch=strtok(line, " \t");
        v->id=atoi(pch);
        strtok(NULL," \t");
        CliqueValue & val = v->value;
        while((pch=strtok(NULL, " ")) != NULL)
        {
            val.push_back(atoi(pch));
        }
        return v;
    }

    virtual void task_spawn(VertexT * v, vector<TaskT> & tcollector)
	{
    	CliqueAgg* agg = get_aggregator();
		VSet Qmax;
		agg->finishPartial(Qmax);//cannot directly use agg->Qmax without rdlock it first
		if(Qmax.size() >= 1 + v->value.size()) return; //==========> pruning with Qmax right at spawning
		//------
		TaskT t;
		tcollector.push_back(t);
		TaskT & task = tcollector.back();
		task.context.round = 1;
		task.context.stop_loc = 0;
		task.context.context.push_back(v->id); //====> this is Q = {v}
		for(int i=0; i<v->value.size(); i++)
		{
			VertexID nb = v->value[i];
			task.pull(nb);
		}
		task.context.start = 0;
		task.context.end = task.to_pull.size();
		task.context.to_split = false;
	}
};

int main(int argc, char* argv[])
{
    init_worker(&argc, &argv);
    WorkerParams param;
    if(argc != 7){
    	cout<<"arg1 = input path in HDFS, arg2 = number of threads, arg3 = split-time threshold, arg4 = BIGTASK_THRESHOLD,"
    			<<"arg5 = SPLIT_SIZE, arg6 = VNUM_ALLOWED_BEFORE_SPLIT"<<endl;
    	return -1;
    }
    param.input_path = argv[1];  //input path in HDFS
    int thread_num = atoi(argv[2]);  //number of threads per process
    TIME_THRESHOLD = atof(argv[3]);  //split-time threshold
    BIGTASK_THRESHOLD = atoi(argv[4]);
    SPLIT_SIZE = atoi(argv[5]);
    VNUM_ALLOWED_BEFORE_SPLIT = atoi(argv[6]);
    param.force_write=true;
    param.native_dispatcher=false;
    //------
    MCQTrimmer trimmer;
    CliqueAgg aggregator;
    CliqueWorker worker(thread_num);
    worker.setTrimmer(&trimmer);
    worker.setAggregator(&aggregator);
    worker.run(param);
    //report ------
    //job finished, no need to lock "agg_rwlock"
    if(_my_rank == MASTER_RANK)
    {
    	VSet * agg = (VSet *)global_agg;
    	cout<<"Qmax: ";
    	for(auto it = agg->begin(); it != agg->end(); it++) cout<<*it<<" ";
    	cout<<endl;
    	cout<<"|Qmax| = "<<agg->size()<<endl;
    }
    //-------------
    worker_finalize();
    return 0;
}
