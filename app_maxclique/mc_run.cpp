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

#include "maxclique.h"

int VNUM_ALLOWED_BEFORE_SPLIT;//### note for split ###: edge number is counted once for (u, v) and (v, u)

typedef Task<CliqueVertex, CliqueValue> CliqueTask;
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

	//check whether task is bigtask
	virtual bool is_bigtask(CliqueTask * task){
		if(task->subG.vertexes.size() > BIGTASK_THRESHOLD
				|| task->to_pull.size() > BIGTASK_THRESHOLD)
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
		t->context.push_back(v->id); //====> this is Q = {v}
		for(int i=0; i<v->value.size(); i++)
		{
			VertexID nb = v->value[i];
			t->pull(nb);
		}
		bool result = is_bigtask(t);
		add_task(t);
		return result;
	}

    virtual bool compute(SubgraphT & g, ContextT & context, vector<VertexT *> & frontier)
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
        CliqueAgg* agg = get_aggregator();
		VSet Qmax;
		agg->finishPartial(Qmax);//cannot directly use agg->Qmax without rdlock it first
		//------
        if(vertices.size() > VNUM_ALLOWED_BEFORE_SPLIT)
        {//split
        	cout<<_my_rank<<">> split: V = "<<vertices.size()<<endl;//@@@@@@@@@@@@@
        	for(int i=0; i<vertices.size(); i++)
			{
        		CliqueTask * t = new CliqueTask;
				VertexID newRoot = vertices[i].id;
				t->context = context; //inherits old context
				t->context.push_back(newRoot); //### note for split ###: context expands for one more level
				for(int j=0; j<vertices[i].value.size(); j++)
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
					for (int k = 0; k < fval.size(); k++) {
						CliqueVertex *v1 = t->subG.getVertex(fval[k]);
						if (v1 != NULL) {
							v_nb->value.push_back(fval[k]);
							v1->value.push_back(nb);
						}
					}
				}
				if(t->context.size() + t->subG.vertexes.size() > Qmax.size())
				{
					 add_task(t); //### note for split ###: add only if bigger than Qmax
				}
				else delete t; //if you do not give the system to manage (and delete), you need to delete yourself to avoid memory leak
			}
        }
        else
        {
        	//run single-threaded mining code
			//--- init Qmax as current max for best pruning: max {prev_max, new_max_to_sync}
			if(vertices.size() + context.size() <= Qmax.size()) return false;//===========================> add a pruning
			//====> below: to construct an arbitrary set of size |Qmax| - 1, for pruning purpose when calling MCQ(.):
			int remain = Qmax.size() - context.size();
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
				MCQ(g, Qmax_dummy);
			//cout<<"Qmax_dummy.size() = "<<Qmax_dummy.size()<<endl;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			if(Qmax_dummy.size() > old_size)
			{
				Qmax_dummy.insert(context.begin(), context.end()); //root unions maxclique in subgraph g
				agg->aggregate(Qmax_dummy);
			}
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

    virtual void task_spawn(VertexT * v, vector<TaskT*> & tcollector)
	{
    	CliqueAgg* agg = get_aggregator();
		VSet Qmax;
		agg->finishPartial(Qmax);//cannot directly use agg->Qmax without rdlock it first
		if(Qmax.size() >= 1 + v->value.size()) return; //==========> pruning with Qmax right at spawning
		//------
		TaskT* task = new TaskT;
		task->context.push_back(v->id); //====> this is Q = {v}
		for(int i=0; i<v->value.size(); i++)
		{
			VertexID nb = v->value[i];
			task->pull(nb);
		}
		tcollector.push_back(task);
	}
};

int main(int argc, char* argv[])
{
    init_worker(&argc, &argv);
    WorkerParams param;
    if(argc != 5){
		cout<<"arg1 = input path in HDFS, arg2 = number of threads, arg3 = VNUM_ALLOWED_BEFORE_SPLIT, arg4 = BIGTASK_THRESHOLD"<<endl;
		return -1;
	}
    param.input_path = argv[1];  //input path in HDFS
    int thread_num = atoi(argv[2]);  //number of threads per process
    VNUM_ALLOWED_BEFORE_SPLIT = atoi(argv[3]);
    BIGTASK_THRESHOLD = atoi(argv[4]);
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
