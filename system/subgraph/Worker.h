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

#ifndef WORKER_H_
#define WORKER_H_

#include <iostream>

#include "util/global.h"
#include "util/ydhdfs.h"
#include "util/communication.h"
#include "Trimmer.h"
#include "adjCache.h"
#include "ReqServer.h"
#include "RespServer.h"
#include "Comper.h"
#include "GC.h"
#include "AggSync.h"
#include "Profiler.h"
#include <unistd.h> //sleep(sec)
#include <queue> //for std::priority_queue
#include <numeric> //for std::accumulate

using namespace std;

struct steal_plan
{
	int id; // steal from/to which server. negative number means receive
	int num; // steal batch size
	steal_plan(){}
	steal_plan(int id, int num):id(id), num(num){}
};

obinstream & operator>>(obinstream & m, steal_plan & s)
{
	m >> s.id;
	m >> s.num;
	return m;
}

ibinstream & operator<<(ibinstream & m, const steal_plan & s)
{
	m << s.id;
	m << s.num;
	return m;
}

ofbinstream & operator>>(ofbinstream & m, steal_plan & s)
{
	m >> s.id;
	m >> s.num;
	return m;
}

ifbinstream & operator<<(ifbinstream & m, const steal_plan & s)
{
	m << s.id;
	m << s.num;
	return m;
}

template <class Comper>
class Worker
{
public:
	typedef typename Comper::TaskType TaskT;
    typedef typename Comper::AggregatorType AggregatorT;
	typedef typename Comper::TaskMapT TaskMapT;

    typedef typename TaskT::VertexType VertexT;
    typedef typename TaskT::SubgraphT SubgraphT;
    typedef typename TaskT::ContextType ContextT;

    typedef typename VertexT::KeyType KeyT;
	typedef typename VertexT::ValueType ValueT;
	typedef typename VertexT::HashType HashT;

    typedef vector<VertexT*> VertexVec;

    typedef hash_map<KeyT, VertexT*> VTable;
    typedef typename VTable::iterator TableIter;

    typedef AdjCache<TaskT> CTable;

    typedef typename AggregatorT::PartialType PartialT;
    typedef typename AggregatorT::FinalType FinalT;

    typedef GC<TaskT> GCT;

    typedef Trimmer<VertexT> TrimmerT;

    typedef deque<TaskT*> TaskQueue;

    //=======================================================
    //worker's data structures
    HashT hash;
    VTable local_table; //key-value store of my vertex portion
    CTable * cache_table; //cached remote vertices, it creates ReqQueue for appeding reqs

	VertexVec vertexes;
	typedef stack<VertexT*> VertexStack;
	VertexStack vertex_stack;

    TaskMapT** taskmap_vec;

    bool local_idle; //indicate whether the current worker is idle
    //logic with work stealing:
    //if it tries to steal from all workers but failed to get any job, this should be set
    //after master's idle-condition sync, if job does not terminate, should steal for another round
    //ToDo："local_idle" may need to change to "one per thread", or aggregated over threads

    Comper* compers; //dynamic array of compers

    TaskQueue& q_bigtask(){// get the global big task queue
    	return *(TaskQueue *)big_task_queue;
    }

    TaskMapT& get_big_maptask(){// get the global big map_task
        	return *(TaskMapT *)big_map_task;
    }
    //=======================================================
    //Trimmer
    void setTrimmer(TrimmerT* trimmer)
	{
    	global_trimmer = trimmer;
	}

    //=======================================================
    //constructor & destructor

    Worker(int comper_num, string local_disk_path = "/opt/hadoop/dfs/yan_share/guimu/buffered_tasks", string report_path = "/opt/hadoop/dfs/yan_share/guimu/report")
    {
    	num_compers = comper_num;
    	TASK_DISK_BUFFER_DIR = local_disk_path;
    	_mkdir(TASK_DISK_BUFFER_DIR.c_str());
    	REPORT_DIR = report_path;
    	_mkdir(REPORT_DIR.c_str());
    	//------
    	//create global big task queue, table and buffer
    	big_task_queue = new TaskQueue;
    	big_map_task = new TaskMapT;
    	get_big_maptask().thread_rank = num_compers;
    	get_big_maptask().need_lock = true;

    	global_end_label = false;
    	local_idle = false;
    	global_trimmer = NULL;
    	global_aggregator = NULL;
    	global_agg = NULL;
    	req_counter = new atomic<size_t>[_num_workers];
    	for(int i=0; i<_num_workers; i++) req_counter[i] = 0; //must be before the next line
    	global_vcache = cache_table = new CTable;
    	global_local_table = &local_table;
    	global_vertexes_stack = &vertex_stack;
		idle_set = new atomic<bool>[comper_num];
		for(int i=0; i<comper_num; i++) idle_set[i] = false;
    }

    void setAggregator(AggregatorT* ag)
    {
        global_aggregator = ag;
        global_agg = new FinalT;
        ag -> init();
    }

    AggregatorT* get_aggregator() //get aggregator
    //cannot use the same name as in global.h (will be understood as the local one, recursive definition)
	{
		return (AggregatorT*)global_aggregator;
	}

	virtual ~Worker()
	{
        for(int i=0;i<vertexes.size();i++){
            delete vertexes[i];
        }
		delete[] compers;
		delete[] taskmap_vec;
		delete[] global_tasknum_vec;
		delete[] idle_set;
		delete[] req_counter;
		delete cache_table;
		delete (TaskQueue *)big_task_queue;
		delete (TaskMapT *)big_map_task;
		//ToDo: release aggregator
        if (global_agg != NULL)
            delete (FinalT*)global_agg;
	}

	//=======================================================
	//graph loading:

	//user-defined loading function
	virtual VertexT* toVertex(char* line) = 0;

	void load_graph(const char* inpath, VertexVec & vVec)
	{
		TrimmerT* trimmer = NULL;
		if(global_trimmer != NULL) trimmer = (TrimmerT*)global_trimmer;
		//------
		hdfsFS fs = getHdfsFS();
		hdfsFile in = getRHandle(inpath, fs);
		LineReader reader(fs, in);
		while(true)
		{
			reader.readLine();
			if (!reader.eof())
			{
				VertexT * v = toVertex(reader.getLine());
				if(trimmer) trimmer->trim(*v);
				vVec.push_back(v);
			}
			else
				break;
		}
		hdfsCloseFile(fs, in);
		hdfsDisconnect(fs);
	}

	void sync_graph(VertexVec & vVec)
	{
		//ResetTimer(4);
		//set send buffer
		vector<VertexVec> _loaded_parts(_num_workers);
		for (int i = 0; i < vVec.size(); i++) {
			VertexT* v = vVec[i];
			_loaded_parts[hash(v->id)].push_back(v);
		}
		//exchange vertices to add
		all_to_all(_loaded_parts, GRAPH_LOAD_CHANNEL);

		vVec.clear();
		//collect vertices to add
		for (int i = 0; i < _num_workers; i++) {
			vVec.insert(vVec.end(), _loaded_parts[i].begin(), _loaded_parts[i].end());
		}
		_loaded_parts.clear();
		//StopTimer(4);
		//PrintTimer("Reduce Time",4);
	};

	void set_local_table(VertexVec & vVec)
	{
		for(int i=0; i<vVec.size(); i++)
		{
			VertexT * v = vVec[i];
			local_table[v->id] = v;
		}
	}

	//=======================================================
	void create_compers()
	{
		compers = new Comper[num_compers];
		//set global_taskmap_vec
		taskmap_vec = new TaskMapT*[num_compers];
		global_tasknum_vec = new atomic<size_t>[num_compers+1];
		global_tasknum_vec[num_compers] = 0;
		global_taskmap_vec = taskmap_vec;
		for(int i=0; i<num_compers; i++)
		{
			compers[i].map_task.need_lock = false;
			taskmap_vec[i] = &(compers[i].map_task);
			global_tasknum_vec[i] = 0;
			compers[i].start(i);
		}
	}

	//called by the main worker thread, to sync computation-progress, and aggregator
	void status_sync(bool sth2steal)
	{
		bool worker_idle = (sth2steal == false) && (global_num_idle.load(memory_order_relaxed) == num_compers);
		if(_my_rank != MASTER_RANK)
		{
			send_data(worker_idle, MASTER_RANK, STATUS_CHANNEL);
			bool all_idle = recv_data<bool>(MASTER_RANK, STATUS_CHANNEL);
			if(all_idle) global_end_label = true;
		}
		else
		{
			bool all_idle = worker_idle;
			for(int i=0; i<_num_workers; i++)
			{
				if(i != MASTER_RANK) all_idle = (recv_data<bool>(i, STATUS_CHANNEL) && all_idle);
			}
			if(all_idle) global_end_label = true;
			for(int i=0; i<_num_workers; i++)
			{
				if(i != MASTER_RANK) send_data(all_idle, i, STATUS_CHANNEL);
			}
		}
	}

	//=======================================================
	//task stealing
	size_t get_remaining_task_num()
	//not counting number of active tasks in memory (for simplicity)
	{
		bigtask_que_lock.lock();
		int q_bigtask_size = q_bigtask().size();
		bigtask_que_lock.unlock();
		return q_bigtask_size + global_bigTask_file_num * BIG_TASK_FLUSH_BATCH;
	}


	struct max_heap_entry
	{
		size_t num_remain;
		int rank;

		bool operator<(const max_heap_entry& o) const
		{
			return num_remain < o.num_remain;
		}
	};

	struct min_heap_entry
	{
		size_t num_remain;
		int rank;

		bool operator<(const min_heap_entry& o) const
		{
			return num_remain > o.num_remain;
		}
	};

	//UDF for stea1ing seed tasks
	virtual void task_spawn(VertexT * v, vector<TaskT*> & tvec) = 0;

	/*//=== deprecated, 50 vertices may just spawn 0 task or 2 tasks (etc.), so the quota of 50 is wasted during plan generation
	//get tasks from local-table
	//returns false if local-table is exhausted
	bool locTable2vec(vector<TaskT> & tvec)
	{
		size_t begin, end; //[begin, end) are the assigned vertices (their positions in local-table)
		//note that "end" is exclusive
		int size = local_table.size();
		//======== critical section on "global_vertex_pos"
		global_vertex_pos_lock.lock();
		if(global_vertex_pos < size)
		{
			begin = global_vertex_pos; //starting element
			end = begin + TASK_BATCH_NUM;
			if(end > size) end = size;
			global_vertex_pos = end; //next position to spawn
		}
		else begin = -1; //meaning that local-table is exhausted
		global_vertex_pos_lock.unlock();
		//======== spawn tasks from local-table[begin, end)
		if(begin == -1) return false;
		else
		{
			VertexVec & gb_vertexes = *(VertexVec*) global_vertexes;
			for(int i=begin; i<end; i++)
			{//call UDF to spawn tasks
				task_spawn(gb_vertexes[i], tvec);
			}
			return true;
		}
	}
	*/

	//get tasks from local-table
	//returns false if local-table is exhausted
	/*bool locTable2vec(vector<TaskT> & tvec)
	{
		size_t begin, end; //[begin, end) are the assigned vertices (their positions in local-table)
		//note that "end" is exclusive
		int size = local_table.size();
		//======== critical section on "global_vertex_pos"
		while(tvec.size() < TASK_BATCH_NUM)
		{
			global_vertex_pos_lock.lock();
			if(global_vertex_pos < size)
			{
				begin = global_vertex_pos; //starting element
				end = begin + MINI_BATCH_NUM;
				if(end > size) end = size;
				global_vertex_pos = end; //next position to spawn
			}
			else begin = -1; //meaning that local-table is exhausted
			global_vertex_pos_lock.unlock();
			//======== spawn tasks from local-table[begin, end)
			if(begin == -1) return false;
			else
			{
				VertexVec & gb_vertexes = *(VertexVec*) global_vertexes;
				for(int i=begin; i<end; i++)
				{//call UDF to spawn tasks
					task_spawn(gb_vertexes[i], tvec);
				}
			}
		}
		return true;
	}*/
	//temporary for bigtask revision
//	bool locTable2vec(vector<TaskT> & tvec){
//		return false;
//	}

	//get big tasks from disk files
	//returns false if "global_bigTask_fileList" is empty
	bool bigTask_file2vec(vector<TaskT*> & tvec)
	{
		string file;
		bool succ = global_bigTask_fileList.dequeue(file);
		//@@@@@@@@@@@@@@@@@
		//if(succ) cout<<"!!!!!!!!!!! there is big file !!!!!!!!"<<endl;
		if(!succ) return false; //"global_file_list" is empty
		else
		{
			global_bigTask_file_num --;
			ofbinstream in(file.c_str());
			while(!in.eof())
			{
				TaskT * task = new TaskT;
				in >> *task;
				tvec.push_back(task);
			}
			in.close();
			//------
			if (remove(file.c_str()) != 0) {
				cout<<"Error removing file: "<<file<<endl;
				perror("Error printed by perror");
			}
			return true;
		}
	}

	//=== for handling task streaming on disk ===
	char fname[1000], num[20];
	long long fileSeqNo = 1;
	void set_fname() //will proceed file seq #
	{
		strcpy(fname, TASK_DISK_BUFFER_DIR.c_str());
		sprintf(num, "/%d_", _my_rank);
		strcat(fname, num);
		sprintf(num, "%d_", num_compers); //compers have rank 0, 1, ... comper_num-1; so there's no conflict
		strcat(fname, num);
		sprintf(num, "%lld", fileSeqNo);
		strcat(fname, num);
		fileSeqNo++;
	}

	//set name for big task file
	long long bigFileSeqNo = 1;
	void set_bigTask_fname()
	{
		strcpy(fname, TASK_DISK_BUFFER_DIR.c_str());
		sprintf(num, "/bt_%d_", _my_rank);
		strcat(fname, num);
		sprintf(num, "%d_", num_compers);
		strcat(fname, num);
		sprintf(num, "%lld", bigFileSeqNo);
		strcat(fname, num);
		bigFileSeqNo++;
	}

	void add_bigTask(TaskT * task)
	{
		unique_lock<mutex> lck(bigtask_que_lock);
		TaskQueue& btq = q_bigtask();
		//get the ref of global big task queue
		if(btq.size() == BIG_TASK_QUEUE_CAPACITY){
			//@@@@@
			//cout<<"!!!!!!save the file"<<endl;
			set_bigTask_fname();
			ifbinstream bigTask_out(fname);
			int i = 0;
			while(i < BIG_TASK_FLUSH_BATCH)
			{
				//get task at the tail
				TaskT * t = btq.back();
				btq.pop_back();
				//stream to file
				bigTask_out << t;
				//release from memory
				delete t;
				i++;
			}
			bigTask_out.close();
			global_bigTask_fileList.enqueue(fname);
			global_bigTask_file_num ++;
		}
		btq.push_back(task);
	}

	//load bigtasks from a file (from "global_bigTask_fileList" to the task queue)
	//returns false if "global_bigTask_fileList" is empty
	bool file2bigTask_queue()
	{
		string file;
		bool succ = global_bigTask_fileList.dequeue(file);
		if(!succ) return false; //"global_bigTask_fileList" is empty
		else
		{
			global_bigTask_file_num --;
			TaskQueue& btq = q_bigtask();
			ofbinstream in(file.c_str());
			while(!in.eof())
			{
				TaskT* task;
				in >> task;
				btq.push_back(task);
			}
			in.close();

			if (remove(file.c_str()) != 0) {
				cout<<"Error removing file: "<<file<<endl;
				perror("Error printed by perror");
			}
			return true;
		}
	}

	bool steal_planning() //whether there's something to steal from/to others
	{
		//these are assigned by master, indicating how to steal
		vector<steal_plan> my_single_steal_list; // if my_single_steal_list[i] = 3, I need to steal a tasks from Worker 3's global queue
		vector<int> my_batch_steal_list; // if my_batch_steal_list[i] = 3, I need to steal a task-file from Worker 3's global file list
		//====== set my_steal_list
		if(_my_rank != MASTER_RANK)
		{
			send_data(get_remaining_task_num(), MASTER_RANK, STATUS_CHANNEL);
			recv_data<vector<steal_plan> >(MASTER_RANK, STATUS_CHANNEL, my_single_steal_list);
			recv_data<vector<int> >(MASTER_RANK, STATUS_CHANNEL, my_batch_steal_list);
		}
		else
		{
			//collect remaining workloads
			vector<size_t> remain_vec(_num_workers);
			for(int i=0; i<_num_workers; i++)
			{
				if(i != MASTER_RANK)
					remain_vec[i] = recv_data<size_t>(i, STATUS_CHANNEL);
				else
					remain_vec[i] = get_remaining_task_num();
			}
			//------
			priority_queue<max_heap_entry> max_heap;
			priority_queue<min_heap_entry> min_heap;
			avg_num = ceil((float)accumulate(remain_vec.begin(), remain_vec.end(), 0)/remain_vec.size());
			for(int i=0; i<_num_workers; i++)
			{
				if(remain_vec[i] > avg_num)
				{
					max_heap_entry en;
					en.num_remain = remain_vec[i];
					en.rank = i;
					max_heap.push(en);
				}
				else if(remain_vec[i] < avg_num)
				{
					min_heap_entry en;
					en.num_remain = remain_vec[i];
					en.rank = i;
					min_heap.push(en);
				}
			}
			//------
			//plan generation
			//to track stolen task number, each element should not exceed MAX_STEAL_TASK_NUM
			vector<int> steal_num(_num_workers, 0); // steal_num[i] = Worker i currently will steal how many tasks
			vector<vector<steal_plan> > single_steal_lists(_num_workers); //steal_list[i] = {j{id,num}, ...} means Worker_i will steal a j.num tasks.. from Worker j.id, ...
			vector<vector<int> > batch_steal_lists(_num_workers); //steal_list[i] = {j, k, ...} means Worker_i will steal a task-file from each Worker j, k, ...
			int total_plans_num = 0;
			while(!max_heap.empty() && !min_heap.empty())
			{
				max_heap_entry max = max_heap.top();
				max_heap.pop();
				min_heap_entry min = min_heap.top();
				min_heap.pop();
				if(avg_num - min.num_remain > BIG_TASK_FLUSH_BATCH
						&& max.num_remain - avg_num > BIG_TASK_FLUSH_BATCH)
				{// both has a gap >= task-batchsize, steal file
					max.num_remain -= BIG_TASK_FLUSH_BATCH;
					min.num_remain += BIG_TASK_FLUSH_BATCH;
					steal_num[min.rank] += BIG_TASK_FLUSH_BATCH;
					total_plans_num += BIG_TASK_FLUSH_BATCH; // only for printing

					//a negative tag (-x-1) means receiving
					batch_steal_lists[min.rank].push_back(-max.rank-1);
					batch_steal_lists[max.rank].push_back(min.rank);
					//---
					if(max.num_remain > avg_num) max_heap.push(max);
					if(steal_num[min.rank] < MAX_STEAL_TASK_NUM &&
							min.num_remain < avg_num)
						min_heap.push(min);
				}
				else
				{//steal n task; n < BIG_TASK_FLUSH_BATCH
					int steal_batch = std::min((max.num_remain - avg_num), (avg_num - min.num_remain));
					max.num_remain -= steal_batch;
					min.num_remain += steal_batch;
					steal_num[min.rank] += steal_batch;
					total_plans_num += steal_batch;

					single_steal_lists[min.rank].push_back(steal_plan(-max.rank-1,steal_batch));
					single_steal_lists[max.rank].push_back(steal_plan(min.rank,steal_batch));
					//---
					if(max.num_remain > avg_num) max_heap.push(max);
					if(steal_num[min.rank] < MAX_STEAL_TASK_NUM &&
							min.num_remain < avg_num)
						min_heap.push(min);
				}
			}
			//------
			if(total_plans_num > 0) cout<<total_plans_num<<" stealing plans generated at the master"<<endl;//@@@@@@
			//------
			//distribute the plans to machines
			for(int i=0; i<_num_workers; i++)
			{
				if(i == _my_rank)
				{
					single_steal_lists[i].swap(my_single_steal_list);
					batch_steal_lists[i].swap(my_batch_steal_list);
				}
				else
				{
					send_data(single_steal_lists[i], i, STATUS_CHANNEL);
					send_data(batch_steal_lists[i], i, STATUS_CHANNEL);
				}
			}
		}
		//====== execute my_steal_list
		//return false if there is no stealing plan for this worker
		if(my_single_steal_list.empty() && my_batch_steal_list.empty()) return false;

		//steal tasks from bigTask files
		/*if(!my_batch_steal_list.empty())
		{
			for(int i=0; i<my_batch_steal_list.size(); i++)
			{
				int other = my_batch_steal_list[i];
				if(other < 0)
				{
					char* buffer;
					recv_data<char*>(-other-1, STATUS_CHANNEL, buffer);
					if(strlen(buffer) > 0)
					{
						set_bigTask_fname();
						num_stolen += BIG_TASK_FLUSH_BATCH;
						//------
						//register with "global_file_list"
						ofstream outfile(fname, ios::binary);
						outfile.write(buffer,strlen(buffer));
						if(!outfile)
						{
							cerr<<"open error!"<<endl;
							abort();
						}
						global_bigTask_fileList.enqueue(fname);
						global_bigTask_file_num ++;
						outfile.close();
					}
					delete []buffer;
				}
				else
				{
					char* buffer;
					if(get_remaining_task_num() > avg_num)
					//check this since time has passed, and more tasks may have been processed
					//send empty filebuf if no longer a task heavy-hitter
					{
						string file;
						bool succ = global_file_list.dequeue(file);
						if(!succ)
						{
							send_data(buffer, other, STATUS_CHANNEL); //send even if it's empty
						}
						else
						{
							filebuf *pbuf;
							ifstream filestr;
							long size;
							global_bigTask_file_num --;
							filestr.open (file.c_str(), ios::binary);
							pbuf=filestr.rdbuf();
							size=pbuf->pubseekoff (0,ios::end,ios::in);
							pbuf->pubseekpos (0,ios::in);
							buffer=new char[size];
							pbuf->sgetn (buffer,size);
							send_data(buffer, other, STATUS_CHANNEL);
							filestr.close();
							//------
							if (remove(file.c_str()) != 0) {
								cout<<"Error removing file: "<<file<<endl;
								perror("Error printed by perror");
							}
						}
					}
					else
						send_data(buffer, other, STATUS_CHANNEL); //send even if it's empty
					delete []buffer;
				}
			}
		}*/

		if(!my_batch_steal_list.empty())
		{
			for(int i=0; i<my_batch_steal_list.size(); i++)
			{
				int other = my_batch_steal_list[i];
				if(other < 0)
				{
					vector<TaskT> tvec;
					recv_data<vector<TaskT> >(-other-1, STATUS_CHANNEL, tvec);
					if(!tvec.empty())
					{
						set_bigTask_fname();
						ifbinstream out(fname);
						//------
						for(int i=0; i<tvec.size(); i++)
						{
							//@@@@@@@@@@@@@@@@@
							//cout<<"++++++++++ tvec.size()"<<endl;
							out << tvec[i];
						}
						out.close();
						num_stolen += tvec.size();
						//------
						//register with "global_file_list"
						global_bigTask_fileList.enqueue(fname);
						global_bigTask_file_num ++;
					}
				}
				else
				{
					vector<TaskT*> tvec;
					if(get_remaining_task_num() > avg_num)
					//check this since time has passed, and more tasks may have been processed
					//send empty tvec if no longer a task heavy-hitter
						bigTask_file2vec(tvec);
					send_data(tvec, other, STATUS_CHANNEL); //send even if it's empty
					for(int i=0; i<tvec.size(); i++)
						delete tvec[i];
				}
			}
		}


		//steal tasks from bigTask queue
		if(!my_single_steal_list.empty())
		{
			for(int i=0; i<my_single_steal_list.size(); i++)
			{
				int other = my_single_steal_list[i].id;
				if(other < 0)
				{
					vector<TaskT*> tvec;
					recv_data<vector<TaskT*> >(-other-1, STATUS_CHANNEL, tvec);
					for(int i=0; i<tvec.size(); i++)
						add_bigTask(tvec[i]);
					num_stolen += tvec.size();
				}
				else
				{
					int steal_num = my_single_steal_list[i].num;
					vector<TaskT*> tvec;
					for(int i=0; i<steal_num; i++)
					{
						if(get_remaining_task_num() <= avg_num) break;
						//check this since time has passed, and more tasks may have been processed
						//send empty task-vec if no longer a task heavy-hitter

						TaskT * task = NULL;
						TaskQueue& btq = q_bigtask();
						bigtask_que_lock.lock();
						if(btq.size() <= BIG_TASK_FLUSH_BATCH)
							file2bigTask_queue();

						if(!btq.empty())
						{
							task = btq.front();
							btq.pop_front();
						}
						bigtask_que_lock.unlock();
						if(task != NULL) tvec.push_back(task);
					}
					send_data(tvec, other, STATUS_CHANNEL); //send even if it's empty
					for(int i=0; i<tvec.size(); i++)
						delete tvec[i];
				}
//				else
//				{
//					int steal_num = my_single_steal_list[i].num;
//					vector<TaskT> tvec;
//					for(int i=0; i<steal_num; i++)
//					{
//						if(get_remaining_task_num() <= avg_num) break;
//						//check this since time has passed, and more tasks may have been processed
//						//send empty task-vec if no longer a task heavy-hitter
//
//						TaskT * task = NULL;
//						TaskQueue& btq = q_bigtask();
//						bigtask_que_lock.lock();
//						if(btq.size() <= BIG_TASK_FLUSH_BATCH)
//							file2bigTask_queue();
//
//						if(!btq.empty())
//						{
//							task = btq.front();
//							btq.pop_front();
//						}
//						bigtask_que_lock.unlock();
//						if(task != NULL) tvec.push_back(*task);
//					}
//					send_data(tvec, other, STATUS_CHANNEL); //send even if it's empty
//				}
			}
		}

		return true;
	}

	//=======================================================
	//program entry point
    void run(const WorkerParams& params)
    {
        //check path + init
        if (_my_rank == MASTER_RANK)
        {
            if (dirCheck(params.input_path.c_str()) == -1)
                return;
        }
        init_timers();

		//dispatch splits
		ResetTimer(WORKER_TIMER);
		vector<vector<string> >* arrangement;
		if (_my_rank == MASTER_RANK) {
			arrangement = params.native_dispatcher ? dispatchLocality(params.input_path.c_str()) : dispatchRan(params.input_path.c_str());
			//reportAssignment(arrangement);//DEBUG !!!!!!!!!!
			masterScatter(*arrangement);
			vector<string>& assignedSplits = (*arrangement)[0];
			//reading assigned splits (map)
			for (vector<string>::iterator it = assignedSplits.begin();
				 it != assignedSplits.end(); it++)
				load_graph(it->c_str(), vertexes);
			delete arrangement;
		} else {
			vector<string> assignedSplits;
			slaveScatter(assignedSplits);
			//reading assigned splits (map)
			for (vector<string>::iterator it = assignedSplits.begin();
				 it != assignedSplits.end(); it++)
				load_graph(it->c_str(), vertexes);
		}

		//send vertices according to hash_id (reduce)
		sync_graph(vertexes);

		//use "vertexes" to set local_table
		set_local_table(vertexes);

		//barrier for data loading
		worker_barrier();
		StopTimer(WORKER_TIMER);
		PrintTimer("Load Time", WORKER_TIMER);

		//copy vertex* in vertexes to vertex_stack (tail to head)
		int verVec_len = vertexes.size();
		for(int i=verVec_len-1; i>=0; i--){
			vertex_stack.push(vertexes[i]);
		}

		vertexes.clear();

		//ReqQueue already set, by Worker::cache_table
		//>> by this time, ReqQueue occupies about 0.3% CPU

		//set up ReqServer (containing RespQueue), let it know local_table for responding reqs
		ReqServer<VertexT> server_req(local_table);

		//set up computing threads
		create_compers(); //side effect: set global_comper_vec

		//set up RespServer, let it know cache_table so that it can update it when getting resps
		RespServer<Comper> server_resp(*cache_table); //it would read global_comper_vec

		//set up vcache GC
		GCT gc(*cache_table);

		//set up AggSync
		AggSync<AggregatorT> * agg_thread; //the thread that runs agg_sync()
		if(global_aggregator != NULL) agg_thread = new AggSync<AggregatorT>;

		Profiler* profiler = new Profiler;

		//call status_sync() periodically
		while(global_end_label == false)
		{
			clock_t last_tick = clock();
			bool sth2steal = steal_planning();
            status_sync(sth2steal);
            //------
            //reset idle status of Worker, compers will add back if idle
            mtx_go.lock();
            for(int i=0; i<num_compers; i++) idle_set[i] = false;
            global_num_idle = 0;
            cv_go.notify_all(); //release threads to compute tasks
            mtx_go.unlock();
            usleep(STATUS_SYNC_TIME_GAP);
		}

		if(global_aggregator != NULL) delete agg_thread; //make sure destructor of agg_thread is called to do agg_sync() before exiting run()
		delete profiler;
    }
};

#endif
