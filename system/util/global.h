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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <mpi.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h> //for ease of debug
#include <sys/stat.h>
#include <ext/hash_set>
#include <ext/hash_map>

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#include <atomic>
#include <mutex>
#include "conque.h"

#include "rwlock.h"

//============================
#include <time.h>
#include <stack>
#include <fstream>
#include <chrono>

int POLLING_TIME; //unit: usec, user-configurable, used by sender
//set in init_worker()
static clock_t polling_ticks; // = POLLING_TIME * CLOCKS_PER_SEC / 1000000;

#define SLEEP_PARAM 2000 // POLLING_TIME = SLEEP_PARAM * _num_workers
//this ratio is tested on Azure
//if too small, communication is unbalanced and congestion happens
//if too big, low bandwidth utilization

#define MAX_BATCH_SIZE 1000 //number of bytes sent in a batch

#define WAIT_TIME_WHEN_IDLE 100 //unit: usec, user-configurable, used by recv-er

#define STATUS_SYNC_TIME_GAP 100000 //unit: usec, used by Worker main-thread

#define AGG_SYNC_TIME_GAP 1000000 //unit: usec, used by AggSync main-thread

#define PROGRESS_SYNC_TIME_GAP 1000000 //unit: usec, used by Profiler main-thread

#define TASK_BATCH_NUM 150 //minimal number of tasks processed as a unit
#define TASKMAP_LIMIT 8 * TASK_BATCH_NUM //number of tasks allowed in a task map

#define VCACHE_LIMIT 2000000 //how many vertices allowed in vcache (pull-cache + adj-cache)
#define VCACHE_OVERSIZE_FACTOR 0.2
#define VCACHE_OVERSIZE_LIMIT VCACHE_LIMIT * VCACHE_OVERSIZE_FACTOR

#define MAX_STEAL_TASK_NUM BIG_TASK_FLUSH_BATCH //how many tasks to steal at a time at most

#define MINI_BATCH_NUM 10 //used by spawning from local
#define REQUEST_BOUND 50000 //the maximal number of requests could be sent between each two workers //tuned on GigE

#define GRAPH_LOAD_CHANNEL 200
#define REQ_CHANNEL 201
#define RESP_CHANNEL 202
#define STATUS_CHANNEL 203
#define AGG_CHANNEL 204
#define PROGRESS_CHANNEL 205

void* global_trimmer = NULL;
void* global_taskmap_vec; //set by Worker using its compers, used by RespServer
void* global_vcache;
void* global_local_table;

atomic<int> global_num_idle(0);

conque<string> global_file_list; //tasks buffered on local disk; each element is a file name
atomic<int> global_file_num; //number of files in global_file_list

void* global_vertexes_stack;
mutex global_vertex_stack_lock; //lock for global_vertex_stack

#define TASK_GET_NUM 1
#define TASK_RECV_NUM 1
//try "TASK_GET_NUM" times of fetching and processing tasks
//try "TASK_RECV_NUM" times of inserting processed tasks to task-queue
//============================

#define hash_map __gnu_cxx::hash_map
#define hash_set __gnu_cxx::hash_set

using namespace std;

atomic<bool> global_end_label(false);

//============================
///worker info
#define MASTER_RANK 0

int _my_rank;
int _num_workers;
inline int get_worker_id()
{
    return _my_rank;
}
inline int get_num_workers()
{
    return _num_workers;
}

void init_worker(int * argc, char*** argv)
{
	int provided;
	MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
	if(provided != MPI_THREAD_MULTIPLE)
	{
	    printf("MPI do not Support Multiple thread\n");
	    exit(0);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &_num_workers);
	MPI_Comm_rank(MPI_COMM_WORLD, &_my_rank);
    POLLING_TIME = SLEEP_PARAM * _num_workers;
    polling_ticks = POLLING_TIME * CLOCKS_PER_SEC / 1000000;
}

void worker_finalize()
{
    MPI_Finalize();
}

void worker_barrier()
{
    MPI_Barrier(MPI_COMM_WORLD); //only usable before creating threads
}

//------------------------
// worker parameters

struct WorkerParams {
    string input_path;
    bool force_write;
    bool native_dispatcher; //true if input is the output of a previous blogel job

    WorkerParams()
    {
    	force_write = true;
        native_dispatcher = false;
    }
};

//============================
//big task info
#define BIG_TASK_FLUSH_BATCH 10
#define BIG_TASK_QUEUE_CAPACITY 50
#define BIG_TASKMAP_LIMIT  8 * BIG_TASK_FLUSH_BATCH
conque<string> global_bigTask_fileList;
atomic<int> global_bigTask_file_num; //number of files in global_bigTask_fileList

void* big_task_queue;
mutex bigtask_que_lock;

void* big_map_task;
size_t avg_num;
//============================
//general types
typedef int VertexID;

void* global_aggregator = NULL;

void* global_agg = NULL; //for aggregator, FinalT of previous round
rwlock agg_rwlock;

//============================
string TASK_DISK_BUFFER_DIR;
string REPORT_DIR;

//disk operations
void _mkdir(const char *dir) {//taken from: http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html
	char tmp[256];
	char *p = NULL;
	size_t len;

	snprintf(tmp, sizeof(tmp), "%s", dir);
	len = strlen(tmp);
	if(tmp[len - 1] == '/') tmp[len - 1] = '\0';
	for(p = tmp + 1; *p; p++)
		if(*p == '/') {
				*p = 0;
				mkdir(tmp, S_IRWXU);
				*p = '/';
		}
	mkdir(tmp, S_IRWXU);
}

void _rmdir(string path){
    DIR* dir = opendir(path.c_str());
    struct dirent * file;
    while ((file = readdir(dir)) != NULL) {
        if(strcmp(file->d_name, ".") == 0 || strcmp(file->d_name, "..") == 0)
        	continue;
        string filename = path + "/" + file->d_name;
        remove(filename.c_str());
    }
    if (rmdir(path.c_str()) == -1) {
    	perror ("The following error occurred");
        exit(-1);
    }
}

atomic<bool>* idle_set; //to indicate whether a comper has notified worker of its idleness
mutex mtx_go;
condition_variable cv_go;

//used by profiler
atomic<size_t>* global_tasknum_vec; //set by Worker using its compers, updated by comper, read by profiler
atomic<size_t> num_stolen(0); //number of tasks stolen by the current worker since previous profiling barrier

atomic<size_t>* req_counter; //to count how many requests were sent to each worker

int num_compers;
int BIGTASK_THRESHOLD = 100;

//============= to allow long long to be ID =============
namespace __gnu_cxx {
    template <>
    struct hash<long long> {
        size_t operator()(long long key) const
        {
            return (size_t)key;
        }
    };
}
//====================================================

#endif
