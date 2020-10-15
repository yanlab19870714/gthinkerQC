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

#include "../divide_conquer_version/graph.h"
int COMPACT_SIZE = 10;

struct ContextValue
{
	int round; // the iteration number in compute(.)
	map<VertexID, kc_value> g_map; // for k-core pruning: g_map[v] = (neighbors, v_deleted?)
	Graph split_g; //for time delay split

	//for Expand() input
	int nclique_size;
	int num_of_cands;
	int num_of_tail_vertices;
	VERTEX *pvertices;

	ContextValue(){
		round = -1;
		pvertices = NULL;
	}

	~ContextValue(){
		if(round == 3){
			split_g.DestroySplitGraph();
			if(pvertices != NULL)
				delete []pvertices;
		}
	}
};


obinstream & operator>>(obinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.g_map;
    if(c.round == 3)
    {
    	m >> c.split_g;
		m >> c.nclique_size;
		m >> c.num_of_cands;
		m >> c.num_of_tail_vertices;
		int num_of_vertices = c.nclique_size + c.num_of_cands + c.num_of_tail_vertices;
		c.pvertices = new VERTEX[num_of_vertices];
		for(int i = 0; i < num_of_vertices; i++)
			m >> c.pvertices[i];
    }
    return m;
}

ibinstream & operator<<(ibinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.g_map;
    if(c.round == 3)
    {
    	m << c.split_g;
		m << c.nclique_size;
		m << c.num_of_cands;
		m << c.num_of_tail_vertices;
		int num_of_vertices = c.nclique_size + c.num_of_cands + c.num_of_tail_vertices;
		for(int i = 0; i < num_of_vertices; i++)
			m << c.pvertices[i];
    }
    return m;
}

ofbinstream & operator>>(ofbinstream & m, ContextValue & c)
{
    m >> c.round;
    m >> c.g_map;
    if(c.round == 3)
    {
    	m >> c.split_g;
		m >> c.nclique_size;
		m >> c.num_of_cands;
		m >> c.num_of_tail_vertices;
		int num_of_vertices = c.nclique_size + c.num_of_cands + c.num_of_tail_vertices;
		c.pvertices = new VERTEX[num_of_vertices];
		for(int i = 0; i < num_of_vertices; i++)
			m >> c.pvertices[i];
    }
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const ContextValue & c)
{
    m << c.round;
    m << c.g_map;
    if(c.round == 3)
    {
    	m << c.split_g;
		m << c.nclique_size;
		m << c.num_of_cands;
		m << c.num_of_tail_vertices;
		int num_of_vertices = c.nclique_size + c.num_of_cands + c.num_of_tail_vertices;
		for(int i = 0; i < num_of_vertices; i++)
			m << c.pvertices[i];
    }
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
		if(task->context.round == 3)
		{
			if(task->context.num_of_cands > BIGTASK_THRESHOLD)
				return true;
			else
				return false;
		}
		else
		{
			if(task->to_pull.size() > BIGTASK_THRESHOLD)
				return true;
			else
				return false;
		}
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

	int Expand(VERTEX *pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, FILE *gfpout, Graph& gograph)
	{
		VERTEX *pnew_vertices, *pnew_cands, *pclique;
		int num_of_vertices, num_of_new_cands, i, j, num_of_new_tail_vertices, nmin_deg;
		bool bis_subsumed, blook_succeed, bgen_new_lvl2nbs;
		int nisvalid, nremoved_vertices;
		int nsuperclique_size, nmax_clique_size, nnew_clique_size;
		struct timeb cur_time;
		double drun_time;
		CLQ_STAT one_clq_stat;

		num_of_vertices = nclique_size+num_of_cands+num_of_tail_vertices;
		pnew_vertices = new VERTEX[num_of_vertices];
		pclique = new VERTEX[nclique_size+1];
		nmax_clique_size = 0;

		for(i=nclique_size;i<nclique_size+num_of_cands && pvertices[i].bis_cand && pvertices[i].bto_be_extended;i++) // not iterating covered vertices (segment 3)
		{
			gograph.VerifyVertices(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, false);
			if(i>nclique_size)
				nisvalid = gograph.RemoveCandVertex(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i);
			else
				nisvalid = 1;
			if(nisvalid==-1)
				break;
			else if(nisvalid==1 && nclique_size+1<gnmax_size)
			{
				if(i<nclique_size+num_of_cands-1)
				{
					blook_succeed = gograph.Lookahead(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, pnew_vertices, gfpout);
					if(blook_succeed)
					{
						if(nmax_clique_size<nclique_size+(nclique_size+num_of_cands-i))
							nmax_clique_size = nclique_size+(nclique_size+num_of_cands-i);
						break;
					}
				}

				if(gdmin_deg_ratio==1)
					num_of_new_cands = gograph.AddOneVertex(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, true, pnew_vertices, num_of_new_tail_vertices, &one_clq_stat);
				else
					num_of_new_cands = gograph.AddOneVertex(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, false, pnew_vertices, num_of_new_tail_vertices, &one_clq_stat);
				nnew_clique_size = nclique_size+1;
				if(num_of_new_cands>0)
					gograph.CrtcVtxPrune(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, pclique, &one_clq_stat);

				bis_subsumed = false;
				pnew_cands = &pnew_vertices[nnew_clique_size];
				if(num_of_new_cands>0)
				{
					ftime(&cur_time);
					drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
					if(drun_time < TIME_THRESHOLD)
					{
						bgen_new_lvl2nbs = gograph.GenCondGraph(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices);
						nremoved_vertices = gograph.ReduceCands(pnew_vertices, nnew_clique_size, num_of_new_cands+num_of_new_tail_vertices, bis_subsumed);

						if(num_of_new_cands-nremoved_vertices>0) // there is still some candidate left for expansion
							nsuperclique_size = Expand(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, gfpout, gograph);
						else
							nsuperclique_size = 0;
						if(bgen_new_lvl2nbs)
							gograph.DelCondGraph();
					} else {
//						cout<<"split ..."<<endl;
						//to split
						QCliqueTask * t = new QCliqueTask;// should I move the following into next {}
						t->context.split_g.mnum_of_vertices = gograph.mnum_of_vertices;
						t->context.split_g.mblvl2_flag = gograph.mblvl2_flag;
						t->context.round = 3;
						t->context.pvertices = NULL; //for ~ContextValue()
						gograph.ForceGenCondGraph(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, t->context.split_g);
						nremoved_vertices = gograph.ReduceCands(pnew_vertices, nnew_clique_size, num_of_new_cands+num_of_new_tail_vertices, bis_subsumed);
						if(num_of_new_cands-nremoved_vertices>0) // there is still some candidate left for expansion
						{
							//Expand(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, gfpout, gograph);
							//prepare task
							//"=" operator overloading in Graph
//							t->context.split_g = gograph;
							t->context.split_g.index2id = new int[gograph.mnum_of_vertices];
							memcpy(t->context.split_g.index2id, gograph.index2id, sizeof(int)*gograph.mnum_of_vertices);

							//prepare input of Expand()
							t->context.pvertices = new VERTEX[num_of_vertices];
							memcpy(t->context.pvertices, pnew_vertices, sizeof(VERTEX)*(num_of_vertices));
							t->context.nclique_size = nnew_clique_size;
							t->context.num_of_cands = num_of_new_cands;
							t->context.num_of_tail_vertices = num_of_new_tail_vertices;
							add_task(t);
						} else {
//							Graph& split_g =  t->context.split_g;
//							for(i=0;i<split_g.mnum_of_vertices;i++)
//								delete []split_g.mpplvl2_nbs[i];
//							delete []split_g.mpplvl2_nbs;
//
//							for(i=0;i<split_g.mnum_of_vertices;i++)
//								delete []split_g.mppadj_lists[i];
//							delete []split_g.mppadj_lists;
							delete t;
						}

						//alway need to check current task when split happen
						nsuperclique_size = 0;

//						if(bgen_new_lvl2nbs)
//							gograph.DelCondGraph();
					}
				}
				else
					nsuperclique_size = 0;

				//check whether current set is a QC

				if(nsuperclique_size==0 && !bis_subsumed)
				{
					if(nnew_clique_size>=gnmin_size)
					{
						nmin_deg = gograph.GetMinDeg(nnew_clique_size);
						for(j=0;j<nnew_clique_size;j++)
						{
							if(pnew_vertices[j].nclique_deg<nmin_deg)
								break;
						}
						if(j>=nnew_clique_size)
						{
							if(gdmin_deg_ratio<1)
								num_of_new_tail_vertices = 0;
							else if(num_of_new_tail_vertices==0)
								num_of_new_tail_vertices = gograph.GenTailVertices(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, pnew_vertices, nnew_clique_size);
							gograph.OutputOneClique(pnew_vertices, nnew_clique_size, num_of_new_tail_vertices, gfpout);
							if(nmax_clique_size<nnew_clique_size)
								nmax_clique_size = nnew_clique_size;
						}
						else if(nnew_clique_size>nclique_size+1 && nclique_size+1>=gnmin_size)
						{
							nmin_deg = gograph.GetMinDeg(nclique_size+1);
							for(j=0;j<=nclique_size;j++)
							{
								if(pclique[j].nclique_deg<nmin_deg)
									break;
							}
							if(j>nclique_size)
							{
								memcpy(pnew_vertices, pclique, sizeof(VERTEX)*(nclique_size+1));
								num_of_new_tail_vertices = 0;
								gograph.OutputOneClique(pnew_vertices, nclique_size+1, num_of_new_tail_vertices, gfpout);
								if(nmax_clique_size<nclique_size+1)
									nmax_clique_size = nclique_size+1;
							}
						}
					}
				}
				else if(nsuperclique_size>0)
				{
					if(nmax_clique_size<nsuperclique_size)
						nmax_clique_size = nsuperclique_size;
				}
				else if(nmax_clique_size<nnew_clique_size)
					nmax_clique_size = nnew_clique_size;
			}
		}
		delete []pnew_vertices;
		delete []pclique;

	//	if(num_of_cands>=10 && nmax_clique_size==0 && nisvalid==1)
	//		printf("stop\n");

		return nmax_clique_size;
	}

	int ExpandOnce(VERTEX *pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, FILE *gfpout, Graph& gograph)
	{
		VERTEX *pnew_vertices, *pnew_cands, *pclique;
		int num_of_vertices, num_of_new_cands, i, j, num_of_new_tail_vertices, nmin_deg;
		bool bis_subsumed, blook_succeed, bgen_new_lvl2nbs;
		int nisvalid, nremoved_vertices;
		int nsuperclique_size, nmax_clique_size, nnew_clique_size;
		struct timeb cur_time;
		double drun_time;
		CLQ_STAT one_clq_stat;

		num_of_vertices = nclique_size+num_of_cands+num_of_tail_vertices;
		pnew_vertices = new VERTEX[num_of_vertices];
		pclique = new VERTEX[nclique_size+1];
		nmax_clique_size = 0;

		for(i=nclique_size;i<nclique_size+1 && pvertices[i].bis_cand && pvertices[i].bto_be_extended;i++) // not iterating covered vertices (segment 3)
		{
			gograph.VerifyVertices(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, false);
			if(i>nclique_size)
				nisvalid = gograph.RemoveCandVertex(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i);
			else
				nisvalid = 1;
			if(nisvalid==-1)
				break;
			else if(nisvalid==1 && nclique_size+1<gnmax_size)
			{
				if(i<nclique_size+num_of_cands-1)
				{
					blook_succeed = gograph.Lookahead(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, pnew_vertices, gfpout);
					if(blook_succeed)
					{
						if(nmax_clique_size<nclique_size+(nclique_size+num_of_cands-i))
							nmax_clique_size = nclique_size+(nclique_size+num_of_cands-i);
						break;
					}
				}

				if(gdmin_deg_ratio==1)
					num_of_new_cands = gograph.AddOneVertex(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, true, pnew_vertices, num_of_new_tail_vertices, &one_clq_stat);
				else
					num_of_new_cands = gograph.AddOneVertex(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, false, pnew_vertices, num_of_new_tail_vertices, &one_clq_stat);
				nnew_clique_size = nclique_size+1;
				if(num_of_new_cands>0)
					gograph.CrtcVtxPrune(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, pclique, &one_clq_stat);

				bis_subsumed = false;
				pnew_cands = &pnew_vertices[nnew_clique_size];
				if(num_of_new_cands>0)
				{
					ftime(&cur_time);
					drun_time = cur_time.time-gograph.gtime_start.time+(double)(cur_time.millitm-gograph.gtime_start.millitm)/1000;
					if(drun_time < TIME_THRESHOLD)
					{
						bgen_new_lvl2nbs = gograph.GenCondGraph(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices);
						nremoved_vertices = gograph.ReduceCands(pnew_vertices, nnew_clique_size, num_of_new_cands+num_of_new_tail_vertices, bis_subsumed);

						if(num_of_new_cands-nremoved_vertices>0) // there is still some candidate left for expansion
							nsuperclique_size = Expand(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, gfpout, gograph);
						else
							nsuperclique_size = 0;
						if(bgen_new_lvl2nbs)
							gograph.DelCondGraph();
					} else {
//						cout<<"split ..."<<endl;
						//to split
						QCliqueTask * t = new QCliqueTask;// should I move the following into next {}
						t->context.split_g.mnum_of_vertices = gograph.mnum_of_vertices;
						t->context.split_g.mblvl2_flag = gograph.mblvl2_flag;
						t->context.round = 3;
						t->context.pvertices = NULL; //for ~ContextValue()
						gograph.ForceGenCondGraph(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, t->context.split_g);
						nremoved_vertices = gograph.ReduceCands(pnew_vertices, nnew_clique_size, num_of_new_cands+num_of_new_tail_vertices, bis_subsumed);
						if(num_of_new_cands-nremoved_vertices>0) // there is still some candidate left for expansion
						{
							//Expand(pnew_vertices, nnew_clique_size, num_of_new_cands, num_of_new_tail_vertices, gfpout, gograph);
							//prepare task
							//"=" operator overloading in Graph
//							t->context.split_g = gograph;
							t->context.split_g.index2id = new int[gograph.mnum_of_vertices];
							memcpy(t->context.split_g.index2id, gograph.index2id, sizeof(int)*gograph.mnum_of_vertices);

							//prepare input of Expand()
							t->context.pvertices = new VERTEX[num_of_vertices];
							memcpy(t->context.pvertices, pnew_vertices, sizeof(VERTEX)*(num_of_vertices));
							t->context.nclique_size = nnew_clique_size;
							t->context.num_of_cands = num_of_new_cands;
							t->context.num_of_tail_vertices = num_of_new_tail_vertices;
							add_task(t);
						} else {
//							Graph& split_g =  t->context.split_g;
//							for(i=0;i<split_g.mnum_of_vertices;i++)
//								delete []split_g.mpplvl2_nbs[i];
//							delete []split_g.mpplvl2_nbs;
//
//							for(i=0;i<split_g.mnum_of_vertices;i++)
//								delete []split_g.mppadj_lists[i];
//							delete []split_g.mppadj_lists;
							delete t;
						}

						//alway need to check current task when split happen
						nsuperclique_size = 0;

//						if(bgen_new_lvl2nbs)
//							gograph.DelCondGraph();
					}
				}
				else
					nsuperclique_size = 0;

				if(nsuperclique_size==0 && !bis_subsumed)
				{
					if(nnew_clique_size>=gnmin_size)
					{
						nmin_deg = gograph.GetMinDeg(nnew_clique_size);
						for(j=0;j<nnew_clique_size;j++)
						{
							if(pnew_vertices[j].nclique_deg<nmin_deg)
								break;
						}
						if(j>=nnew_clique_size)
						{
							if(gdmin_deg_ratio<1)
								num_of_new_tail_vertices = 0;
							else if(num_of_new_tail_vertices==0)
								num_of_new_tail_vertices = gograph.GenTailVertices(pvertices, nclique_size, num_of_cands, num_of_tail_vertices, i, pnew_vertices, nnew_clique_size);
							gograph.OutputOneClique(pnew_vertices, nnew_clique_size, num_of_new_tail_vertices, gfpout);
							if(nmax_clique_size<nnew_clique_size)
								nmax_clique_size = nnew_clique_size;
						}
						else if(nnew_clique_size>nclique_size+1 && nclique_size+1>=gnmin_size)
						{
							nmin_deg = gograph.GetMinDeg(nclique_size+1);
							for(j=0;j<=nclique_size;j++)
							{
								if(pclique[j].nclique_deg<nmin_deg)
									break;
							}
							if(j>nclique_size)
							{
								memcpy(pnew_vertices, pclique, sizeof(VERTEX)*(nclique_size+1));
								num_of_new_tail_vertices = 0;
								gograph.OutputOneClique(pnew_vertices, nclique_size+1, num_of_new_tail_vertices, gfpout);
								if(nmax_clique_size<nclique_size+1)
									nmax_clique_size = nclique_size+1;
							}
						}
					}
				}
				else if(nsuperclique_size>0)
				{
					if(nmax_clique_size<nsuperclique_size)
						nmax_clique_size = nsuperclique_size;
				}
				else if(nmax_clique_size<nnew_clique_size)
					nmax_clique_size = nnew_clique_size;
			}
		}
		delete []pnew_vertices;
		delete []pclique;

	//	if(num_of_cands>=10 && nmax_clique_size==0 && nisvalid==1)
	//		printf("stop\n");

		return nmax_clique_size;
	}


	int Cliques(map<VertexID, kc_value> & gmap, FILE *gfpout, Graph& gograph)
	{
		int i, j, num_of_vertices, num_of_cands, nmax_deg, nmaxdeg_vertex;
		int nrm_start, nrm_end, nvertex_no, norder, num_of_freq_cands;
		struct timeb start;
		VERTEX *pvertices, onevertex;

		ftime(&start);
		gograph.gtime_start = start;

		//	gfpout = fopen(szoutput_filename, "wt");

		if(gfpout==NULL)
		{
	//		printf("Error: cannot open file %s for write\n", szoutput_filename);
			printf("Error: cannot open file for write\n");
			return 0;
		}


	//	VertexID root = 0;
		vector<vector<VertexID> > adj_list(gmap.size());
		gograph.TranGraph(gmap, adj_list);
		num_of_vertices = gograph.LoadGraph(adj_list);
		gograph.gpmin_degs = new int[num_of_vertices+1];
		gograph.CalcMinDegs(gograph.gpmin_degs, num_of_vertices); // precompute min-deg threshold for different clique size
		gograph.gpremoved_vertices = new int[num_of_vertices]; // ---1---: those vertices that fail the check: degree >= r * (min_size - 1)  AND  |vertices_within-2hops| > min_size - 1
		gograph.gpremain_vertices = new int[num_of_vertices]; // ---2---: those vertices that pass the check: degree >= r * (min_size - 1)  AND  |vertices_within-2hops| > min_size - 1
		gograph.gpvertex_order_map  = new int[num_of_vertices]; // map[vertex_no] = position in pvertices[.], this is needed as vertices change their positions in pvertices[.]
		memset(gograph.gpvertex_order_map, -1, sizeof(int)*num_of_vertices);
		gograph.gptemp_array = new int[num_of_vertices];

		pvertices = new VERTEX[num_of_vertices];
		num_of_cands = num_of_vertices;
		num_of_freq_cands = 0;
		nrm_start = 0;
		nrm_end = 0;
		for(i=0;i<num_of_vertices;i++)
		{
			pvertices[i].nvertex_no = i;
			pvertices[i].nclique_deg = 0;
			pvertices[i].ncand_deg = gograph.mppadj_lists[i][0];
			if(gograph.mblvl2_flag)
				pvertices[i].nlvl2_nbs = gograph.mpplvl2_nbs[i][0];
			else
				pvertices[i].nlvl2_nbs = 0;
			if(gograph.mppadj_lists[i][0]>=gnmin_deg && (!gograph.mblvl2_flag || gograph.mblvl2_flag && gograph.mpplvl2_nbs[i][0]>=gnmin_size-1)) // degree >= r * (min_size - 1)  AND  |vertices_within-2hops| > min_size - 1
			{
				//todo return if root id is pruned here
				pvertices[i].bis_cand = true;
				pvertices[i].bto_be_extended = true;
				gograph.gpremain_vertices[num_of_freq_cands++] = i;
			}
			else
			{
				pvertices[i].bis_cand = false;
				pvertices[i].bto_be_extended = false;
				gograph.gpremoved_vertices[nrm_end++] = i;
			}
		}
		//set ncand_deg based on qualified v: if this v is qualified, then it's neighbor's ncand_deg ++
		while(num_of_freq_cands<num_of_cands/2)
		{//?? how to set the degree here given the S and ext_S
			num_of_cands = num_of_freq_cands;
			for(i=0;i<num_of_cands;i++)
				pvertices[gograph.gpremain_vertices[i]].ncand_deg = 0;
			for(i=0;i<num_of_cands;i++)
			{
				nvertex_no = gograph.gpremain_vertices[i];
				for(j=1;j<=gograph.mppadj_lists[nvertex_no][0];j++)
				{
					norder = gograph.mppadj_lists[nvertex_no][j];
					if(pvertices[norder].bis_cand)
						pvertices[norder].ncand_deg++;
				}
			}
			num_of_freq_cands = 0;
			nrm_start = 0;
			nrm_end = 0;
			for(i=0;i<num_of_cands;i++)
			{
				norder = gograph.gpremain_vertices[i];
				if(pvertices[norder].ncand_deg>=gnmin_deg) // degree updated: degree >= r * (min_size - 1)
					gograph.gpremain_vertices[num_of_freq_cands++] = norder;
				else
				{
					pvertices[norder].bis_cand = false;
					pvertices[norder].bto_be_extended = false;
					gograph.gpremoved_vertices[nrm_end++] = norder;
				}
			}
		}
		//todo return if root id is removed in above check
		num_of_cands = num_of_freq_cands;
	//update the degree based on removed v. I will go to neigbor's neigbor, just like our kcore-prune
		while(nrm_end>nrm_start)
		{
			nvertex_no = gograph.gpremoved_vertices[nrm_start];
			for(i=1;i<=gograph.mppadj_lists[nvertex_no][0];i++)
			{
				norder = gograph.mppadj_lists[nvertex_no][i];
				if(pvertices[norder].bis_cand)
				{
					pvertices[norder].ncand_deg--;
					if(pvertices[norder].ncand_deg<gnmin_deg) // degree updated: degree < r * min_size
					{
						pvertices[norder].bis_cand = false;
						pvertices[norder].bto_be_extended = false;
						num_of_cands--;
						gograph.gpremoved_vertices[nrm_end++] = norder;
					}
				}
			}
			nrm_start++;
		}
		//todo return if root id is removed in above check

		qsort(&pvertices[1], num_of_vertices-1, sizeof(VERTEX), comp_vertex_freq);

		gograph.VerifyVertices(pvertices, 0, num_of_cands, 0, false);
	//should revise as our S will not be empty.
		if(num_of_cands>=gnmin_size)
		{
			gograph.mnmax_cond_graphs = MAX_COND_GRAPH;
			gograph.mpcond_graphs = new COND_GRAPH[gograph.mnmax_cond_graphs];
			gograph.mnum_of_cond_graphs = 0;
			gograph.gocondgraph_buf.phead = new INT_PAGE;
			gograph.gocondgraph_buf.phead->pnext = NULL;
			gograph.gocondgraph_buf.pcur_page = gograph.gocondgraph_buf.phead;
			gograph.gocondgraph_buf.ncur_pos = 0;
			gograph.gocondgraph_buf.ntotal_pages = 1;
			gograph.gpcondgraph_listpt_buf = new int*[gograph.mnum_of_vertices*2*MAX_COND_GRAPH];
			gograph.gpcounters = new int[num_of_cands];

			ExpandOnce(pvertices, 0, num_of_cands, 0, gfpout, gograph);

			delete []gograph.gpcounters;
			gograph.DelCGIntBuf();
			delete []gograph.gpcondgraph_listpt_buf;
			delete []gograph.mpcond_graphs;
			if(gograph.mnum_of_cond_graphs!=0)
				printf("Error: the number of conditiaonal database should be 0\n");
		}

		delete []pvertices;
		delete []gograph.gpremoved_vertices;
		delete []gograph.gpremain_vertices;
		delete []gograph.gpvertex_order_map;
		delete []gograph.gptemp_array;
		delete []gograph.gpmin_degs;
		delete []gograph.index2id;
		gograph.DestroyGraph();

		return num_of_vertices;
	}

    virtual bool compute(SubgraphT & g, ContextT & context, vector<VertexT *> & frontier)
    {
    	if(context.round == 1)
    	{
    		VertexID rootID = g.vertexes[0].id;
    		map<VertexID, kc_value>& g_map = context.g_map; //used for k-core pruning
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
    		VertexID rootID = g.vertexes[0].id;
    		map<VertexID, kc_value>& g_map = context.g_map; //used for k-core pruning
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
			int num_of_vertices = Cliques(g_map, gfpout, gograph);
			//no time delay version only need 2 round
			return false;
			//### since there is no vertex requested, step 3 will run directly
			//### so we only need round 3 to call QCQ to simplify code
    	}else{
    		Graph split_g = context.split_g;
    		split_g.SetupGraph(context.nclique_size, context.num_of_cands, context.num_of_tail_vertices);
    		ftime(&split_g.gtime_start);
    		Expand(context.pvertices, context.nclique_size, context.num_of_cands, context.num_of_tail_vertices, gfpout, split_g);

    		split_g.ClearGraph(); //delete graph's variable
//    		delete []context.pvertices;
    		//pvertices and graph data will be cleared by DestroySplitGraph() in ~ContextValue()
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
    if(argc != 7){
    	cout<<"arg1 = input path in HDFS, arg2 = number of threads"
    			<<", arg3 = degree ratio, arg4 = min_size, arg5 = time delay threshold"
				<<", arg6 = BIGTASK_THRESHOLD"<<endl;
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
    TIME_THRESHOLD = atof(argv[5]);
    BIGTASK_THRESHOLD = atoi(argv[6]);
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
