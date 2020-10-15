#pragma once

#pragma warning (disable:4996)

struct VERTEX // variables start with 'b' is boolean, 'n' is integer
{
	int nvertex_no; // vertex ID
	int nclique_deg; // degree in the context of X
	int ncand_deg; // degree in the context of cand(X)
	int nlvl2_nbs;
	bool bis_cand;
	bool bto_be_extended;
};

#define CG_PAGE_SIZE (1<<12)
#define MAX_COND_GRAPH 5

struct INT_PAGE
{
	int parray[CG_PAGE_SIZE];
	INT_PAGE *pnext;
};

struct COND_GRAPH_BUF
{
	INT_PAGE *phead;
	INT_PAGE *pcur_page;
	int ncur_pos;
	int ntotal_pages;
};

struct COND_GRAPH
{
	int **ppadj_lists;
	int **pplvl2_nbs;
	int num_of_vertices;
	INT_PAGE *pcur_page;
	int ncur_pos;
};

struct CLQ_STAT //lower and upper bounds: see Lines 602-607 of graph.cpp
{
	int nmin_cands; //L_X
	int nmax_cands; //U_X
	int nmin_ext_deg;
};

class Graph
{
// prefix 'm' means it's member of Graph, 'p' is a pointer, 'pp' is a pointer (array) of pointers
	int **mppadj_lists, *mpadj_list_buf, **mpplvl2_nbs; // mppadj_lists[i] = node i's 1-hop neighbors, or more accurately, pointer to the adjlist segment in mpadj_list_buf;    mpplvl2_nbs[i] = node i's 2-hop neighbor array
	int mnum_of_vertices; // total vertex number
	COND_GRAPH *mpcond_graphs; // conditional graphs
	int mnum_of_cond_graphs, mnmax_cond_graphs;
	bool mblvl2_flag; // whether need to check 2-hop neighbors, set by Graph::LoadGraph(.)

public:
	int LoadGraph(char* szgraph_file);
	void GenLevel2NBs();
	void DestroyGraph();
	void OutputLvl2Graph(char* szoutput_filename);

	int Cliques(char *szgraph_file, char* szoutput_filename);
	int Expand(VERTEX *pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices);

	int RemoveCandVertex(VERTEX *pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, int ncur_pos);
	bool Lookahead(VERTEX* pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, int ncur_pos, VERTEX *pnew_vertices);
	int AddOneVertex(VERTEX *pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, int ncur_pos, bool bhas_tail, VERTEX *pnew_vertices, int &num_of_new_tails, CLQ_STAT *pclq_stat);
	void CrtcVtxPrune(VERTEX *pvertices, int &nclique_size, int &num_of_cands, int &num_of_tail_vertices, VERTEX *pclique, CLQ_STAT *pclq_stat);
	bool GenCondGraph(VERTEX* pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices);
	void DelCondGraph();
	int ReduceCands(VERTEX *pvertices, int nclique_size, int num_of_cands, bool &bis_subsumed);
	int GenTailVertices(VERTEX* pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, int ncur_pos, VERTEX *pnew_vertices, int nnew_clique_size); // cover vertex pruning
	int ReduceTailVertices(VERTEX* pvertices, int nclique_size, int num_of_tail_vertices, int **ppadj_lists);
	void OutputOneClique(VERTEX *pvertices, int nclique_size, int num_of_tail_vertices);
	bool CheckMaximal(VERTEX* pvertices, int nclique_size, int num_of_exts);
	

	void VerifyVertices(VERTEX* pvertices, int nclique_size, int num_of_cands, int num_of_tail_vertices, bool bonly_tobeextended);
};

//------------------------------------------------------------------------
extern COND_GRAPH_BUF gocondgraph_buf;
extern int** gpcondgraph_listpt_buf;

inline int* NewCGIntArray(int nlen)
{
	INT_PAGE *pnew_page;
	int *parray;

	if(nlen>CG_PAGE_SIZE)
		return new int[nlen];
	else if(gocondgraph_buf.ncur_pos+nlen>CG_PAGE_SIZE)
	{
		if(gocondgraph_buf.pcur_page->pnext==NULL)
		{
			pnew_page = new INT_PAGE;
			pnew_page->pnext = NULL;
			gocondgraph_buf.pcur_page->pnext = pnew_page;
			gocondgraph_buf.pcur_page = pnew_page;
			gocondgraph_buf.ntotal_pages++;
		}
		else
			gocondgraph_buf.pcur_page = gocondgraph_buf.pcur_page->pnext;

		gocondgraph_buf.ncur_pos = 0;
	}

	parray = &gocondgraph_buf.pcur_page->parray[gocondgraph_buf.ncur_pos];
	gocondgraph_buf.ncur_pos += nlen;

	return parray;
}

void DelCGIntBuf();



//===============================================================
int comp_int(const void *e1, const void *e2);
int comp_int_des(const void* e1, const void *e2);
int comp_vertex_freq(const void *e1, const void *e2);
int comp_vertex_clqdeg(const void *e1, const void *e2);

extern double gdmin_deg_ratio;
extern int gnmin_size;
extern int gnmax_size;
extern int gnmin_deg;
extern int *gpmin_degs;

extern int gntotal_cliques;
extern int gnmax_clique_size;
extern int gntotal_calls;
extern int gndepth;
extern int gnmax_depth;
extern int gnmaxcheck_calls;
extern int gntotal_maxcheck_calls;
extern int gnmaxcheck_depth;
extern int gnmax_maxcheck_depth;
extern int gnum_of_condgraphs;
extern int gnlookahead_suceeds;
extern double gdmining_time;
extern double gdrm_nonmax_time;

extern FILE *gfpout; // output file for keeping QCs


inline int CalcMaxExts(int nmin_clq_totaldeg, int nclique_size, int num_of_cands)
{
	return (int)(nmin_clq_totaldeg/gdmin_deg_ratio)-nclique_size;
//	return num_of_cands;
}


inline void CalcMinDegs(int *pmin_degs, int num_of_vertices)
{
	int i, nmin_deg;

	pmin_degs[0] = 0;
	for(i=1;i<=num_of_vertices;i++)
	{
		nmin_deg = (int)(gdmin_deg_ratio*(i-1));
		if(nmin_deg<gdmin_deg_ratio*(i-1))
			nmin_deg++;
		pmin_degs[i] = nmin_deg; //^^^^^ precomputed here
	}
}

inline int GetMinDeg(int nclique_size)
{
	if(nclique_size<=gnmin_size)
		return gnmin_deg; // Q cannot be less the min-size, so not point to have a smaller degree than r * min_size
	else
		return gpmin_degs[nclique_size]; //^^^^^ use precomputed mindeg for the particular clique-size
}

inline bool IsValidCand(VERTEX* pvertex, int nclique_size, CLQ_STAT *pclq_stat)
{
	if(pvertex->nclique_deg+pvertex->ncand_deg<pclq_stat->nmin_ext_deg)
		return false;
	else if(pvertex->nclique_deg+pvertex->ncand_deg<GetMinDeg(nclique_size+pvertex->ncand_deg+1))
		return false;
	else if(pvertex->nclique_deg+pclq_stat->nmax_cands-1<gpmin_degs[nclique_size+pclq_stat->nmax_cands])
		return false;
	else if(pvertex->nclique_deg+pvertex->ncand_deg<GetMinDeg(nclique_size+pclq_stat->nmin_cands))
		return false;
	else
		return true;
}

inline void Output1Clique(VERTEX *pclique, int nclique_size)
{
	gntotal_cliques++;
	if(gnmax_clique_size<nclique_size)
		gnmax_clique_size = nclique_size;

	int i;
	fprintf(gfpout, "%d ", nclique_size);
	for(i=0; i<nclique_size; i++) 
		fprintf(gfpout, "%d ", pclique[i].nvertex_no);
	fprintf(gfpout, "\n");

}




