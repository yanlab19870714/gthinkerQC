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

//test whether the candidate order matter.

#include "kcore.h"
#define hash_set __gnu_cxx::hash_set

using namespace std;
typedef int VertexID;

using get_time = chrono::high_resolution_clock;
using ms = chrono::microseconds;

bool COVER_VERTEX_PRUNE = true;
bool LOOKAHEAD_PRUNE = true;
bool UPPER_BOUND_PRUNE = true;
bool LOWER_BOUND_PRUNE = true;
bool CRITICAL_VERTEX_PRUNE = true;

int COMPACT_SIZE = 10;

using get_time = chrono::high_resolution_clock;
using ms = chrono::microseconds;

double _gamma;
int min_size;
int expand_count = 0;

void load(char* fname, map<VertexID, kc_value>& g_map) {
	ifstream in(fname);
	VertexID vid;
	int temp, num;
	while(in >> vid)
	{
		g_map[vid].del = false;
		in >> num;
		QCValue nbs;
		for(int i=0; i<num; i++)
		{
			in >> temp;
			g_map[vid].nbs.insert(temp);
		}
	}
    in.close();
}

//add vertex from cand_exts to X
void add_vertex(QCSubgraph & new_g, QCVertex * v) {
	QCValue & adj = v->value;
	VertexID v_id = v->id;
	QCVertex new_v;
	new_v.id = v_id;
	int adj_size = adj.size();

	for (int j = 0; j < adj_size; j++) {
		QCVertex * temp = new_g.getVertex(adj[j]);
		if (temp != NULL) {
			temp->value.push_back(v_id);
			new_v.value.push_back(adj[j]);
		}
	}
	new_g.addVertex(new_v);
}

//output the quasi-clique
void output_qcq(QCSubgraph & gs, ofstream & fout) {
	vector<QCVertex> & v_vec = gs.vertexes;
	int gs_size = v_vec.size();
	for (int j = 0; j < gs_size; j++)
		fout << v_vec[j].id << " ";
	fout << "\n";
}

//whether the subgraph is a quasi-clique
bool is_QC(QCSubgraph & gs) {
	vector<QCVertex> & v_vec = gs.vertexes;
	int gs_size = v_vec.size();
	float min_deg = _gamma * (gs_size - 1);
	for (int j = 0; j < gs_size; j++) {
		if (v_vec[j].value.size() < min_deg)
			return false;
	}
	return true;
}

//get the intersection of cand_exts and Nk(v)
void filter_by_2hop(vector<QCVertex*> & new_cand, QCVertex * v,
		QCSubgraph & g, QCHopSubgraph& g_2hop) {
	//get the v's 2 hop neighbors
	set<VertexID>& k_nbs_set = g_2hop.getVertex(v->id)->value;

	//get the intersection of the 2-hop neighbors and new_cand
	//use a temp QCValue to hold the intersection
	int cand_size = new_cand.size();
	vector<QCVertex*> cand_intersect;
	for (int j = 0; j < cand_size; j++) {
		if (k_nbs_set.find(new_cand[j]->id) != k_nbs_set.end())
			cand_intersect.push_back(new_cand[j]);
	}
	//swap the temp(intersection) with new_cand. new_cand is the output here
	new_cand.swap(cand_intersect);
}

//get the indeg of all the candidate.
//get the exdeg of all the vertice in gs.
//!!!!!! make sure all the values in the vectors have been reset to 0
void get_cand_indeg(vector<int> & cand_indeg, vector<int>& X_exdeg,
		QCSubgraph& X_g, vector<QCVertex*>& cand) {
	int n_cand = cand.size();
	hash_map<VertexID, int> & Xg_map = X_g.vmap;
	for (int j = 0; j < n_cand; j++) {
		QCValue & cand_j_adj = cand[j]->value;
		int adj_size = cand_j_adj.size();
		for (int k = 0; k < adj_size; k++) {
			VertexID nb = cand_j_adj[k];
			auto it = Xg_map.find(nb);
			if (it != Xg_map.end()) {
				cand_indeg[j]++;
				X_exdeg[it->second]++;
			}
		}
	}
}

int get_lower_bound(int new_gs_size, int new_cand_size, vector<int> & subg_indeg,
		vector<int> & cand_indeg_sorted) {
	//get the deg_min
	auto it = std::min_element(subg_indeg.begin(), subg_indeg.end());
	int indeg_min = *it;

	//get the loose lower bound L_min and tight lower bound L
	int t = 0; //!!!!!! L_min = t
	while (t <= new_cand_size) {
		if (indeg_min + t >= ceil(_gamma * (new_gs_size + t - 1)))
			break;
		t++;
	}
	//#### ???????????????
	/*if (t < new_cand_size + 1)
		return t;*/
	if (t == new_cand_size + 1)
		return t;
	double l_left = accumulate(subg_indeg.begin(), subg_indeg.end(), 0) +
		accumulate(cand_indeg_sorted.begin(), cand_indeg_sorted.begin() + t, 0);
	double l_right = new_gs_size * ceil(_gamma * (new_gs_size + t - 1));
	//t = L
	//if new_cand_size == 0 and t == 0, there will be a problem for "cand_indeg_sorted[t]".
	//Make sure candidate is not empty before lb calculation.
	while (t < new_cand_size) {
		if (l_left >= l_right) break;
		l_left += cand_indeg_sorted[t];
		l_right = new_gs_size * ceil(_gamma * (new_gs_size + t - 1));
		t++;
	}
	return t;
}

int get_upper_bound(int new_gs_size, int new_cand_size, vector<int> & subg_indeg,
		vector<int> & subg_exdeg, vector<int> & cand_indeg_sorted) {
	//get the deg_min
	int deg_min = INT_MAX;
	for (int j = 0; j < new_gs_size; j++) {
		int subg_deg = subg_indeg[j] + subg_exdeg[j];
		if (deg_min > subg_deg)
			deg_min = subg_deg;
	}

	//get the loose upper bound U_min and tight upper bound U
	int U_min = floor(deg_min / _gamma) + 1 - new_gs_size;
    if(U_min < 1) return 0; // no need to execute on; note from Line 138 that lb >= 1
	int U = U_min < new_cand_size ? U_min : new_cand_size;
	double u_left = accumulate(subg_indeg.begin(), subg_indeg.end(), 0) +
		accumulate(cand_indeg_sorted.begin(), cand_indeg_sorted.begin() + U, 0);
	while (U >= 1) {
		double u_right = new_gs_size * ceil(_gamma * (new_gs_size + U - 1));
		if (u_left >= u_right) break;
		u_left -= cand_indeg_sorted[U - 1];
		U--;
	}
	return U;
}


//pruning based on cover vertices
//see C_X(u) def on paper Page 9, 2nd to last paragraph, Line 2
int cover_prune(QCSubgraph & gs, QCSubgraph & g, vector<QCVertex*> & cand_exts){
	//====== current max
	int max_cover_size = 0;
	vector<VertexID> max_cover_set;
	//======
	int cand_size = cand_exts.size(); // |cand(X)|
	vector<QCVertex> sub_vertices = gs.vertexes; // X
	int sub_size = sub_vertices.size(); // |X|
	//====== create cand(X), to be used for intersecting N(u)
	vector<VertexID> cand_set;
	for(int i = 0; i < cand_size; i++) cand_set.push_back(cand_exts[i]->id);
	sort(cand_set.begin(), cand_set.end());
	//====== create X
	vector<VertexID> sub_set;
	for(int i = 0; i < sub_size; i++) sub_set.push_back(sub_vertices[i].id);
	sort(sub_set.begin(), sub_set.end());
	set<int> X_set(sub_set.begin(), sub_set.end());
	int indeg_threshold = ceil(sub_set.size()*_gamma);

	//====== go through every u in cand(X) to check its cover set
	for (int i = 0; i < cand_size; i++) {
		QCVertex* u = cand_exts[i];
//		if(u->id != 1) continue;
		QCValue & cand_adj = u->value;

		//revise 2(1): quit this u's cover prune check if its indegX is not qualified.
		int indegX_u = 0;
		for(auto v: cand_adj)
			if(X_set.find(v) != X_set.end())
				indegX_u++;
		if(indegX_u < indeg_threshold)
			continue;

		int adj_size = cand_adj.size();
		//------ compute cover_set = intersect(cand_X, adj_u)
		vector<VertexID> cover_set(cand_size); //u's cover set
		auto it = set_intersection(cand_set.begin(), cand_set.end(),
				cand_adj.begin(), cand_adj.end(), cover_set.begin());
		//cand_adj should've already been sorted when loading the graph
		//todo https://en.cppreference.com/w/cpp/algorithm/set_intersection
		cover_set.resize(it - cover_set.begin());
		//------
		if(cover_set.size() <= max_cover_size) continue;
		//------ compute pruned_X = X - adj_u
//		set<VertexID> u_X_set; // intersection of u's nbs and X
//		set_intersection(sub_set.begin(), sub_set.end(),
//				cand_adj.begin(), cand_adj.end(), inserter(u_X_set,u_X_set.begin()));
//
//		vector<VertexID> v_vec;
//		for(auto v:sub_set)
//			if(u_X_set.find(v) == u_X_set.end())
//				v_vec.push_back(v);

		//revise 1: fix the first bug. v should have no edge with u.
		set<VertexID> u_X_set(cand_adj.begin(), cand_adj.end());
		vector<VertexID> v_vec;
		for(auto v:sub_set)
			if(u_X_set.find(v) == u_X_set.end())
				v_vec.push_back(v);

		//revise 2(2): fix a bug. break if v's indeg is not qualified.
		//------ use pruned_X to filter cover_set
		bool v_flag = true;
		for(auto it1 = v_vec.begin(); it1 != v_vec.end(); it1++){
			QCVertex* v = g.getVertex(*it1);
			QCValue & v_adj = v->value;
			int indegX_v = 0;
			for(auto v: v_adj)
				if(X_set.find(v) != X_set.end())
					indegX_v++;
			if(indegX_v < indeg_threshold)
			{
				cover_set.clear();
				v_flag = false;
				break;
			}
		}

		if(v_flag)
		{
			for(auto it1 = v_vec.begin(); it1 != v_vec.end(); it1++)
			{
				QCVertex* v = g.getVertex(*it1);
				QCValue & v_adj = v->value;

				//intersect v_adj and cover_set
				vector<VertexID> new_cover_set(cover_set.size()); //u's cover set
				auto it2 = set_intersection(v_adj.begin(), v_adj.end(),
						cover_set.begin(), cover_set.end(), new_cover_set.begin());
				//v_adj should've already been sorted when loading the graph
				//todo https://en.cppreference.com/w/cpp/algorithm/set_intersection
				new_cover_set.resize(it2 - new_cover_set.begin());
				cover_set.swap(new_cover_set);
			}
		}

		if(cover_set.size() > max_cover_size){
			max_cover_size = cover_set.size();
			max_cover_set.swap(cover_set);
		}
	}
	//====== adjust cand(X) to move cover_set(u) to its end
	if(max_cover_size > 0){
		vector<QCVertex*> cover_vec, uncover_vec;
		for(int i = 0; i < cand_size; i++){
			if(binary_search(max_cover_set.begin(), max_cover_set.end(), cand_exts[i]->id)){
				cover_vec.push_back(cand_exts[i]);
			}else{
				uncover_vec.push_back(cand_exts[i]);
			}
		}
		uncover_vec.insert(uncover_vec.end(), cover_vec.begin(), cover_vec.end());
		cand_exts.swap(uncover_vec);
	}
	return max_cover_size;
}

//true iff the case of extending S (excluding S itself) is pruned;
//ext(S) is passed as a reference,
//and some elements may be pruned when the function returns
bool iterative_bounding(vector<QCVertex*>& new_cand, QCSubgraph& new_gs, ofstream & fout, QCSubgraph& g){
	//until Ly > Uy or Z is empty or cand_y is empty
	 int lb, ub;
	vector<QCVertex>& subg_vertices = new_gs.vertexes;
	int new_gs_size = subg_vertices.size();
	int new_cand_size = new_cand.size();
do {
		//-------I. get the critical vertex

		vector<int> subg_indeg(new_gs_size, 0), subg_exdeg(new_gs_size, 0),
			cand_indeg(new_cand_size, 0);
		//1. get the indeg_x
		for (int j = 0; j < new_gs_size; j++)
			subg_indeg[j] = subg_vertices[j].value.size();
		//2. get the exdeg_x & indeg_cand
		get_cand_indeg(cand_indeg, subg_exdeg, new_gs, new_cand);

		//#### do not sort X-candidate twice
		//presort cand_indeg for lb and ub calculation
		vector<int> cand_indeg_sorted = cand_indeg;
		sort(cand_indeg_sorted.begin(), cand_indeg_sorted.end(), greater<int>());

		//3. get the lower bound for critical vertex
		lb = get_lower_bound(new_gs_size, new_cand_size, subg_indeg, cand_indeg_sorted);
		if (lb > new_cand_size)
			return true;

		//4. get the upper bound for critical vertex
		ub = get_upper_bound(new_gs_size, new_cand_size, subg_indeg, subg_exdeg, cand_indeg_sorted);
		if (ub < 1){
			if(new_gs_size >= min_size && is_QC(new_gs)){
				output_qcq(new_gs, fout);
				return true;
			}
		}
		if(ub < lb)
			return true;

		bool critical_hit;
		do{
			critical_hit = false;
			if(CRITICAL_VERTEX_PRUNE){
				//1. get candidate map
				map<VertexID, QCVertex*> cand_map;
				for (int j = 0; j < new_cand_size; j++)
				{
					QCVertex* v = new_cand[j];
					cand_map[v->id] = v;
				}
				//2. find the critical vertex
				int c_right = ceil(_gamma * (new_gs_size + lb - 1));
				set<VertexID> crit_nbs_set;
				for (int j = 0; j < new_gs_size; j++) {
					if (subg_indeg[j] + subg_exdeg[j] == c_right) {
						//#### check G(S) as in Lines~23--24 before expanding S if find the critical vertex v
						if(new_gs_size >= min_size && is_QC(new_gs))
							output_qcq(new_gs, fout);

						QCValue & critical_nbs = g.getVertex(subg_vertices[j].id)->value;;
						crit_nbs_set.insert(critical_nbs.begin(), critical_nbs.end());
					}
				}

				//add the intersection of the critical vertex's adjacent list and candidate
				//add critical nbs into S; delete critical nbs from candidate
				if(!crit_nbs_set.empty())
				{
					set<VertexID> intersection;
					for (VertexID c_nbs_id: crit_nbs_set) {
						auto it = cand_map.find(c_nbs_id);
						if (it != cand_map.end()) {
							QCVertex* v = it->second;
							//add the vertex to vertex set
							add_vertex(new_gs, v);
							intersection.insert(v->id);
						}
					}
					if(!intersection.empty())
					{
						critical_hit = true;
						// remove it from candidate
						vector<QCVertex*> new_cand_temp;
						for (int k = 0; k < new_cand_size; k++) {
							QCVertex* v = new_cand[k];
							if (intersection.find(v->id) == intersection.end())
								new_cand_temp.push_back(v);
						}
						new_cand.swap(new_cand_temp);
					}

					if(new_cand.empty())
					{
						if(new_gs.vertexes.size() >= min_size && is_QC(new_gs))
							output_qcq(new_gs, fout);
						return true;
					}
				}
			}

			//#### If critical vertex prune never hit, do not need to recalculate the indeg_x, exdeg_x, indeg_cand
			//------II. upper bound and lower bound prune
			//update g_size and candidate size only if critical vertex pruning was executed.
			if(critical_hit){
				new_gs_size = new_gs.vertexes.size();
				new_cand_size = new_cand.size();
				//#### check the candidate size every time it is pruned. It will cause lb calculation error if not check here
				if(new_cand_size == 0) break;

				cand_indeg.assign(new_cand_size, 0);
				subg_indeg.assign(new_gs_size, 0);
				subg_exdeg.assign(new_gs_size, 0);

				//1. update the in_degree of vertices in subgraph
				for (int j = 0; j < new_gs_size; j++)
					subg_indeg[j] = subg_vertices[j].value.size();

				//2. update the in_degree of every vertex in candidate
				//and the ex_degree of every vertex in subgraph
				get_cand_indeg(cand_indeg, subg_exdeg, new_gs, new_cand);

				//#### do not sort X-candidate twice
				//presort cand_indeg for lb and ub calculation
				vector<int> cand_indeg_sorted = cand_indeg;
				sort(cand_indeg_sorted.begin(), cand_indeg_sorted.end(), greater<int>());

				//3. update the update lower bound
				lb = LOWER_BOUND_PRUNE ? get_lower_bound(new_gs_size, new_cand_size, subg_indeg,
						cand_indeg_sorted) : 0;
				//#### return true if can't find valid t
				if (lb > new_cand_size)
					return true;

				//4. update the update upper bound
				ub = UPPER_BOUND_PRUNE ? get_upper_bound(new_gs_size, new_cand_size, subg_indeg,
								subg_exdeg, cand_indeg_sorted) : new_cand_size;
				//#### check G(S) and return true if can't find valid t
				if (ub < 1){
					if(new_gs_size >= min_size && is_QC(new_gs)){
						output_qcq(new_gs, fout);
						return true;
					}
				}
				if(ub < lb)
					return true;
			}
		} while(critical_hit);

		//--------III check vertex in vertex set (Theorems 4,6 and 8)
		bool cond;
		bool deg_cond1_hit = false;
		for (int j = 0; j < new_gs_size; j++) {
			//Theorems 4 (Type II Degree Pruning Condition(i))
			if(!deg_cond1_hit){
				cond = subg_indeg[j] < ceil(_gamma * new_gs_size) && subg_exdeg[j] == 0;
				if(cond) deg_cond1_hit = true;
			}

			//#### Check the S if Theorems 4 (i) hit.
			//Theorems 4 (Type II Degree Pruning Condition(ii))
			cond = subg_indeg[j] + subg_exdeg[j]
					< ceil(_gamma * (new_gs_size + subg_exdeg[j] - 1));
			if (cond) return true;

			//Theorems 6 (Type II Upper Bound Pruning)
			if(UPPER_BOUND_PRUNE){
				cond = subg_indeg[j] + ub
						< ceil(_gamma * (new_gs_size + ub - 1));
				if (cond) return true;
			}
			//Theorems 8 (Type II Lower Bound Pruning)
			if(LOWER_BOUND_PRUNE){
				cond = subg_indeg[j] + subg_exdeg[j]
						< ceil(_gamma * (new_gs_size + lb - 1));
				if (cond) return true;
			}
		}
		if(deg_cond1_hit){
			if(new_gs_size >= min_size && is_QC(new_gs)){
				output_qcq(new_gs, fout);
				return true;
			}
		}
		//#### Move cand_exdeg calculate before part 3
		//calculate the ex_degree of the vertices in candidate
		vector<int> cand_exdeg(new_cand_size, 0);
		//put all candidate in a set
		set<VertexID> new_cand_set;
		for (int j = 0; j < new_cand_size; j++)
			new_cand_set.insert(new_cand[j]->id);
		//check whether each adjacent vertex of candidate vertex is in set
		for (int j = 0; j < new_cand_size; j++) {
			QCValue & temp_adj = new_cand[j]->value;
			int temp_size = temp_adj.size();
			for (int k = 0; k < temp_size; k++) {
				if (new_cand_set.find(temp_adj[k])
						!= new_cand_set.end())
					cand_exdeg[j]++;
			}
		}

		//--------IV check vertex in candidate
		vector<QCVertex*> new_cand_temp;
		for (int j = 0; j < new_cand_size; j++) {
			//Theorems 3 (Type I Degree Pruning)
			cond = cand_indeg[j] + cand_exdeg[j]
					< ceil(_gamma * (new_gs_size + cand_exdeg[j]));
			if(cond) continue;
			//Theorems 5 (Type I Upper Bound Pruning)
			if(UPPER_BOUND_PRUNE){
				cond = cand_indeg[j] + ub - 1
						< ceil(_gamma * (new_gs_size + ub - 1));
				if(cond) continue;
			}
			//Theorems 7 (Type I Lower Bound Pruning)
			if(LOWER_BOUND_PRUNE){
				cond = cand_indeg[j] + cand_exdeg[j]
						< ceil(_gamma * (new_gs_size + lb - 1));
				if(cond) continue;
			}
			new_cand_temp.push_back(new_cand[j]);
		}
		//if no vertex in ext(S) was type-I-pruned, then break
		if(new_cand.size() == new_cand_temp.size()) break;
		else new_cand.swap(new_cand_temp);
		new_cand_size = new_cand.size();
	} while (new_cand_size > 0);
//#### Add QC check for the ext(S) = O condition
if(new_cand_size == 0){
	if(new_gs_size >= min_size && is_QC(new_gs)){
			output_qcq(new_gs, fout);
		}
	return true;
	}
return false;
}


struct deg_sorter{
	QCVertex* v;
	int indeg = 0;
	int exdeg = 0;
	int hop2 = 0;
};
bool comp(const deg_sorter &x, const deg_sorter &y)
{
	if(x.indeg == y.indeg)
	{
		if(x.exdeg == y.exdeg)
		{
			if(x.hop2 == y.hop2)
				return x.v->id < y.v->id;
			else
				return x.hop2 < y.hop2;
		}
		else
			return x.exdeg < y.exdeg;
	}
	else
		return x.indeg < y.indeg;
}

void get_cand_deg(vector<QCVertex*>& cand, QCSubgraph& X_g, QCHopSubgraph& g_2hop, vector<deg_sorter>& cand_deg)
{
	set<VertexID> cand_set;
	for(auto v:cand)
		cand_set.insert(v->id);

	int n_cand = cand.size();
	hash_map<VertexID, int> & Xg_map = X_g.vmap;
//	cand_deg[0].v = cand[0];
	for (int j = 0; j < n_cand; j++)
//	for (int j = 1; j < n_cand; j++)
	{
		cand_deg[j].v = cand[j];
		cand_deg[j].hop2 = g_2hop.getVertex(cand[j]->id)->value.size();
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

void get_cand_deg_1st(vector<QCVertex*>& cand, QCHopSubgraph& g_2hop, vector<deg_sorter>& cand_deg)
{
	int n_cand = cand.size();
//	cand_deg[0].v = cand[0];
//	for (int j = 1; j < n_cand; j++)
	for (int j = 0; j < n_cand; j++)
	{
		cand_deg[j].v = cand[j];
		cand_deg[j].hop2 = g_2hop.getVertex(cand[j]->id)->value.size();
		cand_deg[j].exdeg = cand[j]->value.size();
	}
}

void sort_deg(vector<QCVertex*> & cand_exts, QCSubgraph & gs, QCHopSubgraph& g_2hop){
	vector<deg_sorter> cand_deg(cand_exts.size());
	if(!gs.vertexes.empty())
		get_cand_deg(cand_exts, gs, g_2hop, cand_deg);
	else
		get_cand_deg_1st(cand_exts, g_2hop, cand_deg);
	sort(cand_deg.begin(), cand_deg.end(),comp);
//	sort(cand_deg.begin()+1, cand_deg.end(),comp);
	cand_exts.clear();
	for(auto c: cand_deg)
		cand_exts.push_back(c.v);
}

//revise 5: shrink g
void compact_graph(QCSubgraph & gs, vector<QCVertex*> & cand, QCSubgraph & g,
		QCSubgraph & new_g, QCHopSubgraph & new_g_2hop){
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
		QC2Hop v2h;
		v.id = v2h.id = it->id;
		QCValue temp_nbs;
		QC2HopValue temp_nbs_2hop;
		for(VertexID adj : it->value)
		{
			if(s_cand_set.find(adj) != s_cand_set.end())
			{
				temp_nbs.push_back(adj);
				temp_nbs_2hop.insert(adj);
				QCVertex* nb_it = g.getVertex(adj);
				for(VertexID nb_2hop : nb_it->value)
					if(s_cand_set.find(nb_2hop) != s_cand_set.end())
						temp_nbs_2hop.insert(nb_2hop);
			}
		}
		v.value.swap(temp_nbs);
		new_g.addVertex(v);

		v2h.value.swap(temp_nbs_2hop);
		new_g_2hop.addVertex(v2h);
	}

	for(auto cv : cand)
	{
		QCVertex v;
		QC2Hop v2h;
		v.id = v2h.id = cv->id;

		QCValue temp_nbs;
		QC2HopValue temp_nbs_2hop;
		for(VertexID adj: cv->value)
		{
			if(s_cand_set.find(adj) != s_cand_set.end())
			{
				temp_nbs.push_back(adj);
				temp_nbs_2hop.insert(adj);
				QCVertex* nb_it = g.getVertex(adj);
				for(VertexID nb_2hop : nb_it->value)
					if(s_cand_set.find(nb_2hop) != s_cand_set.end())
						temp_nbs_2hop.insert(nb_2hop);
			}
		}
		v.value.swap(temp_nbs);
		new_g.addVertex(v);

		v2h.value.swap(temp_nbs_2hop);
		new_g_2hop.addVertex(v2h);
	}

	vector<QCVertex*> new_cand;
	for(auto v : cand)
		new_cand.push_back(new_g.getVertex(v->id));
	cand.swap(new_cand);
}

//put the vertex whose degree is max at the beginning
void max_deg_ahead(vector<QCVertex*> & cand_exts)
{
	//first find the max degree index
	//using std swap
	int max_index = 0;
	int max_value = 0;
	for(int i = 0; i<cand_exts.size();i++)
	{
		int v_deg = cand_exts[i]->value.size();
		if(v_deg > max_value)
		{
			max_index = i;
			max_value = v_deg;
		}
	}
	iter_swap(cand_exts.begin(), cand_exts.begin()+max_index);
}

//QCQ
bool QCQ(QCSubgraph & gs, QCSubgraph & g, vector<QCVertex*> & cand_exts,
		ofstream & fout, QCHopSubgraph& g_2hop) {
	//for debug
	expand_count++;
	//put max deg ahead
	//max_deg_ahead(cand_exts);
	//revise 4: sort cand increasingly based on 1)indeg 2)exdeg.
	sort_deg(cand_exts, gs, g_2hop);
	int cand_size = cand_exts.size();
	int gs_size = gs.vertexes.size();
	bool bhas_qclq = false;
	int cover_size = COVER_VERTEX_PRUNE ? cover_prune(gs, g, cand_exts) : 0;
	//cover_prune() already sorted the candidates
	for (int i = 0; i < cand_size - cover_size; i++) {
		if ((gs_size + cand_size - i) < min_size)
			return bhas_qclq;

		//whether the union of subgraph with candidate is QCQ
		if(LOOKAHEAD_PRUNE){
			QCSubgraph union_q = gs;
			for (int k = i; k < cand_size; k++)
				add_vertex(union_q, cand_exts[k]);

			if (is_QC(union_q)) {
				output_qcq(union_q, fout);
				return true;
			}
		}

		//construct new subgraph
		QCVertex * v = cand_exts[i];
		QCSubgraph new_gs = gs;
		add_vertex(new_gs, v);

		//construct new candidate

		vector<QCVertex*> new_cand;
		new_cand.insert(new_cand.begin(), cand_exts.begin() + i + 1,
				cand_exts.end());

		filter_by_2hop(new_cand, v, g, g_2hop);
		int new_gs_size = new_gs.vertexes.size();
		int new_cand_size = new_cand.size();

		if (new_cand_size == 0) {
			if (new_gs_size >= min_size && is_QC(new_gs)) {
				bhas_qclq = true;
				output_qcq(new_gs, fout);
			}

		} else {
			//#### no need to pass lb and ub
			bool ext_prune = iterative_bounding(new_cand, new_gs, fout, g);
			new_gs_size = new_gs.vertexes.size();
			new_cand_size = new_cand.size();
			//-----------------------------------
			if(!ext_prune && new_cand_size + new_gs_size >= min_size){
				bool bhas_super_qclq;
				if((new_cand_size+new_gs_size)*2 < g.vertexes.size()
						&& g.vertexes.size() > COMPACT_SIZE)
				{
					QCSubgraph new_g;
					QCHopSubgraph new_g_2hop;
					compact_graph(new_gs, new_cand, g, new_g, new_g_2hop);
					bhas_super_qclq = QCQ(new_gs, new_g, new_cand, fout, new_g_2hop);
				} else {
					bhas_super_qclq = QCQ(new_gs, g, new_cand, fout, g_2hop);
				}

				bhas_qclq = bhas_qclq || bhas_super_qclq;
				if (!bhas_super_qclq && new_gs_size >= min_size && is_QC(new_gs)) {
					bhas_qclq = true;
					output_qcq(new_gs, fout);
				}
			}
		}
	}
	return bhas_qclq;
}

/*
//get 2hop graph from one hop g_map
void get_2hop_graph(map<VertexID, kc_value>& g_map, QCHopSubgraph& g_2hop)
{
	for(auto it = g_map.begin(); it != g_map.end(); ++it)
	{
		VertexID vid = it->first;
		kc_value & val = it->second;
		set<VertexID>& nbs = val.nbs;
		QC2Hop v;
		v.id = vid;
		for(auto nb : nbs)
		{
			auto it2 = g_map.find(nb);
			if(it2 != NULL)
			{
				v.value.insert(nb);
				set<VertexID>& nbs2hop = it2->second.nbs;
				v.value.insert(nbs2hop.begin(),nbs2hop.end());
			}
		}
		g_2hop.addVertex(v);
	}

}
*/

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
			nbs_2hop.insert(nb_nbs.begin(), nb_nbs.end());
		}
	}
}


//bool deg_comp (QCVertex* i,QCVertex* j)
//{
//	if(i->value.size() == j->value.size())
//		return i->id < j->id;
//	else
//		return i->value.size() < j->value.size();
//}

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        cout<<"arg1 = input file name, arg2 = degree ratio, arg3 = min_size"<<endl;
        return -1;
    }
    string fn = argv[1];
    _gamma = atof(argv[2]);
    min_size = atoi(argv[3]);
    int min_deg = ceil(_gamma * (min_size - 1));
    QCSubgraph g;
    QCHopSubgraph g_2hop;
    map<VertexID, kc_value> g_map;
	load(argv[1], g_map);
    cout<<"loaded"<<endl;

    set_2hop(g_map);

    vector<VertexID> to_del; //output of k-core pruning
	for(auto it = g_map.begin(); it != g_map.end(); ++it)
	{
		kc_value & val = it->second;
		if(val.nbs.size() < min_deg || val.nbs2hop.size() < min_size - 1)
			if(!val.del) prune_2hop(it->first, val, g_map, min_deg, min_size, to_del);
	}
	//now, g_map is degree-pruned
	//do the actual deletes
	for(auto it = to_del.begin(); it != to_del.end(); ++it) g_map.erase(*it);
	//check all adj-items in g_map
	//-- if item is not in 1-hop, move to to_pull
	//-- we remove the item from g_map, since it may be pruned in round 2
	for(auto it = g_map.begin(); it != g_map.end(); ++it){
		QCVertex v;
		QC2Hop v2hop;
		v.id = v2hop.id = it->first;
		set<VertexID> & nbs = it->second.nbs;
		QCValue adj_temp;
		for(auto itr = nbs.begin(); itr != nbs.end(); ++itr){
			adj_temp.push_back(*itr);
		}
		v.value.swap(adj_temp);
		v2hop.value.swap(it->second.nbs2hop);
		g.addVertex(v);
		g_2hop.addVertex(v2hop);
	}

	vector<QCVertex*> cand_exts;
	vector<QCVertex> & g_vertex = g.vertexes;
	for(int i=0; i < g_vertex.size(); i++){
		cand_exts.push_back(&g_vertex[i]);
	}
	//put max degree vertex ahead
	max_deg_ahead(cand_exts);

	//revise 3: sort cand degree increasingly before it start.
//	sort(cand_exts.begin()+1, cand_exts.end(), deg_comp);
	QCSubgraph gs;
	sort_deg(cand_exts, gs, g_2hop);

	string report_path = fn + "_output";
	ofstream fout;
	fout.open(report_path);

    auto Start = get_time::now();

//    int cover_size = COVER_VERTEX_PRUNE ? cover_prune(gs, g, cand_exts) : 0;
    int cover_size = 0;
//    for (int i = 0; i <= cand_exts.size() - min_size; i++)
    int exclude_size = cover_size > min_size - 1 ? cover_size : min_size - 1;
    QCValue& max_u_nbs = cand_exts[0]->value;
    set<VertexID> u_nbs_set(max_u_nbs.begin(),max_u_nbs.end());
    vector<QCVertex*> cand_temp_in, cand_temp_out;
    for(auto cv: cand_exts)
    	if(u_nbs_set.find(cv->id) != u_nbs_set.end())
    		cand_temp_out.push_back(cv);
    	else
    		cand_temp_in.push_back(cv);
    cand_temp_in.insert(cand_temp_in.end(), cand_temp_out.begin(), cand_temp_out.end());

    cand_exts.swap(cand_temp_in);

	for (int i = 0; i < cand_exts.size() - exclude_size; i++)
	{
		QCVertex * v = cand_exts[i];
		QCSubgraph gs;
		add_vertex(gs, v);

		vector<QCVertex*> new_cand;
		new_cand.insert(new_cand.begin(), cand_exts.begin() + i + 1,
				cand_exts.end());
		filter_by_2hop(new_cand, v, g, g_2hop);

		if (new_cand.size() < min_size - 1) continue;
		QCSubgraph new_g;
		QCHopSubgraph new_g_2hop;
		compact_graph(gs, new_cand, g, new_g, new_g_2hop);
		QCQ(gs, new_g, new_cand, fout, new_g_2hop);
	}
    auto  End = get_time::now();
    auto elapsedTime = chrono::duration_cast<ms>(End - Start);
    cout<<"time: "<<elapsedTime.count()<<endl;
    cout<<"expand count: "<<expand_count<<endl;
    fout.close();
    return 0;
}
