//#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>

#include "graph.h"
#include "RmNonMax.h"

Graph gograph;

double gdmin_deg_ratio; // minimum degree ratio threshold
int gnmin_size; // minimum size threhsold
int gnmax_size; // (maximum size threhsold)
int gnmin_deg; // (min_size - 1) * gamma

void PrintSum(char* szgraph_file, int num_of_vertices);

int main(int argc, char *argv[])
{
	int narg_no, num_of_vertices;

	if(argc!=3)
	{
		printf("Usage:\n");
		printf("\t%s graph_file output_file\n", argv[0]);
		return 0;
	}

	RemoveNonMax(argv[1], argv[2]); // postprocessing to remove non-maximal qcqs

//	PrintSum(argv[1], num_of_vertices);

//	_CrtDumpMemoryLeaks();

	return gntotal_max_cliques;
}


void PrintSum(char* szgraph_file, int num_of_vertices)
{
	FILE *fpsum;
	
	fpsum = fopen("qcliques-new.sum", "a+");
	if(fpsum==NULL)
	{
		printf("Error: cannot open file cliques.sum for write\n");
		return;
	}

	fprintf(fpsum, "%s %d %.3f %d\t", szgraph_file, num_of_vertices, gdmin_deg_ratio, gnmin_size);
	fprintf(fpsum, "%d %d %d\t", gnmax_clique_size, gntotal_cliques, gntotal_max_cliques);
	fprintf(fpsum, "%.3f %.3f\t", gdmining_time, gdrm_nonmax_time);
	fprintf(fpsum, "%d %d %d %d\t", gntotal_calls, gnmaxcheck_calls, gnum_of_condgraphs, gocondgraph_buf.ntotal_pages);
	fprintf(fpsum, "%d\t", gnlookahead_suceeds);
	fprintf(fpsum, "%d %d\t", gntotal_maxcheck_calls, gnmax_maxcheck_depth);
	fprintf(fpsum, "\n");

	fclose(fpsum);
}

