#include <stdlib.h>
#include <stdio.h>


void RemoveNonMax(char* szset_filename, char* szoutput_filename);
extern int gntotal_max_cliques;


struct TREE_NODE
{
	int nid;
	TREE_NODE *pchild;
	TREE_NODE *pright_sib;
	bool bis_max;
};


#define TNODE_PAGE_SIZE (1<<10)

struct TNODE_PAGE
{
	TREE_NODE ptree_nodes[TNODE_PAGE_SIZE];
	TNODE_PAGE *pnext;
};

struct TNODE_BUF
{
	TNODE_PAGE *phead;
	TNODE_PAGE *pcur_page; 
	int ncur_pos;
	int ntotal_pages;
};

extern TNODE_BUF gotreenode_buf;

inline TREE_NODE* NewTreeNode()
{
	TREE_NODE *ptnode;
	TNODE_PAGE *pnew_page;

	if(gotreenode_buf.ncur_pos==TNODE_PAGE_SIZE)
	{
		if(gotreenode_buf.pcur_page->pnext==NULL)
		{
			pnew_page = new TNODE_PAGE;
			pnew_page->pnext = NULL;
			gotreenode_buf.pcur_page->pnext = pnew_page;
			gotreenode_buf.pcur_page = pnew_page;
			gotreenode_buf.ntotal_pages++;
		}
		else 
			gotreenode_buf.pcur_page = gotreenode_buf.pcur_page->pnext;
		gotreenode_buf.ncur_pos = 0;
	}
	
	ptnode = &(gotreenode_buf.pcur_page->ptree_nodes[gotreenode_buf.ncur_pos]);
	gotreenode_buf.ncur_pos++;

	ptnode->bis_max = true;
	
	return ptnode;
}

inline void OutputOneSet(FILE *fp, int *pset, int nlen)
{
	int i;

	gntotal_max_cliques++;

	fprintf(fp, "%d ", nlen);
	for(i=0;i<nlen;i++)
		fprintf(fp, "%d ", pset[i]);
	fprintf(fp, "\n");

}



