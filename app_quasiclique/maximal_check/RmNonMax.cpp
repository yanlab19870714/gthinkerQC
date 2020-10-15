#include <stdio.h>
#include <time.h>
#include <sys/timeb.h>

#include "graph.h"
#include "RmNonMax.h"

int gntotal_max_cliques;

TNODE_BUF gotreenode_buf;

void DelTNodeBuf()
{
	TNODE_PAGE *ppage;

	ppage = gotreenode_buf.phead;
	while(ppage!=NULL)
	{
		gotreenode_buf.phead = gotreenode_buf.phead->pnext;
		delete ppage;
		gotreenode_buf.ntotal_pages--;
		ppage = gotreenode_buf.phead;
	}
	if(gotreenode_buf.ntotal_pages!=0)
		printf("Error: inconsistent number of pages\n");
}

void InsertOneSet(int* pset, int nlen, TREE_NODE *&proot)
{
	TREE_NODE *pnode, *pparent, *pleftsib, *pnew_node;
	int i, j;

	i = 0;
	pparent = NULL;
	pnode = proot;
	pleftsib = NULL;

	while(i<nlen)
	{
		while(pnode!=NULL && pnode->nid<pset[i])
		{
			pleftsib = pnode;
			pnode = pnode->pright_sib;
		}

		if(pnode==NULL || pnode->nid>pset[i])
		{
			pnew_node = NewTreeNode();
			pnew_node->nid = pset[i];
			pnew_node->pchild = NULL;
			pnew_node->pright_sib = pnode;
			if(pleftsib!=NULL)
				pleftsib->pright_sib = pnew_node;
			else if(pparent!=NULL)
				pparent->pchild = pnew_node;
			if(i==0 && pleftsib==NULL)
				proot = pnew_node;
			pparent = pnew_node;
			for(j=i+1;j<nlen;j++)
			{
				pnew_node = NewTreeNode();
				pnew_node->nid = pset[j];
				pnew_node->pchild = NULL;
				pnew_node->pright_sib = NULL;
				pparent->pchild = pnew_node;
				pparent = pnew_node;
			}
			break;
		}
		else 
		{
			pparent = pnode;
			pnode = pnode->pchild;
			pleftsib = NULL;
		}
		i++;
	}
}

int BuildTree(char* szset_filename, TREE_NODE* &proot)
{
	FILE *fp;
	int nlen, *pset, nset_size, i, nmax_len, num_of_sets;

	fp = fopen(szset_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szset_filename);
		return 0;
	}

	gotreenode_buf.phead = new TNODE_PAGE;
	gotreenode_buf.phead->pnext = NULL;
	gotreenode_buf.pcur_page = gotreenode_buf.phead;
	gotreenode_buf.ntotal_pages = 1;
	gotreenode_buf.ncur_pos = 0;

	proot = NULL;

	num_of_sets = 0;

	nset_size = 100;
	pset = new int[nset_size];

	nmax_len = 0;
	fscanf(fp, "%d", &nlen);
	while(!feof(fp))
	{
		if(nmax_len<nlen)
			nmax_len = nlen;
		if(nlen>nset_size)
		{
			delete []pset;
			nset_size *= 2;
			if(nset_size<nlen)
				nset_size = nlen;
			pset = new int[nset_size];
		}
		for(i=0;i<nlen;i++)
			fscanf(fp, "%d", &pset[i]);
		qsort(pset, nlen, sizeof(int), comp_int);
		InsertOneSet(pset, nlen, proot);

		num_of_sets++;
		fscanf(fp, "%d", &nlen);
	}
	fclose(fp);

	delete []pset;

	return nmax_len;
}

void SearchSubset(int *pset, int nset_len, TREE_NODE *proot, TREE_NODE** pstack, int *ppos)
{
	TREE_NODE *pnode;
	int ntop, npos;

	if(proot==NULL)
		return;
	ntop = 0;
	npos =  0;
	pnode = proot;

	while(ntop>=0)
	{
		while(pnode!=NULL && npos<nset_len && pnode->nid!=pset[npos])
		{
			if(pnode->nid<pset[npos])
				pnode = pnode->pright_sib;
			else 
				npos++;
		}
		if(pnode!=NULL && npos<nset_len)
		{
			if(pnode->pchild==NULL && pnode->bis_max)
				pnode->bis_max = false;
			pstack[ntop] = pnode;
			ppos[ntop] = npos;
			ntop++;
			pnode = pnode->pchild;
			npos++;
		}
		else 
		{
			ntop--;
			if(ntop>=0)
			{
				pnode = pstack[ntop]->pright_sib;
				npos = ppos[ntop]+1;
			}
		}
	}

}

void RmNonMax(TREE_NODE *proot, int nmax_len)
{
	TREE_NODE *pnode, **pstack, **psearch_stack;
	int *pset, ntop, i, *ppos;

	pset = new int[nmax_len];
	pstack = new TREE_NODE*[nmax_len];
	psearch_stack = new TREE_NODE*[nmax_len];
	ppos = new int[nmax_len];

	pstack[0] = proot;
	pset[0] = proot->nid;
	ntop = 1;
	pnode = proot;

	while(ntop>0)
	{
		if(pnode->pchild!=NULL)
		{
			pnode = pnode->pchild;
			pstack[ntop] = pnode;
			pset[ntop] = pnode->nid;
			ntop++;
		}
		else
		{
			if(ntop>=2 && pnode->bis_max)
			{
				for(i=ntop-1;i>=1;i--)
				{
					if(pstack[i-1]->pright_sib!=NULL)
						SearchSubset(&pset[i], ntop-i, pstack[i-1]->pright_sib, psearch_stack, ppos);
				}
			}

			while(ntop>0 && pnode->pright_sib==NULL)
			{
				ntop--;
				if(ntop>0)
					pnode = pstack[ntop-1];
			}
			if(ntop==0)
				break;
			else //if(pnode->pright_sib!=NULL)
			{
				pnode = pnode->pright_sib;
				pstack[ntop-1] = pnode;
				pset[ntop-1] = pnode->nid;
			}
		}
	}

	delete []pset;
	delete []pstack;
	delete []psearch_stack;
	delete []ppos;
}

void OutputMaxSet(TREE_NODE *proot, int nmax_len, char* szoutput_filename)
{
	FILE *fp;
	TREE_NODE **pstack, *pnode;
	int *pset, ntop;

	fp = fopen(szoutput_filename, "wt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for write\n", szoutput_filename);
		return;
	}

	pstack = new TREE_NODE*[nmax_len];
	pset = new int[nmax_len];

	pstack[0] = proot;
	pset[0] = proot->nid;
	ntop = 1;
	pnode = proot;

	while(ntop>0)
	{
		if(pnode->pchild!=NULL)
		{
			pnode = pnode->pchild;
			pstack[ntop] = pnode;
			pset[ntop] = pnode->nid;
			ntop++;
		}
		else
		{
			if(pnode->bis_max)
				OutputOneSet(fp, pset, ntop);

			while(ntop>0 && pnode->pright_sib==NULL)
			{
				ntop--;
				if(ntop>0)
					pnode = pstack[ntop-1];
			}
			if(ntop==0)
				break;
			else //if(pnode->pright_sib!=NULL)
			{
				pnode = pnode->pright_sib;
				pstack[ntop-1] = pnode;
				pset[ntop-1] = pnode->nid;
			}
		}
	}

	delete []pstack;
	delete []pset;

	fclose(fp);
}

void RemoveNonMax(char* szset_filename, char* szoutput_filename)
{
	TREE_NODE *proot;
	int nmax_len;
	struct timeb start, end;

	ftime(&start);

	gntotal_max_cliques = 0;

	nmax_len = BuildTree(szset_filename, proot);
	RmNonMax(proot, nmax_len);
	OutputMaxSet(proot, nmax_len, szoutput_filename);

	DelTNodeBuf();

	ftime(&end);
	gdrm_nonmax_time = end.time-start.time+(double)(end.millitm-start.millitm)/1000;

	printf("#maximal cliques: %d\n", gntotal_max_cliques);
	printf("Maximality checking time: %.3f\n", gdrm_nonmax_time);
}


