#include <stdio.h>
#include "bmatch.h"

int CCbmat_matchingedges (CCbmat_graph *G, long *pval, int *pmcount,
        int **pmelist, int **pmelen)
{
    int i, rval = 0, mcount = 0;
    int *melist = (int *) NULL, *melen = (int *) NULL;
    CCbmat_edge *e;
    long val = 0;

    for (i = 0, e = G->edgelist; i < G->ecount; i++, e++) {
        if (e->x) {
            val += e->weight;
            mcount++;
        }
    }
    for (e = G->newedges; e; e = e->next) {
        if (e->x) {
            val += e->weight;
            mcount++;
        }
    }

    if (mcount && (pmelist != (int **) NULL || pmelen != (int **) NULL)) {
        if (pmelist != (int **) NULL) {
            melist = CC_SAFE_MALLOC (2*mcount, int);
            CCcheck_NULL (melist, "out of memory for melist");
        }
        if (pmelen != (int **) NULL) {
            melen = CC_SAFE_MALLOC (mcount, int);
            CCcheck_NULL (melen, "out of memory for melen");
        }
        mcount = 0;
        for (i = 0, e = G->edgelist; i < G->ecount; i++, e++) {
            if (e->x) {
                if (melist) {
                    melist[2*mcount]   = e->ends[0] - G->nodelist;
                    melist[2*mcount+1] = e->ends[1] - G->nodelist;
                }
                if (melen) {
                    melen[mcount] = e->weight;
                }
                mcount++;
            }
        }
        for (e = G->newedges; e; e = e->next) {
            if (e->x) {
                if (melist) {
                    melist[2*mcount]   = e->ends[0] - G->nodelist;
                    melist[2*mcount+1] = e->ends[1] - G->nodelist;
                }
                if (melen) {
                    melen[mcount] = e->weight;
                }
                mcount++;
            }
        }
    }
    if (pval) *pval = val;
    if (pmcount) *pmcount = mcount;
    if (pmelist) *pmelist = melist;
    if (pmelen)  *pmelen  = melen;

CLEANUP:
    return rval;
}

static CCbmat_edge *findedge_blossom (CCbmat_node *p, CCbmat_node *n1,
        CCbmat_node *n2)
{
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;

    if (!p) return (CCbmat_edge *) NULL;
    for (pee = p->edgelist; pee; pee = pee->next) {
        pe = CCbmat_THISEDGE(pee);
        if ((pe->ends[0] == n2 && pe->ends[1] == n1) ||
            (pe->ends[1] == n2 && pe->ends[0] == n1)) {
            return pe;
        }
    }
    return findedge_blossom (p->nest.parent, n1, n2);
}

CCbmat_edge *findedge (CCbmat_node *n1, CCbmat_node *n2)
{
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;

    for (pee = n1->edgelist; pee; pee = pee->next) {
        pe = CCbmat_THISEDGE(pee);
        if (pe->ends[0] == n2 || pe->ends[1] == n2) {
            return pe;
        }
    }
    return findedge_blossom (n1->nest.parent, n1, n2);
}

CCbmat_node *smallest (CCbmat_graph *G, CCbmat_node *p)
{
    CCbmat_node *p1, *p2, *psmall = G->nodelist + G->ncount;

    if (!p->status.pseudo)
        return p;
    for (p1 = p->nest.child; p1; p1 = p1->nest.sibling) {
        if ((p2 = smallest(G, p1)) < psmall) {
            psmall = p2;
        }
    }
    return psmall;
}

/* find the second outermost pseudonode containing n */
CCbmat_node *surfsecond (CCbmat_node *n)
{
    CCbmat_node *n1 = n->nest.parent, *n2;

    assert (n1);
    while ((n2 = n1->nest.parent)) {
        n = n1;
        n1 = n2;
    }
    return n;
}

int CCbmat_surfcheck (CCbmat_node *n, CCbmat_node *p)
{
    if (n == p) return 1;
    while ((n = n->nest.parent))
        if (n == p) return 1;
    return 0;
}

CCbmat_node *surfbelow (CCbmat_node *n, CCbmat_node *p)
{
    CCbmat_node *n1;

    while ((n1 = n->nest.parent)) {
        if (n1 == p)
            return n;
        n = n1;
    }
    return NULL;
}

