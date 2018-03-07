#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "edgegen.h"
#include "bmatch.h"
#include "kdtree.h"

int CCbmat_christofides (int ncount, CCdatagroup *dat, int *val, int *tour,
    CCrandstate *rstate);

typedef struct edge {
    int ends[2];
    int mark;
} edge;

typedef struct node {
    int degree;
    edge **elist;
} node;

#define OTHEREND(e,n) (e->ends[0] == n ? e->ends[1] : e->ends[0])

static int euler_tour (int ncount, int ecount, int *elist, int *etour);
static void  find_euler (int v, int ncount, node *nodelist, edge *edgelist,
    int *tourcnt, int *etour);

int CCbmat_christofides (int ncount, CCdatagroup *dat, int *val, int *outtour,
        CCrandstate *rstate)
{
    int i, j, nodd, norm, tval, rval = 0;
    int *tree_edges = (int *) NULL;
    int *degrees = (int *) NULL, *odd_nodes = (int *) NULL;
    long mval;
    int ecount, mcount;
    int *melist = (int *) NULL, *melen = (int *) NULL;
    int *elist = (int *) NULL, *etour = (int *) NULL;
    int *marks = (int *) NULL, *tour = (int *) NULL;
    
    double treeval;
    CCdatagroup odd_dat;

    printf ("CCbmat_christofides ...\n"); fflush (stdout);

    CCutil_init_datagroup (&odd_dat);

    tree_edges = CC_SAFE_MALLOC (2*(ncount-1), int);
    CCcheck_NULL (tree_edges, "out of memory for tree_edges");
    degrees = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (degrees, "out of memory for degrees");

    rval = CCkdtree_prim_spanningtree ((CCkdtree * ) NULL, ncount, dat,
                  (double *) NULL, tree_edges, &treeval, rstate);
    CCcheck_rval (rval, "CCkdtree_prim_spanningtree failed");

    printf ("Spanning Tree: %.0f\n", treeval);
    fflush (stdout);

    for (i = 0; i < ncount; i++) degrees[i] = 0;
    for (i = 0; i < ncount-1; i++) {
        degrees[tree_edges[2*i]]++;
        degrees[tree_edges[2*i+1]]++;
    }

    nodd = 0;
    for (i = 0; i < ncount; i++) {
        if (degrees[i] % 2) nodd++;
    }
    printf ("Odd-degree nodes: %d\n", nodd); fflush (stdout);

    odd_nodes = CC_SAFE_MALLOC (nodd, int);
    CCcheck_NULL (odd_nodes, "out of memory for odd_nodes");
    nodd = 0;
    for (i = 0; i < ncount; i++) {
        if (degrees[i] % 2) {
            odd_nodes[nodd++] = i;
        }
    }

    CCutil_dat_getnorm (dat, &norm);
    CCutil_dat_setnorm (&odd_dat, norm);

    odd_dat.x = CC_SAFE_MALLOC (nodd, double);
    CCcheck_NULL (odd_dat.x, "out of memory for odd_dat.x");
    odd_dat.y = CC_SAFE_MALLOC (nodd, double);
    CCcheck_NULL (odd_dat.y, "out of memory for odd_dat.y");
    for (i = 0; i < nodd; i++) {
        odd_dat.x[i] = dat->x[odd_nodes[i]];
        odd_dat.y[i] = dat->y[odd_nodes[i]];
    }

    rval = CCbmat_onematch_geom (nodd, &odd_dat, &mval, &mcount, &melist,
                                 &melen, rstate);
    CCcheck_rval (rval, "CCbat_onematch_geom failed");
    printf ("Matching: %ld\n", mval); fflush (stdout);
    printf ("Tree + Odd Matching: %d\n", (int) treeval + (int) mval);
    fflush (stdout);

    ecount = (ncount-1) + mcount;
    elist = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_rval (rval, "out of memory for elist");
    for (i = 0; i < ncount-1; i++) {
        elist[2*i]   = tree_edges[2*i];
        elist[2*i+1] = tree_edges[2*i+1];
    }
    for (i = 0; i < mcount; i++) {
        elist[2*(ncount-1+i)]   = odd_nodes[melist[2*i]];
        elist[2*(ncount-1+i)+1] = odd_nodes[melist[2*i+1]];
    }

    etour = CC_SAFE_MALLOC (ecount+1, int);
    CCcheck_NULL (etour, "out of memory for etour");

    rval = euler_tour (ncount, ecount, elist, etour);
    CCcheck_rval (rval, "euler_tour failed");

    marks = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (tour, "out of memory for tour");

    for (i = 0, j = 0; j < ecount+1; j++) {
        if (marks[etour[j]] == 0) {
            tour[i++] = etour[j];
            marks[etour[j]] = 1;
        }
    }
    if (i != ncount) {
        printf ("j = %d, Euler tour did not reach all nodes\n", i);
        rval = 1; goto CLEANUP;
    }

    tval = 0;
    for (i = 1; i < ncount; i++) {
        tval += CCutil_dat_edgelen (tour[i-1], tour[i], dat);
    }
    tval += CCutil_dat_edgelen (tour[0], tour[ncount-1], dat);
    printf ("Christofides Tour: %d\n", tval);

    if (outtour) {
        for (i = 0; i < ncount; i++) outtour[i] = tour[i];
    }
    if (val) *val = tval;


CLEANUP:
    CCutil_freedatagroup (&odd_dat);
    CC_IFFREE (tree_edges, int);
    CC_IFFREE (degrees, int);
    CC_IFFREE (odd_nodes, int);
    CC_IFFREE (melist, int);
    CC_IFFREE (melen, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (etour, int);
    CC_IFFREE (marks, int);
    CC_IFFREE (tour, int);
    return rval;
}

static int euler_tour (int ncount, int ecount, int *elist, int *etour)
{
    int i, tourcnt, rval = 0;
    node *nodelist = (node *) NULL;
    edge *edgelist = (edge *) NULL;
    node *n0, *n1;
 
    edgelist = CC_SAFE_MALLOC (ecount, edge);
    CCcheck_NULL (edgelist, "out of memory for edgelist");
    for (i = 0; i < ecount; i++) {
        edgelist[i].ends[0] = elist[2*i];
        edgelist[i].ends[1] = elist[2*i+1];
        edgelist[i].mark = 0;
    }

    nodelist = CC_SAFE_MALLOC (ncount, node);
    CCcheck_NULL (nodelist, "out of memory for nodelist");
    for (i = 0; i < ncount; i++) {
        nodelist[i].degree = 0;
        nodelist[i].elist = (edge **) NULL;
    }
    for (i = 0; i < ecount; i++) {
        nodelist[edgelist[i].ends[0]].degree++;
        nodelist[edgelist[i].ends[1]].degree++;
    }
    for (i = 0; i < ncount; i++) {
        nodelist[i].elist = CC_SAFE_MALLOC (nodelist[i].degree, edge *);
        CCcheck_NULL (nodelist[i].elist, "out of memory for nodelist[i]");
        nodelist[i].degree = 0;
    }
    for (i = 0; i < ecount; i++) {
        n0 = &nodelist[edgelist[i].ends[0]];
        n1 = &nodelist[edgelist[i].ends[1]];
        n0->elist[n0->degree++] = &edgelist[i];
        n1->elist[n1->degree++] = &edgelist[i];
    }


    tourcnt = 0;
    find_euler (0, ncount, nodelist, edgelist, &tourcnt, etour);
    if (tourcnt != ecount+1) {
        printf ("Euler tour does not use all edges\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    if (nodelist) {
        for (i = 0; i < ncount; i++) {
            CC_IFFREE (nodelist[i].elist, edge *);
        }
        CC_IFFREE (nodelist, node);
    }
    CC_IFFREE (edgelist, edge);
    return rval;
}


static void  find_euler (int v, int ncount, node *nodelist, edge *edgelist,
        int *tourcnt, int *etour)
{
    int i, u;
    node *nv = &nodelist[v];
    edge *e;

    for (i = 0; i < nv->degree; i++) {
        e = nv->elist[i];
        if (e->mark == 0) {
            e->mark = 1;
            u = OTHEREND (e,v);
            find_euler (u, ncount, nodelist, edgelist, tourcnt, etour);
        }
    }
    etour[(*tourcnt)++] = v;
}

