#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "edgegen.h"
#include "bmatch.h"

static int call_bmatch (int ncount, int ecount, int *elist, int *elen,
    int bdegree, int *degrees, long *val, char *pedgefname,
    CCdatagroup *dat, int *mcount, int **melist, int **melen);
static int call_bmatch_geom (int ncount, CCdatagroup *dat, int bdegree,
    int *degrees, long *val, int *mcount, int **melist, int **melen,
    CCrandstate *rstate);
static int bmatch (CCbmat_graph *G, int bdegree, int *degrees, long *val,
    CCdatagroup *dat, char *pedgefname, CCbmat_node **pseudo_freelist,
    CCbmat_edgeptr **edge_supply);
static int initmatch (CCbmat_graph *G, CCbmat_heaps *H);
static void clearaugment (CCbmat_node *n);
static long cleanup (CCbmat_graph *G, int bdegree, int *degrees);
static void fixmark (CCbmat_node *n);
static long cleanupnode (CCbmat_node *n, int bdegree);
static int countteeth (CCbmat_node *n);
static void mark_node (CCbmat_node *p);
static int build_graph (CCbmat_graph *G, int ncount, int ecount, int *elist,
    int *elen);
static void free_graph (CCbmat_graph *G);
static void doubleweights (CCbmat_graph *G);
static long minedge (CCbmat_edgeptr *pee);
static void setupmatch (CCbmat_graph *G, int bdegree, int *degrees);
static void jumpsetupmatch (CCbmat_graph *G, int bdegree, int *degrees);
static int flipout (CCbmat_edge *e);
static void free_nodelist (CCbmat_node *nhead);
static int get_init_edges (CCdatagroup *dat, int ncount, int *pecount,
    int **pelist, int **pelen, int bdegree, int *degrees, CCrandstate *rstate);

static int ddoo = 0;

int CCbmat_onematch_geom (int ncount, CCdatagroup *dat, long *val, int *mcount,
        int **melist, int **melen, CCrandstate *rstate)
{
    int rval = 0;
    rval = call_bmatch_geom (ncount, dat, 1, (int *) NULL, val, mcount,
                             melist, melen, rstate);
    CCcheck_rval (rval, "call_bmatch_geom failed");

CLEANUP:
    return rval;
}

int CCbmat_twomatch_geom (int ncount, CCdatagroup *dat, long *val, int *mcount,
        int **melist, int **melen, CCrandstate *rstate)
{
    int rval = 0;
    rval = call_bmatch_geom (ncount, dat, 2, (int *) NULL, val, mcount,
                             melist, melen, rstate);
    CCcheck_rval (rval, "call_bmatch_geom failed");

CLEANUP:
    return rval;
}

int CCbmat_bmatch_geom (int ncount, CCdatagroup *dat, int *degrees, long *val,
        int *mcount, int **melist, int **melen, CCrandstate *rstate)
{
    int rval = 0;
    rval = call_bmatch_geom (ncount, dat, 0, degrees, val, mcount,
                             melist, melen, rstate);
    CCcheck_rval (rval, "call_bmatch_geom failed");

CLEANUP:
    return rval;
}

int CCbmat_onematch (int ncount, int ecount, int *elist, int *elen,
        long *val, CCdatagroup *dat, char *pedgefname, int *mcount,
        int **melist, int **melen)
{
    int rval = 0;
    rval = call_bmatch (ncount, ecount, elist, elen, 1, (int *) NULL,
                        val, pedgefname, dat, mcount, melist, melen);
    CCcheck_rval (rval, "call_bmatch failed");

CLEANUP:
    return rval;
}

int CCbmat_twomatch (int ncount, int ecount, int *elist, int *elen,
        long *val, CCdatagroup *dat, char *pedgefname, int *mcount,
        int **melist, int **melen)
{
    int rval = 0;
    rval = call_bmatch (ncount, ecount, elist, elen, 2, (int *) NULL,
                        val, pedgefname, dat, mcount, melist, melen);
    CCcheck_rval (rval, "call_bmatch failed");

CLEANUP:
    return rval;
}

int CCbmat_bmatch (int ncount, int ecount, int *elist, int *elen,
        int *degrees, long *val, CCdatagroup *dat, char *pedgefname,
        int *mcount, int **melist, int **melen)
{
    int rval = 0;
    rval = call_bmatch (ncount, ecount, elist, elen, 0, degrees,
                        val, pedgefname, dat, mcount, melist, melen);
    CCcheck_rval (rval, "call_bmatch failed");

CLEANUP:
    return rval;
}

static int call_bmatch_geom (int ncount, CCdatagroup *dat, int bdegree,
        int *degrees, long *val, int *mcount, int **melist, int **melen,
        CCrandstate *rstate)
{
    int rval = 0;
    int ecount, *elist = (int *) NULL, *elen = (int *) NULL;

    rval = get_init_edges (dat, ncount, &ecount, &elist, &elen, bdegree,
                      degrees, rstate);
    CCcheck_rval (rval, "get_init_edges failed");

    printf ("Nodes: %d, Edges %d\n", ncount, ecount);
    fflush (stdout);
    
    rval = call_bmatch (ncount, ecount, elist, elen, bdegree, degrees, val,
                 (char *) NULL, dat, mcount, melist, melen);
    CCcheck_rval (rval, "call_bmatch failed");

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    return rval;
}

static int call_bmatch (int ncount, int ecount, int *elist, int *elen,
        int bdegree, int *degrees, long *val, char *pedgefname,
        CCdatagroup *dat, int *mcount, int **melist, int **melen)
{
    int rval = 0;
    int jumpstart = 1;
    long pval, dval;
    CCbmat_node *pseudo_freelist = (CCbmat_node *) NULL;
    CCbmat_edgeptr *edge_supply = (CCbmat_edgeptr *) NULL, *e, *enext;
    CCbmat_graph G;

    G.nodelist = (CCbmat_node *) NULL;
    G.edgelist = (CCbmat_edge *) NULL;
    G.newedges = (CCbmat_edge *) NULL;

    rval = build_graph (&G, ncount, ecount, elist, elen);
    CCcheck_rval (rval, "build_graph failed");

    doubleweights (&G);
    if (jumpstart) {
        rval = CCbmat_jump (&G, bdegree, degrees);
        CCcheck_rval (rval, "jump failed");
        jumpsetupmatch (&G, bdegree, degrees);
    } else {
        setupmatch (&G, bdegree, degrees);
    }

    rval = bmatch (&G, bdegree, degrees, &dval, dat, pedgefname, 
                   &pseudo_freelist, &edge_supply);
    CCcheck_rval (rval, "bmatch failed")

    rval = CCbmat_matchingedges (&G, &pval, mcount, melist, melen);
    CCcheck_rval (rval, "CCbmat_matchinedges failed");

    if (pval <= dval-1 || dval <= pval-1) {
        fprintf (stderr,"Error, dual (%ld) != primal (%ld)\n",dval,pval);
        goto CLEANUP;
    }

    if (val) *val = pval / 2;

    if (degrees) printf ("\nOptimum b-Matching:  ");
    else         printf ("\nOptimum %d-Matching:  ", bdegree);
    printf ("%ld\n", pval / 2);
    fflush (stdout);

    rval = CCbmat_freecombs (&G, &pseudo_freelist);
    CCcheck_rval (rval, "CCbmat_freecombs failed");

    printf ("ddoo = %d\n", ddoo);

CLEANUP:
    /* free edge_supply */
    e = edge_supply;
    while (e) {
        enext = e->next;
        CC_IFFREE (e->this, CCbmat_edge);
        CC_FREE (e, CCbmat_edgeptr);
        e = enext;
    }

    free_nodelist (pseudo_freelist);
    free_graph (&G);
    return rval;
}

static int bmatch (CCbmat_graph *G, int bdegree, int *degrees, long *val,
        CCdatagroup *dat, char *pedgefname, CCbmat_node **pseudo_freelist,
        CCbmat_edgeptr **edge_supply)
{
    int i, rval = 0;
    CCbmat_node *p;
    CCbmat_heaps H;

    CCbmat_init_heaps_struct (&H);

    rval = initmatch (G, &H);
    CCcheck_rval (rval, "initmatch failed");

    for (i=0, p = G->nodelist; i<G->ncount; i++, p++) {
        while (p->deficiency) {
            printf ("."); fflush (stdout);
            CCbmat_augment2 (&H, p, CCbmat_MAXIWEIGHT, pseudo_freelist);
        }
    }
    if (dat) {
        rval = CCbmat_pricedats (dat, G, &H, pseudo_freelist,
                                edge_supply);
        CCcheck_rval (rval, "CCbmat_pricedats failed");
    }
    if (pedgefname) {
        rval = CCbmat_priceedges (pedgefname, G, &H, pseudo_freelist,
                                 edge_supply);
        CCcheck_rval (rval, "CCbmat_priceedges failed");
    }
    i = cleanup (G, bdegree, degrees);
    if (val) *val = (long) i;

CLEANUP:
    CCbmat_free_heaps (&H);
    return rval;
}

static int initmatch (CCbmat_graph *G, CCbmat_heaps *H)
{
    int i, rval = 0;
    CCbmat_edge *e;
    CCbmat_node *n;

    for (i=0, e=G->edgelist; i<G->ecount; i++, e++) {
        e->status.needschange = CCbmat_FALSE;
        e->count = 0;
        e->istooth[0] = e->istooth[1] = 0;
        e->surfends[0] = e->ends[0];
        e->surfends[1] = e->ends[1];
        e->heap.type = CCbmat_NOTHEAPED;
    }
    for (i=0, n=G->nodelist; i<G->ncount; i++,n++) {
        n->parentedge = NULL;
        n->tree.parent = n->tree.child = n->tree.sibling = NULL;
        n->nest.parent = n->nest.child = n->nest.sibling = NULL;
        n->status.pseudo = CCbmat_FALSE;
        n->status.inpath = CCbmat_FALSE;
        n->label = CCbmat_UNLABELED;
        n->surf = n;
        n->heap.type = CCbmat_NOTHEAPED;
    }

    rval = CCbmat_init_heaps (H);
    CCcheck_rval (rval, "CCbmat_init_heaps failed");

CLEANUP:
    return rval;
}

void CCbmat_augment2 (CCbmat_heaps *H, CCbmat_node *n, int ybound,
        CCbmat_node **pseudo_freelist)
{
    /* stop when n->ymat == ybound */

    assert (n->ymat <= ybound);

    n->label = CCbmat_PLUSMAT;
    n->tree.parent = n->tree.child = n->tree.sibling = NULL;
    n->deficiency--;

    CCbmat_heap_label_update (H, n);

    while (!CCbmat_heap_dispatch (H, ybound, pseudo_freelist)) ddoo++;
    CCbmat_clear_heaps (H);
    clearaugment (n->surf);
}

static void clearaugment (CCbmat_node *n)
{
    CCbmat_node *p, *p1;

    /* mark base nodes as changed, for pricing */
    if (n->label == CCbmat_PLUSMAT) mark_node (n);

    for (p = n->tree.child; p; p = p1) {
        p1 = p->tree.sibling;
        clearaugment (p);
    }
    n->tree.parent = n->tree.child = n->tree.sibling = NULL;
    n->parentedge = NULL;
    n->label = CCbmat_UNLABELED;
}

static long cleanup (CCbmat_graph *G, int bdegree, int *degrees)
{
    long dualval = 0;
    int i, k;
    CCbmat_node *n;
    CCbmat_edge *e;

    for (i=0, n=G->nodelist; i<G->ncount; i++, n++) {
        k = (degrees ? degrees[i] : bdegree);
        dualval += cleanupnode (n, k);
    }
    for (i=0, n=G->nodelist; i<G->ncount; i++, n++) {
        fixmark (n);
    }
    for (i=0, e=G->edgelist; i<G->ecount; i++, e++) {
        if (e->slack < 0) dualval += e->slack;
    }
    for (e = G->newedges; e; e = e->next) {
        if (e->slack < 0) dualval += e->slack;
    }
    return dualval;
}

static void fixmark (CCbmat_node *n)
{
    if (n && n->status.inpath) {
        fixmark (n->nest.parent);
        n->status.inpath = CCbmat_FALSE;
    }
}

static long cleanupnode (CCbmat_node *n, int bdegree)
{
    long dval;

    if (!n || n->status.inpath)
        return 0;
    
    dval = cleanupnode (n->nest.parent, bdegree);
    dval += n->ymat * (n->status.pseudo ? 1 - countteeth(n) : bdegree);
    CCbmat_fixmatching (n);
    n->status.inpath = CCbmat_TRUE;
    return dval;
}

static int countteeth (CCbmat_node *n)
{
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;
    int cnt = 0;

    for (pee = n->edgelist; pee; pee = pee->next) {
        pe = CCbmat_THISEDGE(pee);
        cnt += CCbmat_surfcheck (pe->ends[0], n) ?
                  pe->istooth[1] : pe->istooth[0];
    }
    return cnt;
}

static void mark_node (CCbmat_node *p)
{
    CCbmat_node *q;

    if (p->status.pseudo) {
        for (q = p->nest.child; q; q = q->nest.sibling) {
            mark_node (q);
        }
    } else {
        p->status.marked = 1;
    }
}

static int build_graph (CCbmat_graph *G, int ncount, int ecount, int *elist,
        int *elen)
{
    int i, end1, end2, w, rval = 0;

    G->ncount = ncount;
    G->ecount = ecount;
    G->nodelist = CC_SAFE_MALLOC (ncount, CCbmat_node);
    CCcheck_NULL (G->nodelist, "out of memory for nodelist");
    G->edgelist = CC_SAFE_MALLOC (ecount, CCbmat_edge);
    CCcheck_NULL (G->edgelist, "out of memory for edgelist");

    for (i=0; i<ncount; i++) {
        G->nodelist[i].edgelist = (CCbmat_edgeptr *) NULL;
    }

    for (i=0; i<ecount; i++) {
        end1 = elist[2*i];
        end2 = elist[2*i+1];
        w = elen[i];

        CCbmat_SETTHIS (&G->edgelist[i]);
        G->edgelist[i].ends[0] = G->nodelist + end1;
        G->edgelist[i].ends[1] = G->nodelist + end2;
        G->edgelist[i].weight = w;
        G->edgelist[i].ptrs[0].next = G->nodelist[end1].edgelist;
        G->nodelist[end1].edgelist = &G->edgelist[i].ptrs[0];
        G->edgelist[i].ptrs[1].next = G->nodelist[end2].edgelist;
        G->nodelist[end2].edgelist = &G->edgelist[i].ptrs[1];
    }

CLEANUP:
    return rval;
}

static void free_graph (CCbmat_graph *G)
{
    if (G) {
        CC_IFFREE (G->edgelist, CCbmat_edge);
        CC_IFFREE (G->nodelist, CCbmat_node);
    }
}

/****** from old matchmain.c ******/


static void doubleweights (CCbmat_graph *G)
{
    int i;
    CCbmat_edge *pe;

    for (i=0, pe=G->edgelist; i<G->ecount; i++, pe++) pe->weight *= 2;
}

static long minedge (CCbmat_edgeptr *pee)
{
    long m = CCbmat_MAXIWEIGHT;
    CCbmat_edge *e;

    for (; pee; pee = pee->next) {
        e = CCbmat_THISEDGE(pee);
        if (e->weight < m)
            m = e->weight;
    }
    return m;
}

static void setupmatch (CCbmat_graph *G, int bdegree, int *degrees)
{
    int i;
    CCbmat_node *p;
    CCbmat_edge *e;

    for (i=0, p = G->nodelist; i<G->ncount; i++, p++) {
        p->tree.parent = p->tree.child = p->tree.sibling =
        p->nest.parent = p->nest.child = p->nest.sibling = NULL;
        p->parentedge = NULL;
        p->deficiency = (degrees ? degrees[i]: bdegree);
        p->status.pseudo = p->status.inpath = CCbmat_FALSE;
        p->ymat = minedge (p->edgelist) / 2;  
        p->label = CCbmat_UNLABELED;
    }
    for (i = 0, e = G->edgelist; i < G->ecount; i++, e++) {
        e->slack = e->weight - e->ends[0]->ymat - e->ends[1]->ymat;
        e->x = 0;
        e->istooth[0] = e->istooth[1] = CCbmat_FALSE;
        e->status.needschange = CCbmat_FALSE;
    }
}

static void jumpsetupmatch (CCbmat_graph *G, int bdegree, int *degrees)
{
    int i, k = 0;
    CCbmat_node *p;
    CCbmat_edge *e;
    int l;

    for (i=0, p = G->nodelist; i<G->ncount; i++, p++) {
        p->tree.parent = p->tree.child = p->tree.sibling =
        p->nest.parent = p->nest.child = p->nest.sibling = NULL;
        p->parentedge = NULL;
        p->deficiency = (degrees ? degrees[i] : bdegree);
        p->status.pseudo = p->status.inpath = CCbmat_FALSE;
        p->label = CCbmat_UNLABELED;
    }
    /* flip cycles */
    for (i = 0, l = 0, e = G->edgelist; i < G->ecount; i++, e++) {
        if (e->x == CCbmat_ONE) l += flipout (e);
    }

    for (i = 0, e = G->edgelist; i < G->ecount; i++, e++) {
        e->slack = e->weight - e->ends[0]->ymat - e->ends[1]->ymat;
        if (e->x == CCbmat_TWO) {
            e->ends[0]->deficiency--;
            e->ends[1]->deficiency--;
            e->x = 1;
            k++;
        } else {
            e->x = 0;
        }
        e->istooth[0] = e->istooth[1] = CCbmat_FALSE;
        e->status.needschange = CCbmat_FALSE;
    }
    printf ("%d node graph, grabbed %d edges in fmatch\n", G->ncount, k);
}

static int flipout (CCbmat_edge *e)
{
    int val = CCbmat_ZERO;
    CCbmat_node *oldnode = e->ends[0];
    CCbmat_edge *newe;
    CCbmat_node *n;
    CCbmat_edgeptr *pe;
    int done = 0;

    newe = e;
    do {
        newe->x = val;
        if (val == CCbmat_TWO) {
            done++;
            val = CCbmat_ZERO;
        } else {
            val = CCbmat_TWO;
        }
        n = newe->ends[0];
        if (n == oldnode) n = newe->ends[1];
        oldnode = n;
        for (pe = n->edgelist; pe && (CCbmat_THISEDGE(pe)->x != CCbmat_ONE);
             pe = pe->next);
        if (pe == (CCbmat_edgeptr *) NULL) {   /* Bico 120330 */
            newe = (CCbmat_edge *) NULL; 
        } else {
            newe = CCbmat_THISEDGE(pe);
        }
    } while (newe);
    return done;
}

static void free_nodelist (CCbmat_node *nhead)
{
    CCbmat_node *n = nhead, *nnext;

    while (n) {
        nnext = n->next;
        CC_IFFREE (n, CCbmat_node);
        n = nnext;
    }
}

static int get_init_edges (CCdatagroup *dat, int ncount, int *pecount,
        int **pelist, int **pelen, int bdegree, int *degrees,
        CCrandstate *rstate)
{
    int i, k, ecount, end1, end2, w, rval = 0;
    CCedgegengroup plan;
    int *elist = (int *) NULL, *elen = (int *) NULL;

    CCedgegen_init_edgegengroup (&plan);

    if (degrees) {
        k = 0;
        for (i = 0; i < ncount; i++) {
            if (degrees[i] > k) k = degrees[i];
        }
    } else {
        k = bdegree;
    }

    if (k <= 3) {
        plan.f2match_nearest.number = 5*k;
    } else {
        plan.f2match_nearest.number = 2*k;
    }
    plan.f2match_nearest.basic = 0;
    plan.f2match_nearest.priced = 1;
    plan.tour.greedy = 1;
 
    printf ("Generating %d fractional 2-matching neighbors\n", 5*k);
    fflush (stdout);

    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &ecount,
                            &elist, 0, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");
    
    elen = CC_SAFE_MALLOC (ecount, int);
    CCcheck_NULL (elen, "out of memory for elen");
    for (i = 0; i < ecount; i++) {
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], dat);
    }

    *pecount = ecount;
    *pelist = elist;
    *pelen = elen;

CLEANUP:
    if (rval) {
        CC_IFFREE (elist, int);
        CC_IFFREE (elen, int);
    }
    return rval;
}

