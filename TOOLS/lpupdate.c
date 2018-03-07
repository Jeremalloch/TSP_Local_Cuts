/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                         REPLACE EDGE SET OF AN LP                        */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: February 21, 2015                                                 */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

typedef struct node {
    int *neighbors;
    int deg;
} node;

typedef struct edge {
    int ends[2];
    int len;
} edge;

typedef struct graph {
    int ncount;
    int ecount;
    node *nodelist;
    int *neighborspace;
} graph;

static char *rootfname    = (char *) NULL;
static char *masterfname = (char *) NULL;
static char *edgefname = (char *) NULL;
static char *fixfname = (char *) NULL;
static int seed = 0;

int main (int ac, char **av);
static int update_edges (CCtsp_lp *rootlp, int ecount, int *elist, int fcount,
    int *flist);
static int update_edges_full (CCtsp_lp *lp, int ecount, int *elist, graph *E,
    int fcount, int *flist, graph *F);
static int find_edge_full (CCtsp_lp *lp, int from, int to);
static void init_graph (graph *G);
static void free_graph (graph *G);
static int edge_in_graph (int a, int b, graph *G);
static int buildgraph (graph *G, int ncount, int ecount, int *elist);
static int lp_value (CCtsp_lp *lp, double *val);
static int parseargs (int ac, char **av);
static void usage (char *fname);

int main (int ac, char **av)
{
    int  i, ncount, ecount = 0, fcount = 0, rval = 0;
    int infeasible = 0;
    int *ptour = (int *) NULL;
    int *elist = (int *) NULL, *elen = (int *) NULL, *pelist = (int *) NULL;
    int *flist = (int *) NULL, *flen = (int *) NULL, *pflist = (int *) NULL;
    int *invperm = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dominopool = (CCtsp_lpcuts *) NULL;
    double val, szeit;
    graph E, F, *pE = (graph *) NULL, *pF = (graph *) NULL;
    
    CCutil_init_datagroup (&dat);
    init_graph (&E);
    init_graph (&F);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname) {
        fprintf (stderr, "Must specify a master file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (!edgefname && !fixfname) {
        fprintf (stderr, "Must specify an edge file and/or fixed file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &dominopool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");

    rval = CCtsp_init_lp (&rootlp, (char *) NULL, -1, rootfname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               pool, dominopool, 0, &rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    if (!edgefname && !rootlp->full_edges_valid) {
        fprintf (stderr, "Need edge file or LP with valid full edge set\n");
        rval = 1; goto CLEANUP;
    }

    rval = lp_value (rootlp, &val);
    CCcheck_rval (rval, "lp_value failed");
    printf ("LP Value: %f\n", val); fflush (stdout);

    CC_MALLOC (invperm, ncount, int);
    for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;

    if (edgefname) {
        rval = CCutil_getedgelist(ncount, edgefname, &ecount, &pelist,
                                  &elen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");

        /* permute names to match masterfile */
        elist = CC_SAFE_MALLOC (2*ecount, int);
        CCcheck_NULL (elist, "out of memory for elist");
        for (i = 0; i < ecount; i++) {
            elist[2*i]   = invperm[pelist[2*i]];
            elist[2*i+1] = invperm[pelist[2*i+1]];
        }
        CC_IFFREE (pelist, int);
        buildgraph (&E, ncount, ecount, elist);
        pE = &E;
    }

    if (fixfname) {
        rval = CCutil_getedgelist(ncount, fixfname, &fcount, &pflist, &flen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");

        flist = CC_SAFE_MALLOC (2*fcount, int);
        CCcheck_NULL (flist, "out of memory for flist");
        for (i = 0; i < fcount; i++) {
            flist[2*i]   = invperm[pflist[2*i]];
            flist[2*i+1] = invperm[pflist[2*i+1]];
        }
        CC_IFFREE (pflist, int);
        buildgraph (&F, ncount, fcount, flist);
        pF = &F;
    }

    if (rootlp->full_edges_valid) {
        rval = update_edges (rootlp, ecount, elist, fcount, flist);
        CCcheck_rval (rval, "update_edges failed\n");
    } else {
        rval = update_edges_full (rootlp, ecount, elist, pE, fcount, flist, pF);
        CCcheck_rval (rval, "update_edges failed full\n");
    }

    printf ("Run pricing loop ...\n"); fflush (stdout);
    rval = CCtsp_pricing_loop (rootlp, &val, 0, &rstate);
    CCcheck_rval (rval, "CCtsp_pricing_loop failed");

    printf ("Write updated LP ...\n"); fflush (stdout);
    rval = CCtsp_write_probfile_sav (rootlp);
    CCcheck_rval (rval, "CCtsp_write_probfile_save failed");

CLEANUP:
    if (pool) { CCtsp_free_cutpool (&pool); }
    if (dominopool) { CCtsp_free_cutpool (&dominopool); }
    CC_IFFREE (ptour, int);
    CC_IFFREE (pelist, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (pflist, int);
    CC_IFFREE (flist, int);
    CC_IFFREE (flen, int);
    CC_IFFREE (invperm, int);
    CCutil_freedatagroup (&dat);
    free_graph (&E);
    free_graph (&F);
    if (rootlp) CCtsp_free_tsp_lp_struct (&rootlp);

    return rval;
}

static int update_edges_full (CCtsp_lp *lp, int ecount, int *elist, graph *E,
        int fcount, int *flist, graph *F)
{
    CCtsp_lpgraph *g = &lp->graph;
    int rval = 0, ncount = g->ncount, i, end0, end1, ek, itmp, exfix = 0;
    int *elen = (int *) NULL;
    CCtsp_predge *prlist = (CCtsp_predge *) NULL;

    if (lp->nfixededges) {
        printf ("Not set up to for full update when LP already has fixed\n");
        rval = 1;  goto CLEANUP;
    }

    if (E && F) {
        for (i = 0; i < fcount; i++) {
            if (!edge_in_graph (flist[2*i], flist[2*i+1], E)) {
                fprintf (stderr, "ERROR: fixed edge not in edge list\n");
                rval = 1; goto CLEANUP;
            }
        }
    }

    /* add any missing fixed edge to the LP */

    if (fcount) {
        CC_MALLOC (prlist, fcount, CCtsp_predge);
        for (i = 0; i < fcount; i++) {
            end0 = flist[2*i];  end1 = flist[2*i+1];
            if (end0 > end1) { CC_SWAP (end0, end1, itmp); }
            ek = CCtsp_find_edge (g, end0, end1);
            if (ek == -1) {
                printf ("LP missing fixed edge (%d,%d)\n", end0, end1);
                fflush (stdout);
                prlist[exfix].ends[0] = end0;
                prlist[exfix].ends[1] = end1;
                prlist[exfix].len = CCutil_dat_edgelen (end0, end1, lp->dat);
                exfix++;
            }
        }
        if (exfix) {
            printf ("Adding %d fixed edges to the LP\n", exfix);
            fflush (stdout);
            rval = CCtsp_add_vars_to_lp (lp, prlist, exfix);
            CCcheck_rval (rval, "CCtsp_add_vars_to_lp failed");
        }
    }

    /* build the adj structure for the new edge list */

    CC_MALLOC (elen, ecount, int);
    for (i = 0; i < ecount; i++) {
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], lp->dat);
    }

    rval = CCtsp_edgelist_to_genadj (ncount, ecount, elist, elen,
                           &(lp->fulladj), &(lp->fulladjspace));
    CCcheck_rval (rval, "CCtsp_edgelist_to_genadj failed");
    lp->fullcount = ecount;
    lp->full_edges_valid = 1;

    /* set fixed edges */

    if (fcount) {
        CC_MALLOC (lp->fixededges, 2*fcount, int);
        for (i = 0; i < fcount; i++) {
            lp->fixededges[2*i]   = flist[2*i];
            lp->fixededges[2*i+1] = flist[2*i+1];
        }
        lp->nfixededges = fcount;
        printf ("Fixed %d edges to 1\n", fcount); fflush (stdout);
    }

    printf ("Call CCtsp_eliminate_rebuild ...\n"); fflush (stdout);
    rval = CCtsp_eliminate_rebuild (lp, 0);
    CCcheck_rval (rval, "CCtsp_eliminate_rebuild failed");


CLEANUP:
    CC_IFFREE (prlist, CCtsp_predge);
    CC_IFFREE (elen, int);
    return rval;
}

static int update_edges (CCtsp_lp *lp, int ecount, int *elist, int fcount,
        int *flist)
{
    CCtsp_lpgraph *g = &lp->graph;
    int rval = 0, ncount = g->ncount, i, end0, end1, ek, itmp;
    edge *newlist = (edge *) NULL;
    int *newfix = (int *) NULL, nfixed = 0, nremain = 0;
    CCtsp_genadj *adj = (CCtsp_genadj *) NULL;
    CCtsp_genadjobj *adjspace = (CCtsp_genadjobj *) NULL, *pa;

    printf ("Starting LP Edge Count: %d\n", g->ecount);
    printf ("Starting LP Full Count: %d\n", lp->fullcount);
    fflush (stdout);

    if (ecount) {
        CC_MALLOC (newlist, ecount, edge);
        CC_MALLOC (adj, ncount, CCtsp_genadj);
        for (i = 0; i < ncount; i++) adj[i].deg = 0;

        for (i = 0; i < ecount; i++) {
            end0 = elist[2*i];  end1 = elist[2*i+1];
            if (end0 > end1) { CC_SWAP (end0, end1, itmp); }
            if (find_edge_full (lp, end0, end1)) {
                newlist[nremain].ends[0] = end0;
                newlist[nremain].ends[1] = end1;
                newlist[nremain].len = CCutil_dat_edgelen (end0, end1, lp->dat);
                adj[end0].deg++;
                nremain++;
            }
        }
        printf ("Remaining edges: %d\n", nremain); fflush (stdout);
    }

    if (fcount) {
        CC_MALLOC (newfix, 2*fcount, int);
        for (i = 0; i < fcount; i++) {
            end0 = flist[2*i];  end1 = flist[2*i+1];
            if (end0 > end1) { CC_SWAP (end0, end1, itmp); }
            ek = CCtsp_find_edge (g, end0, end1);
            if (ek == -1) {
                printf ("Warning: fixed edge (%d,%d) not in lp\n", end0, end1);
                fflush (stdout);
            } else if (lp->graph.edges[ek].fixed == 0) {
                newfix[2*nfixed]   = end0;
                newfix[2*nfixed+1] = end1;
                nfixed++;
            }
        }
        printf ("Set %d new fixed edges\n", nfixed); fflush (stdout);
    }

    if (nremain) {
        CC_MALLOC (adjspace, nremain, CCtsp_genadjobj);
        pa = adjspace;
        for (i = 0; i < ncount; i++) {
            adj[i].list = pa;
            pa += adj[i].deg;
            adj[i].deg = 0;
        }

         for (i = 0; i < nremain; i++) {
            end0 = newlist[i].ends[0];
            adj[end0].list[adj[end0].deg].end = newlist[i].ends[1];
            adj[end0].list[adj[end0].deg].len = newlist[i].len;
            adj[end0].deg++;
        }

        CC_IFFREE (lp->fulladjspace, CCtsp_genadjobj);
        CC_IFFREE (lp->fulladj, CCtsp_genadj);
        lp->fullcount = nremain;
        lp->fulladjspace = adjspace;
        lp->fulladj = adj;
    }

    if (nfixed) {
        void *tmp_ptr = (void *) lp->fixededges;
        rval = CCutil_reallocrus_count (&tmp_ptr,
                  2 * (lp->nfixededges + nfixed), sizeof (int));
        CCcheck_rval (rval, "out of memory for fixed edges");
        lp->fixededges = (int *) tmp_ptr;

        for (i = 0; i < nfixed; i++) {
            lp->fixededges[2*lp->nfixededges]   = newfix[2*i];
            lp->fixededges[2*lp->nfixededges+1] = newfix[2*i+1];
            lp->nfixededges++;
        }
    }

    printf ("Call CCtsp_eliminate_rebuild ...\n"); fflush (stdout);
    rval = CCtsp_eliminate_rebuild (lp, 0);
    CCcheck_rval (rval, "CCtsp_eliminate_rebuild failed");
        
CLEANUP:
    CC_IFFREE (newlist, edge);
    CC_IFFREE (newfix, int);
    return rval;
}

static int find_edge_full (CCtsp_lp *lp, int from, int to)
{
    int i;
    CCtsp_genadjobj *a;

    if (from > to) { CC_SWAP (from, to, i); }
    a = lp->fulladj[from].list;
    for (i = lp->fulladj[from].deg-1; i >= 0; i--) {
        if (a[i].end == to) {
            return 1;
        }
    }
    return 0;
}

static void init_graph (graph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->nodelist = (node *) NULL;
        G->neighborspace = (int *) NULL;
    }
}

static void free_graph (graph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->neighborspace, int);
        init_graph (G);
    }
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist)
{
    int rval = 0, i, *p;
    node *n0, *n1;

    init_graph (G);
    G->ncount = ncount;
    G->ecount = ecount;
    CC_MALLOC (G->nodelist, ncount, node);
    CC_MALLOC (G->neighborspace, 2*ecount, int);
    for (i = 0; i < ncount; i++) G->nodelist[i].deg = 0;

    for (i = 0; i < ecount; i++) {
        G->nodelist[elist[2*i]].deg++;
        G->nodelist[elist[2*i+1]].deg++;
    }

    p = G->neighborspace;
    for (i = 0; i < ncount; i++) {
        G->nodelist[i].neighbors = p;
        p += G->nodelist[i].deg;
        G->nodelist[i].deg = 0;
    }

    for (i = 0; i < ecount; i++) {
        n0 = &G->nodelist[elist[2*i]];
        n1 = &G->nodelist[elist[2*i+1]];
        n0->neighbors[n0->deg++] = elist[2*i+1];
        n1->neighbors[n1->deg++] = elist[2*i];
    }

CLEANUP:
    return rval;
}

static int edge_in_graph (int a, int b, graph *G)
{
    int i;

    for (i = 0; i < G->nodelist[a].deg; i++) {
        if (G->nodelist[a].neighbors[i] == b) return 1;
    }
    return 0;
}

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval = 0, infeasible = 0;

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed");
    if (infeasible) {
        fprintf (stderr, "Infeasible LP\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "e:f:M:s:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'e':
            edgefname = boptarg;
            break;
        case 'f':
            fixfname = boptarg;
            break;
        case 'M':
            masterfname = boptarg;
            break;
        case 's':
            seed = atoi(boptarg);
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        rootfname = av[boptind++];
    } else{
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] savfile_root\n", fname);
    fprintf (stderr, "   -e f  specify an edge file\n");
    fprintf (stderr, "   -f f  specify an fixed edge file\n");
    fprintf (stderr, "   -M f  specify a master file (required)\n");
    fprintf (stderr, "   -s #  seed\n");
}


