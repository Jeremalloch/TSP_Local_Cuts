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
/*             ROUTINES TO COMPUTE and EXTRACT ZONES AND MOATS              */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 19, 2014 (taken from iPhone code)                           */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_geom_dual (int ncount, CCdatagroup *indat,                    */
/*      CCtsp_moatlist *moats, int boundtype, CCrandstate *rstate)          */
/*    COMPUTES zones and moats for a geometric dual  dgelist.               */
/*     -ncount is number of nodes                                           */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "fmatch.h"
#include "edgegen.h"
#include "linkern.h"
#include "tsp.h"
#include "lp.h"
#include "bigguy.h"
#include "cut.h"
#include "pq.h"
#include "cuttree.h"
#include "verify.h"

/* Code for Moats and Zones */

typedef struct CCtsp_moat {
    double pi;
    int count;
    int *nodes;
} CCtsp_moat;

typedef struct CCtsp_moatlist {
    int count;
    double *zones;
    CCtsp_moat *moats;
    double bnd;
    int *tour;
    double tourlen;
    int ready;
} CCtsp_moatlist;

void CCtsp_moatlist_init (CCtsp_moatlist *m);
void CCtsp_moatlist_free (CCtsp_moatlist *m);
int CCtsp_geom_dual (int ncount, CCdatagroup *indat, CCtsp_moatlist *moats,
    int boundtype, CCrandstate *rstate);

static void perm_bound (double *bound, int count, CCdatagroup *dat);
static int build_edges (int ncount, CCdatagroup *dat, int *ecount, int **elist,
    int **elen, CCrandstate *rstate, int just_subtour);
static int grab_plan_edges (int ncount, CCdatagroup *dat,
    CCedgegengroup *plan, int *ecount, int **elist, int **elen,
    CCrandstate *rstate);
static int add_edge_subtours (CCtsp_lp *lp, CCrandstate *rstate);
static int get_objval (CCtsp_lp *lp, double *val);
static int collect_the_moats (CCtsp_lp *lp, CCtsp_moatlist *moats);

void CCtsp_moatlist_init (CCtsp_moatlist *m)
{
    if (m) {
        m->count = 0;
        m->zones = (double *) NULL;
        m->moats = (CCtsp_moat *) NULL;
        m->bnd = 0.0;
        m->tour = (int *) NULL;
        m->tourlen = 0.0;
        m->ready = 0;
    }
}

void CCtsp_moatlist_free (CCtsp_moatlist *m)
{
    int i;

    if (m) {
        CC_IFFREE (m->zones, double);
        for (i = 0; i < m->count; i++) {
            CC_IFFREE (m->moats[i].nodes, int);
        }
        m->count = 0;
        CC_IFFREE (m->moats, CCtsp_moat);
        m->bnd = 0.0;
        CC_IFFREE (m->tour, int);
        m->tourlen = 0.0;
        m->ready = 0;
    }
}

static void perm_bound (double *bound, int ncount, CCdatagroup *dat)
{
    double bnd;
    int i;

    bnd = CCutil_dat_edgelen (ncount - 1, 0, dat);
    for (i = 1; i < ncount; i++) {
        bnd += CCutil_dat_edgelen (i-1, i, dat);
    }
    *bound = bnd;
}

static int build_edges (int ncount, CCdatagroup *dat, int *ecount, int **elist,
        int **elen, CCrandstate *rstate, int just_subtour)
{
    int rval = 0;
    CCedgegengroup plan;
    int norm;

    *ecount = 0;
    *elist  = (int *) NULL;
    *elen   = (int *) NULL;

    CCutil_dat_getnorm (dat, &norm);

    if (norm == CC_SPARSE) {
        if (dat->sparse_ecount <=  2 * ncount) {
            printf ("Use entire sparse graph as initial edge set\n");
            fflush (stdout);
            rval = CCutil_get_sparse_dat_edges (ncount, dat, ecount, elist,
                                                elen);
            if (rval) {
                fprintf (stderr, "CCutil_get_sparse_dat_edges failed\n");
                goto CLEANUP;
            }
        } else {
            int tecount;
            int *telist = (int *) NULL;
            int *telen  = (int *) NULL;

            CCedgegen_init_edgegengroup (&plan);
            plan.nearest = 4;
            rval = grab_plan_edges (ncount, dat, &plan, &tecount, &telist,
                                    &telen, rstate);
            if (rval) {
                fprintf (stderr, "grab_plan_edges failed\n"); goto CLEANUP;
            }

            rval = CCutil_sparse_strip_edges (dat, tecount, telist, telen,
                                              ecount, elist, elen);
            if (rval) {
                fprintf (stderr, "CCutil_sparse_strip_edges failed\n");
                CC_IFFREE (telist, int);
                CC_IFFREE (telen, int);
                goto CLEANUP;
            }

            CC_IFFREE (telist, int);
            CC_IFFREE (telen, int);
        }
    } else {
        CCedgegen_init_edgegengroup (&plan);
        if (ncount <= 100) {
           plan.nearest = ncount-1;
        } if (just_subtour) {
           plan.tour.greedy = 1;
           plan.f2match_nearest.number = 20 /* 4 */;
        } else {
            plan.linkern.count = 10;
            plan.linkern.quadnearest = 2;
            plan.linkern.greedy_start = 0;
            plan.linkern.nkicks = (ncount / 100) + 1;
        }

        rval = grab_plan_edges (ncount, dat, &plan, ecount, elist, elen,
                                rstate);
        if (rval) {
            fprintf (stderr, "grab_plan_edges failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }

    return rval;
}

static int grab_plan_edges (int ncount, CCdatagroup *dat,
        CCedgegengroup *plan, int *ecount, int **elist, int **elen,
        CCrandstate *rstate)
{
    int i, rval = 0;

    rval = CCedgegen_edges (plan, ncount, dat, (double *) NULL, ecount,
                            elist, 1, rstate);
    if (rval) {
        fprintf (stderr, "CCedgegen_edges failed\n"); goto CLEANUP;
    }

    *elen = CC_SAFE_MALLOC (*ecount, int);
    if (!(*elen)) {
        fprintf (stderr, "out of memory in grab_plan_edges\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i],
                                         (*elist)[(2*i) + 1], dat);
    }

CLEANUP:
    return rval;
}

#if 0
static int grab_edge_moats (CCtsp_lp *lp, double *zones, CCtsp_moatlist *moats,
    int *perm);
#endif

int CCtsp_geom_dual (int ncount, CCdatagroup *indat, CCtsp_moatlist *moats,
        int boundtype, CCrandstate *rstate)
{
    /* boundtype: 0 zones, 1 zonespenalty, 2 moats */
    int rval = 0, silent = 1, i, ecount, infeasible = 0;
    int *perm = (int *) NULL, *elist  = (int *) NULL, *elen = (int *) NULL;
    double upperbound, val, val2;
    double *node_pi = (double *) NULL, *cut_pi = (double *) NULL;
    CCdatagroup dat;
    char pname[1024];
    CCtsp_lp *lp = (CCtsp_lp *) NULL;

    /* boundtype = 0; */

    CCutil_init_datagroup (&dat);
    sprintf (pname, "noname");

    if (!moats) {
        fprintf (stderr, "CCtsp_geom_dual called without moats struct\n");
        rval = 1; goto CLEANUP;
    }

    moats->zones = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (moats->zones, "out of memory for zones");
    for (i = 0; i < ncount; i++) moats->zones[i] = 0.0;

    rval = CCutil_copy_datagroup (ncount, indat, &dat);
    CCcheck_rval (rval, "CCutil_copy_datagroup failed");

    perm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (perm, "out of memory for perm");

    rval = CCtsp_call_linkern (ncount, &dat, perm, &val, 0.0, 1, rstate);
    CCcheck_rval (rval, "CCtsp_call_linkern failed");
    moats->tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (moats->tour, "out of memory for tour");
    for (i = 0; i < ncount; i++) moats->tour[i] = perm[i];

    rval = CCutil_datagroup_perm (ncount, &dat, perm);
    CCcheck_rval (rval, "CCutil_datgroup_perm failed");

    rval = build_edges (ncount, &dat, &ecount, &elist, &elen, rstate, 1);
    CCcheck_rval (rval, "build_edges failed");

    perm_bound (&upperbound, ncount, &dat);
    moats->tourlen = upperbound;

    rval = CCtsp_init_lp (&lp, pname, -1, (char *) NULL, ncount, &dat,
               ecount, elist, elen, 0, (int *) NULL, (int *) NULL, 0,
               perm, upperbound, (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL,
               silent, rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        printf ("Initial infeasible LP\n"); rval = 1; goto CLEANUP;
    }

    rval = get_objval (lp, &val);
    CCcheck_rval (rval, "get_objval failed");

    if (boundtype == 0) {
        for (i = 0; i < ecount; i++) {
            rval = CClp_setbnd (lp->lp, i, 'U', 1000.0);
            CCcheck_rval (rval, "CClp_setbnd failed");
        }
    }

    for (i = 0; i < ncount; i++) {
        rval = CClp_change_sense (lp->lp, i, 'G');
        CCcheck_rval (rval, "CClp_change_sense failed");
    }

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed");
    if (infeasible) {
        printf ("Created infeasible LP\n"); rval = 1; goto CLEANUP;
    }
    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed");

    rval = CCtsp_get_lp_result (lp, &val2, (double *) NULL,
                     (int *) NULL, (int **) NULL, (double **) NULL,
                     (double **) NULL, &node_pi, &cut_pi);
    CCcheck_rval (rval, "CCtsp_lp_result failed");

    if (boundtype == 0) {
        printf ("LP without bounds: %f\n", val2); fflush (stdout);
    }

    if (boundtype == 0 || boundtype == 1) {
        moats->bnd = val2;
        for (i = 0; i < ncount; i++) {
            if (node_pi[i] < 0.0) {
                moats->zones[perm[i]] = 0.0;
                printf ("Negative zone: %f", node_pi[i]); fflush (stdout);
            } else {
                moats->zones[perm[i]] = node_pi[i];
            }
        }
        moats->ready = 1;
        goto CLEANUP;
    }

    printf ("Computing moats ...\n");  fflush (stdout);

    rval = CCtsp_subtour_loop (lp, 1, 0.01, rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_subtour_loop failed");
    if (infeasible) {
        printf ("Created infeasible LP\n"); rval = 1; goto CLEANUP;
    }

    rval = get_objval (lp, &val);
    CCcheck_rval (rval, "get_objval failed");
    printf ("LP with subtours: %f\n", val); fflush (stdout);
    moats->bnd = val;

    rval = add_edge_subtours (lp, rstate);
    CCcheck_rval (rval, "add_edge_subtours failed");
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    if (infeasible) {
        printf ("Created infeasible LP\n"); rval = 1; goto CLEANUP;
    }
    CCcheck_rval (rval, "CClp_opt failed");
    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed");

    rval = collect_the_moats (lp, moats);
    CCcheck_rval (rval, "collect_the_moats failed");
    moats->ready = 1;

/* HACK for edge subtours 
    rval = get_objval (lp, &val);
    CCcheck_rval (rval, "get_objval failed");
    printf ("LP with edge subtours: %f\n", val); fflush (stdout);
    moats->bnd = val;
    rval = collect_the_moats (lp, moats);
    CCcheck_rval (rval, "collect_the_moats failed");
    moats->ready = 1;
    END HACK for edge subtours */


CLEANUP:
    CCtsp_free_tsp_lp_struct (&lp);
    CCutil_freedatagroup (&dat);
    CC_IFFREE (perm, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (node_pi, double);
    CC_IFFREE (cut_pi, double);
    return rval;
}

static int add_edge_subtours (CCtsp_lp *lp, CCrandstate *rstate)
{
    int rval = 0, nadd = 0, i, ar[2], infeasible = 0;
    int ecount = lp->graph.ecount, ncount = lp->graph.ncount;
    double eps = -1000.0;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;

    for (i = 0; i < ecount; i++) {
        ar[0] = lp->graph.edges[i].ends[0];
        ar[1] = lp->graph.edges[i].ends[1];
        rval = CCtsp_array_to_subtour (&c, ar, 2, ncount);
        CCcheck_rval (rval, "CCtsp_array_to_subtour failed");
        c->next = cuts;
        cuts = c;
    }

    CCtsp_add_cuts_to_queue (lp, &cuts);
    rval = CCtsp_process_cuts (lp, &nadd, 0, 1, rstate, &eps, &infeasible);
    CCcheck_rval (rval, "CCtsp_process_cuts failed");
    if (infeasible) {
        printf ("Created infeasible LP\n"); rval = 1; goto CLEANUP;
    }

    if (nadd != ecount) {
        fprintf (stderr, "edge-cut not be added");
        rval = 1; goto CLEANUP;
    }
    
    for (i = 0; i < ecount; i++) {
        rval = CClp_setbnd (lp->lp, i, 'U', 1000.0);
        CCcheck_rval (rval, "CClp_setbnd failed");
    }

CLEANUP:
    return rval;
}

#if 0
static int grab_edge_moats (CCtsp_lp *lp, double *zones, CCtsp_moatlist *moats,
        int *perm)
{
    int i, k, j, n0, n1, rval = 0;
    int ncount = lp->graph.ncount;
    int ncuts = lp->cuts.cutcount;
    int nmoats = 0;
    double *node_pi = (double *) NULL;
    double *cut_pi = (double *) NULL;
    double cx;
    int *hit = (int *) NULL;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
                     (int *) NULL, (int **) NULL, (double **) NULL,
                     (double **) NULL, &node_pi, &cut_pi);
    CCcheck_rval (rval, "CCtsp_lp_result failed");

    for (i = 0; i < ncuts; i++) {
        cx = cut_pi[i];
        for (k = 0; k < lp->cuts.cuts[i].modcount; k++) {
            node_pi[lp->cuts.cuts[i].mods[k].node] += cx *
                (((int) lp->cuts.cuts[i].mods[k].mult) - 128);
        }
    }

    for (i = 0; i < ncount; i++) {
        zones[perm[i]] = node_pi[i];
    }

    for (i = 0; i < ncuts; i++) {
        if (cut_pi[i] > 0.01) nmoats++;
    }

    hit = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (hit, "out of memory for hit");
    for (i = 0; i < ncount; i++) hit[i] = 0;

    moats->moats = CC_SAFE_MALLOC(nmoats, CCtsp_moat);
    CCcheck_NULL (moats->moats, "out of memory for moats");
    nmoats = 0;
    for (i = 0; i < ncuts; i++) {
        j = ncuts-i-1;
        n0 = perm[lp->graph.edges[j].ends[0]];
        n1 = perm[lp->graph.edges[j].ends[1]];
        if (cut_pi[i] > 0.0 && !hit[n0] && !hit[n1]) {
            moats->moats[nmoats].pi = cut_pi[i];
            moats->moats[nmoats].count = 2;
            moats->moats[nmoats].nodes = CC_SAFE_MALLOC(2, int);
            CCcheck_NULL (moats->moats[nmoats].nodes,
                          "out of memory for moat nodes");
            moats->moats[nmoats].nodes[0] = n0;
            moats->moats[nmoats].nodes[1] = n1;
            hit[n0] = 1;
            hit[n1] = 1;
            nmoats++;
        }
    }
    moats->count = nmoats;
    printf ("have %d moats\n", nmoats);  fflush (stdout);

CLEANUP:
    CC_IFFREE (node_pi, double);
    CC_IFFREE (cut_pi, double);
    CC_IFFREE (hit, int);
    return rval;
}
#endif

static int get_objval (CCtsp_lp *lp, double *val)
{
    return CCtsp_get_lp_result (lp, val, (double *) NULL,
                     (int *) NULL, (int **) NULL, (double **) NULL,
                     (double **) NULL, (double **) NULL, (double **) NULL);
}

static int subby_cross (CCtsp_moat *s, CCtsp_moat *t, int *hit);
static int subby_equal (CCtsp_moat *s, CCtsp_moat *t, int *hit);
static int subby_uncross (int is, int it, int *scount, CCtsp_moat **slist,
        int *hit, int ncount, double *pi, int *inv);
static int subby_complement (CCtsp_moat *s, int ncount);

#define MAXMOATS 100000

static int collect_the_moats (CCtsp_lp *lp, CCtsp_moatlist *moats)
{
    int rval = 0;
    double *pi = (double *) NULL;
    double *x  = (double *) NULL;
    double *cutpi;
    double cx;
    int ncount = lp->graph.ncount;
    int ncuts = lp->cuts.cutcount;
    int nrows = ncount + ncuts;
    int i, j, k, acount, icnt, ucnt = 0, ncols;
    int *inv = (int *) NULL;
    int *hit = (int *) NULL;
    int *ar = (int *) NULL;
    CCtsp_lpcut_in c;
    CCtsp_moat **slist = (CCtsp_moat **) NULL;
    CCtsp_moat *ssupply = (CCtsp_moat *) NULL;
    int *sperm = (int *) NULL;
    int *sval = (int *) NULL;
    int smax = MAXMOATS, scount = 0;
    double len;
    int n0, n1;

    ncols = CClp_ncols (lp->lp);

    slist = CC_SAFE_MALLOC (smax, CCtsp_moat *);
    CCcheck_NULL (slist, "out of memory in print_the_subs");
    ssupply = CC_SAFE_MALLOC (smax, CCtsp_moat);
    CCcheck_NULL (ssupply, "out of memory in print_the_subs");
    for (i = 0; i < smax; i++) slist[i] = &ssupply[i];

    pi = CC_SAFE_MALLOC (nrows, double);
    CCcheck_NULL (pi, "out of memory in print_the_subs");
    inv = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (inv, "out of memory in print_the_subs");
    hit = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (hit, "out of memory in print_the_subs");
    x = CC_SAFE_MALLOC (ncols, double);
    CCcheck_NULL (x, "out of memory in print_the_subs");

    for (i = 0; i < ncount; i++)  {
        hit[i] = 0;
        inv[lp->perm[i]] = i;
    }

    rval = CClp_x (lp->lp, x);
    CCcheck_rval (rval, "CClp_x failed");

    rval = CClp_pi (lp->lp, pi);
    CCcheck_rval (rval, "CClp_pi failed");
    cutpi = pi + ncount;

    for (i = 0; i < ncuts; i++) {
        cx = cutpi[i];
        for (k = 0; k < lp->cuts.cuts[i].modcount; k++) {
            pi[lp->cuts.cuts[i].mods[k].node] += cx *
                (((int) lp->cuts.cuts[i].mods[k].mult) - 128);
        }
    }

    for (i = 0; i < ncuts; i++) {
        if (cutpi[i] >= 0.1 && lp->cuts.cuts[i].cliquecount == 1) {
            rval = CCtsp_lpcut_to_lpcut_in (&lp->cuts, &lp->cuts.cuts[i], &c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            rval = CCtsp_clique_to_array (&c.cliques[0], &ar, &acount);
            CCcheck_rval (rval, "CCtsp_clique_to_array failed");

            for (j = 0; j < acount; j++) { ar[j] = lp->perm[ar[j]]; }
            CCutil_int_array_quicksort (ar, acount);

            if (scount == smax) {
                fprintf (stderr, "Abort: more than %d moats", smax);
                rval = 1; goto CLEANUP;
            }

            slist[scount]->pi = cutpi[i];
            slist[scount]->count = acount;
            slist[scount++]->nodes = ar;
            CCtsp_free_lpcut_in (&c);

            if (acount > ncount/2) {
                rval = subby_complement (slist[scount-1], ncount);
                CCcheck_rval (rval, "subby_complement failed");
            }
        }
    }

    printf ("Check for crossings\n"); fflush (stdout);

DOGGY:

    ucnt++;
    if (ucnt > MAXMOATS) {
        fprintf (stderr, "Abort: uncrossed %d moats", MAXMOATS);
        rval = 1; goto CLEANUP;
    }

    icnt = 0;
    for (i = 0; i < scount; i++) {
        for (j = i+1; j < scount; j++) {
            if (subby_cross (slist[i], slist[j], hit)) {
                printf ("Cross X[%d,%d]\n", i, j);
                if (scount == smax) {
                    fprintf (stderr, "Abort: more than %d moats", smax);
                    rval = 1; goto CLEANUP;
                }
                rval = subby_uncross (i, j, &scount, slist, hit, ncount, pi,
                                      inv);
                CCcheck_rval (rval, "subby_uncross failed");
                goto DOGGY;
                icnt++;
            }
        }
    }
    if (icnt) printf ("\n");
    printf ("Number of intersections: %d\n", icnt); fflush (stdout);

{
    int t, z = 0, q = 0;
    double *rc = (double *) NULL;

    rc = CC_SAFE_MALLOC (ncols, double);
    CCcheck_NULL (rc, "out of memory for rc");

    rval = CClp_rc (lp->lp, rc);
    CCcheck_rval (rval, "CClp_rc failed");

    for (i = 0; i < ncols; i++) {
        if (rc[i] < -0.01) q++;
        n0 = lp->graph.edges[i].ends[0];
        n1 = lp->graph.edges[i].ends[1];
        len = (double) CCutil_dat_edgelen (n0, n1, lp->dat);
        len -= (pi[n0] + pi[n1]);
        n0 = lp->perm[n0];
        n1 = lp->perm[n1];
        for (j = 0; j < scount; j++) {
            t = 0; 
            for (k = 0; k < slist[j]->count; k++) {
                if (slist[j]->nodes[k] == n0) t++;
                if (slist[j]->nodes[k] == n1) t++;
            }
            if (t == 1) len -= (slist[j]->pi);
        }
        if (len < -1.0) {
            z++;
            printf ("Negative edge: %f\n", len); fflush (stdout);
        }
    }
    printf ("Now have %d negative edges\n", z); fflush (stdout);
    printf ("But have %d negative edges from rc\n", q); fflush (stdout);
}

    sperm = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (sperm, "out of memory in print_the_subs");
    sval = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (sval, "out of memory in print_the_subs");
    for (i = 0; i < scount; i++) {
        sperm[i] = i;
        sval[i] = slist[i]->count;
    }
    CCutil_int_perm_quicksort (sperm, sval, scount);

    for (i = 0; i < ncount; i++) {
        moats->zones[i] = pi[inv[i]];
    }

    moats->moats = CC_SAFE_MALLOC(scount, CCtsp_moat);
    CCcheck_NULL (moats->moats, "out of memory for moats");

    for (i = 0; i < scount; i++) {
        k = sperm[i];
        moats->moats[i].pi = slist[k]->pi;
        moats->moats[i].count = slist[k]->count;
        moats->moats[i].nodes = CC_SAFE_MALLOC(slist[k]->count, int);
        CCcheck_NULL (moats->moats[i].nodes, "out of memory for moats");
        for (j = 0; j < slist[k]->count; j++) {
            moats->moats[i].nodes[j] = slist[k]->nodes[j];
        }
    }
    moats->count = scount;

CLEANUP:
    
    CC_IFFREE (slist, CCtsp_moat *);
    CC_IFFREE (ssupply, CCtsp_moat);
    CC_IFFREE (pi, double);
    CC_IFFREE (x, double);
    CC_IFFREE (inv, int);
    return rval;
}

static int subby_complement (CCtsp_moat *s, int ncount)
{
    int i, k = 0,  rval = 0;
    char *hit = (char *) NULL;

    hit = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (hit, "out of memory in subby_complement");

    for (i = 0; i < ncount; i++) hit[i] = 0;
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
  
    CC_IFFREE (s->nodes, int);
    s->nodes = CC_SAFE_MALLOC (ncount - s->count, int);
    CCcheck_NULL (s->nodes, "out of memory in subby_complement");

    for (i = 0; i < ncount; i++) {
        if (!hit[i]) {
            s->nodes[k++] = i;
        }
    }
    if (k != ncount - s->count) {
        fprintf (stderr, "lost a node\n");
        rval = 1;  goto CLEANUP;
    }
    s->count = k;

CLEANUP:

    CC_IFFREE (hit, char);
    return rval;
}

static int subby_cross (CCtsp_moat *s, CCtsp_moat *t, int *hit)
{
    int i, in = 0, out = 0;
    CCtsp_moat *tmp;

    if (s->count < t->count) {
        CC_SWAP (s, t, tmp);
    }

    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
    for (i = 0; i < t->count; i++) {
        if (hit[t->nodes[i]]) {
            in = 1;
        } else {
            out = 1;
        }
    }
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 0;

    if (in && out) return 1;
    else           return 0;
}

static int subby_equal (CCtsp_moat *s, CCtsp_moat *t, int *hit)
{
    int i, val = 1;

    if (s->count != t->count) return 0;

    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
    for (i = 0; i < t->count; i++) {
        if (!hit[t->nodes[i]]) {
            val = 0;  break;
        }
    }
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 0;

    return val;
}

static int subby_uncross (int is, int it, int *scount, CCtsp_moat **slist,
       int *hit, int ncount, double *pi, int *inv)
{
    int rval = 0, i, k, iu, iv, no_u = 0, no_v = 0;
    CCtsp_moat *s = slist[is];
    CCtsp_moat *t = slist[it];
    CCtsp_moat *tmp, *u, *v;
    int icount = 0, ucount = 0;

    printf ("uncross subtour %d (%d, %f) and subtour %d (%d, %f)\n",
              is, s->count, s->pi, it, t->count, t->pi);
    fflush (stdout);

    if (s->pi < t->pi) {
        CC_SWAP (s, t, tmp);
        CC_SWAP (is, it, i);
    }

    if (s->count == 2 && t->count == 2) printf ("Edge crossing\n");

    s->pi -= t->pi;
    
    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 1;
    ucount = s->count;
    for (i = 0; i < t->count; i++) {
        if (hit[t->nodes[i]]) icount++;
        else                  ucount++;
    }

    if (ucount == ncount) {
       fprintf (stderr, "union of subtours is entire node set\n");
       rval = 1;  goto CLEANUP;
    }

    iu = *scount;
    iv = *scount + 1;
    u = slist[iu];
    v = slist[iv];

    u->pi = t->pi;
    u->count = ucount;
    u->nodes = CC_SAFE_MALLOC (u->count, int);
    CCcheck_NULL (u->nodes, "out of memory in subby_uncross");

    v->pi = t->pi;
    v->count = icount;
    v->nodes = CC_SAFE_MALLOC (v->count, int);
    CCcheck_NULL (v->nodes, "out of memory in subby_uncross");


    icount = ucount = 0;
    for (i = 0; i < s->count; i++) {
        u->nodes[ucount++] = s->nodes[i];
    }

    for (i = 0; i < t->count; i++) {
        if (hit[t->nodes[i]]) {
            v->nodes[icount++] = t->nodes[i];
        } else {
            u->nodes[ucount++] = t->nodes[i];
        }
    }
    if (ucount != u->count || icount != v->count) {
        fprintf (stderr, "bad counts in uncross\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < s->count; i++) hit[s->nodes[i]] = 0;

    for (i = 0; i < *scount; i++) {
        if (subby_equal (slist[i], u, hit)) {
            slist[i]->pi += u->pi;
            no_u = 1;
            break;
        }
    }

    if (v->count == 1) {
        k = inv[v->nodes[0]];
        pi[k] += v->pi; 
        no_v = 1;
    } else {
        for (i = 0; i < *scount; i++) {
            if (subby_equal (slist[i], v, hit)) {
                slist[i]->pi += v->pi;
                no_v = 1;
                break;
            }
        }
    }

    CC_IFFREE (t->nodes, int);
    if (no_u) {
       CC_IFFREE (u->nodes, int);
    }
    if (no_v) {
       CC_IFFREE (v->nodes, int);
    }

    if (no_u && no_v) {
        CC_SWAP (slist[it], slist[*scount-1], tmp);
        (*scount)--;
    } else if (no_u) {
        CC_SWAP (slist[it], slist[iv], tmp);
    } else if (no_v) {
        CC_SWAP (slist[it], slist[iu], tmp);
    } else {
        CC_SWAP (slist[it], slist[iv], tmp);
        (*scount)++;
    }

    if (s->pi < 0.1) {
        CC_SWAP (slist[is], slist[*scount-1], tmp);
        (*scount)--;
    }

CLEANUP:

    for (i = 0; i < ncount; i++) hit[i] = 0;

    return rval;
}

