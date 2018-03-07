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
/*                     ELIMINATE EDGES BY DEPTH-1 BRANCHING                 */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: January 23, 2011                                                  */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

static char *rootfname    = (char *) NULL;
static char *masterfname = (char *) NULL;
static int seed = 0;
static double initial_ub = CCtsp_LP_MAXDOUBLE;

#define EDGEBATCH 10

int
    main (int ac, char **av);

static int
    add_edges (CCtsp_lp *lp, int acount, int *alist, CCdatagroup *dat,
        double dtarget, CCbigguy target, int *results),
    set_edge_to_one (CCtsp_lp *lp, int col, int n1, int n2, double dtarget,
        CCbigguy target, int *fixit),
    build_working_lp (char *fname, CCtsp_lp **lp, CCdatagroup *dat,
        int *ptour, CCrandstate *rstate),
    lpedges_sanity_check (CCtsp_lp *lp, int fullcount, int *fullelist),
    edge_in_lpgraph (CCtsp_lpgraph *g, int n0, int n1, int *ind),
    lp_value (CCtsp_lp *lp, double *val),
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, k, ncount, cnt, rval = 0, infeasible = 0;
    int *ptour = (int *) NULL;
    int *inv_ptour = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    CCtsp_lp *worklp = (CCtsp_lp *) NULL;
    CCbigguy bound, target;
    double val, szeit, dtarget;
    CClp_warmstart *warmstart = (CClp_warmstart *) NULL;
    int fullcount = 0, nwinners = 0;
    int *fullelist = (int *) NULL, *fullperm = (int *) NULL;
    int *testlist = (int *) NULL, *results = (int *) NULL;
    double *rc = (double *) NULL;

    CCutil_init_datagroup (&dat);
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

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    inv_ptour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (inv_ptour, "out of memory for inv_ptour");

    for (i = 0; i < ncount; i++) {
        inv_ptour[ptour[i]] = i;
    }

    rval = CCtsp_init_lp (&rootlp, (char *) NULL, -1, rootfname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 0, &rstate,
               &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "init LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Upper bound from LP struct: %.0f\n", rootlp->upperbound);
    fflush (stdout);

    if (initial_ub < rootlp->upperbound) {
        printf ("Resetting upperbound: %.0f\n", initial_ub);
        fflush (stdout);
        rootlp->upperbound = initial_ub;
    }

    rval = CCtsp_exact_price (rootlp, &bound, 0, 0, 0, 0, 0);
    CCcheck_rval (rval, "CCtsp_exact_price failed");
    printf ("Exact lower bound: %.6f\n", CCbigguy_bigguytod (bound));
    fflush (stdout);

    target = CCbigguy_dtobigguy (rootlp->upperbound);
    CCbigguy_sub (&target, CCbigguy_ONE);
    printf ("Elimination Target: %f\n", CCbigguy_bigguytod (target));
    fflush (stdout);

    dtarget = rootlp->upperbound - 0.9999;
    printf ("Exact Test for each LP Value Above: %f\n", dtarget);
    fflush (stdout);

    rval = lp_value (rootlp, &val);
    CCcheck_rval (rval, "lp_value failed");
    printf ("LP Value: %f\n", val); fflush (stdout);

    rval = CClp_get_warmstart (rootlp->lp, &warmstart);
    CCcheck_rval (rval, "CClp_get_warmstart failed");

    rval = CCtsp_reduced_cost_all (rootlp, &fullcount, &fullelist, &rc,
                                   (int *) NULL);
    CCcheck_rval (rval, "CCtsp_reduced_cost_all failed");

    fullperm = CC_SAFE_MALLOC (fullcount, int);
    CCcheck_NULL (fullperm, "out of memory for fullperm");
    for (i = 0; i < fullcount; i++) {
        fullperm[i] = i;
        rc[i] = -rc[i];
    }
    CCutil_double_perm_quicksort (fullperm, rc, fullcount);
    for (i = 0; i < fullcount; i++) rc[i] = -rc[i];

#if 0
    {
        /* Try a random ordering rather than RC-ordering */
        int temp;

        for (i = fullcount; i > 1; i--) {
            k = CCutil_lprand (&rstate) % i;
            CC_SWAP (fullperm[i - 1], fullperm[k], temp);
        }
    }
#endif

    rval = lpedges_sanity_check (rootlp, fullcount, fullelist);
    CCcheck_rval (rval, "lpedges_sanity_check failed");

    rval = build_working_lp (rootfname, &worklp, &dat, ptour, &rstate);
    CCcheck_rval (rval, "build_working_lp failed");

    testlist = CC_SAFE_MALLOC (2*EDGEBATCH, int);
    CCcheck_NULL (testlist, "out of memory for testlist");
    results = CC_SAFE_MALLOC (EDGEBATCH, int);
    CCcheck_NULL (results, "out of memory for results");

    k = 0;
    while (k < fullcount) {
        for (cnt = 0; cnt < EDGEBATCH && k < fullcount; cnt++, k++) {
             testlist[2*cnt]   = fullelist[2*fullperm[k]];
             testlist[2*cnt+1] = fullelist[2*fullperm[k]+1];
        }
        rval = add_edges (worklp, cnt, testlist, &dat, dtarget, target,
                          results);
        CCcheck_rval (rval, "add_edges failed");
        for (i = 0; i < cnt; i++) {
            if (results[i] == 1) {
                nwinners++;
                if(nwinners % 50 == 0) printf ("\n");
            }
        }
    }
    printf ("\n");
    printf ("Final Working LP Edge Count = %d,  Column Count = %d\n",
                worklp->graph.ecount, CClp_ncols (worklp->lp));
    printf ("Eliminated Edges: %d\n", nwinners);
    printf ("Remaining Edges:  %d\n", fullcount - nwinners);

CLEANUP:

    CC_IFFREE (ptour, int);
    CC_IFFREE (inv_ptour, int);
    CC_IFFREE (fullelist, int);
    CC_IFFREE (testlist, int);
    CC_IFFREE (rc, double);
    CCutil_freedatagroup (&dat);
    if (warmstart) CClp_free_warmstart (&warmstart);
    if (rootlp) CCtsp_free_tsp_lp_struct (&rootlp);
    if (worklp) CCtsp_free_tsp_lp_struct (&worklp);

    return rval;
}

static int add_edges (CCtsp_lp *lp, int acount, int *alist, CCdatagroup *dat,
        double dtarget, CCbigguy target, int *results)
{
    int i, k, n0, n1, outcnt, ind, fixit, rval = 0;
    int ecount = lp->graph.ecount;
    CCtsp_predge *prlist = (CCtsp_predge *) NULL;
    int *delmarks = (int *) NULL, *lindex = (int *) NULL;

    prlist = CC_SAFE_MALLOC (acount, CCtsp_predge);
    CCcheck_NULL (prlist, "out of memory for prlist");
    lindex = CC_SAFE_MALLOC (acount, int);
    CCcheck_NULL (prlist, "out of memory for lindex");

    for (outcnt = 0, i = 0; i < acount; i++) {
        n0 = alist[2*i];
        n1 = alist[2*i+1];
        if (edge_in_lpgraph (&lp->graph, n0, n1, &ind) == 0) {
            prlist[outcnt].ends[0] = n0;
            prlist[outcnt].ends[1] = n1;
            prlist[outcnt].len = (double) CCutil_dat_edgelen (n0, n1, dat);
            lindex[i] = ecount+outcnt;
            outcnt++;
        } else {
            lindex[i] = ind;
        }
    }
    if (outcnt) {
        rval = CCtsp_add_vars_to_lp (lp, prlist, outcnt);
        CCcheck_rval (rval, "CCtsp_add_vars_to_lp failed");
    }

    for (i = 0; i < acount; i++) {
        rval = set_edge_to_one (lp, lindex[i], alist[2*i], alist[2*i+1],
                                dtarget, target, &fixit);
        CCcheck_rval (rval, "set_edge_to_one failed");
        results[i] = fixit;
    }
  
    if (outcnt) {
        delmarks = CC_SAFE_MALLOC (ecount + EDGEBATCH, int);
        CCcheck_NULL (delmarks, "out of memory for delmarks");
        for (i = 0; i < ecount + EDGEBATCH; i++) delmarks[i] = 0;

        for (k = 0; k < outcnt; k++) delmarks[ecount+k] = 1;
        rval = CClp_delete_set_of_columns (lp->lp, delmarks);
        CCcheck_rval (rval, "CClp_delete_set_of_columns failed");

        lp->graph.ecount = ecount;
        rval = CCtsp_build_lpadj (&lp->graph, 0, lp->graph.ecount);
        CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
    }

CLEANUP:

    CC_IFFREE (prlist, CCtsp_predge );
    CC_IFFREE (delmarks, int);
    return rval;
}

static int set_edge_to_one (CCtsp_lp *lp, int col, int n0, int n1,
        double dtarget, CCbigguy target, int *fixit)
{
    int rval = 0, infeasible = 0;
    double val, szeit;
    CCbigguy bound;

    *fixit = 0;
    if (lp->graph.edges[col].fixed) goto CLEANUP;

    rval = CClp_setbnd (lp->lp, col, 'L', 1.0);
    CCcheck_rval (rval, "CClp_setbnd failed");

    szeit = CCutil_zeit();
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed");
    if (infeasible) {
        rval = CClp_setbnd (lp->lp, col, 'L', 0.0);
        CCcheck_rval (rval, "CClp_setbnd failed");
        goto CLEANUP;
    }
    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed");
    CCtsp_free_bigdual (&lp->exact_dual);

    rval = lp_value (lp, &val);
    CCcheck_rval (rval, "lp_value failed");

    printf ("LP Solve Time: %.2f   Val = %f  Gap = %f\n",
           CCutil_zeit() - szeit, val, dtarget - val);
    fflush (stdout);

    if (val > dtarget) {
        lp->graph.edges[col].fixed = 1;
        lp->fixededges[2*lp->nfixededges] = n0;
        lp->fixededges[2*lp->nfixededges+1] = n1; 
        lp->nfixededges++;
        rval = CCtsp_exact_price (lp, &bound, 0, 0, 0, 0, 0);
        CCcheck_rval (rval, "CCtsp_exact_price failed");
        lp->graph.edges[col].fixed = 0;
        lp->nfixededges--;
        printf ("Val = %f  Exact = %.6f\n", val, CCbigguy_bigguytod (bound));
        fflush (stdout);
        if (CCbigguy_cmp (bound, target) > 0) {
            *fixit = 1;
            printf ("."); fflush (stdout);
        }
    }

    rval = CClp_setbnd (lp->lp, col, 'L', 0.0);
    CCcheck_rval (rval, "CClp_setbnd failed");

CLEANUP:
    return rval;
}

static int build_working_lp (char *fname, CCtsp_lp **lp, CCdatagroup *dat,
        int *ptour, CCrandstate *rstate)
{
    int rval = 0, silent = 1, infeasible = 0;
    void *tmp_ptr;

    rval = CCtsp_init_lp (lp, (char *) NULL, -1, fname, 0,
               dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, silent, rstate,
               &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    /* Create space to add an edge to fixded-edge list */

    tmp_ptr = (void *) (*lp)->fixededges;
    rval = CCutil_reallocrus_count (&tmp_ptr,
                  2 * ((*lp)->nfixededges + 1), sizeof (int));
    CCcheck_rval (rval, "failed to realloc for fixededges");
    (*lp)->fixededges = (int *) tmp_ptr;

CLEANUP:

    return rval;
}

static int lpedges_sanity_check (CCtsp_lp *lp, int fullcount, int *fullelist)
{
    CCtsp_lpgraph *g = &lp->graph;
    CCtsp_lpedge *e;
    int i, k, n, n0, n1, got = 0;
    int rval = 0;

    printf ("Full Edges Valid: %d (%d edges)\n", lp->full_edges_valid,
             lp->fullcount);
    fflush (stdout);

    for (i = 0; i < g->ecount; i++) {
        e = &(g->edges[i]);
        if (e->ends[0] >= e->ends[1]) {
            printf ("End in wrong order\n");
            rval = 1;  goto CLEANUP;
        }
        n = e->ends[0];
        for (k = 0; k < lp->fulladj[n].deg; k++) {
            if (lp->fulladj[n].list[k].end == e->ends[1]) break;
        }
        if (k == lp->fulladj[n].deg) {
            printf ("Missing LP edge in full adj list\n");
            rval = 1;  goto CLEANUP;
        }
    }

    printf ("All %d LP edges, present and accounted for Sir.\n", i);
    fflush (stdout);

    for (i = 0; i < fullcount; i++) {
        n0 = fullelist[2*i];
        n1 = fullelist[2*i+1];
        if (n0 >= n1) {
            printf ("End in wrong order\n");
            rval = 1;  goto CLEANUP;
        }
        got += edge_in_lpgraph (g, n0, n1, (int *) NULL);
    }

    printf ("Found %d of the LP edges Sir.\n", got);
    fflush (stdout);

CLEANUP:

    return rval;
}

static int edge_in_lpgraph (CCtsp_lpgraph *g, int n0, int n1, int *ind)
{
    int k, yesno = 0;

    for (k = 0; k < g->nodes[n0].deg; k++) {
        if (g->nodes[n0].adj[k].to == n1) {
            if (ind) *ind = g->nodes[n0].adj[k].edge;
            yesno = 1; break;
        }
    }
    return yesno;
}

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval;

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
    return rval;
}


static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "M:s:u:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'M':
            masterfname  = boptarg;
            break;
        case 's':
            seed = atoi(boptarg);
            break;
        case 'u':
            initial_ub = atof (boptarg);
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
    fprintf (stderr, "   -M f  specify a master file (required)\n");
    fprintf (stderr, "   -s #  seed\n");
    fprintf (stderr, "   -u #  upperbound on optimal tour length\n");
}


