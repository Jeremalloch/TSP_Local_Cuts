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
/*  Date: November 5, 2013                                                  */
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
static char *edgefname = (char *) NULL;
static char *fixfname = (char *) NULL;
static char *extrafname = (char *) NULL;
static int seed = 0;

int
    main (int ac, char **av);

static int
    replace_all_edges (CCtsp_lp *oldlp, CCtsp_lp **newlp,
        int ecount, int *elist, int *elen, CCrandstate *rstate, int silent,
        CCtsp_lpcuts *pool, CCtsp_lpcuts *dominopool),
    add_extra_cuts (CCtsp_lp *lp, CCtsp_lp *extralp, CCrandstate *rstate,
        int silent),
    lpedges_sanity_check (CCtsp_lp *lp, int fullcount, int *fullelist),
    edge_in_lpgraph (CCtsp_lpgraph *g, int n0, int n1, int *ind),
    lp_value (CCtsp_lp *lp, double *val),
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, k, ncount, ecount = 0, fcount = 0, rval = 0, silent = 1;
    int infeasible = 0;
    int *ptour = (int *) NULL, *elist = (int *) NULL, *elen = (int *) NULL;
    int *flist = (int *) NULL, *flen = (int *) NULL;
    int *invperm = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    CCtsp_lp *worklp = (CCtsp_lp *) NULL, *extralp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dominopool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *epool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *edominopool = (CCtsp_lpcuts *) NULL;
    double val, szeit;
    
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

    if (!edgefname && !fixfname && !extrafname) {
        fprintf (stderr, "Must specify an edge, fixed, or extra file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    if (edgefname) {
        rval = CCutil_getedgelist(ncount, edgefname, &ecount, &elist, &elen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
    } else if (fixfname) {
        rval = CCutil_getedgelist(ncount, fixfname, &fcount, &flist, &flen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
    } else {
        rval = CCtsp_init_lp (&extralp, (char *) NULL, -1, extrafname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 0, &rstate,
               &infeasible);
        CCcheck_rval (rval, "CCtsp_init_lp failed");
        if (infeasible) {
            fprintf (stderr, "initial LP is infeasible\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &epool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &edominopool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");

    rval = CCtsp_init_lp (&rootlp, (char *) NULL, -1, rootfname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               epool, edominopool, 0, &rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    rval = lp_value (rootlp, &val);
    CCcheck_rval (rval, "lp_value failed");
    printf ("LP Value: %f\n", val); fflush (stdout);

    if (ecount > 0) {
        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed");
        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &dominopool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");

        rval = replace_all_edges (rootlp, &worklp, ecount, elist, elen,
                                 &rstate, silent, pool, dominopool);
        CCcheck_rval (rval, "replace_all_edge failed");

        rval = lpedges_sanity_check (worklp, ecount, elist);
        CCcheck_rval (rval, "lpedges_sanity_check failed");

        printf ("Optimizing LP with new edge set ...\n"); fflush (stdout);
        rval = lp_value (worklp, &val);
        CCcheck_rval (rval, "lp_value failed");
        printf ("Copied LP Value: %f\n", val); fflush (stdout);

        rval = CCtsp_write_probfile_sav (worklp);
        CCcheck_rval (rval, "CCtsp_write_probfile_sav\n");
    }

    if (fcount > 0) {
        if (rootlp->nfixededges > 0) {
            printf ("TERMINATE: Already have fixed edges in the LP\n");
            rval = 1; goto CLEANUP;
        }
  
        printf ("FIX %d edges to value 1\n", fcount);
        fflush (stdout);

        CC_MALLOC (invperm, ncount, int);
        for (i = 0; i < ncount; i++) invperm[rootlp->perm[i]] = i;

        rootlp->fixededges = CC_SAFE_MALLOC (2*fcount, int);
        CCcheck_NULL (rootlp->fixededges,
                     "out of memory for rootlp->fixededges");
        for (i = 0; i < fcount; i++) {
            rootlp->fixededges[2*i]   = invperm[flist[2*i]];
            rootlp->fixededges[2*i+1] = invperm[flist[2*i+1]];
        }
        rootlp->nfixededges = fcount;

        for (i = 0; i < rootlp->nfixededges; i++) {
            k = CCtsp_find_edge (&(rootlp->graph), rootlp->fixededges[2*i],
                                                   rootlp->fixededges[2*i+1]);
            if (k != -1) {
                rval = CClp_setbnd (rootlp->lp, k, 'L', 1.0);
                rootlp->graph.edges[k].fixed = 1;
            } else {
                printf ("ERROR: Fixed edge (%d, %d) is not in LP\n",
                         rootlp->fixededges[2*i], rootlp->fixededges[2*i+1]);
                fflush (stdout);
                rval = 1; goto CLEANUP;
            }
        }

        printf ("Optimizing LP with fixed edges ...\n"); fflush (stdout);
        rval = lp_value (rootlp, &val);
        CCcheck_rval (rval, "lp_value failed");
        printf ("Copied LP Value: %f\n", val); fflush (stdout);

        rval = CCtsp_write_probfile_sav (rootlp);
        CCcheck_rval (rval, "CCtsp_write_probfile_sav\n");
    }

    if (extralp) {
        printf ("Add cuts from extra LP file\n"); fflush (stdout);
        rval = add_extra_cuts (rootlp, extralp, &rstate, silent);
        CCcheck_rval (rval, "add_extra_cuts");

        rval = lp_value (rootlp, &val);
        CCcheck_rval (rval, "lp_value failed");
        printf ("Final LP Value: %f\n", val); fflush (stdout);
        rval = CCtsp_write_probfile_sav (rootlp);
        CCcheck_rval (rval, "CCtsp_write_probfile_sav\n");
    }

CLEANUP:
    if (pool) { CCtsp_free_cutpool (&pool); }
    if (dominopool) { CCtsp_free_cutpool (&dominopool); }
    if (epool) { CCtsp_free_cutpool (&pool); }
    if (edominopool) { CCtsp_free_cutpool (&dominopool); }
    CC_IFFREE (ptour, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (flist, int);
    CC_IFFREE (flen, int);
    CC_IFFREE (invperm, int);
    CCutil_freedatagroup (&dat);
    if (rootlp) CCtsp_free_tsp_lp_struct (&rootlp);
    if (worklp) CCtsp_free_tsp_lp_struct (&worklp);
    if (extralp) CCtsp_free_tsp_lp_struct (&extralp);

    return rval;
}

static int replace_all_edges (CCtsp_lp *oldlp, CCtsp_lp **newlp,
        int ecount, int *elist, int *elen, CCrandstate *rstate, int silent,
        CCtsp_lpcuts *pool, CCtsp_lpcuts *dominopool)
{
    int rval = 0, i, j, ncount = oldlp->graph.ncount, cutadded = 0;
    int tighten = 0, *perm = oldlp->perm, *permuted_elist = (int *) NULL;
    int *invperm = (int *) NULL, infeasible = 0;
    double val;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    char *probname = (char *) NULL;
    CCtsp_lpcuts *cuts = &oldlp->cuts;
    CCtsp_lpcut_in *c, *newcuts = (CCtsp_lpcut_in *) NULL;

    *newlp = (CCtsp_lp *) NULL;

    probname = oldlp->problabel;

    CC_MALLOC (invperm, ncount, int);
    for (i = 0; i < ncount; i++) invperm[perm[i]] = i;

    permuted_elist = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (permuted_elist, "out of memory for permuted_elist");
    for (i = 0; i < ecount; i++) {
        permuted_elist[2*i]   = invperm[elist[2*i]];
        permuted_elist[2*i+1] = invperm[elist[2*i+1]];
    }

    rval = CCtsp_init_lp (&lp, probname, -1, (char *) NULL, ncount,
              oldlp->dat, ecount, permuted_elist, elen, ecount,
              permuted_elist, elen, 1, perm, oldlp->upperbound,
              pool, dominopool, silent, rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "Initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    printf ("New LP is initialized\n"); fflush (stdout);

    for (j = 0; j < 10; j++) {
        printf ("Starting loop %d\n", j); fflush (stdout);
        for (i = 0; i < cuts->cutcount; i++) {
            CC_MALLOC (c, 1, CCtsp_lpcut_in);
            CCtsp_init_lpcut_in (c);
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            c->next = newcuts;
            newcuts = c;
        }
        CCtsp_add_cuts_to_queue (lp, &newcuts);
        rval = CCtsp_process_cuts (lp, &cutadded, tighten, silent, rstate,
                                   (double *) NULL, &infeasible);
        CCcheck_rval (rval, "CCtsp_process_cuts failed");
        if (infeasible) {
            fprintf (stderr, "Infeasible LP after cuts added\n");
            rval = 1; goto CLEANUP;
        }
        printf ("Added %d of %d cuts\n", cutadded, cuts->cutcount);
        fflush (stdout);
        rval = CCtsp_get_lp_result (lp, &val, (double *) NULL, (int *) NULL,
                     (int **) NULL, (double **) NULL, (double **) NULL,
                     (double **) NULL, (double **) NULL);
        if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
        printf ("Round %d LP Value: %f\n", j, val); fflush (stdout);
    }

    *newlp = lp;

CLEANUP:
    CC_IFFREE (invperm, int);
    CC_IFFREE (permuted_elist, int);
    return rval;
}

static int add_extra_cuts (CCtsp_lp *lp, CCtsp_lp *extralp,
        CCrandstate *rstate, int silent)
{
    int rval = 0, i, j, tighten = 0, cutadded = 0, infeasible = 0;
    CCtsp_lpcuts *cuts = &extralp->cuts;
    CCtsp_lpcut_in *c, *newcuts = (CCtsp_lpcut_in *) NULL;
    double val;

    for (j = 0; j < 10; j++) {
        printf ("Starting extra, loop %d\n", j); fflush (stdout);
        for (i = 0; i < cuts->cutcount; i++) {
            CC_MALLOC (c, 1, CCtsp_lpcut_in);
            CCtsp_init_lpcut_in (c);
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            c->next = newcuts;
            newcuts = c;
        }
        CCtsp_add_cuts_to_queue (lp, &newcuts);
        rval = CCtsp_process_cuts (lp, &cutadded, tighten, silent, rstate,
                                   (double *) NULL, &infeasible);
        CCcheck_rval (rval, "CCtsp_process_cuts failed");
        if (infeasible) {
            fprintf (stderr, "Infeasible LP after cuts added\n");
            rval = 1; goto CLEANUP;
        }
        printf ("Added %d of %d cuts\n", cutadded, cuts->cutcount);
        fflush (stdout);
        rval = CCtsp_get_lp_result (lp, &val, (double *) NULL, (int *) NULL,
                     (int **) NULL, (double **) NULL, (double **) NULL,
                     (double **) NULL, (double **) NULL);
        if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
        printf ("Round %d LP Value: %f\n", j, val); fflush (stdout);
    }

CLEANUP:
    return rval;
}

static int lpedges_sanity_check (CCtsp_lp *lp, int fullcount, int *fullelist)
{
    CCtsp_lpgraph *g = &lp->graph;
    CCtsp_lpedge *e;
    int i, k, n, n0, n1, tmp, got = 0, rval = 0;
    int *invperm = (int *) NULL;

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

    CC_MALLOC (invperm, lp->graph.ncount, int);
    for (i = 0; i < lp->graph.ncount; i++) invperm[lp->perm[i]] = i;

    for (i = 0; i < fullcount; i++) {
        n0 = invperm[fullelist[2*i]];
        n1 = invperm[fullelist[2*i+1]];
        if (n0 >= n1) {
            CC_SWAP (n0, n1, tmp);
        }
        got += edge_in_lpgraph (g, n0, n1, (int *) NULL);
    }

    printf ("Found %d of the LP edges Sir.\n", got);
    fflush (stdout);

CLEANUP:
    CC_IFFREE (invperm, int);
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

    while ((c = CCutil_bix_getopt (ac, av, "e:f:M:s:x:", &boptind, &boptarg)) != EOF) { 
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
        case 'x':
            extrafname = boptarg;
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
    fprintf (stderr, "   -x f  second savfile for extra cuts\n");
}


