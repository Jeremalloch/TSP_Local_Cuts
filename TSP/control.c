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
/*                  THE CONTROLLER FOR CUTTING PLANES                       */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: June 27, 1997                                                     */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCtsp_init_cutselect (CCtsp_cutselect *s)                          */
/*    INITIALIZES the cut selections                                        */
/*                                                                          */
/*  void CCtsp_cutselect_dominos (CCtsp_cutselect *s, int domsel)           */
/*    SETS the domino field according to value of domsel                    */
/*                                                                          */
/*  void CCtsp_cutselect_tighten (CCtsp_cutselect *s, int tighten)          */
/*    SETS the usetighten field according to value of tighten               */
/*                                                                          */
/*  void CCtsp_cutselect_chunksize (CCtsp_cutselect *s, int chunksize)      */
/*    SETS the maxchunksize filed according to value of chunksize           */
/*                                                                          */
/*  void CCtsp_cutselect_filecuts (CCtsp_cutselect *s, char *fname)         */
/*    SETS the cutselector to read cuts from file fname                     */
/*                                                                          */
/*  void CCtsp_cutselect_remotepool (CCtsp_cutselect *s, char *cutbossname) */
/*    SETS the cutselector to use the specified cutboss                     */
/*                                                                          */
/*  void CCtsp_cutselect_remotedompool (CCtsp_cutselect *s,                 */
/*       char *domcutbossname)                                              */
/*    SETS the cutselector to use the specified cutboss                     */
/*                                                                          */
/*  void CCtsp_cutselect_domboss (CCtsp_cutselect *s, char *dombossname)    */
/*    SETS the cutselector to use the specified domino boss                 */
/*                                                                          */
/*  void CCtsp_cutselect_keep_cutting (CCtsp_cutselect *s)                  */
/*    SETS the cutselector to use keep cutting the root LP                  */
/*                                                                          */
/*  void CCtsp_init_tentative_cutselect (CCtsp_cutselect *s)                */
/*    INITIALIZES the cut selections for tenative branching                 */
/*                                                                          */
/*  int CCtsp_cutselect_set_tols (CCtsp_cutselect *s, CCtsp_lp *lp,         */
/*      int level, int silent)                                              */
/*    SETS the tolerances for the cut selections                            */
/*     -level should be set to 0 for tentative cutting                      */
/*    NOTES: The lp should be solved before this call.                      */
/*                                                                          */
/*  int CCtsp_cutting_multiple_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,    */
/*      int *snowtour, int savelp, int maxlocal, int update_tol,            */
/*      int silent, CCrandstate *rstate, int *infeasible, int *snowtour)    */
/*    CALLS CCtsp_cutting_loop multiple times, incrementing maxchunksize    */
/*      from 16 up to maxlocal, going up by 4 each time.                    */
/*     -infeasible set to 1 if LP becomes infeasible                        */
/*                                                                          */
/*  int CCtsp_cutting_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,             */
/*      int *snowtour, int savelp, int silent, CCrandstate *rstate,         */
/*      int *infeasible)                                                    */
/*    CALLS the cutting plane and pricing routines.                         */
/*     -sel should be set with the desired cut selection.                   */
/*     -snowtour (if not NULL) should input a tour in permuation format to  */
/*      be used to generate closex vectors for snowcuts                     */
/*     -savelp should be set to a nonzero value to write the lps to after   */
/*      rounds of cuts                                                      */
/*     -silent turns off most output if set to a nonzero value              */
/*     -infeasible set to 1 if LP becomes infeasible                        */
/*                                                                          */
/*  int CCtsp_subtour_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate,  */
/*      int *infeasible)                                                    */
/*    CALLS the cutting and pricing to optimize over the subtour polytope.  */
/*     -infeasible set to 1 if LP becomes infeasible                        */
/*                                                                          */
/*  int CCtsp_blossom_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate   */
/*      int *infeasible)                                                    */
/*    CALLS the cutting and pricing to optimize over the blossom polytope.  */
/*     -infeasible set to 1 if LP becomes infeasible                        */
/*                                                                          */
/*  int CCtsp_subtour_and_blossom_loop (CCtsp_lp *lp, int silent,           */
/*      CCrandstate *rstate)                                                */
/*    CALLS the cutting and princing to optimize over subtours and          */
/*     trivial blossoms.                                                    */
/*     -infeasible set to 1 if LP becomes infeasible                        */
/*                                                                          */
/*  int CCtsp_pricing_loop (CCtsp_lp *lp, double *bnd, int silent,          */
/*      CCrandstate *rstate)                                                */
/*    ADDS negative reduced costs edges to lp and returns the current       */
/*     lowerbound.                                                          */
/*     -bnd can be NULL                                                     */
/*    NOTES: The LP must have full_edges_valid.                             */
/*                                                                          */
/*  int CCtsp_call_x_heuristic (CCtsp_lp *lp, double *val, int *outcyc,     */
/*      int silent, CCrandstate *rstate)                                    */
/*    CALLS the x-greedy LK heuristic with the current LP solution.         */
/*     -val returns the length of the tour.                                 */
/*     -outcyc will return the tour in node-node-node format if the         */
/*      length of the tour is less than lp->upperbound; the array should    */
/*      at least of length ncount (it can be NULL)                          */
/*                                                                          */
/*  int CCtsp_bb_cutting (char *probname, int probnum, int prob_newnum,     */
/*      int ncount, CCdatagroup *dat, int *ptour, double *upbound,          */
/*      CCtsp_lpcuts *pool, CCtsp_lpcuts *dominopool, CCtsp_cutselect *sel, */
/*      double *val, int *prune, int *foundtour, int *besttour, int level,  */
/*      int silent, CCrandstate *rstate)                                    */
/*    CALLS the cutting loop after reading the lp; writes the result to     */
/*     prob file prob_newnum; using exact price to verify pruned runs       */
/*     -upbound should be passed in as the current bound; if a better       */
/*      tour is found then upbound will be updated                          */
/*     -val returns the lp bound; it is CCtsp_LP_MAXDOUBLE if infeasible    */
/*     -prune is set to 1 if bbnode can be pruned                           */
/*     -foundtour is set to 1 if a better tour is found.                    */
/*     -besttour (if not NULL) will return a better tour if one is found.   */
/*     -level should be set to 0 for tentative cutting                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "bigguy.h"
#include "cut.h"
#include "pq.h"
#include "cuttree.h"
#include "consec1.h"
#include "necklace.h"
#include "localcut.h"
#if 0
#include "verify.h"
#endif

static int
    call_add_cuts (CCtsp_lp *lp, CCtsp_lpcut_in **cuts, int *cut_added,
        int *xcount, int **xlist, double **x, double *val, int tighten,
        int *istour, int silent, CCrandstate *rstate, int *infeasible),
    lp_value (CCtsp_lp *lp, double *val),
    lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x),
    lp_upperbound (CCtsp_lp *lp, double *ub),
    full_edge_check (CCtsp_lp *lp, int *nadded, double *bnd, int silent,
        CCrandstate *rstate),
    sparse_edge_check (CCtsp_lp *lp, CCtsp_edgegenerator *eg, int *nadded,
        double *bnd, int silent, CCrandstate *rstate),
    bb_cutting_work (CCtsp_lp **lp, char *probname, int probnum,
        int ncount, CCdatagroup *dat, int *ptour, double upbound,
        CCtsp_lpcuts *pool, CCtsp_lpcuts *dominopool,
        CCtsp_cutselect *sel, double *val, int level,
        int silent, CCrandstate *rstate),
    grab_local_x (int ncount, int *tour, int xcount, int *xlist,
        double *x, int *newcount, int **newlist, double **newx, double mult),
    no_tighten (int ncount, int xcount, int *xlist, double *x, int *test,
        double tol);

static void
    zero_cutselect (CCtsp_cutselect *s),
    cut_msg (int silent, const char *msg, int count, double sec,
        double maxviol);

static int subtours_to_general (CCtsp_lpcut_in *cuts, int ncount);

static void zero_cutselect (CCtsp_cutselect *s)
{
    s->cutpool             = 0;
    s->remotepool          = 0;
    s->remotehost          = (char *) NULL;
    s->remoteport          = 0;
    s->remotedompool       = 0;
    s->remotedomhost       = (char *) NULL;
    s->remotedomport       = 0;
    s->dominopool          = 0;
    s->domboss             = 0;
    s->dombosshost         = (char *) NULL;
    s->connect             = 0;
    s->segments            = 0;
    s->blockcombs          = 0;
    s->growcombs           = 0;
    s->prclique            = 0;
    s->exactsubtour        = 0;
    s->exactblossom        = 0;
    s->tighten_lp          = 0;
    s->xtighten_lp         = 0;
    s->teething_lp         = 0;
    s->tighten_pool        = 0;
    s->teething_pool       = 0;
    s->cliquetree_lp       = 0;
    s->fastblossom         = 0;
    s->ghfastblossom       = 0;
    s->consecutiveones     = 0;
    s->necklace            = 0;
    s->usetighten          = 0;
    s->extra_connect       = 0;
    s->decker_lp           = 0;
    s->decker_pool         = 0;
    s->star_lp             = 0;
    s->star_pool           = 0;
    s->handling_lp         = 0;
    s->handling_pool       = 0;
    s->maxchunksize        = 0; 
    s->filecuts            = 0;
    s->filecutname         = (char *) NULL;
    s->nexttol             = 0.0;
    s->roundtol            = 0.0;
    s->domino_tighten_lp   = 0;
    s->domino_tighten_pool = 0;
    s->DP_cuts             = 0;
    s->TP_cuts             = 0;
    s->domino_safe_shrink  = 0;
    s->fastcuts            = 0;
    s->keep_cutting        = 0;
}

void CCtsp_init_cutselect (CCtsp_cutselect *s)
{
    s->cutpool             = 1;  /* 1 */
    s->remotepool          = 0;
    s->remotehost          = (char *) NULL;
    s->remoteport          = 0;
    s->remotedompool       = 0;
    s->remotedomhost       = (char *) NULL;
    s->remotedomport       = 0;
    s->dominopool          = 1;  /* 1 */
    s->domboss             = 0;
    s->dombosshost         = (char *) NULL;
    s->connect             = 1;  /* 1 */
    s->segments            = 1;  /* 1 */
    s->blockcombs          = 1;  /* 1 */
    s->growcombs           = 0;  /* 0 */
    s->prclique            = 0;  /* 0 */
    s->exactsubtour        = 1;  /* 1 */
    s->exactblossom        = 0;  /* 0 */
    s->tighten_lp          = 1;  /* 1 */
    s->xtighten_lp         = 0;  /* 0 */
    s->teething_lp         = 1;  /* 1 */
    s->tighten_pool        = 0;  /* 1 */   /* Slow when pool is large */
    s->teething_pool       = 0;  /* 0 */
    s->cliquetree_lp       = 0;  /* 0 */
    s->fastblossom         = 1;  /* 1 */
    s->ghfastblossom       = 1;  /* 1 */
    s->consecutiveones     = 0;  /* 1 */   /* Slow on large instances */
    s->necklace            = 0;  /* 0 */   /* Uses a big chunk of memory */
    s->usetighten          = 0;  /* 0 */
    s->extra_connect       = 0;  /* 0 */
    s->decker_lp           = 1;  /* 1 */
    s->decker_pool         = 0;  /* 0 */
    s->star_lp             = 0;  /* 0 */   /* Slow on very large instances */
    s->star_pool           = 0;  /* 0 */   /* Slow on most instances */
    s->handling_lp         = 1;  /* 1 */
    s->handling_pool       = 0;  /* 0 */
    s->maxchunksize        = 16; /* 16 */
    s->filecuts            = 0;  /* 0 */
    s->filecutname         = (char *) NULL;
    s->nexttol             = 0.0;
    s->roundtol            = 0.0;
    s->domino_tighten_lp   = 1;
    s->domino_tighten_pool = 1;   /* 0 */
    s->DP_cuts             = 0;
    s->TP_cuts              = 0;
    s->domino_safe_shrink  = 1;
    s->fastcuts            = 0;           /* Keep this at 0 (changes tols) */
    s->keep_cutting        = 0;
}

void CCtsp_cutselect_dominos (CCtsp_cutselect *s, int domsel, int safeshrink)
{
    if (domsel == 1 || domsel == 3) s->DP_cuts = 1;
    if (domsel == 2 || domsel == 3) s->TP_cuts = 1;
    s->domino_safe_shrink = safeshrink;
}

void CCtsp_cutselect_tighten (CCtsp_cutselect *s, int tighten)
{
    s->usetighten = tighten;
}

void CCtsp_cutselect_xtighten (CCtsp_cutselect *s, int tighten)
{
    s->xtighten_lp = tighten;
}

void CCtsp_cutselect_chunksize (CCtsp_cutselect *s, int chunksize)
{
    s->maxchunksize = chunksize;
}

void CCtsp_cutselect_filecuts (CCtsp_cutselect *s, char *fname)
{
    s->filecuts = 1;
    s->filecutname = fname;
}

void CCtsp_cutselect_remotepool (CCtsp_cutselect *s, char *cutbossname)
{
    s->remotepool = 1;
    s->remotehost = cutbossname;
    s->remoteport = CCtsp_CUT_PORT;
}

void CCtsp_cutselect_remotedompool (CCtsp_cutselect *s, char *domcutbossname)
{
    s->remotedompool = 1;
    s->remotedomhost = domcutbossname;
    s->remotedomport = CCtsp_DOMCUT_PORT;
}

void CCtsp_cutselect_domboss (CCtsp_cutselect *s, char *dombossname)
{
    s->domboss = 1;
    s->dombosshost = dombossname;
}

void CCtsp_cutselect_keep_cutting (CCtsp_cutselect *s)
{
    s->keep_cutting = 1;
}

void CCtsp_init_tentative_cutselect (CCtsp_cutselect *s)
{
    zero_cutselect (s);
    s->cutpool             = 1;
    s->connect             = 1;
    s->segments            = 1;
    s->blockcombs          = 1;
    s->exactsubtour        = 1;
    s->tighten_lp          = 1;
    s->teething_lp         = 1;
    s->fastblossom         = 1;
    s->ghfastblossom       = 1;
    s->decker_lp           = 1;
    s->handling_lp         = 1;
    s->domino_safe_shrink  = 1;
}

void CCtsp_init_simple_cutselect (CCtsp_cutselect *s)
{
    zero_cutselect (s);
    s->connect             = 1;
    s->segments            = 1;
    s->tighten_lp          = 1;
    s->decker_lp           = 1;
    s->domino_safe_shrink  = 1;
}

void CCtsp_init_fast_cutselect (CCtsp_cutselect *s)
{
    zero_cutselect (s);
    s->connect             = 1;
    s->segments            = 1;
    s->blockcombs          = 1;
    s->exactsubtour        = 1;
    s->tighten_lp          = 1;
    s->fastblossom         = 1;
    s->ghfastblossom       = 1;
    s->decker_lp           = 1;
    s->domino_safe_shrink  = 1;
    s->fastcuts            = 1;
}

int CCtsp_cutselect_set_tols (CCtsp_cutselect *s, CCtsp_lp *lp, int level,
        int silent)
{
    int rval = 0;
    double ub, lb, beta;

    rval = CCtsp_get_lp_result (lp, &lb, &ub, (int *) NULL,
              (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

#ifdef CCtsp_CUTS_DELTA
    if (ub - lb < 1.0) {
        beta = 1.0;
    } else {
        if (ub > 2*lb) beta = lb;
        else           beta = ub - lb;
    }
#else
    if (ub > 2*lb) beta = lb;
    else           beta = ub;
#endif

    if (level == -1) {
        s->nexttol  = 10.0 * CCtsp_TENTATIVE_CUTS_NEXT_TOL * beta;
        s->roundtol = 10.0 * CCtsp_TENTATIVE_CUTS_NEXT_ROUND * beta;
    } else if (level == 0) {
        s->nexttol  = /* 0.001 * */ CCtsp_TENTATIVE_CUTS_NEXT_TOL * beta;
        s->roundtol = /* 0.001 * */ CCtsp_TENTATIVE_CUTS_NEXT_ROUND * beta;
    } else if (level == 1){
        s->nexttol  = /* 0.001 * */ CCtsp_CUTS_NEXT_TOL * beta; 
        s->roundtol = /* 0.001 * */ CCtsp_CUTS_NEXT_ROUND * beta;
    } else {
        s->nexttol  = 0.0001  * lb; 
        s->roundtol = 0.00001 * lb;
    }

    /* Simple Branching 
    s->nexttol =  2 * CCtsp_CUTS_NEXT_TOL * beta;
    s->roundtol = 2 * CCtsp_CUTS_NEXT_ROUND * beta;
    */

    /* For fast runs in US50K
    s->nexttol =  100.0;
    s->roundtol = 25.0;
    */

    /* For fast runs in UK50
    s->nexttol =  100.0;
    s->roundtol = 10.0;
    */

    if (!silent) {
        printf ("Setting tolerances: next cuts %.4f next round %.4f\n",
                s->nexttol, s->roundtol);
        fflush (stdout);
    }

CLEANUP:
    return rval;
}

int CCtsp_cutting_multiple_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,
        int *snowtour, int savelp, int maxlocal, int update_tol, int silent,
        CCrandstate *rstate, int *infeasible)
{
    int rval = 0, k;

    *infeasible = 0;

    if (maxlocal < 16) {
        rval = CCtsp_cutting_loop (lp, sel, snowtour, savelp, silent, rstate,
                                   infeasible);
        CCcheck_rval (rval, "CCtsp_cutting_loop failed");
    } else {
        for (k = 16; k <= maxlocal; k += 4) {
            sel->maxchunksize = k;
            printf ("SETTING MAXCHUNKSIZE = %d\n", k); fflush (stdout);
            if (update_tol) {
                rval = CCtsp_cutselect_set_tols (sel, lp, 1, 0);
                CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
            }
            rval = CCtsp_cutting_loop (lp, sel, snowtour, savelp, silent,
                                       rstate, infeasible);
            CCcheck_rval (rval, "CCtsp_cutting_loop failed");
        }
        if (maxlocal % 4 != 0) {
            sel->maxchunksize = maxlocal;
            printf ("SETTING MAXCHUNKSIZE = %d\n", maxlocal); fflush (stdout);
            if (update_tol) {
                rval = CCtsp_cutselect_set_tols (sel, lp, 1, 0);
                CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
            }
            rval = CCtsp_cutting_loop (lp, sel, snowtour, savelp, silent,
                                       rstate, infeasible);
            CCcheck_rval (rval, "CCtsp_cutting_loop failed");
        }
    }

CLEANUP:
    return rval;
}


#define LOOP_FULL (25)      /* to force a full price after 25 inner loops */
#define CC_NO_NEAREST (50)  /* the initial sparse graph for pricing       */
#define KEEP_MAX (100)      /* keep_cutting loops without improvement     */

#define PROCESS_CUTS(xmsg) {                                                \
    cut_msg (silent, xmsg, cutcount, z, maxviol);                           \
    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount, &xlist,           \
            &x, &newval, sel->usetighten, &istour, silent, rstate,          \
            infeasible);                                                    \
    CCcheck_rval (rval, "call_add_cuts failed");                            \
    totalcuts += cut_added;                                                 \
    if (*infeasible) goto CLEANUP;                                          \
    if (istour) goto OUT_LOOP;                                              \
}

#define PROCESS_CUTS_0(xmsg) {                                              \
    cut_msg (silent, xmsg, cutcount, z, 0.0);                               \
    rval = call_add_cuts (lp, &cuts, &cut_added, &xcount, &xlist,           \
            &x, &newval, sel->usetighten, &istour, silent, rstate,          \
            infeasible);                                                    \
    CCcheck_rval (rval, "call_add_cuts failed");                            \
    totalcuts += cut_added;                                                 \
    if (*infeasible) goto CLEANUP;                                          \
    if (istour) goto OUT_LOOP;                                              \
}

int CCtsp_cutting_loop (CCtsp_lp *lp, CCtsp_cutselect *sel, int *snowtour, 
        int savelp, int silent, CCrandstate *rstate, int *infeasible)
{
    int rval = 0, tval, ttest = 0, outside = 0, qqq, kloops = 0;
    int xcount, cutcount, cutcount_connect, cut_added, edge_added, closecount;
    int loopcount = 0, istour = 0, ncount = lp->graph.ncount, totalcuts = 0;
    int *xlist = (int *) NULL, *closelist = (int *) NULL;
#ifdef CCtsp_USE_DOMINO_CUTS
    int domstarter = 0;
#endif
    double newval, oldval, ub, priceval, maxviol, z, dval, qval;
    double *x = (double *) NULL, *closex = (double *) NULL;
    double save_nexttol  = sel->nexttol, save_roundtol = sel->roundtol;
    double lcimprove = CCtsp_LP_MAXDOUBLE, otherimprove = 0.0;
    double szeit = CCutil_zeit ();
    char buf[1024];
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCchunk_localcut_timer lc_timer;
    CCtsp_edgegenerator eginside;

    CCutil_start_timer (&lp->stats.cutting_loop);
    CCchunk_init_localcut_timer (&lc_timer);
    
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed\n");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, ncount, lp->dat,
                 (CCtsp_genadj *) NULL, CC_NO_NEAREST, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator (sparse) failed\n");
    }

    rval = lp_upperbound (lp, &ub);
    CCcheck_rval (rval, "lp_upperbound failed");

    do {
        loopcount = 0;
        cutcount_connect = 0;
        totalcuts = 0;
        do {
            CCutil_start_timer (&lp->stats.cutting_inner_loop);

            /* STAR ROTATE: Turn off to get same LPs as non-rotate */
            rval = CCtsp_resparsify_lp (lp, rstate, silent);
            CCcheck_rval (rval, "CCtsp_resparsify_lp failed");

            cut_added = 0;
            rval = lp_value (lp, &oldval);
            CCcheck_rval (rval, "lp_value failed");
            newval = oldval;

            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

            rval = CCtsp_check_integral (lp, &dval, (int **) NULL, &istour,
                                         silent);
            CCcheck_rval (rval, "CCtsp_check_integral failed");
            if (istour) goto OUT_LOOP;

            if (sel->filecuts) {
                rval = CCtsp_file_cuts (sel->filecutname, &cuts, &cutcount,
                                        ncount, lp->perm, &z);
                CCcheck_rval (rval, "CCtsp_file_cuts failed");
                PROCESS_CUTS_0("file cuts")
            }

            if (sel->cutpool) {
                qqq = 10;  /* 10 */
                do {
                    rval = CCtsp_search_cutpool (lp->pool, &cuts, &cutcount,
                       &maxviol, ncount, xcount, xlist, x, 0, rstate, &z);
                    CCcheck_rval (rval, "CCtsp_search_cutpool failed");
                    PROCESS_CUTS("pool cuts")
                } while (cutcount && --qqq);
            }

            if (sel->connect) {
                rval = CCtsp_connect_cuts (&cuts, &cutcount, ncount, xcount,
                                           xlist, x, &z);
                CCcheck_rval (rval, "CCtsp_connect_cuts failed");
                PROCESS_CUTS_0("connect cuts")
            }

            if (sel->segments) {
                rval = CCtsp_segment_cuts (&cuts, &cutcount, ncount, xcount,
                                           xlist, x, &z);
                CCcheck_rval (rval, "CCtsp_segment_cuts failed");
                PROCESS_CUTS_0("segment cuts")
            }

            if (sel->fastblossom) {
                rval = CCtsp_fastblossom (&cuts, &cutcount, ncount, xcount,
                                          xlist, x, &z);
                CCcheck_rval (rval, "CCtsp_fastblossom failed");
                PROCESS_CUTS_0("fast blossoms cuts")
            }

            if (sel->ghfastblossom) {
                rval = CCtsp_ghfastblossom (&cuts, &cutcount, ncount, xcount,
                                            xlist, x, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_ghfastblossom failed");
                PROCESS_CUTS("Groetschel-Holland blossoms")
            }

            if (sel->remotepool) {
                qqq = 1;  /* 25 */
                do {
                    rval = CCtsp_search_remotepool (sel->remotehost,
                       sel->remoteport, &cuts, &cutcount, &maxviol, ncount,
                       xcount, xlist, x, &z);
                    CCcheck_rval (rval, "CCtsp_search_remotepool failed");
                    PROCESS_CUTS("remote pool cuts")
                } while (cutcount && --qqq);
            }

            if (0 && sel->remotedompool) {  /* US50K: Turn off - too slow */
                printf ("DOM CUTPOOL: %s\n", sel->remotedomhost);
                fflush (stdout);
                qqq = 1;  /* 1 */
                do {
                    rval = CCtsp_search_remotepool (sel->remotedomhost,
                       sel->remotedomport, &cuts, &cutcount, &maxviol,
                       ncount, xcount, xlist, x, &z);
                    CCcheck_rval (rval, "CCtsp_search_remotedompool failed");
                    PROCESS_CUTS("remote dompool cuts")
                } while (cutcount && --qqq);
            }

            if (sel->blockcombs) {
                rval = CCtsp_block_combs (&cuts, &cutcount, ncount, xcount,
                   xlist, x, silent, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_block_combs failed");
                PROCESS_CUTS("block combs")
            }

            if (sel->growcombs) {
                rval = CCtsp_edge_comb_grower(&cuts, &cutcount, ncount, xcount,
                   xlist, x, &lp->stats.extra_tighten_stats, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_block_combs failed");
                PROCESS_CUTS("grown combs")
            }

            if (sel->prclique) {
                rval = CCtsp_pr_cliquetree (&cuts, &cutcount, ncount, xcount,
                   xlist, x, &lp->stats.extra_tighten_stats, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_pr_cliquetree failed");
                PROCESS_CUTS("PR cliquetrees")
            }

            if (sel->exactsubtour) {
                rval = CCtsp_exact_subtours (&cuts, &cutcount,
                   ncount, xcount, xlist, x, (double *) NULL, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_exact_subtours failed");
                PROCESS_CUTS("exact subtours")
            }

#ifdef CCtsp_USE_DOMINO_CUTS
#if 0
            {
                int save_silent = silent;
                silent = 0;
                rval = CCtsp_comb2dp_lp(lp->pool, &cuts, &cutcount,
                   ncount, xcount, xlist, x, 2.0, &maxviol, &z);
                CCcheck_rval (rval, "CCtsp_comb2dp_lp failed");
                PROCESS_CUTS("comb-to-dp")
                silent = save_silent;
            }

#endif

            if (sel->dominopool && lp->dominopool) {
                qqq = 10;  /* 5 */
                do {
                    rval = CCtsp_search_dominopool (lp->dominopool, &cuts,
                       &cutcount, &maxviol, ncount, xcount, xlist, x,
                       &domstarter, 25000 /* 10000 */, 500, rstate, &z);
                    CCcheck_rval (rval, "CCtsp_search_dominopool failed");
                    PROCESS_CUTS("dompinopool cuts")
                } while (cutcount && --qqq);
            }

#if 1  /* STAR TSP turn off tighten domino */
            if (sel->domino_tighten_lp) {
                rval = CCtsp_domino_tighten_lp (&lp->cuts, &cuts, &cutcount,
                   ncount, xcount, xlist, x, 0.01, 500, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_domino_tighten_lp failed");
                PROCESS_CUTS("domino tighten cuts")
            }

            if (sel->domino_tighten_pool && lp->dominopool) {
                rval = CCtsp_domino_tighten_pool (lp->dominopool, &cuts,
                   &cutcount, ncount, xcount, xlist, x, 0.01, 1000,
                   &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_domino_tighten_lp failedd");
                PROCESS_CUTS("dominopool tighten cuts")
            }
#endif

            if (sel->DP_cuts) {
                rval = lp_value (lp, &qval);
                CCcheck_rval (rval, "lp_value failed");

                /* DP and TP separators need to be nearly in subtour polytope */
                qqq = 25;  /* 25 */
                do {
                    cut_added = 0;
                    rval = CCtsp_exact_subtours (&cuts, &cutcount, ncount,
                       xcount, xlist, x, (double *) NULL, &z, &maxviol);
                    CCcheck_rval (rval, "CCtsp_exact_subtours failed");
                    PROCESS_CUTS("exact subtour")

                    if (newval > qval + (sel->roundtol / 10.0)) {
                        qval = newval;  /* If small improvement, keep going */
                        qqq = 26;
                    }
                } while (cut_added > 0 && --qqq);

                /* and definitely connected */
                qqq = 100;
                while (cut_added && --qqq) {
                    cut_added = 0;
                    rval = CCtsp_connect_cuts (&cuts, &cutcount, ncount,
                                               xcount, xlist, x, &z);
                    CCcheck_rval (rval, "CCtsp_connect_cuts failed");
                    PROCESS_CUTS_0("connect subtours")
                }
   

                rval = CCtsp_DP_cuts (&cuts, &cutcount, ncount, xcount,
                         xlist, x, sel->domino_safe_shrink, sel->dombosshost,
                         &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_DP_cuts failed");
                PROCESS_CUTS("DP cuts")
            }

            if (sel->TP_cuts) {
                int isilent = silent;
                rval = lp_value (lp, &qval);
                CCcheck_rval (rval, "lp_value failed");

                qqq = 25;  /* 25 */
                do {
                    cut_added = 0;
                    rval = CCtsp_exact_subtours (&cuts, &cutcount, ncount,
                       xcount, xlist, x, (double *) NULL, &z, &maxviol);
                    CCcheck_rval (rval, "CCtsp_exact_subtours failed");
                    PROCESS_CUTS("exact subtours")
                    if (newval > qval + (sel->roundtol / 10.0)) {
                        qval = newval;  /* If small improvement, keep going */
                        qqq = 26;
                    }
                } while (cut_added > 0 && --qqq);

                qqq = 100;
                while (cut_added && --qqq) {
                    cut_added = 0;
                    rval = CCtsp_connect_cuts (&cuts, &cutcount, ncount,
                                               xcount, xlist, x, &z);
                    CCcheck_rval (rval, "CCtsp_connect_cuts failed");
                    PROCESS_CUTS_0("connect subtours")
                }
   
                if (cut_added == 0) {
                    rval = CCtsp_TP_cuts (&cuts, &cutcount, ncount,
                                    xcount, xlist, x, &z, &maxviol);
                    CCcheck_rval (rval, "CCtsp_TP_cuts failed");
                    silent = 0;
                    PROCESS_CUTS("TP cuts")
                    silent = isilent;
                }
            }
#endif  /* CCtsp_USE_DOMINO_CUTS */

            if (sel->exactblossom) {
                rval = CCtsp_exactblossom (&cuts, &cutcount, ncount, xcount,
                                           xlist, x, rstate, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_exactblossom failed");
                PROCESS_CUTS("exact blossoms")
            }

            if (sel->tighten_lp) {
                rval = CCtsp_tighten_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   xcount, xlist, x, 0.5, 500, &maxviol, 0, rstate, &z);
                CCcheck_rval (rval, "CCtsp_tighten_lp failed");
                PROCESS_CUTS("tighten_lp cuts")

                rval = grab_local_x (ncount, snowtour, xcount, xlist, x,
                                     &closecount, &closelist, &closex, 0.5);
                CCcheck_rval (rval, "grab_local_x failed");
                rval = CCtsp_tighten_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   closecount, closelist, closex, 0.5, 500, &maxviol, 0,
                   rstate, &z);
                CCcheck_rval (rval, "CCtsp_tighten_lp failed");
                PROCESS_CUTS("CLOSE tighten_lp cuts")
            }

            if (sel->xtighten_lp) {
                rval = CCtsp_tighten_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   xcount, xlist, x, 0.00005, 2000, &maxviol, 1, rstate, &z);
                CCcheck_rval (rval, "CCtsp_xtighten_lp failed");
                PROCESS_CUTS("xtighten_lp cuts")
            }

            if (sel->decker_lp) {
                rval = CCtsp_double_decker_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   xcount, xlist, x, 10.0, 500, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_double_decker_lp failed");
                PROCESS_CUTS("double deckers")

                rval = grab_local_x (ncount, snowtour, xcount, xlist, x,
                                     &closecount, &closelist, &closex, 0.5);
                CCcheck_rval (rval, "grab_local_x failed");
                rval = CCtsp_double_decker_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   closecount, closelist, closex, 10.0, 500, &maxviol, rstate,
                   &z);
                CCcheck_rval (rval, "CCtsp_double_decker_lp failed");
                PROCESS_CUTS("CLOSE double deckers")
            } 

            if (sel->star_lp) {
                rval = CCtsp_star_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   xcount, xlist, x, 10.0, 500, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_star_lp failed");
                PROCESS_CUTS("star inequalities")
            }

            if (sel->handling_lp) {
                rval = no_tighten (ncount, xcount, xlist, x, &ttest, 0.20);
                CCcheck_rval (rval, "no_tighten failed");
                if (!ttest) {
                    rval = CCtsp_handling_lp (&lp->cuts,
                       &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                       ncount, xcount, xlist, x, 10.0, 500, &maxviol, rstate,
                       &z);
                    CCcheck_rval (rval, "CCtsp_handling_lp failed");
                    PROCESS_CUTS("handling inequalities")
                }
            }

            if (sel->cliquetree_lp) {
                rval = CCtsp_cliquetree_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                   ncount, xcount, xlist, x, 10.0, 500, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_cliquetree_lp failed");
                PROCESS_CUTS("comb cliquetrees")
            }

            if (sel->teething_lp) {
                rval = no_tighten (ncount, xcount, xlist, x, &ttest, 0.20);
                CCcheck_rval (rval, "no_tighten failed");
                if (!ttest) {
                    rval = CCtsp_teething_lp (&lp->cuts,
                       &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                       xcount, xlist, x, 10.0, 500, &maxviol, rstate, &z);
                    CCcheck_rval (rval, "CCtsp_teething_lp failed");
                    PROCESS_CUTS("teethed combs")
                }
            }

            if (sel->tighten_pool && newval < oldval+sel->nexttol) {
                rval = no_tighten (ncount, xcount, xlist, x, &ttest, 0.2);
                CCcheck_rval (rval, "no_tighten failed");
                if (!ttest) {
                    qqq = 10;  /* 25 */
                    do {
                        rval = CCtsp_tighten_lp (lp->pool,
                           &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                           ncount, xcount, xlist, x, 0.1, 1000, &maxviol,
                           0, rstate, &z);
                        CCcheck_rval (rval, "CCtsp_tighten_lp failed");
                        PROCESS_CUTS("tighten pool cuts")
                    } while (cutcount && --qqq);
                }
            }

            if (sel->decker_pool && newval < oldval+sel->nexttol) {
                printf ("DDpool ....\n"); fflush (stdout);
                rval = CCtsp_double_decker_lp (lp->pool,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   xcount, xlist, x, 2.0, 1000, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_double_decker_lp failed");
                PROCESS_CUTS("pool double deckers")
            }

            if (sel->star_pool && newval < oldval+sel->nexttol) {
                rval = CCtsp_star_lp (lp->pool, &lp->stats.extra_tighten_stats,
                   &cuts, &cutcount, ncount, xcount, xlist, x, 2.0, 1000,
                   &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_star_lp failed");
                PROCESS_CUTS("pool star inequalities")
            }

            if (sel->handling_pool && newval < oldval+sel->nexttol) {
                rval = CCtsp_handling_lp (lp->pool,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount, ncount,
                   xcount, xlist, x, 2.0, 1000, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_handling_lp failed");
                PROCESS_CUTS("pool handling inequalities")
            }

            if (sel->teething_pool && newval < oldval+sel->nexttol) {
                rval = CCtsp_teething_lp (lp->pool,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                   ncount, xcount, xlist, x, 0.5, 1000, &maxviol, rstate, &z);
                CCcheck_rval (rval, "CCtsp_teething_lp failed");
                PROCESS_CUTS("pool teething combs")
            }

            if (sel->consecutiveones && newval < oldval+sel->nexttol) {
                rval = CCpq_cuttree_improve_quick (&lp->tightcuts, lp->pool,
                            xcount, xlist, x);
                CCcheck_rval (rval, "CCpq_cuttree_improve_quick failed");
                rval = CCpq_consecutiveones (&cuts, &cutcount, &lp->tightcuts,
                            lp->pool, xcount, xlist, x, &z, &maxviol);
                CCcheck_rval (rval, "CCpq_consecutiveones failed");
                PROCESS_CUTS("consecutive ones cuts")
            }

            if (sel->necklace && newval < oldval+sel->nexttol) {
                rval = CCpq_cuttree_improve_quick (&lp->tightcuts, lp->pool,
                            xcount, xlist, x);
                CCcheck_rval (rval, "CCpq_cuttree_improve_quick failed");
                rval = CCpq_necklaces (&cuts, &cutcount, &lp->tightcuts,
                            xcount, xlist, x, rstate, &z, &maxviol);
                CCcheck_rval (rval, "CCpq_necklaces failed");
                PROCESS_CUTS("necklace cuts")
            }

            otherimprove = newval - oldval;
            if (sel->maxchunksize > 0 && newval < oldval+sel->nexttol &&
                otherimprove <  0.5 * lcimprove) {
                int  maxchunksize, firstsize;
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 0;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                else                        firstsize = 8;
                for (maxchunksize = firstsize;
                     maxchunksize <= sel->maxchunksize; maxchunksize++) {
                    flags.maxchunksize = maxchunksize;
                    flags.spheresize   = maxchunksize - 2;

                    rval = CCchunk_localcuts (&cuts, &cutcount, ncount, xcount,
                       xlist, x, 0.0, flags, &lc_timer, silent, rstate, &z,
                       &maxviol);
                    CCcheck_rval (rval, "CCchunk_localcuts failed");
                    PROCESS_CUTS("localcuts")
                    if (newval >= oldval+sel->nexttol)  break;
                }
            }

            if (sel->maxchunksize > 0 && newval < oldval+sel->nexttol) {
                int  maxchunksize, firstsize;
                CCchunk_flag flags;
                double beforeval = newval;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 1;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                if (sel->maxchunksize < 8)  firstsize = sel->maxchunksize;
                else                        firstsize = 8;
                for (maxchunksize = firstsize;
                     maxchunksize <= sel->maxchunksize; maxchunksize++) {
                    flags.maxchunksize = maxchunksize;
                    flags.spheresize = maxchunksize - 2;

                    rval = CCchunk_localcuts (&cuts, &cutcount, ncount, xcount,
                       xlist, x, 0.0, flags, &lc_timer, silent, rstate, &z,
                       &maxviol);
                    CCcheck_rval (rval, "CCchunk_localcuts failed");
                    PROCESS_CUTS("localcuts")
                    if (newval >= oldval + sel->nexttol)  break;

#if 0
                    rval = grab_local_x (ncount, snowtour, xcount, xlist, x,
                                &closecount, &closelist, &closex, 0.50);
                    CCcheck_rval (rval, "grab_local_x failed");
                    rval = CCchunk_localcuts (&cuts, &cutcount, ncount,
                       closecount, closelist, closex, 0.0, flags, &lc_timer,
                       silent, rstate, &z, &maxviol);
                    CCcheck_rval (rval, "CCchunk_localcuts failed");
                    silent = 0;
                    PROCESS_CUTS("CLOSE .50 localcuts")
                    silent = 1;

                    rval = grab_local_x (ncount, snowtour, xcount, xlist, x,
                                &closecount, &closelist, &closex, 0.25);
                    CCcheck_rval (rval, "grab_local_x failed");
                    rval = CCchunk_localcuts (&cuts, &cutcount, ncount,
                       closecount, closelist, closex, 0.0, flags, &lc_timer,
                       silent, rstate, &z, &maxviol);
                    CCcheck_rval (rval, "CCchunk_localcuts failed");
                    silent = 0;
                    PROCESS_CUTS("CLOSE .25 localcuts")
                    silent = 1;
                    if (newval >= oldval + sel->nexttol)  break;

                    rval = grab_local_x (ncount, snowtour, xcount, xlist, x,
                                &closecount, &closelist, &closex, 0.1);
                    CCcheck_rval (rval, "grab_local_x failed");
                    rval = CCchunk_localcuts (&cuts, &cutcount, ncount,
                       closecount, closelist, closex, 0.0, flags, &lc_timer,
                       silent, rstate, &z, &maxviol);
                    CCcheck_rval (rval, "CCchunk_localcuts failed");
                    silent = 0;
                    PROCESS_CUTS("CLOSE .1 localcuts")
                    silent = 1;
                    if (newval >= oldval + sel->nexttol)  break;
#endif

                }
                lcimprove = newval - beforeval;
            }

OUT_LOOP:
            CC_IFFREE (xlist, int);
            CC_IFFREE (x, double);

            CCutil_start_timer (&lp->stats.sparse_edge_check);
            rval = sparse_edge_check (lp, &eginside, &edge_added,
                                      (double *) NULL, silent, rstate);
            CCcheck_rval (rval, "sparse_edge_check failed");
            CCutil_stop_timer (&lp->stats.sparse_edge_check, silent);
            
            if (savelp) {
                rval = CCtsp_write_probfile_sav (lp);
                CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
            }
            if (lp->pool && savelp) {
                if (!silent) {
                    printf ("Write Pool: %d cuts\n", lp->pool->cutcount);
                    fflush (stdout);
                }
                sprintf (buf, "%s.pul", lp->problabel);
                rval = CCtsp_write_cutpool (ncount, buf, lp->pool);
                CCcheck_rval (rval, "CCtsp_wrte_cutpool failed");
            }
            if (lp->pool && sel->remotepool &&
                lp->pool->cutcount > lp->pool->savecount) {
                tval = CCtsp_send_newcuts (ncount, lp->pool, sel->remotehost,
                                           sel->remoteport);
                if (tval) fprintf (stderr, "CCtsp_send_newcuts failed\n");
            }
            if (lp->dominopool && sel->remotedompool &&
                lp->dominopool->cutcount > lp->dominopool->savecount) {
                tval = CCtsp_send_newcuts (ncount, lp->dominopool,
                        sel->remotedomhost, sel->remotedomport);
                if (tval) fprintf (stderr, "CCtsp_send_newcuts failed\n");
            }
            if (lp->dominopool && lp->dominopool->cutcount && savelp) {
                if (!silent) {
                    printf ("Write Domino Pool: %d cuts\n",
                                 lp->dominopool->cutcount);
                    fflush (stdout);
                }
                sprintf (buf, "%s.dompul", lp->problabel);
                rval = CCtsp_write_cutpool (ncount, buf, lp->dominopool);
                CCcheck_rval (rval, "CCtsp_wrte_cutpool failed");
            }
            CCutil_stop_timer (&lp->stats.cutting_inner_loop, silent);
            rval = lp_value (lp, &priceval);
            CCcheck_rval (rval, "lp_value failed");

            if (lp->lowerbound >= lp->upperbound - 0.9) {
                if (!silent) {
                    printf ("Stop cutting, LP within 0.9 of upperbound\n");
                    fflush (stdout);
                }
                goto CLEANUP;
            }
            loopcount++;
            if (silent < 2 && !lp->full_edges_valid) {
                printf ("  LP Value %2d: %f  (%.2f seconds)\n", loopcount,
                     priceval, CCutil_zeit () - szeit);
                fflush (stdout);
            }
        } while ((newval > oldval + sel->roundtol ||
                priceval < newval - sel->roundtol) &&
               loopcount < LOOP_FULL &&
               (lp->full_edges_valid || priceval < lp->upperbound));

        /* STAR Skip full edge check, bb_cutting calls CCtsp_pricing_loop
        goto CLEANUP;
        */

        rval = full_edge_check (lp, &edge_added, (double *) NULL, silent,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        if (savelp) {
            rval = CCtsp_write_probfile_sav (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
        }
        rval = lp_value (lp, &priceval);
        CCcheck_rval (rval, "lp_value failed");

        if (sel->extra_connect && priceval >= newval-sel->roundtol &&
                                 loopcount != LOOP_FULL) {
            if (!silent) {
                printf ("Check connectivity before exiting cutting_loop\n");
                fflush (stdout);
            }

            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");
            rval = CCtsp_connect_cuts (&cuts, &cutcount_connect, ncount,
                                       xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_connect_cuts failed");
            cut_msg (silent, "extra connect cuts", cutcount_connect, z, 0.0);

            rval = call_add_cuts (lp, &cuts, &cut_added, &xcount, &xlist, &x,
               &newval, sel->usetighten, (int *) NULL, silent, rstate,
               infeasible);
            CCcheck_rval (rval, "call_add_cuts failed");
            if (*infeasible) goto CLEANUP;

            CC_FREE (xlist, int);
            CC_FREE (x, double);
        }
        outside++;

        /* This last bit will cause a second pass for large probs, with */
        /* an updated tolerance.                                        */

        if (!sel->fastcuts && lp->full_edges_valid == 0 && outside == 1 &&
             ncount >= 400 && lp->lowerbound < lp->upperbound - 0.9) {
            rval = CCtsp_cutselect_set_tols (sel, lp, 1, silent);
            CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
            loopcount = LOOP_FULL;  /* to run again */
        }
        kloops++;
    } while (priceval < newval - sel->roundtol || loopcount == LOOP_FULL ||
             cutcount_connect ||
            (sel->keep_cutting && totalcuts > 0 && kloops < KEEP_MAX));

CLEANUP:
    if (!rval && *infeasible) {
        if (silent < 2) {
            printf ("LP is infeasible in cutting_loop\n"); fflush (stdout);
        }
    }
    CCutil_stop_timer (&lp->stats.cutting_loop, silent);
    if (!silent) {
        printf ("Number of outside rounds: %d\n", outside); fflush (stdout);
    }

    if (eginside.ncount) CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    CC_IFFREE (closelist, int);
    CC_IFFREE (closex, double);
    sel->nexttol = save_nexttol;
    sel->roundtol = save_roundtol;
    return rval;
}

static void cut_msg (int silent, const char *msg, int count, double sec,
        double maxviol)
{
    if (!silent && count > 0) {
        if (maxviol != 0.0) {
            printf ("Found %2d %s (max viol %.4f) in %.2f seconds\n",
                                               count, msg, maxviol, sec);
        } else {
            printf ("Found %2d %s in %.2f seconds\n", count, msg, sec);
        }
        fflush (stdout);
    }
}

int CCtsp_test_general_cuts (CCtsp_lp *lp, CCrandstate *rstate);

int CCtsp_test_general_cuts (CCtsp_lp *lp, CCrandstate *rstate)
{
    int rval = 0, silent = 0, ncount = lp->graph.ncount, infeasible = 0;
    int xcount, cutcount, cut_added, *xlist = (int *) NULL, edge_added;
    double *x = (double *) NULL, newval, z, ztol = 0.001, maxviol;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;

    printf ("Testing the general form of cutting planes\n");
    fflush (stdout);

    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, ncount, lp->dat,
                (CCtsp_genadj *) NULL, 100, 0, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    }

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
              (int *) NULL, (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

    do {
        edge_added = 0;
        do {
            cut_added = 0;
            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

            rval = CCtsp_connect_cuts (&cuts, &cutcount, ncount, xcount, xlist,
                                       x, &z);
            CCcheck_rval (rval, "CCtsp_connect_cuts failed");
            printf ("Found %2d connect cuts in %.2f seconds\n", cutcount, z);
            fflush (stdout);

            if (cutcount) {
                rval = subtours_to_general (cuts, ncount);
                CCcheck_rval (rval, "subtours_to_general failed");
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                         &xlist, &x, &newval, 0, (int *) NULL, silent, rstate,
                         &infeasible);
                CCcheck_rval (rval, "call_add_cuts failed");
                if (infeasible) goto CLEANUP;
            } else {
                rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
                                     xcount, xlist, x, &ztol, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_exact_subtours failed");
                printf ("Found %2d exact in %.2f seconds\n", cutcount, z);

                rval = subtours_to_general (cuts, ncount);
                CCcheck_rval (rval, "subtours_to_general failed");
                rval = call_add_cuts (lp, &cuts, &cut_added, &xcount,
                         &xlist, &x, &newval, 0, (int *) NULL, silent, rstate,
                         &infeasible);
                CCcheck_rval (rval, "call_add_cuts failed");
                if (infeasible) goto CLEANUP;
            }

            /* CClp_dump_lp (lp->lp, "dump.lp"); */
        } while (cut_added);
        rval = full_edge_check (lp, &edge_added, (double *) NULL, silent,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
    } while (edge_added);

CLEANUP:
    if (eginside.ncount) CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

static int subtours_to_general (CCtsp_lpcut_in *cuts, int ncount)
{
    int rval = 0;
    CCtsp_lpcut_in *c;

    c = cuts;
    while (c) {
        int *a, acount, *mark = (int *) NULL, i, j, *ar = (int *) NULL;
        int arcount, nsemi, s[1], t[1];
        CCtsp_lpdomino *d;

        if (!CCtsp_hypergraph_type_in (c)) {
            printf ("Not a hypergraph!\n"); rval = 1; goto CLEANUP;
        }
        if (c->cliquecount != 1) {
            printf ("Not a subtour cut!\n"); rval = 1; goto CLEANUP;
        }

#if 0
        /* Convert subtour to mults */
        CC_MALLOC (c->cliquemult, 1, int);
        c->cliquemult[0] = 2;
        c->rhs = 4;
#endif

#if 0
        /* convert to outside form with semicuts */
        rval = CCtsp_clique_to_array (&c->cliques[0], &a, &acount);
        CCcheck_rval (rval, "CCtsp_clique_to_array failed");

        CC_MALLOC (mark, ncount, int);
        for (i = 0; i < ncount; i++) mark[i] = 0;
        for (i = 0; i < acount; i++) mark[a[i]] = 1;

        for (i = 0; i < acount; i++) {
            for (j = 0; j < ncount; j++) {
                if (mark[j] == 0) nsemi++;
            }
        }
        CC_MALLOC (d, nsemi, CCtsp_lpdomino);
        nsemi = 0;
        for (i = 0; i < acount; i++) {
            for (j = 0; j < ncount; j++) {
                if (mark[j] == 0) {
                    s[0] = a[i];
                    t[0] = j;
                    rval = CCtsp_array_to_lpclique (s,1,&d[nsemi].sets[0]);
                    CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                    rval = CCtsp_array_to_lpclique (t,1,&d[nsemi].sets[1]);
                    CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                    nsemi++;
                }
            }
        }

        c->semicount = nsemi;
        c->semicuts = d;
        c->rhs = 2;
        c->cliques = (CCtsp_lpclique *) NULL;
        c->cliquecount = 0;

        CC_FREE (mark, int); 
        CC_FREE (a, int); 
#endif

#if 1
        /* convert to inside form with negative semicuts */
        rval = CCtsp_clique_to_array (&c->cliques[0], &a, &acount);
        CCcheck_rval (rval, "CCtsp_clique_to_array failed");
        if (acount > (ncount+1)/2) {
            CC_MALLOC (mark, ncount, int);
            for (i = 0; i < ncount; i++) mark[i] = 0;
            for (i = 0; i < acount; i++) mark[a[i]] = 1;
            CC_MALLOC (ar, ncount-acount, int);
            arcount = 0;
            for (i = 0; i < ncount; i++){
                if (mark[i] == 0) {
                    ar[arcount++] = i;
                }
            }
            CC_FREE (a, int);
            a = ar;
            acount = arcount;
            CC_FREE (mark, int);
        }

/*
        printf ("acount = %d\n", acount); fflush (stdout);
        for (i = 0; i < acount; i++) printf ("%d ", a[i]);
        printf ("\n");
*/

        nsemi = (acount * (acount-1)) / 2;
        CC_MALLOC (d, nsemi, CCtsp_lpdomino);
        nsemi = 0;
        for (i = 0; i < acount; i++) {
            for (j = i+1; j < acount; j++) {
                s[0] = a[i];
                t[0] = a[j];
                rval = CCtsp_array_to_lpclique (s, 1, &d[nsemi].sets[0]);
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                rval = CCtsp_array_to_lpclique (t, 1, &d[nsemi].sets[1]);
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                nsemi++;
            }
        }

        c->semicount = nsemi;
        c->semicuts = d;
        CC_MALLOC (c->semimult, nsemi, int);
        for (i = 0; i < nsemi; i++) c->semimult[i] = -1;
        c->rhs = -(acount-1);
        c->cliques = (CCtsp_lpclique *) NULL;
        c->cliquecount = 0;

        CC_FREE (a, int); 
#endif
        c = c->next;
    }

CLEANUP:
   return rval;
}


#define CC_NO_NEAREST_SUBTOUR 50
#define CC_SUBTOUR_ROUNDS     5

#define PROCESS_LOOP_CUTS(xmsg) {                                           \
    cut_msg (silent, xmsg, cutcount, z, maxviol);                           \
    if (cutcount) {                                                         \
        rval = call_add_cuts (lp, &cuts, &tcut_added, &xcount, &xlist,      \
           &x, &newval, tighten, (int *) NULL, silent, rstate, infeasible); \
        CCcheck_rval (rval, "call_add_cuts failed");                        \
        if (*infeasible) goto CLEANUP;                                      \
        cut_added += tcut_added;                                            \
    }                                                                       \
}

#define PROCESS_LOOP_CUTS_0(xmsg) {                                         \
    cut_msg (silent, xmsg, cutcount, z, 0.0);                               \
    if (cutcount) {                                                         \
        rval = call_add_cuts (lp, &cuts, &tcut_added, &xcount, &xlist,      \
           &x, &newval, tighten, (int *) NULL, silent, rstate, infeasible); \
        CCcheck_rval (rval, "call_add_cuts failed");                        \
        if (*infeasible) goto CLEANUP;                                      \
        cut_added += tcut_added;                                            \
    }                                                                       \
}

int CCtsp_subtour_loop (CCtsp_lp *lp, int silent, double ztol,
        CCrandstate *rstate, int *infeasible)
{
    int rval = 0, outside = 0, inside = 0, tighten = 0;
    int xcount, cutcount, cut_added, tcut_added, edge_added, savelp = 1;
    int *xlist = (int *) NULL;
    double newval, priceval, z, maxviol, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;

    *infeasible = 0;
    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    }

    do {
        do {
            cut_added = 0;
            CCutil_start_timer (&lp->stats.cutting_inner_loop);
            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

            /**** Connect Cuts ****/

            rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                       xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_connect_cuts failed");
            PROCESS_LOOP_CUTS_0("connect cuts")

            /**** Shrink Cuts ****/

            rval = CCtsp_shrink_subtours (&cuts, &cutcount, lp->graph.ncount,
                                xcount, xlist, x, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_shrink_subtours failed");
            PROCESS_LOOP_CUTS("shrink subtours")

            /**** Linear Cuts ****/

            rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                      xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_segment_cuts failed");
            PROCESS_LOOP_CUTS_0("segment cuts")

            /**** Exact Cuts ****/

            rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
               xcount, xlist, x, &ztol, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_exact_subtours failed");
            PROCESS_LOOP_CUTS("exact subtours")

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                CCcheck_rval (rval, "sparse_edge_check failed");
                CCutil_stop_timer (&lp->stats.sparse_edge_check, silent);
            }
            CCutil_stop_timer (&lp->stats.cutting_inner_loop, silent);
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                CCcheck_rval (rval, "lp_value failed");
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);

        rval = full_edge_check (lp, &edge_added, (double *) NULL, silent,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        outside++;
        if (savelp) {
            rval = CCtsp_write_probfile_sav (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
        }
    } while (edge_added);

CLEANUP:

    if (*infeasible) {
        printf ("LP is infeasible in subtour_loop\n"); fflush (stdout);
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    printf ("Time in cutting routine: %.2f\n", z);
    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount) CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

int CCtsp_subtour_loop_rc (CCtsp_lp *lp, int silent, double ztol,
        CCrandstate *rstate, int *infeasible)
{
    int rval = 0, inside = 0, tighten = 0;
    int xcount, cutcount, cut_added, tcut_added;
    int *xlist = (int *) NULL;
    double newval, z, maxviol, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCbigguy bound;

    *infeasible = 0;
    CCutil_start_timer (&lp->stats.cutting_loop);

    do {
        cut_added = 0;
        rval = lp_x (lp, &xcount, &xlist, &x);
        CCcheck_rval (rval, "lp_x failed");

        /**** Connect Cuts ****/

        rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                   xcount, xlist, x, &z);
        CCcheck_rval (rval, "CCtsp_connect_cuts failed");
        PROCESS_LOOP_CUTS_0("connect cuts")

        /**** Shrink Cuts ****/

        rval = CCtsp_shrink_subtours (&cuts, &cutcount, lp->graph.ncount,
                            xcount, xlist, x, &z, &maxviol);
        CCcheck_rval (rval, "CCtsp_shrink_subtours failed");
        PROCESS_LOOP_CUTS("shrink subtours")

        /**** Linear Cuts ****/

        rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, &z);
        CCcheck_rval (rval, "CCtsp_segment_cuts failed");
        PROCESS_LOOP_CUTS_0("segment cuts")

        /**** Exact Cuts ****/

        rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
           xcount, xlist, x, &ztol, &z, &maxviol);
        CCcheck_rval (rval, "CCtsp_exact_subtours failed");
        PROCESS_LOOP_CUTS("exact subtours")

        CC_FREE (xlist, int);
        CC_FREE (x, double);

        inside++;
        if (silent) {
            rval = lp_value (lp, &newval);
            CCcheck_rval (rval, "lp_value failed");
            printf ("  LP Value %2d: %f\n", inside, newval);
            fflush (stdout);
        }
    } while (cut_added);

    printf ("Call exact pricing ...\n"); fflush (stdout);

    rval = CCtsp_exact_price (lp, &bound, 0, 0, 0, 0, silent);
    CCcheck_rval (rval, "CCtsp_exact_price failed");
    printf ("Exact bound: %.6f\n", CCbigguy_bigguytod (bound));

CLEANUP:

    if (*infeasible) {
        printf ("LP is infeasible in subtour_loop\n"); fflush (stdout);
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    printf ("Time in cutting routine: %.2f\n", z);
    printf ("Number of rounds: %d\n", inside); fflush (stdout);

    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

int CCtsp_blossom_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate,
        int *infeasible)
{
    int rval = 0, outside = 0, inside = 0, tighten = 0, *xlist = (int *) NULL;
    int xcount, cutcount, cut_added, tcut_added, edge_added;
    double z, maxviol, newval, priceval, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;

    *infeasible = 0;
    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    }

    do {
        do {
            cut_added = 0;
            CCutil_start_timer (&lp->stats.cutting_inner_loop);
            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

            /****  Fast Blossoms ****/

            rval = CCtsp_fastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_fastblossom failed");
            PROCESS_LOOP_CUTS_0("fast blossoms")

            /****  Groetschel-Holland Fast Blossoms ****/
 
            rval = CCtsp_ghfastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_ghfastblossom failed");
            PROCESS_LOOP_CUTS("Groetschel-Holland blossoms")

            /**** Exact Blossoms ****/

            rval = CCtsp_exactblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, rstate, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_exact_blossom failed");
            PROCESS_LOOP_CUTS("exact blossoms")

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                CCcheck_rval (rval, "sparse_edge_check failed");
                CCutil_stop_timer (&lp->stats.sparse_edge_check, silent);
            }
            CCutil_stop_timer (&lp->stats.cutting_inner_loop, silent);
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                CCcheck_rval (rval, "lp_value failed");
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);

        rval = full_edge_check (lp, &edge_added, (double *) NULL, silent,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        outside++;
    } while (edge_added);

CLEANUP:
    if (*infeasible) {
        printf ("LP is infeasible in blossom_loop\n"); fflush (stdout);
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    printf ("Time in cutting routine: %.2f\n", z);
    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount) CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

int CCtsp_subtour_and_blossom_loop (CCtsp_lp *lp, int silent,
        CCrandstate *rstate, char *filecutname, int *infeasible)
{
    int rval, outside = 0, inside = 0, tighten = 0, *xlist = (int *) NULL;
    int xcount, cutcount, cut_added, tcut_added, edge_added;
    double z, maxviol, oldval, newval, priceval, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;

    if (filecutname) {
        printf ("Taking cuts from file: %s\n", filecutname); fflush (stdout);
    }

    *infeasible = 0;
    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator (sparse) failed");
    }

    do {
        do {
            CCutil_start_timer (&lp->stats.cutting_inner_loop);
            cut_added = 0;

            rval = lp_value (lp, &oldval);
            CCcheck_rval (rval, "lp_value failed");
            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

            /**** Connect Cuts ****/

            rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                       xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_connect_cuts failed");
            PROCESS_LOOP_CUTS_0("connect cuts");

            /**** Linear Cuts ****/

            rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                      xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_segment_cuts failed");
            PROCESS_LOOP_CUTS_0("segment cuts");

            /**** Exact Cuts ****/

            rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
                              xcount, xlist, x, (double *) NULL, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_exact_subtours failed");
            PROCESS_LOOP_CUTS("exact subtours");

            /****  Fast Blossoms ****/

            rval = CCtsp_fastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_fastblossom failed");
            PROCESS_LOOP_CUTS_0("fast blossoms");

#if 1 /* Start use exact blossoms */
            /****  Groetschel-Holland Fast Blossoms ****/

            rval = CCtsp_ghfastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_ghfastblossom failed");
            PROCESS_LOOP_CUTS("Groetschel-Holland blossoms");

            /**** Exact Blossoms ****/

            rval = CCtsp_exactblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, rstate, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_exactblossom failed");
            PROCESS_LOOP_CUTS("exact blossoms");
#endif  /* End use exact blossoms */

#if 1   /* Start use of combs */
            rval = CCtsp_block_combs (&cuts, &cutcount, lp->graph.ncount,
               xcount, xlist, x, silent, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_block_combs failed");
            PROCESS_LOOP_CUTS("block combs");

            rval = CCtsp_edge_comb_grower (&cuts, &cutcount, lp->graph.ncount,
               xcount, xlist, x, &lp->stats.extra_tighten_stats, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_edge_comb_grower failed");
            PROCESS_LOOP_CUTS("grow combs");

            rval = CCpq_cuttree_improve_quick (&lp->tightcuts, lp->pool,
                        xcount, xlist, x);
            CCcheck_rval (rval, "CCpq_cuttree_improve_quick failed");
            rval = CCpq_consecutiveones (&cuts, &cutcount, &lp->tightcuts,
                        lp->pool, xcount, xlist, x, &z, &maxviol);
            CCcheck_rval (rval, "CCpq_consecutiveones failed");
            PROCESS_LOOP_CUTS("consec one combs");

            rval = CCpq_cuttree_improve_quick (&lp->tightcuts, lp->pool,
                        xcount, xlist, x);
            CCcheck_rval (rval, "CCpq_cuttree_improve_quick failed");

            rval = CCpq_necklaces (&cuts, &cutcount, &lp->tightcuts,
                        xcount, xlist, x, rstate, &z, &maxviol);
            CCcheck_rval (rval, "CCpq_necklaces failed");
            PROCESS_LOOP_CUTS("necklace combs");

            if (filecutname) {
                rval = CCtsp_file_cuts (filecutname, &cuts, &cutcount,
                           lp->graph.ncount, lp->perm, &z);
                CCcheck_rval (rval, "CCtsp_file_cuts failed");
                PROCESS_LOOP_CUTS_0("file cuts");
            }
#endif  /* End use of combs */

#if 0 /* Start use of pool */
            if (lp->pool) {
                int qqq = 10, ncount = lp->graph.ncount;
                do {
                    rval = CCtsp_search_cutpool (lp->pool, &cuts, &cutcount,
                                &maxviol, ncount, xcount, xlist, x, 0,
                                rstate, &z);
                    CCcheck_rval (rval, "CCtsp_search_cutpool failed");
                    PROCESS_LOOP_CUTS("pool cuts");
                } while (cutcount && --qqq);
            }

            if (lp->dominopool) {
                int qqq = 10, ncount = lp->graph.ncount;
                int domstarter = 0;
                do {
                    rval = CCtsp_search_dominopool (lp->dominopool, &cuts,
                       &cutcount, &maxviol, ncount, xcount, xlist, x,
                       &domstarter, 10000, 500, rstate, &z);
                    CCcheck_rval (rval, "CCtsp_search_dominopool failed");
                    PROCESS_LOOP_CUTS("domino pool cuts");
                } while (cutcount && --qqq);

                qqq = 1;
                do {
                    rval = CCtsp_domino_tighten_pool (lp->dominopool, &cuts,
                       &cutcount, ncount, xcount, xlist, x, 0.01, 1000,
                       &maxviol, rstate, &z);
                    CCcheck_rval (rval, "CCtsp_domino_tighten_pool failed");
                    PROCESS_LOOP_CUTS("domino tighten pool cuts");
                } while (cutcount && --qqq);
            }


            /****  Tighten ****/

            {
                printf ("Tighten cuts ...\n"); fflush (stdout);
                rval = CCtsp_tighten_lp (&lp->cuts,
                   &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                   lp->graph.ncount, xcount, xlist, x, 0.5, 2000, &maxviol,
                   0, rstate, &z);
                CCcheck_rval (rval, "CCtsp_tighten_lp");
                PROCESS_LOOP_CUTS("tighten cuts");

                printf ("Tighten pool cuts ...\n"); fflush (stdout);
                rval = CCtsp_tighten_lp (lp->pool,
                        &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                        lp->graph.ncount, xcount, xlist, x, 0.1, 2000, &maxviol,
                        0, rstate, &z);
                CCcheck_rval (rval, "CCtsp_tighten_lp");
                PROCESS_LOOP_CUTS("tighten pool cuts");
            }
#endif  /* End use pool */

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                CCcheck_rval (rval, "sparse_edge_check failed");
                CCutil_stop_timer (&lp->stats.sparse_edge_check, silent);
            }
            CCutil_stop_timer (&lp->stats.cutting_inner_loop, silent);
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                CCcheck_rval (rval, "lp_value failed");
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);
#if 0   /* STOP ACTION for Kannan 60 meeting and save LP */
            rval = CCtsp_write_probfile_sav (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
            rval = lp_value (lp, &newval);
            CCcheck_rval (rval, "lp_value failed");
        } while (edge_added || newval > oldval + .1);
#endif  /* END STOP ACTION for Kannan 60 meeting*/

        rval = full_edge_check (lp, &edge_added, (double *) NULL, silent,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        outside++;
    } while (edge_added);

    rval = sparse_edge_check (lp, &eginside, &edge_added, (double *) NULL,
                              silent, rstate);
    CCcheck_rval (rval, "sparse_edge_check failed");
    rval = CCtsp_write_probfile_sav (lp);
    CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");

CLEANUP:
    if (*infeasible) {
        printf ("LP is infeasible in subtour_and_blossom_loop\n");
        fflush (stdout);
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    printf ("Time in cutting routine: %.2f\n", z);
    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount)
        CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

int CCtsp_domino_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate, 
       int *infeasible)
{
    int rval, outside = 0, inside = 0, tighten = 0, *xlist = (int *) NULL;
    int xcount, cutcount, cut_added, tcut_added = 0, edge_added, qqq = 0;
    int ncount = lp->graph.ncount, domstarter = 0;
    double z, maxviol, oldval, newval, priceval, qval, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;

    printf ("CCtsp_domino_loop ...\n"); fflush (stdout);

    *infeasible = 0;
    CCutil_start_timer (&lp->stats.cutting_loop);
    eginside.ncount = 0;

#ifndef CCtsp_USE_DOMINO_CUTS
    printf ("Code has not been compiled for DP/TP cuts\n"); fflush (stdout);
    rval = 1;  goto CLEANUP;
#endif

    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator (sparse) failed");
    }

    do {
        do {
            CCutil_start_timer (&lp->stats.cutting_inner_loop);
            cut_added = 0;
            rval = lp_value (lp, &oldval);
            CCcheck_rval (rval, "lp_value failed");
            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

#ifdef CCtsp_USE_DOMINO_CUTS

            /* DP and TP separators need to be nearly in subtour polytope */
            qval = oldval;
            qqq = 25;
            do {
                rval = CCtsp_exact_subtours (&cuts, &cutcount, ncount,
                   xcount, xlist, x, (double *) NULL, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_exact_subtours failed");
                PROCESS_LOOP_CUTS("exact subtour")
                if (newval > qval + 0.5) { qval = newval;  qqq = 26; }
            } while (tcut_added > 0 && --qqq);

            /* DP and TP separators definitely need connected solutions  */
            qqq = 100;
            do {
                rval = CCtsp_connect_cuts (&cuts, &cutcount, ncount,
                                           xcount, xlist, x, &z);
                CCcheck_rval (rval, "CCtsp_connect_cuts failed");
                PROCESS_LOOP_CUTS_0("connect subtours")
            } while (tcut_added > 0 && --qqq);

            rval = CCtsp_DP_cuts (&cuts, &cutcount, ncount, xcount,
                     xlist, x, 1, (char *) NULL, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_DP_cuts failed");
            PROCESS_LOOP_CUTS("DP cuts")

            if (tcut_added == 0) {
                rval = CCtsp_TP_cuts (&cuts, &cutcount, ncount, xcount, xlist,
                                      x, &z, &maxviol);
                CCcheck_rval (rval, "CCtsp_TP_cuts failed");
                PROCESS_LOOP_CUTS("TP cuts")
            }

            rval = CCtsp_search_dominopool (lp->dominopool, &cuts,
                       &cutcount, &maxviol, ncount, xcount, xlist, x,
                       &domstarter, 10000, 500, rstate, &z);
            CCcheck_rval (rval, "CCtsp_search_dominopool failed");
            PROCESS_LOOP_CUTS("domino pool cuts");

            rval = CCtsp_domino_tighten_pool (lp->dominopool, &cuts,
                       &cutcount, ncount, xcount, xlist, x, 0.01, 1000,
                       &maxviol, rstate, &z);
            CCcheck_rval (rval, "CCtsp_domino_tighten_pool failed");
            PROCESS_LOOP_CUTS("domino tighten pool cuts");

#endif /* CCtsp_USE_DOMINO_CUTS */

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % CC_SUBTOUR_ROUNDS) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                CCcheck_rval (rval, "sparse_edge_check failed");
                CCutil_stop_timer (&lp->stats.sparse_edge_check, silent);
            }
            CCutil_stop_timer (&lp->stats.cutting_inner_loop, silent);
            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                CCcheck_rval (rval, "lp_value failed");
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);

        rval = full_edge_check (lp, &edge_added, (double *) NULL, silent,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        outside++;
    } while (edge_added);

    rval = CCtsp_write_probfile_sav (lp);
    CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");

CLEANUP:
    if (*infeasible) {
        printf ("LP is infeasible in domino_loop\n"); fflush (stdout);
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    printf ("Time in cutting routine: %.2f\n", z);
    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount) CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

static int call_add_cuts (CCtsp_lp *lp, CCtsp_lpcut_in **cuts, int *cut_added,
        int *xcount, int **xlist, double **x, double *val, int tighten,
        int *istour, int silent, CCrandstate *rstate, int *infeasible)
{
    int rval = 0;
    double dval, szeit = CCutil_zeit ();

    if (istour) *istour = 0;
    if (cut_added) *cut_added = 0;
    if (*cuts == (CCtsp_lpcut_in *) NULL) goto CLEANUP;

    CC_IFFREE (*xlist, int);
    CC_IFFREE (*x, double);

    CCtsp_add_cuts_to_queue (lp, cuts);
    rval = CCtsp_process_cuts (lp, cut_added, tighten, silent, rstate,
                              (double *) NULL, infeasible);
    CCcheck_rval (rval, "CCtsp_process_cuts failed");
    if (*infeasible) {
        printf ("Infeasible LP after processing cuts\n");  fflush (stdout);
        goto CLEANUP;
    }

    rval = lp_value (lp, val);
    CCcheck_rval (rval, "lp_value failed");
    if (!silent) {
        printf ("  Add %2d cuts (Total %d), LP: %f (%.2f seconds)\n",
           *cut_added, lp->cuts.cutcount, *val, CCutil_zeit () - szeit);
        fflush (stdout);
    }

    rval = lp_x (lp, xcount, xlist, x);
    CCcheck_rval (rval, "lp_x failed");

    if (istour) {
        rval = CCtsp_check_integral (lp, &dval, (int **) NULL, istour, silent);
        CCcheck_rval (rval, "CCtsp_check_integral failed");
    } 

CLEANUP:
    return rval;
}

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval = 0;

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}

static int lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x)
{
    int rval = 0;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, xcount,
                     xlist, x, (double **) NULL, (double **) NULL,
                     (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}

static int lp_upperbound (CCtsp_lp *lp, double *ub)
{
    int rval = 0;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, ub, (int *) NULL,
              (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}

int CCtsp_pricing_loop (CCtsp_lp *lp, double *bnd, int silent,
        CCrandstate *rstate)
{
    int rval = 0, nadded;
    CCtsp_edgegenerator eg;

    eg.ncount = 0;

    if (!lp->full_edges_valid) {
        if (!lp->dat) {
            fprintf (stderr, "CCtsp_pricing_loop called without datagroup\n");
            rval = 1; goto CLEANUP;
        }
        rval = full_edge_check (lp, &nadded, bnd, silent, rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        goto CLEANUP;
    }

    rval = CCtsp_init_edgegenerator (&eg, lp->graph.ncount, lp->dat,
                                     lp->fulladj, 0, silent, rstate);
    CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");
    rval = sparse_edge_check (lp, &eg, &nadded, bnd, silent, rstate);
    CCcheck_rval (rval, "sparse_edge_check failed");

CLEANUP:
    if (eg.ncount) CCtsp_free_edgegenerator (&eg);
    return rval;
}

static int full_edge_check (CCtsp_lp *lp, int *nadded, double *bnd, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    double val, penalty;
    CCtsp_edgegenerator eg;

    eg.ncount = 0;
    *nadded = 0;
    if (bnd) *bnd = -CCtsp_LP_MAXDOUBLE;

    if (lp->full_edges_valid) goto CLEANUP;

    if (!lp->dat) {
        fprintf (stderr, "no datagroup available\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_start_timer (&lp->stats.full_edge_check);

    rval = CCtsp_init_edgegenerator (&eg, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CCtsp_PRICE_COMPLETE_GRAPH, silent,
                rstate);
    CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");

    rval = CCtsp_addbad_variables (lp, &eg, &penalty, nadded,
                  CCtsp_PRICE_RCTHRESH, CCtsp_PRICE_MAXPENALTY, 0,
                  (int *) NULL, silent, rstate);
    CCcheck_rval (rval, "CCtsp_addbad_variables failed");

    if (!silent) {
        printf ("%d edges added, penalty %f\n", *nadded, penalty);
        fflush (stdout);
    }

    rval = lp_value (lp, &val);
    CCcheck_rval (rval, "lp_value failed");
    if (bnd) *bnd = val + penalty;

    if (val + penalty > lp->lowerbound) {
        if (silent < 2) {
            printf ("New lower bound: %f\n", val+ penalty); fflush (stdout);
        }
        lp->lowerbound = val + penalty;
    }

    CCutil_stop_timer (&lp->stats.full_edge_check, silent);

CLEANUP:
    if (eg.ncount) CCtsp_free_edgegenerator (&eg);
    return rval; 
}

static int sparse_edge_check (CCtsp_lp *lp, CCtsp_edgegenerator *eg,
        int *nadded, double *bnd, int silent, CCrandstate *rstate)
{
    int rval = 0;
    double val, penalty;

    if (nadded) *nadded = 0;
    if (bnd) *bnd = -CCtsp_LP_MAXDOUBLE;

    if (eg->ncount > 0) {
        rval = CCtsp_addbad_variables (lp, eg, &penalty, nadded,
                  CCtsp_PRICE_RCTHRESH, CCtsp_PRICE_MAXPENALTY, 0,
                  (int *) NULL, silent, rstate);
        CCcheck_rval (rval, "CCtsp_addbad_variables failed");

        rval = lp_value (lp, &val);
        CCcheck_rval (rval, "lp_value failed");

        if (!silent) {
            printf ("(SPARSE) %d edges added, penalty %f, val %f\n",
                      *nadded, penalty, val);
            fflush (stdout);
        }

        if (lp->full_edges_valid) {
            if (val + penalty > lp->lowerbound) {
                if (!silent) {
                    printf ("New (node) lower bound: %f\n", val + penalty);
                    fflush (stdout);
                }
                lp->lowerbound = val + penalty;
            }
            if (bnd) *bnd = val + penalty;
        }
    }

CLEANUP:
    return rval;
}

int CCtsp_bb_cutting (char *probname, int probnum, int prob_newnum, int ncount,
        CCdatagroup *dat, int *ptour, double *upbound, CCtsp_lpcuts *pool,
        CCtsp_lpcuts *dominopool, CCtsp_cutselect *sel, double *val,
        int *prune, int *foundtour, int *besttour, int level, int silent,
        CCrandstate *rstate)
{
    int rval = 0, test;
    double cval, tourval;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;

    /* Note: level == 0 is tentative cutting, 1 is normal, -1 fast */

    *val = 0.0;
    *prune = 0;
    *foundtour = 0;

    rval = bb_cutting_work (&lp, probname, probnum, ncount, dat, ptour,
       *upbound, pool, dominopool, sel, &cval, level, silent, rstate);
    CCcheck_rval (rval, "bb_cutting_work failed");

    if (lp != (CCtsp_lp *) NULL) { lp->id = prob_newnum; }

    if (cval == CCtsp_LP_MAXDOUBLE) {
        rval = CCtsp_verify_infeasible_lp (lp, &test, silent);
        CCcheck_rval (rval, "CCtsp_verify_infeasible_lp failed");
        if (test) {
            printf ("verified infeasible LP\n"); fflush (stdout);
            *val = CCtsp_LP_MAXDOUBLE;
            *prune = 1;
            rval = CCtsp_write_probleaf_id (lp);
            CCcheck_rval (rval, "CCtsp_write_probleaf_id failed");
        } else {
            fprintf (stderr, "did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        rval = CCtsp_pricing_loop (lp, val, silent, rstate);
        CCcheck_rval (rval, "CCtsp_pricing_loop failed");
        lp->lowerbound = *val;
        if (lp->upperbound < *upbound) *upbound = lp->upperbound;

        /* US50K turn off the heuristic */
        if (lp->lowerbound < lp->upperbound - 0.9) {
            CCutil_start_timer (&lp->stats.linkern);
            rval = CCtsp_call_x_heuristic (lp, &tourval, besttour, silent,
                                           rstate);
            CCcheck_rval (rval, "CCtsp_call_x_heuristic failed");
            CCutil_stop_timer (&lp->stats.linkern, silent);
            if (tourval < lp->upperbound) {
                printf ("New upperbound from x-heuristic: %.2f\n", tourval);
                lp->upperbound = tourval;
                *upbound = tourval;
                *foundtour = 1;
            }
        }

        if (lp->lowerbound >= lp->upperbound - 0.9) {
            rval = CCtsp_verify_lp_prune (lp, &test,  silent);
            CCcheck_rval (rval, "CCtsp_verify_lp_prune failed");
            if (test) {
                if (!silent) {
                    printf ("verified that LP can be pruned\n");
                    fflush (stdout);
                }
                *prune = 1;
                rval = CCtsp_write_probleaf_id (lp);
                CCcheck_rval (rval, "CCtsp_write_probleaf_id failed");
            } else {
                fprintf (stderr, "exact pricing could not prune the search\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            rval = CCtsp_write_probfile_id (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_id failed");
        }
    }

CLEANUP:
    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

int CCtsp_call_x_heuristic (CCtsp_lp *lp, double *val, int *outcyc,
        int silent, CCrandstate *rstate)
{
    int rval = 0, ncount = lp->graph.ncount, xcount, i;
    int *cyc = (int *) NULL, *xlist = (int *) NULL;
    double *x = (double *) NULL;

    *val = CCtsp_LP_MAXDOUBLE;

    if (!lp->dat) goto CLEANUP;

    CC_MALLOC (cyc, ncount, int);
    rval = lp_x (lp, &xcount, &xlist, &x);
    CCcheck_rval (rval, "lp_x failed");
    
    rval = CCtsp_x_greedy_tour_lk (lp->dat, ncount, xcount, xlist, x,
                   cyc, val, silent, rstate);
    CCcheck_rval (rval, "CCtsp_x_greedy_tour failed");
    if (!silent) {
        printf ("x-heuristic lk  gives: %.2f\n", *val); fflush (stdout);
    }
    if (*val < lp->upperbound) {
        if (outcyc) {
            for (i = 0; i < ncount; i++) { outcyc[i] = cyc[i]; }
        }
    }

CLEANUP:
    CC_IFFREE (cyc, int);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

static int bb_cutting_work (CCtsp_lp **lp, char *probname, int probnum,
        int ncount, CCdatagroup *dat, int *ptour, double initial_ub,
        CCtsp_lpcuts *pool, CCtsp_lpcuts *dominopool, CCtsp_cutselect *sel,
        double *val, int level, int silent, CCrandstate *rstate)
{
    int rval = 0, infeasible = 0;

    *lp = (CCtsp_lp *) NULL;
    *val = 0.0;

    rval = CCtsp_bb_init_lp (lp, probname, probnum, ncount, dat, ptour,
               initial_ub, pool, dominopool, silent, rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_bb_init_lp failed");
    if (infeasible) {
        printf ("LP is reported to be infeasible\n"); fflush (stdout);
        *val = CCtsp_LP_MAXDOUBLE;
        rval = 0; goto CLEANUP;
    } 
    CCutil_start_timer (&(*lp)->stats.total);

    if ((*lp)->lowerbound >= (*lp)->upperbound - 0.9) {
        printf ("Do not cut, the lp is within 1.0 of the upperbound\n");
        fflush (stdout);
        *val = (*lp)->lowerbound;
        goto CLEANUP;
    } else {

#if 0 /* US50K ELIM: Call edge elimination before cutting task */
    {
        CCbigguy bound;
        rval = CCtsp_exact_price (*lp, &bound, 0, 0, 1, 0, silent);
        CCcheck_rval (rval, "CCtsp_exact_price failed");
        printf ("Exact bound after elim: %.6f\n", CCbigguy_bigguytod (bound));
    }
#endif /* END US50K */

        rval = CCtsp_cutselect_set_tols (sel, *lp, level, silent);
        CCcheck_rval (rval, "CCtsp_cutselect_tols failed");
        rval = CCtsp_cutting_loop (*lp, sel, NULL, 0, silent, rstate,
                                   &infeasible);
        CCcheck_rval (rval, "CCtsp_cutting_loop failed");
        if (infeasible) {
            printf ("Cut LP is reported to be infeasible\n"); fflush (stdout);
            *val = CCtsp_LP_MAXDOUBLE;
            rval = 0;
        } else {
            *val = (*lp)->lowerbound;
        }
    }

CLEANUP:
    CCutil_stop_timer (&(*lp)->stats.total, silent);
    if (!silent) {
        printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                CClp_nrows ((*lp)->lp), CClp_ncols ((*lp)->lp),
                CClp_nnonzeros ((*lp)->lp));
        fflush (stdout);
    }
    return rval;
}

static int grab_local_x (int ncount, int *tour, int xcount, int *xlist,
        double *x, int *newcount, int **newlist, double **newx, double mult)
{
    int rval = 0, havehash = 0, i, cnt = 0, a, b, c;
    CCutil_edgehash H;

    CC_IFFREE (*newlist, int);
    CC_IFFREE (*newx, double);
    CC_MALLOC (*newx, xcount+ncount, double);
    CC_MALLOC (*newlist, 2*(xcount+ncount), int);

    rval = CCutil_edgehash_init (&H, 4*(xcount + ncount));
    CCcheck_rval (rval, "CCutil_edgehash_init failed");
    havehash = 1;

    for (i = 0; i < xcount; i++) {
        a = xlist[2*i];  b = xlist[2*i+1];
        rval = CCutil_edgehash_set (&H, a, b, cnt);
        CCcheck_rval (rval, "CCutil_edgehash_set failed");
        (*newx)[cnt] = mult*x[i];
        (*newlist)[2*cnt] = a;
        (*newlist)[2*cnt+1] = b;
        cnt++;
    }

    for (i = 0; i < ncount; i++) {
        if (tour) { a = tour[i];  b = tour[(i+1)%ncount]; }
        else      { a = i;        b = (i+1) % ncount;     }
        if (CCutil_edgehash_find (&H, a, b, &c) != -1) {
            (*newx)[c] += 1.0 - mult;
        } else {
            (*newx)[cnt] = 1.0 - mult;
            (*newlist)[2*cnt] = a;
            (*newlist)[2*cnt+1] = b;
            cnt++;
        }
    }
    (*newcount) = cnt;

CLEANUP:
    if (havehash) CCutil_edgehash_free (&H);
    return rval;
}

static int no_tighten (int ncount, int xcount, int *xlist, double *x, int *test,
        double tol)
{
    int rval = 0, k;
    CC_SRKgraph G;

    *test = 0;
    CCcut_SRK_init_graph (&G);

    rval = CCcut_SRK_buildgraph (&G, ncount, xcount, xlist, x);
    CCcheck_rval (rval, "CCcut_SRK_buildgraph failed");
    CCcut_SRK_increment_marker (&G);

    rval = CCcut_SRK_defluff (&G);
    CCcheck_rval (rval, "CCcut_SRK_defluff failed");

    CCcut_SRK_identify_paths_to_edges (&G, &k, 0);
    if (k < (tol * ncount)) { *test = 1; }

CLEANUP:
    CCcut_SRK_free_graph (&G);
    return rval;
}


/*************************************************************************/
/**************    Experimental Code for large World TSP  ****************/
/*************************************************************************/

#define PROCESS_WORLD_CUTS(xmsg) {                                          \
    cut_msg (silent, xmsg, cutcount, z, maxviol);                           \
    if (cutcount) {                                                         \
        rval = call_add_cuts (lp, &cuts, &cut_added, &xcount, &xlist,       \
           &x, &newval, tighten, (int *) NULL, silent, rstate, infeasible); \
        CCcheck_rval (rval, "call_add_cuts failed");                        \
        if (*infeasible) goto CLEANUP;                                      \
    }                                                                       \
}

int CCtsp_world_loop (CCtsp_lp *lp, int silent, int tmax, CCrandstate *rstate,
        int *infeasible)
{
    int rval = 0, outside = 0, inside = 0, tighten = 0;
    int xcount, cutcount, cut_added, edge_added, qqq, *xlist = (int *) NULL;
    double newval, priceval, maxviol, z, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_edgegenerator eginside;
    CCchunk_localcut_timer lc_timer;

    printf ("CCtsp_world_loop ...\n"); fflush (stdout);

    *infeasible = 0;
    CCutil_start_timer (&lp->stats.cutting_loop);
    CCchunk_init_localcut_timer (&lc_timer);
    eginside.ncount = 0;
    if (lp->fulladj) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                                         lp->fulladj, 0, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator (sparse) failed");
    } else if (lp->dat) {
        rval = CCtsp_init_edgegenerator (&eginside, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CC_NO_NEAREST_SUBTOUR, silent, rstate);
        CCcheck_rval (rval, "CCtsp_init_edgegenerator (sparse) failed");
    }

    printf ("We have the intial edges ...\n"); fflush (stdout);

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
              (int *) NULL, (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

#if 0   /* For pricing */
    { 
        int psav = lp->full_edges_valid;

        printf ("Try to price the complete graph\n"); fflush (stdout);

        lp->full_edges_valid = 0;
        rval = full_edge_check (lp, &edge_added, silent, (double *) NULL,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");
        lp->full_edges_valid = psav;

        rval = CCtsp_write_probfile_sav (lp);
        CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
        goto CLEANUP;
    }
#endif  /* End Pricing */

    do {
        do {
            CCutil_start_timer (&lp->stats.cutting_inner_loop);
            cut_added = 0;

            rval = lp_x (lp, &xcount, &xlist, &x);
            CCcheck_rval (rval, "lp_x failed");

            /**** Connect Cuts ****/

            printf ("Connect cuts ...\n"); fflush (stdout);
            rval = CCtsp_connect_cuts (&cuts, &cutcount, lp->graph.ncount,
                                       xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_connect_cuts failed");
            maxviol = 0.0;
            PROCESS_WORLD_CUTS("connect cuts");

            /**** Linear Cuts ****/

            printf ("Linear cuts ...\n"); fflush (stdout);
            rval = CCtsp_segment_cuts (&cuts, &cutcount, lp->graph.ncount,
                                      xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_segment_cuts failed");
            maxviol = 0.0;
            PROCESS_WORLD_CUTS("segment cuts");

#if 0
            /**** Exact Cuts ****/

            printf ("Exact subtours ...\n"); fflush (stdout);
            rval = CCtsp_exact_subtours (&cuts, &cutcount, lp->graph.ncount,
                        xcount, xlist, x, (double *) NULL, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_exact_subtours failed");
            PROCESS_WORLD_CUTS("exact subtours");
#endif

            /****  Fast Blossoms ****/

            printf ("Fast blossoms ...\n"); fflush (stdout);
            rval = CCtsp_fastblossom (&cuts, &cutcount, lp->graph.ncount,
                                  xcount, xlist, x, &z);
            CCcheck_rval (rval, "CCtsp_fastblossom failed");
            maxviol = 0.0;
            PROCESS_WORLD_CUTS("fast blossoms");

            /****  Cut Pool ****/

            printf ("Cut pool ...\n"); fflush (stdout);
            qqq = 2;
            do {
                rval = CCtsp_search_cutpool (lp->pool, &cuts, &cutcount,
                      &maxviol, lp->graph.ncount, xcount, xlist, x, 0,
                      rstate, &z);
                CCcheck_rval (rval, "CCtsp_search_cutpool failed");
                PROCESS_WORLD_CUTS("pool cuts");
            } while (cutcount && --qqq);


            /****  Tighten ****/

            printf ("Tighten cuts ...\n"); fflush (stdout);
            rval = CCtsp_tighten_lp (&lp->cuts,
                    &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                    lp->graph.ncount, xcount, xlist, x, 0.5, 2000, &maxviol,
                    0, rstate, &z);
            CCcheck_rval (rval, "CCtsp_tighten_lp");
            PROCESS_WORLD_CUTS("tighten cuts");

            /****  Tighten Pool  ****/

            printf ("Tighten pool cuts ...\n"); fflush (stdout);
            rval = CCtsp_tighten_lp (lp->pool,
                    &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                    lp->graph.ncount, xcount, xlist, x, 0.1, 2000, &maxviol,
                    0, rstate, &z);
            CCcheck_rval (rval, "CCtsp_tighten_lp");
            PROCESS_WORLD_CUTS("tighten pool cuts");


            /****  Double Deckers ****/

            printf ("Double deckers ...\n"); fflush (stdout);
            rval = CCtsp_double_decker_lp (&lp->cuts,
                    &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                    lp->graph.ncount, xcount, xlist, x, 10.0, 2000,
                    &maxviol, rstate, &z);
            CCcheck_rval (rval, "CCtsp_double_decker_lp failed");
            PROCESS_WORLD_CUTS("double deckers");

            /****  Teething ****/

            printf ("Teething ...\n"); fflush (stdout);
            rval = CCtsp_teething_lp (&lp->cuts,
                    &lp->stats.extra_tighten_stats, &cuts, &cutcount,
                    lp->graph.ncount, xcount, xlist, x, 10.0, 2000,
                    &maxviol, rstate, &z);
            CCcheck_rval (rval, "CCtsp_teething_lp failed");
            PROCESS_WORLD_CUTS("teethed combs");

            /****  Domino Pool ****/

            if (lp->dominopool) {
                printf ("Domino pool ...\n"); fflush (stdout);
                int domstarter = 0;
                qqq = 5;
                do {
                    rval = CCtsp_search_dominopool (lp->dominopool, &cuts,
                         &cutcount, &maxviol, lp->graph.ncount, xcount,
                         xlist, x, &domstarter, 10000, 500, rstate, &z);
                    CCcheck_rval (rval, "CCtsp_search_dominopool failed");
                    PROCESS_WORLD_CUTS("domino pool cuts");
                } while (cutcount && --qqq);
            }

            /****  Local Cuts ****/

            printf ("Local cuts (%d) ...\n", tmax); fflush (stdout);
            {
                int  maxchunksize = tmax;
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 0;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                flags.maxchunksize = maxchunksize;
                flags.spheresize   = maxchunksize - 2;

                rval = CCchunk_localcuts (&cuts, &cutcount,
                     lp->graph.ncount, xcount, xlist, x, 0.0, flags,
                     &lc_timer, silent, rstate, &z, &maxviol);
                CCcheck_rval (rval, "CCchunk_localcuts failed");
                PROCESS_WORLD_CUTS("local cuts");
            }

            {
                int  maxchunksize = tmax;
                CCchunk_flag flags;

                flags.dummy = 0;
                flags.permute = 0;
                flags.weighted = 0;
                flags.spheres = 1;
                flags.uncivilized = 0;
                flags.noshrink = 0;
                flags.nolift = 0;

                flags.maxchunksize = maxchunksize;
                flags.spheresize = maxchunksize - 2;

                rval = CCchunk_localcuts (&cuts, &cutcount,
                         lp->graph.ncount, xcount, xlist, x, 0.0, flags,
                         &lc_timer, silent, rstate, &z, &maxviol);
                CCcheck_rval (rval, "CCchunk_localcuts failed");
                PROCESS_WORLD_CUTS("local cuts");
            }

            printf ("End loop stuff ...\n"); fflush (stdout);

            /****  End Loop Stuff ****/

            CC_FREE (xlist, int);
            CC_FREE (x, double);

            if (!cut_added || (inside % 5) == 0) {
                CCutil_start_timer (&lp->stats.sparse_edge_check);
                rval = sparse_edge_check (lp, &eginside, &edge_added,
                                          (double *) NULL, silent, rstate);
                CCcheck_rval (rval, "sparse_edge_check failed");
                CCutil_stop_timer (&lp->stats.sparse_edge_check, silent);
            }
            CCutil_stop_timer (&lp->stats.cutting_inner_loop, silent);

            rval = CCtsp_write_probfile_sav (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");

            if (lp->pool) {
                char buf[1024];
                if (!silent) {
                    printf ("Write Pool: %d cuts\n", lp->pool->cutcount);
                    fflush (stdout);
                }
                sprintf (buf, "%s.pul", lp->problabel);
                rval = CCtsp_write_cutpool (lp->graph.ncount, buf, lp->pool);
                CCcheck_rval (rval, "CCtsp_write_cutpool failed");
            }

            inside++;
            if (silent) {
                rval = lp_value (lp, &priceval);
                if (rval) {rval = 1; goto CLEANUP;}
                printf ("  LP Value %2d: %f\n", inside, priceval);
                fflush (stdout);
            }
        } while (edge_added || cut_added);

/*
        rval = full_edge_check (lp, &edge_added, silent, (double *) NULL,
                                rstate);
        CCcheck_rval (rval, "full_edge_check failed");

        rval = CCtsp_write_probfile_sav (lp);
        CCcheck_rval (rval, "CCtsp_write_profile_sav failed");
*/
        outside++;
    } while (edge_added);

CLEANUP:
    if (*infeasible) {
        printf ("LP is infeasible in world_loop\n"); fflush (stdout);
    }
    z = CCutil_stop_timer (&lp->stats.cutting_loop, 0);
    printf ("Number of outside rounds: %d (%d inside)\n", outside, inside);
    fflush (stdout);

    if (eginside.ncount)
        CCtsp_free_edgegenerator (&eginside);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}
