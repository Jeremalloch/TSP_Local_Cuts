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
/*                    Interface to the Cutters                              */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 17, 1997                                                 */
/*                                                                          */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int CCtsp_connect_cuts (CCtsp_lpcut_in **cuts, int *cutcount,           */
/*      int ncount, int ecount, int *elist, double *x, double *fullzeit)    */
/*    FINDS violated subtour inequalities via connectivity.                 */
/*     -cuts will return any new cuts found (they will be added to the      */
/*      head of the linked list)                                            */
/*     -cutcount will return the number of new cuts added                   */
/*     -ncount is the number of nodes                                       */
/*     -ecount is the number of edges                                       */
/*     -elist contains the LP edges in node node format                     */
/*     -x is an LP solution                                                 */
/*     -fullzeit if not NULL returns the total run time in seconds          */
/*                                                                          */
/*  int CCtsp_segment_cuts (CCtsp_lpcut_in **cuts, int *cutcount,           */
/*      int ncount, int ecount, int *elist, double *x, double *fullzeit)    */
/*    FINDS violated subtour inequalities via linsub.                       */
/*                                                                          */
/*  int CCtsp_exact_subtours (CCtsp_lpcut_in **cuts, int *cutcount,         */
/*      int ncount, int ecount, int *elist, double *x, double *eps,         */
/*      double *fullzeit, double *viol)                                     */
/*    FINDS violated subtour inequalities via a mincut algorithm.           */
/*     -fullzeit returns run time in seconds                                */
/*     -viol is set to 0.0 (place-holder for max violation)                 */
/*                                                                          */
/*  int CCtsp_tighten_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,    */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, int extreme, CCrandstate *rstate,        */
/*      double *fullzeit)                                                   */
/*    CALLS tighten for each cut in the cuts.                               */
/*     -stats contains some running statistics of tighten                   */
/*     -cutsout returns the tightened cuts that are violated (they are      */
/*      added to the tail of the linked list)                               */
/*     -cutcount is the number of cuts in cutsout                           */
/*     -testtol is a tolerance for calling tighten (call only when the      */
/*      cut has slack value within testtol)                                 */
/*     -maxcuts is a bound on the number of cuts to be returned             */
/*     -extreme in non-zero then a flow-based tighten in called             */
/*     -fullzeit returns run time in seconds                                */
/*                                                                          */
/*  int CCtsp_double_decker_lp (CCtsp_lpcuts *cuts,                         */
/*      CCtsp_tighten_info *stats, CCtsp_lpcut_in **cutsout,                */
/*      int *cutcount, int ncount, int ecount, int *elist, double *x,       */
/*      double testtol, int maxcuts, double *viol, CCrandstate *rstate,     */
/*      double *fullzeit)                                                   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_cliquetree_lp (CCtsp_lpcuts *cuts,                            */
/*      CCtsp_tighten_info *stats, CCtsp_lpcut_in **cutsout,                */
/*      int *cutcount, int ncount, int ecount, int *elist, double *x,       */
/*      double testtol, int maxcuts, double *viol, CCrandstate *rstate,     */
/*      double *fullzeit)                                                   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_star_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,       */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate, double *fullzeit)   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_handling_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,   */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate, double *fullzeit)   */
/*    CALLS CCtsp_comb_handling for each comb in cuts.                      */
/*     -agruments as in CCtsp_tighten_lp.                                   */
/*                                                                          */
/*  int CCtsp_teething_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,   */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate, double *fullzeit)   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_file_cuts (char *cutfile, CCtsp_lpcut_in **cuts,              */
/*      int *cutcount, int ncount, int *tour, double *fullzeit)             */
/*    READS a set of cuts from a file; the format of the cuts can be        */
/*     found by examining the code                                          */
/*     -cutfile is an asci file with a list of cuts                         */
/*     -cuts will return any new cuts found (they will be added to the      */
/*      tail of the linked list)                                            */
/*     -cutcount with return the number of new cuts added                   */
/*     -ncount is the number of nodes                                       */
/*     -tour the permutation tour (used to map the incoming nodes)          */
/*     -fullzeit if not NULL returns the total run time in seconds          */
/*                                                                          */
/*  int CCtsp_file_cuts_write (const char *cutfile, CCtsp_lpcuts *cuts,     */
/*      int *tour)                                                          */
/*    WRITES a set of cuts in a text file that can be read by               */
/*     tsp_file_cuts                                                        */
/*     -cutfile is the name of the file to be written                       */
/*     -cuts is the set of cuts to be written                               */
/*     -tour is a permutation tour (used to map the outgoing nodes)         */
/*                                                                          */
/*  int CCtsp_test_pure_comb (int ncount, CCtsp_lpcut_in *c, int *yes_no,   */
/*      int *handle)                                                        */
/*    TEST if the cut is a comb (without flipped teeth or intersections)    */
/*     -ncount is the number of nodes in the TSP                            */
/*     -yes_no will be set to either 0 or 1, with 1 meaning yes             */
/*     -handle with return the index of the handle if the cut is a comb     */
/*      (handle can be NULL)                                                */
/*                                                                          */
/*  int CCtsp_test_pseudocomb (int ncount, CCtsp_lpcut_in *c, int handle,   */
/*      int *yes_no)                                                        */
/*    TEST if the cut is a pseudocomb.                                      */
/*     -handle gives the index of the handle of the pseudocomb              */
/*                                                                          */
/*  int CCtsp_test_teeth_disjoint (int ncount, CCtsp_lpcut_in *c,           */
/*      int handle, int *yes_no)                                            */
/*    TEST if the cliques other than handle are pairwise disjoint.          */
/*     -yes_no is 1 if disjoint and 0 otherwise.                            */
/*                                                                          */
/*  int CCtsp_find_pure_handle (int ncount, CCtsp_lpcut_in *c,              */
/*      int *handle)                                                        */
/*    FINDS a clique that is c's handle if c is a comb; the search          */
/*     assumes that the teeth are disjoint, so if the comb has              */
/*     extra intersections then a tooth may be returned.                    */
/*     -handle returns the potential handle (it will return -1 if no        */
/*      clique is a potential handle)                                       */
/*                                                                          */
/*  int CCtsp_buildcut_begin (CCtsp_cutinfo *cuts, int init_cliquecount)    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_buildcut_addclique (CCtsp_cutinfo *cuts, *arr, int size)      */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_buildcut_finish (CCtsp_cutinfo *cuts, int rhs)                */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_buildcut_abort (CCtsp_cutinfo *cuts)                         */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_truncate_cutlist (CCtsp_lpcut_in **cuts, int ncount,          */
/*      int ecount, int *elist, double *x, int maxcuts,                     */
/*      CCrandstate *rstate)                                                */
/*    RETURNS the maxcuts most violated cuts in the linked list, the        */
/*     remaining cuts are freed.                                            */
/*                                                                          */
/*  int CCtsp_create_cut (int ncount, CCtsp_lpcut_in **cuts, int nsets,     */
/*      int *setsize,  int **sets, int *setmultipliers, int nsemicuts,      */
/*      int **semisize, int ***semicuts, int *semimultipliers, int rhs)     */
/*    CREATES an CCtsp_lpcut_in struct containing the inequality described  */
/*      by the parameters.  The inequality is >= form.                      */
/*     -ncount: number of nodes in TSP                                      */
/*     -cuts: passed in as a NULL-terminated linked list of cuts.  The      */
/*      new cut will be added to the tail of the list. To obtain a single   */
/*      cut, set a local variable CCtsp_lpcut_in *c = NULL and pass in &c.  */
/*     -nsets: the number of sets S_i for \delta(S_i) edge sets             */
/*     -setsize: array of length nsets giving size of each set S_i          */
/*     -sets: array of length nsets with sets[i] an array of the setsize[i] */
/*      members of set S_i                                                  */
/*     -setmultipliers: should be NULL if all multipliers are 1, otherwise  */
/*      an array of length nsets where setmultipliers[i] is the multiplier  */
/*      for \delta(S_i).  The multipliers can be negative.                  */
/*     -nsemicuts: the number of semicuts E[A_i,B_i]                        */
/*     -semisize: semisize[i][0] is size of set A_i and semisize[i][1]      */
/*      is size of B_i                                                      */
/*     -semicuts: semicuts[i][0] is an array of the members of A_i and      */
/*      semicuts[i][1] is an array of the members of B_i                    */
/*     -semimultipliers: should be NULL if all multipliers are 1, otherwise */
/*      an array of length nsemicuts giving multipliers for each E[A_i,B_i] */
/*     -rhs: the right-hand-side value for the >= inequality                */
/*    NOTE: The most efficient way to model a cut is with only \delta(S_i)  */
/*      sets and with setmultipliers = NULL.                                */
/*                                                                          */
/*  int CCtsp_test_dp_cut (int ncount, CCtsp_lpcut_in *c, int *yes_no,      */
/*      int *large_domino)                                                  */
/*     -large_domino: if not NULL, then set to 1 if cut has a domino        */
/*      containing 3/4*ncount or more nodes                                 */
/*    TEST if cut is a domino-parity cut (yes_no is set to 1 if true and    */
/*     set to 0 if false)                                                   */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "tsp.h"
#include "cut.h"
#include "combs.h"
#include "verify.h"
#include "localcut.h"

#ifdef CCtsp_USE_DOMINO_CUTS
#include "KP.h"
#endif

#define X_FLUFF (1e-10)
#undef  DUMP_BUILDCUT

typedef struct exactsub_param {
    int             nodecount;
    int             cutcount;
    CCtsp_lpcut_in *cuts;
} exactsub_param;

static int
    add_segment (double val, int a, int b, void *pass_param),
    add_exact (double val, int count, int *cutarray, void *pass_param),
    work_on_combs_in_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts, int caller,
        double *viol, CCrandstate *rstate, double *fullzeit),
    grab_nonzero_x (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double **new_x, double tol);

static int build_triomino (int ncount, int nTset, int *Tset, int nH0set,
        int *H0set, int nH1set, int *H1set,
        CCtsp_lpclique *ptset, CCtsp_lpdomino *psemi0, CCtsp_lpdomino *psemi1);

int CCtsp_connect_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, double *fullzeit)
{
    int rval = 0, i, k, ncomp;
    int *comps = (int *) NULL, *compscount = (int *) NULL;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;
    int rzeit = CCutil_zeit ();

    *cutcount = 0;
    if (fullzeit) *fullzeit = 0.0;
    rval = CCcut_connect_components (ncount, ecount, elist, x, &ncomp,
                                     &compscount, &comps);
    CCcheck_rval (rval, "CCcut_connect_components failed");

    for (i = 0, k = 0; i < ncomp - 1; k += compscount[i], i++) {
        rval = CCtsp_array_to_subtour (&c, comps + k, compscount[i], ncount);
        CCcheck_rval (rval, "CCtsp_array_to_subtour failed");
        c->next = *cuts;
        *cuts = c;
        (*cutcount)++;
    }

    if (fullzeit) *fullzeit = CCutil_zeit () - rzeit;

CLEANUP:
    CC_IFFREE (comps, int);
    CC_IFFREE (compscount, int);
    return rval;
}

int CCtsp_segment_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, double *fullzeit)
{
    int rval = 0, i;
    int *endmark = (int *) NULL;
    exactsub_param p;
    double szeit = CCutil_zeit ();

    *cutcount = 0;
    if (fullzeit) *fullzeit = 0.0;

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    CC_MALLOC (endmark, ncount, int);
    for (i=0; i<ncount; i++) { endmark[i] = 0; }

    for (i=0; i<ecount; i++) {
        if (x[i] >= 0.999999) {
            endmark[elist[2*i]]++;
            endmark[elist[2*i+1]]++;
        }
    }
    for (i=0; i<ncount; i++) {
        if (endmark[i] == 2) {
            endmark[i] = CC_LINSUB_NO_END;
        } else {
            endmark[i] = CC_LINSUB_BOTH_END;
        }
    }

    rval = CCcut_linsub (ncount, ecount, endmark, elist, x, 2.0 - 0.0001,
                         (void *) &p, add_segment);
    CCcheck_rval (rval, "CCcut_linsub failed");

    *cutcount = p.cutcount;
    *cuts = p.cuts;

    rval = 0;
    if (fullzeit) *fullzeit = CCutil_zeit () - szeit;

CLEANUP:
    CC_IFFREE (endmark, int);
    return rval;
}

int CCtsp_shrink_subtours (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, double *fullzeit, double *maxviol)
{
    int rval = 0;
    double szeit = CCutil_zeit();
    exactsub_param p;

    if (maxviol) *maxviol = 0.0;
    *cutcount = 0;
/*
    rval = CCtsp_connect_cuts (cuts, cutcount, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCtsp_connect_cuts failed\n"); goto CLEANUP;
    }

    if (*cutcount > 0) {
        rval = 0; goto CLEANUP;
    }
*/

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    rval = CCcut_shrink_cuts (ncount, ecount, elist, x, 2.0 - 0.0001,
                              add_exact, (void *) &p);
    CCcheck_rval (rval, "CCcut_shrink_cuts failed");

    *cutcount = p.cutcount;
    *cuts = p.cuts;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    return rval;
}

int CCtsp_exact_subtours (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
            int ecount, int *elist, double *x, double *eps, double *fullzeit,
            double *viol)
{
    int rval = 0;
    double tol;
    exactsub_param p;
    double szeit = CCutil_zeit ();

    if (fullzeit) *fullzeit = 0.0;
    if (viol) *viol = 0.0;

    *cutcount = 0;
    rval = CCtsp_connect_cuts (cuts, cutcount, ncount, ecount, elist, x,
                               (double *) NULL);
    CCcheck_rval (rval, "CCtsp_connect_cuts failed");

    if (*cutcount > 0) {
        rval = 0; goto CLEANUP;
    }

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    if (eps) tol = *eps;
    else     tol = 0.0001;

#if 0
    if (eps) {
        printf ("Tolerance for Subtours: %lf\n", tol);
        fflush (stdout);
    }
#endif

    rval = CCcut_violated_cuts (ncount, ecount, elist, x, 2.0 - tol,
                                add_exact, (void *) &p);
    CCcheck_rval (rval, "CCcut_violated_cuts failed");

    *cutcount = p.cutcount;
    *cuts = p.cuts;

#if 0
  - this is just to check the values of the exact cuts
    if (*cutcount) {
        CCtsp_lpgraph lg;
        CCtsp_lpcut_in *c;
        double t;

        CCtsp_init_lpgraph_struct (&lg);

        rval = CCtsp_build_lpgraph (&lg, ncount, ecount, elist, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
        }
        rval = CCtsp_build_lpadj (&lg, 0, ecount);
        if (rval) {
            CCtsp_free_lpgraph (&lg);
            fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
        }
        for (c = p.cuts; c; c = c->next) {
            t = CCtsp_cutprice (&lg, c, x);
            printf ("[%f] ", 2.0 + t); fflush (stdout);
        }
        printf ("\n"); fflush (stdout);
        CCtsp_free_lpgraph (&lg);
    }
#endif
    rval = 0;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit () - szeit;
    return rval;
}

static int add_segment (double val, int a, int b, void *pass_param)
{
    int rval = 0;
    exactsub_param *p = (exactsub_param *) pass_param;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    if (val > 2.0) {
        printf ("Warning: Cut of value %f in add_segment\n", val);
        fflush (stdout);
        goto CLEANUP;
    }

    rval = CCtsp_segment_to_subtour (&c, a, b, p->nodecount);
    CCcheck_rval (rval, "CCtsp_segment_to_subtour failed");

    c->next = p->cuts;
    p->cuts = c;
    p->cutcount++;

CLEANUP:
    return rval;
}

static int add_exact (double val, int count, int *cutarray, void *pass_param)
{
    int rval = 0;
    exactsub_param *p = (exactsub_param *) pass_param;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    if (count >= p->nodecount) goto CLEANUP;

    if (val > 2.0) {
        printf ("Warning: Cut of value %f in add_exact\n", val);
        fflush (stdout);
        goto CLEANUP;
    }

    rval = CCtsp_array_to_subtour (&c, cutarray, count, p->nodecount);
    CCcheck_rval (rval, "CCtsp_array_to_subtour failed");

    c->next = p->cuts;
    p->cuts = c;
    p->cutcount++;

CLEANUP:

    return rval;
}

int CCtsp_tighten_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, int extreme, CCrandstate *rstate, double *fullzeit)
{
    int rval = 0, i, ea = 0, clistsize = 0, vlistsize = 0, count = 0;
    int newecount, *newelist = (int *) NULL, *perm = (int *) NULL;
    double improve, maxviol = 0.0;
    double *newx = (double *) NULL;
    double *vlist = (double *) NULL, *cutval = (double *) NULL;
    CCtsp_lpcut_in new, old, *c, **clist = (CCtsp_lpcut_in **) NULL;
    CCtsp_lpgraph lg;
    double szeit = CCutil_zeit();

    *cutcount = 0;
    if (!cuts || !cuts->cutcount) return 0;

    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    CCcheck_rval (rval, "grab_nonzero_x failed");

    cutval = CC_SAFE_MALLOC (cuts->cutcount, double);
    CCcheck_NULL (cutval, "out of memory for cutval");
    rval = CCtsp_price_cuts (cuts, ncount, newecount, newelist, newx, cutval);
    CCcheck_rval (rval, "CCtsp_price_cuts failed");

    CCtsp_init_lpgraph_struct (&lg);

    rval = CCtsp_build_lpgraph (&lg, ncount, newecount, newelist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
    CC_FREE (newelist, int);
    rval = CCtsp_build_lpadj (&lg, 0, newecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    for (i = 0; i < cuts->cutcount; i++) {
        if (cutval[i] < testtol && !cuts->cuts[i].branch
            && (cuts->cuts[i].cliquecount > 1 ||
               (cutval[i] < 0.1*testtol && extreme == 0))
            && CCtsp_hypergraph_type (&cuts->cuts[i])) {
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &old);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            if (extreme == 0) {
                rval = CCtsp_tighten_lpcut_in (&lg, &old, newx, &new, stats,
                                               &improve);
                CCcheck_rval (rval, "CCtsp_tighten_lpcut failed");
            } else {
                ea++;
                rval = CCtsp_xtighten_lpcut_in (&lg, &old, newx, &new,
                                               &improve);
                CCcheck_rval (rval, "CCtsp_tighten_lpcut failed");
                if (ea % 100 == 99) {
                    printf ("Checked %d (%d) of %d cuts\n", i, ea,
                             cuts->cutcount);
                    fflush (stdout);
                }
            }
            CCtsp_free_lpcut_in (&old);

            if (improve - cutval[i] > CCtsp_MIN_VIOL) {
                c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
                CCcheck_NULL (c, "out of memory in CCtsp_tighten_lp");
                *c = new;

                if (count >= clistsize) {
                    void *tmp_ptr = (void *) clist;
                    rval = CCutil_reallocrus_scale (&tmp_ptr, &clistsize,
                                count + 1, 1.3, sizeof (CCtsp_lpcut_in *));
                    CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
                    clist = (CCtsp_lpcut_in **) tmp_ptr;
                }
                if (count >= vlistsize) {
                    void *tmp_ptr = (void *) vlist;
                    rval = CCutil_reallocrus_scale (&tmp_ptr, &vlistsize,
                                count + 1, 1.3, sizeof (double));
                    CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
                    vlist = (double *) tmp_ptr;
                }
                clist[count] = c;
                vlist[count] = cutval[i] - improve;
                count++;
            } else {
                CCtsp_free_lpcut_in (&new);
            }
        }
    }

    if (count) {
        perm = CC_SAFE_MALLOC (count, int);
        CCcheck_NULL (perm, "out of memory for perm");
        for (i = 0; i < count; i++) {
            perm[i] = i;
        }
        if (count > maxcuts) {
            CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
            for (i = maxcuts; i < count; i++) {
                CCtsp_free_lpcut_in (clist[perm[i]]);
            }
            count = maxcuts;
        }
        for (i = 0; i < count; i++) {
            if (vlist[perm[i]] < maxviol)
                maxviol = vlist[perm[i]];
            clist[perm[i]]->next = *cutsout;
            *cutsout = clist[perm[i]];
        }
    }

    *cutcount = count;
    if (viol) *viol = -maxviol;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (cutval, double);
    CCtsp_free_lpgraph (&lg);
    return rval;
}

#define CALL_TEETHING   1
#define CALL_DDECKER    2
#define CALL_CLIQUETREE 3
#define CALL_STAR       4
#define CALL_HANDLING   5
#define CALL_DPTIGHT    6

int CCtsp_double_decker_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate, double *fullzeit)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
       elist, x, testtol, maxcuts, CALL_DDECKER, viol, rstate, fullzeit);
    CCcheck_rval (rval, "work_on_combs_in_lp failed");

CLEANUP:

    return rval;
}

int CCtsp_teething_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts, double *viol,
        CCrandstate *rstate, double *fullzeit)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
       elist, x, testtol, maxcuts, CALL_TEETHING, viol, rstate, fullzeit);
    CCcheck_rval (rval, "work_on_combs_in_lp failed");

CLEANUP:
    return rval;
}

int CCtsp_cliquetree_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate, double *fullzeit)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
       elist, x, testtol, maxcuts, CALL_CLIQUETREE, viol, rstate, fullzeit);
    CCcheck_rval (rval, "work_on_combs_in_lp failed");

CLEANUP:
    return rval;
}

int CCtsp_star_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate, double *fullzeit)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
       elist, x, testtol, maxcuts, CALL_STAR, viol, rstate, fullzeit);
    CCcheck_rval (rval, "work_on_combs_in_lp failed");

CLEANUP:
    return rval;
}

int CCtsp_handling_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate, double *fullzeit)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
       elist, x, testtol, maxcuts, CALL_HANDLING, viol, rstate, fullzeit);
    CCcheck_rval (rval, "work_on_combs_in_lp failed");

CLEANUP:
    return rval;
}

static int work_on_combs_in_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts, int caller,
        double *viol, CCrandstate *rstate, double *fullzeit)
{
    int rval = 0, i, test;
    CCtsp_lpcut_in new, old;
    CCtsp_lpcut_in *c  = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *dd = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *ddnext;
    double improve, newslack, dslack;
    CCtsp_lpgraph lg;
    CC_GCgraph gg;
    double *newx = (double *) NULL;
    int *newelist = (int *) NULL;
    int newecount;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    double *vlist = (double *) NULL;
    double maxviol = 0.0;
    int clistsize = 0, vlistsize = 0, count = 0;
    int *perm = (int *) NULL;
    double *cutval = (double *) NULL;
    double szeit = CCutil_zeit();

    *cutcount = 0;
    if (!cuts || !cuts->cutcount) return 0;

    CCtsp_init_lpgraph_struct (&lg);
    CCcombs_GC_init_graph (&gg);

    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    CCcheck_rval (rval, "grab_nonzero_x failed");

    cutval = CC_SAFE_MALLOC (cuts->cutcount, double);
    CCcheck_NULL (cutval, "out of memory for cutval");

    rval = CCtsp_price_cuts (cuts, ncount, newecount, newelist, newx, cutval);
    CCcheck_rval (rval, "CCtsp_price_cuts failed");

    rval = CCtsp_build_lpgraph (&lg, ncount, newecount, newelist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");

    if (caller == CALL_DDECKER || caller == CALL_CLIQUETREE ||
        caller == CALL_STAR    || caller == CALL_HANDLING) {
        rval = CCcombs_GC_build_graph (&gg, ncount, newecount, newelist, newx);
        CCcheck_rval (rval, "CCcombs_GC_build_graph failed");
    }

    CC_FREE (newelist, int);
    rval = CCtsp_build_lpadj (&lg, 0, newecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    for (i = 0; i < cuts->cutcount; i++) {
        if (cuts->cuts[i].branch || cuts->cuts[i].cliquecount % 2 ||
            cuts->cuts[i].cliquecount < 4 || cutval[i] >= testtol ||
            CCtsp_hypergraph_type (&cuts->cuts[i]) == 0) {
            continue;
        }
        rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &old);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
        rval = CCtsp_test_pure_comb (ncount, &old, &test, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCtsp_test_pure_comb failed\n");
            CCtsp_free_lpcut_in (&old);
            goto CLEANUP;
        }
        if (test == 1) {
            switch (caller) {
            case CALL_TEETHING:
                rval = CCtsp_teething (&lg, newx, &old, &dd);
                CCcheck_rval (rval, "CCtsp_teething failed");
                break;
            case CALL_DDECKER:
                rval = CCtsp_comb_to_double_decker (&lg, &gg, newx, &old, &dd);
                CCcheck_rval (rval, "CCtsp_comb_to_double_decker failed");
                break;
            case CALL_CLIQUETREE:
                rval = CCtsp_comb_to_cliquetree (&lg, &gg, newx, &old, &dd);
                CCcheck_rval (rval, "CCtsp_comb_to_cliquetree failed");
                break;
            case CALL_STAR:
                rval = CCtsp_comb_to_star (&lg, &gg, newx, &old, &dd);
                CCcheck_rval (rval, "CCtsp_comb_to_star failed");
                break;
            case CALL_HANDLING:
                rval = CCtsp_comb_handling (&lg, &gg, newx, &old, &dd);
                CCcheck_rval (rval, "CCtsp_comb_handling failed");
                break;
            default:
                fprintf (stderr, "unknown caller in work_on_combs_in_lp\n");
                rval = 1; goto CLEANUP;
            }

            CCtsp_free_lpcut_in (&old);
            while (dd) {
                ddnext = dd->next;
                dslack = CCtsp_cutprice (&lg, dd, newx);
                if (dslack >= 1.0 /* 1.0 */) {
                    CCtsp_free_lpcut_in (dd);
                    CC_FREE (dd, CCtsp_lpcut_in);
                } else {
                    rval = CCtsp_tighten_lpcut_in (&lg, dd, newx, &new,
                                                   stats, &improve);
                    CCcheck_rval (rval, "CCtsp_tighten_lpcut failed");
                    CCtsp_free_lpcut_in (dd);
                    CC_FREE (dd, CCtsp_lpcut_in);

                    newslack = dslack - 2.0*improve;
                    if (-newslack > CCtsp_MIN_VIOL) {
                        c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
                        CCcheck_NULL (c, "out of memory for c");
                        *c = new;
                        if (count >= clistsize) {
                            void *tmp_ptr = (void *) clist;
                            rval = CCutil_reallocrus_scale (&tmp_ptr,
                                    &clistsize, count + 1, 1.3,
                                    sizeof (CCtsp_lpcut_in *));
                            CCcheck_rval(rval,"CCutil_reallocrus_scale failed");
                            clist = (CCtsp_lpcut_in **) tmp_ptr;
                        }
                        if (count >= vlistsize) {
                            void *tmp_ptr = (void *) vlist;
                            rval = CCutil_reallocrus_scale (&tmp_ptr,
                                     &vlistsize, count + 1, 1.3,
                                     sizeof (double));
                            CCcheck_rval(rval,"CCutil_reallocrus_scale failed");
                            vlist = (double *) tmp_ptr;
                        }
                        clist[count] = c;
                        vlist[count] = newslack;
                        count++;
                    } else {
                        CCtsp_free_lpcut_in (&new);
                    }
                }
                dd = ddnext;
            }
        } else {
            CCtsp_free_lpcut_in (&old);
        }
    }

    if (count) {
        perm = CC_SAFE_MALLOC (count, int);
        CCcheck_NULL (perm, "out of memory for perm");
        for (i = 0; i < count; i++) { perm[i] = i; }
        if (count > maxcuts) {
            CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
            for (i = maxcuts; i < count; i++) {
                CCtsp_free_lpcut_in (clist[perm[i]]);
            }
            count = maxcuts;
        }
        for (i = 0; i < count; i++) {
            if (vlist[perm[i]] < maxviol)
                maxviol = vlist[perm[i]];
            clist[perm[i]]->next = *cutsout;
            *cutsout = clist[perm[i]];
        }
    }

    *cutcount = count;
    if (viol) *viol = -maxviol;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (cutval, double);
    CCtsp_free_lpgraph (&lg);
    CCcombs_GC_free_graph (&gg);
    if (dd) {
        CCtsp_free_lpcut_in (dd);
    }
    return rval;
}

int CCtsp_file_cuts (char *cutfile, CCtsp_lpcut_in **cuts, int *cutcount,
        int ncount, int *tour, double *fullzeit)
{
    FILE *in = (FILE *) NULL;
    int *inv = (int *) NULL;
    CCtsp_lpcut_in *c;
    CCtsp_lpcut_in **clast;
    int i, j, k;
    int ncliques, size;
    int *icliq = (int *) NULL;
    int rval = 0;
    int rzeit = CCutil_zeit ();

    *cutcount = 0;
    if (fullzeit) *fullzeit = 0.0;

    in = fopen (cutfile, "r");
    if  (in == (FILE *) NULL) {
        fprintf (stderr, "unable to open %s for reading\n", cutfile);
        return 0;
    }

    CC_MALLOC (inv, ncount, int);
    for (i = 0; i < ncount; i++) { inv[tour[i]] = i; }

    clast = cuts;
    while ((*clast) != (CCtsp_lpcut_in *) NULL) {
        clast = &((*clast)->next);
    }
    
    while (fscanf (in, "%d", &ncliques) != EOF) {
        CC_MALLOC (c, 1, CCtsp_lpcut_in);
        CCtsp_init_lpcut_in (c);

        c->cliquecount = ncliques;
        CC_MALLOC (c->cliques, ncliques, CCtsp_lpclique);
        for (i = 0; i < ncliques; i++) {
            fscanf (in, "%d", &size);
            CC_MALLOC (icliq, size, int);
            for (j = 0; j < size; j++) {
                fscanf (in, "%d", &k);
                icliq[j] = inv[k];
            }
            rval = CCtsp_array_to_lpclique (icliq, size, &(c->cliques[i]));
            CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
            CC_FREE (icliq, int);
        }
        fscanf (in, "%d", &(c->rhs));
        c->sense = 'G';
        c->branch = 0;
        rval = CCtsp_construct_skeleton (c, ncount);
        CCcheck_rval (rval, "CCtsp_construct_skeleton failed");

        (*clast) = c;
        c->next = (CCtsp_lpcut_in *) NULL;
        clast = &(c->next);
        (*cutcount)++;
        rval = CCverify_cut (c, ncount, CC_TYPE_ALL, &i, 0, (int *) NULL,
                             (char *) NULL);
        if (rval) {
            fprintf (stderr, "Invalid file cut\n");
            CCtsp_print_lpcut_in (c);
            goto CLEANUP;
        } 
        /* printf ("%d: File cut type %d\n", *cutcount, i); */
    }

    if (fullzeit) *fullzeit = CCutil_zeit () - rzeit;

CLEANUP:
    CC_IFFREE (inv, int);
    fclose (in);
    return  rval;
}

int CCtsp_file_cuts_write (const char *cutfile, CCtsp_lpcuts *cuts, int *tour)
{
    FILE *out = (FILE *) NULL;
    int rval = 0, i, j, k, p;
    int cutcount = cuts->cutcount;
    CCtsp_lpcut *c;
    CCtsp_lpclique *cl;
    int isize;

    out = fopen (cutfile, "w");
    if  (out == (FILE *) NULL) {
        fprintf (stderr, "unable to open %s for writing\n", cutfile);
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < cutcount; i++) {
        c = &cuts->cuts[i];
        if (!c->branch && CCtsp_hypergraph_type (c)) {
            fprintf (out, "%d\n", c->cliquecount);
            for (j = 0; j < c->cliquecount; j++) {
                cl = &cuts->cliques[c->cliques[j]];
                for (k = 0, isize = 0; k < cl->segcount; k++) {
                    isize += (cl->nodes[k].hi - cl->nodes[k].lo + 1);
                }
                fprintf (out, "%d  ", isize);
                CC_FOREACH_NODE_IN_CLIQUE (p, *cl, k) {
                    fprintf (out, "%d ", tour[p]);
                }
                fprintf (out, "\n");
            }
            fprintf (out, "%d\n", c->rhs);
        }
    }

CLEANUP:
    if (out) fclose (out);
    return rval;
}

int CCtsp_buildcut_begin (CCtsp_cutinfo *cuts, int init_cliquecount)
{
    cuts->current = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    if (!cuts->current) return -1;
    CCtsp_init_lpcut_in (cuts->current);
    cuts->current->cliques = CC_SAFE_MALLOC (init_cliquecount, CCtsp_lpclique);
    if (!cuts->current->cliques) {
        CC_FREE (cuts->current, CCtsp_lpcut_in);
        return -1;
    }
    return 0;
}

int CCtsp_buildcut_addclique (CCtsp_cutinfo *cuts, int *arr, int size)
{
    int i;
    int *newarr = (int *) NULL;
    int newsize;
    int rval;
    CCtsp_lpcut_in *c = cuts->current;
    void *tmp_ptr;

    if (!c) {
        fprintf (stderr, "Trying to add to nonexistent clique\n");
        return -1;
    }

    rval = CCcut_SRK_expand (&cuts->expand, arr, size, &newarr, &newsize);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_expand failed\n");
        CCtsp_buildcut_abort (cuts);
        return rval;
    }

    tmp_ptr = (void *) c->cliques;
    rval = CCutil_reallocrus_count (&tmp_ptr, c->cliquecount+1,
                             sizeof (c->cliques[0]));
    if (rval) {
        fprintf (stderr, "couldn't realloc cliques\n");
        CC_IFFREE (newarr, int);
        CCtsp_buildcut_abort (cuts);
        return rval;
    }
    c->cliques = (CCtsp_lpclique *) tmp_ptr;
    
    i = c->cliquecount;

    rval = CCtsp_array_to_lpclique (newarr, newsize, &(c->cliques[i]));
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
        CC_IFFREE (newarr, int);
        CCtsp_buildcut_abort (cuts);
        return rval;
    }
    c->cliquecount++;
    CC_IFFREE (newarr, int);
    return 0;
}

void CCtsp_buildcut_abort (CCtsp_cutinfo *cuts)
{
    CCtsp_free_lpcut_in (cuts->current);
    CC_IFFREE (cuts->current, CCtsp_lpcut_in);
}

int CCtsp_buildcut_finish (CCtsp_cutinfo *cuts, int rhs)
{
    CCtsp_lpcut_in *c = cuts->current;
    int rval;

#ifdef DUMP_BUILDCUT
    {
        int i, j, tmp;
        printf ("new buildcut (%d):", c->cliquecount);
        for (i=0; i<c->cliquecount; i++) {
            printf (" (");
            CC_FOREACH_NODE_IN_CLIQUE (j, c->cliques[i], tmp) {
                printf ("%d ",j);
            }
            printf (")");
        }
        printf (" >= %d\n", rhs);
        fflush (stdout);
    }
#endif

    c->rhs = rhs;
    c->sense = 'G';
    c->branch = 0;

    rval = CCtsp_construct_skeleton (c,
            CCcut_SRK_original_ncount (&cuts->expand));
    CCcheck_rval (rval, "CCtsp_construct_skeleton failed");

    c->next = *cuts->clist;
    (*cuts->clist) = c;
    cuts->current = (CCtsp_lpcut_in *) NULL;
    (*cuts->cutcount)++;

    rval = 0;
 CLEANUP:
    if (rval) {
        CCtsp_free_lpcut_in (c);
    }
    return rval;
}

int CCtsp_hypergraph_type (CCtsp_lpcut *c)
{
   return (c && c->dominocount == 0 && c->semicount == 0 &&
                c->cliquemult == 0  && c->coefcount == 0 &&
                c->TP_handles == (int *) NULL);
   /* We construct skeletons for branching cliques, so do not check sense */
}

int CCtsp_hypergraph_type_in (CCtsp_lpcut_in *c)
{
   return (c && c->dominocount == 0 && c->semicount == 0 &&
                c->cliquemult == 0  && c->coefcount == 0 &&
                c->TP_handles == (CCtsp_lpclique *) NULL);
}

int CCtsp_create_cut (int ncount, CCtsp_lpcut_in **cuts, int nsets,
        int *setsize, int **sets, int *setmultipliers, int nsemicuts,
        int **semisize, int ***semicuts, int *semimultipliers, int rhs)   
{
    int rval = 0, i, j;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL, **clast;

    if (nsets == 0 && nsemicuts == 0) {
        fprintf (stderr, "no information to create a cut\n");
        rval = 1; goto CLEANUP;
    }

    clast = cuts;
    while ((*clast) != (CCtsp_lpcut_in *) NULL) {
        clast = &((*clast)->next);
    }

    CC_MALLOC (c, 1, CCtsp_lpcut_in);
    CCtsp_init_lpcut_in (c);

#if 0
{
    int k = 0;
    printf ("CCtsp_create_cut (%d, %d) ...\n", nsets, nsemicuts);
    fflush (stdout);
    for (i = 0; i < nsets; i++) {
        printf ("Set %d: ", i);
        for (j = 0; j < setsize[i]; j++)  printf ("%d ", sets[i][j]);
        printf ("\n"); fflush (stdout);
    }
    for (i = 0; i < nsemicuts; i++) {
        printf ("Semi-cut %d\n", i); fflush (stdout);
        for (k = 0; k < 2; k++) {
            printf ("Side %d (%d cities): ", k, semisize[i][k]); 
            for (j = 0; j < semisize[i][k]; j++) {
                printf ("%d ", semicuts[i][k][j]);
            }
            printf ("\n"); fflush (stdout);
        }
    }
}
#endif

    if (nsets) {
        CC_MALLOC (c->cliques, nsets, CCtsp_lpclique);
        for (i = 0; i < nsets; i++) CCtsp_init_lpclique (&c->cliques[i]);
        for (i = 0; i < nsets; i++) {
            rval = CCtsp_array_to_lpclique (sets[i], setsize[i],
                                            &(c->cliques[i]));
            CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        }
        c->cliquecount = nsets;
    }

    if (nsemicuts) {
        CC_MALLOC (c->semicuts, nsemicuts, CCtsp_lpdomino);
        for (i = 0; i < nsemicuts; i++) CCtsp_init_lpdomino (&c->semicuts[i]);
        for (i = 0; i < nsemicuts; i++) {
            for (j = 0; j < 2; j++) {
                rval = CCtsp_array_to_lpclique (semicuts[i][j], semisize[i][j],
                                       &(c->semicuts[i].sets[j]));
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
            }
        }
        c->semicount = nsemicuts;
    }

    if (setmultipliers) {
        CC_MALLOC (c->cliquemult, nsets, int);
        for (i = 0; i < nsets; i++) c->cliquemult[i] = setmultipliers[i];
    }

    if (semimultipliers) {
        CC_MALLOC (c->semimult, nsemicuts, int);
        for (i = 0; i < nsemicuts; i++) c->semimult[i] = semimultipliers[i];
    }

    c->rhs = rhs;
    c->sense = 'G';
    c->branch = 0;
    if (CCtsp_hypergraph_type_in (c)) {
        rval = CCtsp_construct_skeleton (c, ncount);
        CCcheck_rval (rval, "CCtsp_construct_skeleton failed");
    }

    (*clast) = c;
    c->next = (CCtsp_lpcut_in *) NULL;

CLEANUP:
    if (rval && c) {
        CCtsp_free_lpcut_in (c);
        CC_FREE (c, CCtsp_lpcut_in);
    }
    return rval;
}

static int grab_nonzero_x (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double **new_x, double tol)
{
    int rval = 0, i, count;

    *new_ecount = 0;
    *new_elist = (int *) NULL;
    *new_x = (double *) NULL;

    for (i = 0, count = 0; i < ecount; i++) {
        if (x[i] > tol) count++;
    }

    CC_MALLOC (*new_elist, 2*count, int);
    CC_MALLOC (*new_x, count, double);

    for (i = 0, count = 0; i < ecount; i++) {
        if (x[i] > tol) {
            (*new_elist)[2*count] = elist[2*i];
            (*new_elist)[2*count+1] = elist[2*i+1];
            (*new_x)[count] = x[i];
            count++;
        }
    }
    *new_ecount = count;

CLEANUP:
    return rval;
}

int CCtsp_test_pure_comb (int ncount, CCtsp_lpcut_in *c, int *yes_no,
        int *handle)
{
    int rval = 0, i, marked, ihandle, *marks = (int *) NULL;

    *yes_no = 0;
    if (handle) *handle = -1;

    if (!CCtsp_hypergraph_type_in (c) || c->cliquecount < 4 ||
                                         c->cliquecount % 2) {
        goto CLEANUP;
    }

    rval = CCtsp_find_pure_handle (ncount, c, &ihandle);
    CCcheck_rval (rval, "CCtsp_find_pure_handle failed");
    if (ihandle == -1) goto CLEANUP;

    CC_MALLOC (marks, ncount, int);
    CCtsp_mark_cut (c, marks, 0);

    CCtsp_mark_clique (&c->cliques[ihandle], marks, 1);
    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (!marked) goto CLEANUP;
            CCtsp_is_clique_marked (&c->cliques[i], marks, 0, &marked);
            if (!marked) goto CLEANUP;
        }
    }
    CCtsp_mark_clique (&c->cliques[ihandle], marks, 0);

    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (marked) goto CLEANUP;
            CCtsp_mark_clique (&c->cliques[i], marks, 1);
        }
    }

    *yes_no = 1;
    if (handle) *handle = ihandle;

CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

int CCtsp_test_pseudocomb (int ncount, CCtsp_lpcut_in *c, int handle,
        int *yes_no)
{
    int rval = 0, i, k, marked, *ends = (int *) NULL, *marks = (int *) NULL;

    *yes_no = 0;
    if (c->cliquecount <= 1 || c->cliquecount % 2 || c->sense != 'G') {
        printf ("bad cliquecount or sense in pseudocomb\n"); fflush (stdout);
        goto CLEANUP;
    }

    CC_MALLOC (marks, ncount, int);
    CCtsp_mark_cut (c, marks, 0);

    /* Teeth intersect H and are not contained in H */

    CCtsp_mark_clique (&c->cliques[handle], marks, 1);
    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (!marked) goto CLEANUP;
            CCtsp_is_clique_marked (&c->cliques[i], marks, 0, &marked);
            if (!marked) goto CLEANUP;
        }
    }
    CCtsp_mark_clique (&c->cliques[0], marks, 0);

    /* Big teeth are pairwise disjoint */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&c->cliques[i], &k);
            if (k >= 3) {
                CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
                if (marked) goto CLEANUP;
                CCtsp_mark_clique (&c->cliques[i], marks, 1);
            }
        }
    }
    for (i = 1; i < c->cliquecount; i++) {
        CCtsp_mark_clique (&c->cliques[i], marks, 0);
    }

    /* No small tooth is contained in a big tooth */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&c->cliques[i], &k);
            if (k >= 3) {
                CCtsp_mark_clique (&c->cliques[i], marks, i + 1);
            }
        }
    }
    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&c->cliques[i], &k);
            if (k < 3) {
                rval = CCtsp_clique_to_array (&c->cliques[i], &ends, &k);
                if (rval) {
                    fprintf (stderr, "CCtsp_clique_to_array failed\n");
                    goto CLEANUP;
                }
                if (ends[0] != 0 && ends[0] == ends[1]) goto CLEANUP;
                CC_IFFREE (ends, int);
            }
        }
    }

    *yes_no = 1;

CLEANUP:
    CC_IFFREE (marks, int);
    CC_IFFREE (ends, int);
    return rval;
}

int CCtsp_test_teeth_disjoint (int ncount, CCtsp_lpcut_in *c, int handle,
        int *yes_no)
{
    int rval = 0, i, marked, *marks = (int *) NULL;

    *yes_no = 0;

    CC_MALLOC (marks, ncount, int);
    CCtsp_mark_cut (c, marks, 0);

    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (marked) goto CLEANUP;
            CCtsp_mark_clique (&c->cliques[i], marks, 1);
        }
    }

    *yes_no = 1;

CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

int CCtsp_find_pure_handle (int ncount, CCtsp_lpcut_in *c, int *handle)
{
    int rval = 0, i, test, *marks = (int *) NULL;

    *handle = -1;
    if (!CCtsp_hypergraph_type_in (c) || c->cliquecount % 2 ||
                                        c->cliquecount < 4) {
        goto CLEANUP;
    }

    CC_MALLOC (marks, ncount, int);
    CCtsp_mark_cut (c, marks, 0);

    CCtsp_mark_clique (&c->cliques[0], marks, 1);
    CCtsp_is_clique_marked (&c->cliques[1], marks, 1, &test);
    if (test) {
        CCtsp_is_clique_marked (&c->cliques[2], marks, 1, &test);
        if (test) {
            *handle = 0; goto CLEANUP;
        } else {
            *handle = 1; goto CLEANUP;
        }
    } else {
        for (i = 2; i < c->cliquecount; i++) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &test);
            if (test) {
                *handle = i;
                goto CLEANUP;
            }
        }
    }

CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

int CCtsp_truncate_cutlist (CCtsp_lpcut_in **cuts, int ncount, int ecount,
        int *elist, double *x, int maxcuts, CCrandstate *rstate)
{
    int rval = 0, i, count = 0, *perm = (int *) NULL;
    CCtsp_lpcut_in *c, *cnext;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    double *vlist = (double *) NULL;
    CCtsp_lpgraph lg;

    CCtsp_init_lpgraph_struct (&lg);

    if (maxcuts <= 0) {
        for (c = *cuts; c; c = cnext) {
            cnext = c->next;
            CCtsp_free_lpcut_in (c);
        }
        *cuts = (CCtsp_lpcut_in *) NULL;
        goto CLEANUP;
    }

    for (c = *cuts; c; c = c->next) count++;
    if (count > maxcuts) {
        rval = CCtsp_build_lpgraph (&lg, ncount, ecount, elist, (int *) NULL);
        CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
        rval = CCtsp_build_lpadj (&lg, 0, ecount);
        CCcheck_rval (rval, "CCtsp_build_lpadj failed")

        CC_MALLOC (vlist, count, double);
        CC_MALLOC (clist, count, CCtsp_lpcut_in *);
        CC_MALLOC (perm, count, int);
        for (i = 0, c = *cuts; c; c = c->next, i++) {
            clist[i] = c;
            vlist[i] = CCtsp_cutprice (&lg, c, x);
            perm[i] = i;
        }

        CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
        for (i = maxcuts; i < count; i++) {
            CCtsp_free_lpcut_in (clist[perm[i]]);
        }

        *cuts = (CCtsp_lpcut_in *) NULL;
        for (i = 0; i < maxcuts; i++) {
            clist[perm[i]]->next = *cuts;
            *cuts = clist[perm[i]];
        }
    }

CLEANUP:
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CCtsp_free_lpgraph (&lg);
    return rval;
}

#ifdef CCtsp_USE_DOMINO_CUTS

#define CC_DP_NEIGHBORHOOD 0.96 /* 1.0=full, 0.55=fast; recommended 0.96  */
#define CC_DP_TIME  10.0        /* large instances 100.0, medium 10.0 */

#if 0   /* Two definitions for old DP library. */
int DPseparator (int nnodes, int nedges, int* edges, double* weigh, int* nIneq,
    int** nDominoes, int*** nAset, int*** nBset, int** nHandle, int**** Aset,
    int **** Bset, int*** Handle, const char *boss_name, double percentage,
    double ddp_heuristic_maxtime);
int DPtighten (int n_nodes, int n_edges, int *edges, double *weight,
    int n_dominos, int *n_aset, int *n_bset, int n_handle, int **aset,
    int **bset, int *handle, int **new_n_aset, int **new_n_bset,
    int *new_n_handle, int ***new_aset, int ***new_bset, int **new_handle,
    double *violation);
#endif

static int domino_tighten_lpcut_in (int ncount, int ecount, int *elist,
       double *x, CCtsp_lpcut_in *old, CCtsp_lpcut_in **new, double *slack);
static int comb2dp (int ncount, CCtsp_lpcut_in *comb, CCtsp_lpcut_in **dp);
static int grab_sorted_edges (CC_SRKgraph *G, int *oncount, int *oecount,
        int **polist, double **polen, CC_SRKexpinfo *expand);
static int build_dp_cut_expand (CCtsp_lpcut_in **cut, int ndomino, int *Acount,
        int **A, int *Bcount, int **B, int handlecount, int *handle,
        CC_SRKexpinfo *expand, int ncount, int *comb);
static void free_raw_dominos (int nIneq, int *nHandle, int **Handle,
    int *nDominoes, int **nAset, int ***Aset, int **nBset, int ***Bset);
static int check_raw_domino (int ncount, int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B, int *valid, int *comb,
        int *large_domino);
static void  print_raw_domino (int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B);

int CCtsp_DP_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, int safe_shrink,
        const char *dombossname, double *fullzeit, double *maxviol)
{
    int rval = 0, i, k, comb, nIneq = 0, nCombs = 0;
    int oncount = 0, oecount = 0, emptycount = 0, tcount = 0;
    int newecount, *newelist = (int *) NULL;
    int *nDominoes = (int *) NULL;
    int **nAset = (int **) NULL, ***Aset = (int ***) NULL;
    int **nBset = (int **) NULL, ***Bset = (int ***) NULL;
    int *nHandle = (int *) NULL, **Handle = (int **) NULL;
    int *oelist = (int *) NULL;
    double cval, val, ttotal = 0.0, dzeit, tzeit = 0.0;
    double *olen = (double *) NULL, *newx = (double *) NULL;
    CCtsp_lpcut_in *c;
    CC_SRKexpinfo expand;
    CC_SRKgraph G;
    CCtsp_lpgraph g;
    CCtsp_lpcut_in *new = (CCtsp_lpcut_in *) NULL;
    double szeit = CCutil_zeit();

    *cutcount = 0;
    if (maxviol) *maxviol = 0.0;
    CCtsp_init_lpgraph_struct (&g);
    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    CCcheck_rval (rval, "grab_nonzero_x failed");

    rval = CCtsp_build_lpgraph (&g, ncount, newecount, newelist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");

    rval = CCtsp_build_lpadj (&g, 0, newecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    printf ("Call DPseparator ...\n"); fflush (stdout);
    if (dombossname) { printf ("Use Domino Boss: %s\n", dombossname); }

    CCcut_SRK_init_expinfo (&expand);
    CCcut_SRK_init_graph (&G);
    rval = CCcut_SRK_buildgraph (&G, ncount, ecount, elist, x);
    CCcheck_rval (rval, "CCcut_SRK_buildgraph failed");
    rval = CCcut_SRK_defluff (&G);
    CCcheck_rval (rval, "CCcut_SRK_defluff failed");

    if (safe_shrink) {
        CCcut_SRK_identify_paths_to_edges (&G, &k, 0);
        CCcut_SRK_identify_one_triangles (&G, &k, (CC_SRKnode *) NULL, 0.001,
                                          2.0, 0);
    }

    rval = grab_sorted_edges (&G, &oncount, &oecount, &oelist, &olen, &expand);
    CCcheck_rval (rval, "grab_sorted_edges failed");


    rval = DPseparator (oncount, oecount, oelist, olen, &nIneq, &nDominoes,
                        &nAset, &nBset, &nHandle, &Aset, &Bset, &Handle,
                        dombossname, CC_DP_NEIGHBORHOOD, CC_DP_TIME);
    CCcheck_rval (rval, "DPseparator failed");

    for (i = 0; i < nIneq; i++) {
        if (nHandle[i] > 0 && nHandle[i] < ncount) {
            rval = build_dp_cut_expand (&c, nDominoes[i], nAset[i], Aset[i],
                     nBset[i], Bset[i], nHandle[i], Handle[i], &expand,
                     ncount, &comb);
            CCcheck_rval (rval, "build_dp_cut_expand failed");
            if (c) {
                nCombs += comb;
                dzeit = CCutil_zeit ();
                cval = CCtsp_cutprice (&g, c, newx);
                rval = domino_tighten_lpcut_in (ncount, newecount, newelist,
                           newx, c, &new, &val);
                CCcheck_rval (rval, "domino_tighten_lpcut_in failed");

                if (new) {
                    double tval;

                    tval = CCtsp_cutprice (&g, new, newx);
                    if (tval < cval) {
                        ttotal += (cval - tval);
                        tcount++;

                        CCtsp_free_lpcut_in (c);
                        CC_FREE (c, CCtsp_lpcut_in);
                        c = new;
                        new = (CCtsp_lpcut_in *) NULL;
                    } else {
                        CCtsp_free_lpcut_in (new);
                        CC_FREE (new, CCtsp_lpcut_in);
                    }
                }
                tzeit +=  (CCutil_zeit () - dzeit);
                c->next = *cuts;
                *cuts = c;
                (*cutcount)++;
            }
        } else {
            emptycount++;
        }
    }

    if (emptycount) {
        printf ("Skipped %d empty-handled DP cuts\n", emptycount);
        fflush (stdout);
    }
    printf ("Found %d domino cuts (%d combs) ...\n", nIneq, nCombs);
    fflush (stdout);

    if (tcount) {
        printf ("   Tightened %d cuts (avg %f) in %.2f seconds\n",
                    tcount, ttotal / (double) tcount, tzeit);
    }

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    CC_IFFREE (oelist, int);
    CC_IFFREE (olen, double);
    CCcut_SRK_free_expinfo (&expand);
    CCcut_SRK_free_graph (&G);
    free_raw_dominos (nIneq, nHandle, Handle, nDominoes, nAset, Aset, nBset,
                      Bset);
    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CCtsp_free_lpgraph (&g);
    return rval;
}

int CCtsp_TP_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, double *fullzeit, double *maxviol)
{
    int rval = 0, request = 1024, rhs;
    int i, j, tricount, domcount0, domcount1;
    CCtsp_lpcut_in **clast, *c = (CCtsp_lpcut_in *) NULL;
    D3cfg_t cfg = zeroD3cfg;
    D3tpIneq_t TP=zeroD3tpIneq;

/*  D3cfgDump(&cfg,stderr);  -- prints D3cfg parameters */
    D3tpIneqInit(&TP,request);

    if (maxviol) *maxviol = 0.0;
    *cutcount = 0;

    cfg.verb = 0;
    cfg.saveinput = 0;
    cfg.tpwlim = 1.80*ONE32b;   /* Suggested by Daniel 22 May 2016 */

    clast = cuts;  /* Append triomino ineq to end of cuts linked list */
    while ((*clast) != (CCtsp_lpcut_in *) NULL) {
        clast = &((*clast)->next);
    }

    double szeit = CCutil_zeit();

    printf ("  Call TPseparator (%d, %d)\n", ncount, ecount);
    fflush (stdout);
    rval = TPseparator (&cfg, (char *) NULL, ncount, ecount,
                       (const int*const) elist, x, &TP);
    CCcheck_rval (rval, "TPseparator failed");
    printf ("  TPseparator returns %d cuts, %.2f seconds\n", TP.nIneq,
               CCutil_zeit() - szeit);
    fflush (stdout);

    for (i = 0; i < TP.nIneq; i++) {
        if (TP.n2Dom[i] == 0) { /* Treat TP-cut as a DP-cut */
            if (TP.nH0andle[i] != TP.nH1andle[i]) {
                printf ("TP to DP ERROR: handles of different size\n");
                continue;
            }
            for (j = 0; j < TP.nH0andle[i]; j++) {
                if (TP.H0andle[i][j] != TP.H1andle[i][j]) break;
            }
            if (j != TP.nH0andle[i]) {
                printf ("TP to DP ERROR: handles diff (ordered) node sets\n");
                continue;
            }
            for (j = 0; j < TP.n1Dom[i]; j++) {
                if (TP.H1Dom[i][j] != 2) break;
            }
            if (j != TP.n1Dom[i]) {
                printf ("TP to DP ERROR: domino for only one handle\n");
                continue;
            }
            rval = CCtsp_build_dp_cut (&c, TP.n1Dom[i],
                                       TP.nAset[i], TP.Aset[i],
                                       TP.nBset[i], TP.Bset[i],
                                       TP.nH0andle[i], TP.H0andle[i]);
            CCcheck_rval (rval, "CCtsp_build_dp_cut failed");
            if (c) {
                (*clast) = c;
                c->next = (CCtsp_lpcut_in *) NULL;
                clast = &(c->next);
            }
            continue;
        }

        CC_MALLOC (c, 1, CCtsp_lpcut_in);
        CCtsp_init_lpcut_in (c);

        CC_MALLOC (c->TP_handles, 2, CCtsp_lpclique);
        rval = CCtsp_array_to_lpclique (TP.H0andle[i], TP.nH0andle[i],
                                        &(c->TP_handles[0]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        rval = CCtsp_array_to_lpclique (TP.H1andle[i], TP.nH1andle[i],
                                        &(c->TP_handles[1]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

        domcount0 = domcount1 = 0;
        for (j = 0; j < TP.n1Dom[i]; j++) {
            if (TP.H1Dom[i][j] == 0 || TP.H1Dom[i][j] == 2) domcount0++;
            if (TP.H1Dom[i][j] == 1 || TP.H1Dom[i][j] == 2) domcount1++;
        }
        if (domcount0) {
            CC_MALLOC (c->TP_dominos0, domcount0, CCtsp_lpdomino);
            for (j = 0; j < domcount0; j++) {
                CCtsp_init_lpdomino (&c->TP_dominos0[j]);
            }
        }
        if (domcount1) {
            CC_MALLOC (c->TP_dominos1, domcount1, CCtsp_lpdomino);
            for (j = 0; j < domcount1; j++) {
                CCtsp_init_lpdomino (&c->TP_dominos1[j]);
            }
        }
        domcount0 = domcount1 = 0;
        for (j = 0; j < TP.n1Dom[i]; j++) {
            if (TP.H1Dom[i][j] == 0 || TP.H1Dom[i][j] == 2) {
                rval = CCtsp_array_to_lpclique (TP.Aset[i][j], TP.nAset[i][j],
                                       &(c->TP_dominos0[domcount0].sets[0]));
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                rval = CCtsp_array_to_lpclique (TP.Bset[i][j], TP.nBset[i][j],
                                       &(c->TP_dominos0[domcount0].sets[1]));
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                domcount0++;
            }
            if (TP.H1Dom[i][j] == 1 || TP.H1Dom[i][j] == 2) {
                rval = CCtsp_array_to_lpclique (TP.Aset[i][j], TP.nAset[i][j],
                                       &(c->TP_dominos1[domcount1].sets[0]));
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                rval = CCtsp_array_to_lpclique (TP.Bset[i][j], TP.nBset[i][j],
                                       &(c->TP_dominos1[domcount1].sets[1]));
                CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
                domcount1++;
            }
        }
        c->TP_domcount0 = domcount0;
        c->TP_domcount1 = domcount1;

        tricount = TP.n2Dom[i];
        if (tricount) {
            CC_MALLOC (c->TP_tsets, tricount, CCtsp_lpclique);
            CC_MALLOC (c->TP_semicuts0, tricount, CCtsp_lpdomino);
            CC_MALLOC (c->TP_semicuts1, tricount, CCtsp_lpdomino);
        }
        for (j = 0; j < tricount; j++) {
            CCtsp_init_lpclique (&c->TP_tsets[j]);
            CCtsp_init_lpdomino (&c->TP_semicuts0[j]);
            CCtsp_init_lpdomino (&c->TP_semicuts1[j]);
        }
        for (j = 0; j < tricount; j++) {
            rval = build_triomino (ncount, TP.nTset[i][j],  TP.Tset[i][j],
                                           TP.nH0set[i][j], TP.H0set[i][j],
                                           TP.nH1set[i][j], TP.H1set[i][j],
              &(c->TP_tsets[j]), &(c->TP_semicuts0[j]), &(c->TP_semicuts1[j]));
            CCcheck_rval (rval, "build_triomino failed");

        }
        c->TP_tricount = tricount;

        rhs = 2 + 4*tricount + 3*(domcount0 + domcount1);

        c->rhs = rhs;
        c->sense = 'G';
        c->branch = 0;

        /* printf ("Here is a triomino cut\n"); CCtsp_print_lpcut_in (c); */

        (*clast) = c;
        c->next = (CCtsp_lpcut_in *) NULL;
        clast = &(c->next);
    } 

CLEANUP:
    D3tpIneqClear (&TP);
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    return rval;
}

static int build_triomino (int ncount, int nTset, int *Tset, int nH0set,
        int *H0set, int nH1set, int *H1set,
        CCtsp_lpclique *ptset, CCtsp_lpdomino *psemi0, CCtsp_lpdomino *psemi1)
{
    int rval = 0, i, *marks = (int *) NULL, na, acount, *alist = (int *) NULL;

    if (nTset == 0 || nH0set == 0 || nH1set == 0) {
        fprintf (stderr, "triomino has an empty set\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_array_to_lpclique (Tset, nTset, ptset);
    CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

    CC_MALLOC (alist, nTset, int);
    CC_MALLOC (marks, ncount, int);
    for (i = 0; i < nH0set; i++) marks[H0set[i]] = 0;
    for (i = 0; i < nH1set; i++) marks[H1set[i]] = 0;
    for (i = 0; i < nTset; i++) marks[Tset[i]] = 1;

    na = 0;
    for (i = 0; i < nH0set; i++) { if (marks[H0set[i]]) na++; }

    if (na == 0) {
        rval = CCtsp_array_to_lpclique (H0set, nH0set, &(psemi0->sets[0]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        rval = CCtsp_array_to_lpclique (Tset, nTset, &(psemi0->sets[1]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
    } else if (na < nH0set) {
        fprintf (stderr, "H0set has proper intersection with Tset\n");
        rval = 1; goto CLEANUP;
    } else {
        rval = CCtsp_array_to_lpclique (H0set, nH0set, &(psemi0->sets[0]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        for (i = 0; i < nH0set; i++) marks[H0set[i]] = 2;
        acount = 0;
        for (i = 0; i < nTset; i++) {
            if (marks[Tset[i]] == 1) { alist[acount++] = Tset[i]; }
        }
        if (acount == 0) {
            fprintf (stderr, "Tset/H0set is empty\n");
            rval = 1; goto CLEANUP;
        } else {
            rval = CCtsp_array_to_lpclique (alist, acount, &(psemi0->sets[1]));
            CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        }
        for (i = 0; i < nH0set; i++) marks[H0set[i]] = 1;
    }

    na = 0;
    for (i = 0; i < nH1set; i++) { if (marks[H1set[i]]) na++; }

    if (na == 0) {
        rval = CCtsp_array_to_lpclique (H1set, nH1set, &(psemi1->sets[0]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        rval = CCtsp_array_to_lpclique (Tset, nTset, &(psemi1->sets[1]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
    } else if (na < nH1set) {
        fprintf (stderr, "H1set has proper intersection with Tset\n");
        rval = 1; goto CLEANUP;
    } else {
        rval = CCtsp_array_to_lpclique (H1set, nH1set, &(psemi1->sets[0]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        for (i = 0; i < nH1set; i++) marks[H1set[i]] = 2;
        acount = 0;
        for (i = 0; i < nTset; i++) {
            if (marks[Tset[i]] == 1) { alist[acount++] = Tset[i]; }
        }
        if (acount == 0) {
            fprintf (stderr, "Tset/H1set is empty\n");
            rval = 1; goto CLEANUP;
        } else {
            rval = CCtsp_array_to_lpclique (alist, acount, &(psemi1->sets[1]));
            CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        }
        for (i = 0; i < nH1set; i++) marks[H1set[i]] = 1;
    }
    
CLEANUP:
    CC_IFFREE (marks, int);
    CC_IFFREE (alist, int);
    return rval;
}

int CCtsp_domino_tighten_lp (CCtsp_lpcuts *cuts, CCtsp_lpcut_in **cutsout,
        int *cutcount, int ncount, int ecount, int *elist, double *x,
        double testtol, int maxcuts, double *viol, CCrandstate *rstate,
        double *fullzeit)
{
    int rval = 0, i, clistsize = 0, vlistsize = 0, count = 0;
    int *perm = (int *) NULL, *newelist = (int *) NULL, newecount;
    double *newx = (double *) NULL, cval, val, maxviol = 0.0;
    double *vlist = (double *) NULL;
    CCtsp_lpcut_in old, *new = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    CCtsp_lpgraph g;
    double szeit = CCutil_zeit();

    CCtsp_init_lpgraph_struct (&g);
    if (viol) *viol = -1.0;
    *cutcount = 0;

    if (!cuts || !cuts->cutcount) goto CLEANUP;
    for (i = 0; i < cuts->cutcount; i++) {
        if (!cuts->cuts[i].branch && cuts->cuts[i].dominocount > 0) {
            break;
        }
    }
    if (i == cuts->cutcount) goto CLEANUP;  /* No DP cuts in LP */

    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    CCcheck_rval (rval, "grab_nonzero_x failed");

    rval = CCtsp_build_lpgraph (&g, ncount, newecount, newelist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");

    rval = CCtsp_build_lpadj (&g, 0, newecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    for (i = 0; i < cuts->cutcount; i++) {
        if (!cuts->cuts[i].branch && cuts->cuts[i].dominocount > 0) {
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &old);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            cval = CCtsp_cutprice (&g, &old, newx);
            if (cval < testtol) {
                rval = domino_tighten_lpcut_in (ncount, newecount, newelist,
                                                newx, &old, &new, &val);
                CCcheck_rval (rval, "domino_tighten_lpcut_in failed");
                CCtsp_free_lpcut_in (&old);

                if (new && val < -CCtsp_MIN_VIOL) {
                    double tval;

                    tval = CCtsp_cutprice (&g, new, newx);
                    if (0 && tval < cval) {
                        printf ("Tighten Improvement: %f\n", cval-tval);
                        fflush (stdout);
                    }

                    if (count >= clistsize) {
                        void *tmp_ptr = (void *) clist;
                        rval = CCutil_reallocrus_scale (&tmp_ptr, 
                              &clistsize, count + 1, 1.3,
                              sizeof (CCtsp_lpcut_in *));
                        CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
                        clist = (CCtsp_lpcut_in **) tmp_ptr;
                    }
                    if (count >= vlistsize) {
                        void *tmp_ptr = (void *) vlist;
                        rval = CCutil_reallocrus_scale (&tmp_ptr, 
                            &vlistsize, count + 1, 1.3, sizeof (double));
                        CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
                        vlist = (double *) tmp_ptr;
                    }
                    clist[count] = new;
                    vlist[count] = tval;
                    count++;
                    new = (CCtsp_lpcut_in *) NULL;
                } else {
                    if (new) {
                        CCtsp_free_lpcut_in (new);
                        CC_FREE (new, CCtsp_lpcut_in);
                    }
                }
            } else { 
                CCtsp_free_lpcut_in (&old);
            }
        }
    }

    if (count) {
        CC_MALLOC (perm, count, int);
        for (i = 0; i < count; i++)  perm[i] = i;
        if (count > maxcuts) {
            CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
            for (i = maxcuts; i < count; i++) {
                CCtsp_free_lpcut_in (clist[perm[i]]);
            }
            count = maxcuts;
        }
        for (i = 0; i < count; i++) {
            if (vlist[perm[i]] < maxviol)
                maxviol = vlist[perm[i]];
            clist[perm[i]]->next = *cutsout;
            *cutsout = clist[perm[i]];
        }
    }

    *cutcount = count;
    if (viol) *viol = -maxviol;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CCtsp_free_lpgraph (&g);
    return rval;
}

static int CC_DPhead = 0;
#define CC_DPlim 500

int CCtsp_domino_tighten_pool (CCtsp_lpcuts *cuts, CCtsp_lpcut_in **cutsout,
        int *cutcount, int ncount, int ecount, int *elist, double *x,
        double testtol, int maxcuts, double *viol, CCrandstate *rstate,
        double *fullzeit)
{
    int rval = 0, i, trys, look, clistsize = 0, vlistsize = 0, count = 0;
    int *perm = (int *) NULL, *newelist = (int *) NULL, newecount;
    double cval, val, maxviol = 0.0;
    double *newx = (double *) NULL, *vlist = (double *) NULL;
    CCtsp_lpcut_in old, *new = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    CCtsp_lpgraph g;
    double szeit = CCutil_zeit ();

    CCtsp_init_lpgraph_struct (&g);
    if (viol) *viol = -1.0;
    *cutcount = 0;
    if (!cuts || !cuts->cutcount) goto CLEANUP;

/*
    printf ("CCtsp_domino_tighten_lp (%f) ...\n", testtol); fflush (stdout);
*/

    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    CCcheck_rval (rval, "grab_nonzero_x failed");

    rval = CCtsp_build_lpgraph (&g, ncount, newecount, newelist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");

    rval = CCtsp_build_lpadj (&g, 0, newecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    trys = 0;
    look = 0;
    while (trys < CC_DPlim && look < cuts->cutcount) {
        if (CC_DPhead >= cuts->cutcount) CC_DPhead = 0;
        if (!cuts->cuts[CC_DPhead].branch &&
             cuts->cuts[CC_DPhead].dominocount > 0) {
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[CC_DPhead]),
                                            &old);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            cval = CCtsp_cutprice (&g, &old, newx);
            if (cval < testtol) {
                trys++;
                rval = domino_tighten_lpcut_in (ncount, newecount, newelist,
                           newx, &old, &new, &val);
                CCcheck_rval (rval, "domino_tighten_lpcut_in failed");
                CCtsp_free_lpcut_in (&old);

                if (new && val < -CCtsp_MIN_VIOL) {
                    double tval;

                    tval = CCtsp_cutprice (&g, new, newx);
                    if (0 && tval < cval) {
                        printf ("Tighten Improvement: %f\n", cval - tval);
                        fflush (stdout);
                    }

                    if (count >= clistsize) {
                        void *tmp_ptr = (void *) clist;
                        rval = CCutil_reallocrus_scale (&tmp_ptr, 
                              &clistsize, count + 1, 1.3,
                              sizeof (CCtsp_lpcut_in *));
                        CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
                        clist = (CCtsp_lpcut_in **) tmp_ptr;
                    }
                    if (count >= vlistsize) {
                        void *tmp_ptr = (void *) vlist;
                        rval = CCutil_reallocrus_scale (&tmp_ptr, 
                            &vlistsize, count + 1, 1.3, sizeof (double));
                        CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
                        vlist = (double *) tmp_ptr;
                    }
                    clist[count] = new;
                    vlist[count] = tval;
                    count++;
                    new = (CCtsp_lpcut_in *) NULL;
                } else {
                    if (new) {
                        CCtsp_free_lpcut_in (new);
                        CC_FREE (new, CCtsp_lpcut_in);
                    }
                }
            } else { 
                CCtsp_free_lpcut_in (&old);
            }
        }
        CC_DPhead++;
        look++;
    }

    printf ("DPtight Totals: %d look, %d try, %d winners, %.2f seconds\n",
               look, trys, count, CCutil_zeit () - szeit);
    fflush (stdout);

    if (count) {
        perm = CC_SAFE_MALLOC (count, int);
        CCcheck_NULL (perm, "out of memory for perm");
        for (i = 0; i < count; i++) {
            perm[i] = i;
        }
        if (count > maxcuts) {
            CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
            for (i = maxcuts; i < count; i++) {
                CCtsp_free_lpcut_in (clist[perm[i]]);
            }
            count = maxcuts;
        }
        for (i = 0; i < count; i++) {
            if (vlist[perm[i]] < maxviol)
                maxviol = vlist[perm[i]];
            clist[perm[i]]->next = *cutsout;
            *cutsout = clist[perm[i]];
        }
    }

    *cutcount = count;
    if (viol) *viol = -maxviol;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CCtsp_free_lpgraph (&g);
    return rval;
}

static int domino_tighten_lpcut_in (int ncount, int ecount, int *elist,
       double *x, CCtsp_lpcut_in *old, CCtsp_lpcut_in **new, double *slack)
{
    int rval = 0,i, chkvalid, chkcomb, large_domino = 0;
    int ndomino = 0, handlecount, *handle = (int *) NULL;
    int *Acount = (int *) NULL, **A = (int **) NULL;
    int *Bcount = (int *) NULL, **B = (int **) NULL;
    int new_ndomino = 0, new_handlecount, *new_handle = (int *) NULL;
    int *new_Acount = (int *) NULL, **new_A = (int **) NULL;
    int *new_Bcount = (int *) NULL, **new_B = (int **) NULL;

    *slack = 100.0;
    *new = (CCtsp_lpcut_in *) NULL;

    if (old->branch != 0) {
        fprintf (stderr, "try to domino tighten a branch cut\n");
        rval = 1; goto CLEANUP;
    }
    if (old->sense != 'G') {
        fprintf (stderr, "try to domino tighten a <= cut\n");
        rval = 1; goto CLEANUP;
    }
    if (old->dominocount <= 0) {
        fprintf (stderr, "try to domino tighten a non-DP cut\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_extract_dp_cut (old, &ndomino, &Acount, &A, &Bcount, &B,
                                 &handlecount, &handle);
    CCcheck_rval (rval, "CCtsp_extract_dp_cut failed");

    rval = check_raw_domino (ncount, handlecount, handle, ndomino, 
                  Acount, A, Bcount, B, &chkvalid, &chkcomb, &large_domino);
    CCcheck_rval (rval, "check_raw_domino");

    if (chkvalid == 0) {
        printf ("ERROR: Raw DP-cut is not valid\n"); fflush (stdout);
        rval = 1;  goto CLEANUP;
    }

    rval = DPtighten (ncount, ecount, elist, x, ndomino, Acount, Bcount,
                      handlecount, (const int*const*)A, (const int*const*)B,
                      handle, &new_Acount, &new_Bcount, &new_handlecount,
                      &new_A, &new_B, &new_handle, slack); 
    CCcheck_rval (rval, "DPtighten failed");
    new_ndomino = ndomino;

    if (new_Acount) {
        rval = check_raw_domino (ncount, new_handlecount, new_handle,
                   new_ndomino, new_Acount, new_A, new_Bcount, new_B,
                   &chkvalid, &chkcomb, &large_domino);
        CCcheck_rval (rval, "check_raw_domino");

        if (chkvalid == 0) {
            printf ("ERROR: Returned DP-cut is not valid\n"); fflush (stdout);
            print_raw_domino (new_handlecount, new_handle, new_ndomino,
                              new_Acount, new_A, new_Bcount, new_B);
            rval = 1;  goto CLEANUP;
        }

        rval = CCtsp_build_dp_cut (new, new_ndomino, new_Acount, new_A,
                               new_Bcount, new_B, new_handlecount, new_handle);
        CCcheck_rval (rval, "CCtsp_build_dp_cut failed");
    }

CLEANUP:
    CC_IFFREE (Acount, int);
    CC_IFFREE (Bcount, int);
    if (A) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (A[i], int);
        }
        CC_FREE (A, int *);
    }
    if (B) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (B[i], int);
        }
       CC_FREE (B, int *);
    }
    CC_IFFREE (handle, int);

    CC_IFFREE (new_Acount, int);
    CC_IFFREE (new_Bcount, int);
    if (new_A) {
        for (i = 0; i < new_ndomino; i++) {
            CC_IFFREE (new_A[i], int);
        }
        CC_FREE (new_A, int *);
    }
    if (new_B) {
        for (i = 0; i < new_ndomino; i++) {
            CC_IFFREE (new_B[i], int);
        }
       CC_FREE (new_B, int *);
    }
    CC_IFFREE (new_handle, int);

    return rval;
}

int CCtsp_comb2dp_lp (CCtsp_lpcuts *cuts, CCtsp_lpcut_in **cutsout,
        int *cutcount, int ncount, int ecount, int *elist, double *x,
        double testtol, double *viol, double *fullzeit)
{
    int rval = 0, zcount, i, test;
    int *zlist = (int *) NULL;
    double szeit, maxviol = 0.0, slack;
    double *zx = (double *) NULL, *cutval = (double *) NULL;
    CCtsp_lpcut_in comb, *dp = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *new = (CCtsp_lpcut_in *) NULL;

    szeit = CCutil_zeit();

    *cutcount = 0;
    if (!cuts || !cuts->cutcount) goto CLEANUP;

    rval = grab_nonzero_x (ecount, elist, x, &zcount, &zlist, &zx, X_FLUFF);
    CCcheck_rval (rval, "grab_nonzero_x failed");

    CC_MALLOC (cutval, cuts->cutcount, double);
    rval = CCtsp_price_cuts (cuts, ncount, zcount, zlist, zx, cutval);
    CCcheck_rval (rval, "CCtsp_price_cuts failed");

    for (i = 0; i < cuts->cutcount; i++) {
        if (cuts->cuts[i].branch || cuts->cuts[i].cliquecount % 2 ||
            cuts->cuts[i].cliquecount < 4 || cutval[i] >= testtol ||
            CCtsp_hypergraph_type (&cuts->cuts[i]) == 0) {
            continue;
        }
        rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &comb);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
        rval = CCtsp_test_pure_comb (ncount, &comb, &test, (int *) NULL);
        CCcheck_rval (rval, "CCtsp_test_pure_comb failed");
        if (test == 1) {
            rval = comb2dp (ncount, &comb, &dp);
            CCcheck_rval (rval, "comb2dp failed\n");

            rval = domino_tighten_lpcut_in (ncount, zcount, zlist, zx, dp,
                                            &new, &slack);
            CCcheck_rval (rval, "domino_tighten_lpcut_in failed");

            CCtsp_free_lpcut_in (dp);
            CC_FREE (dp, CCtsp_lpcut_in);

            if (new && slack < -CCtsp_MIN_VIOL) {
                if (new) {
                    new->next = *cutsout;
                    *cutsout = new;
                    new = (CCtsp_lpcut_in *) NULL;
                    if (slack < maxviol) maxviol = slack;
                }
            } else {
                if (new) {
                    CCtsp_free_lpcut_in (new);
                    CC_FREE (new, CCtsp_lpcut_in);
                }
            }
        }

        CCtsp_free_lpcut_in (&comb);
    }

    if (viol) *viol = -maxviol;

CLEANUP:
    if (fullzeit) *fullzeit = CCutil_zeit() - szeit;
    CC_IFFREE (zlist, int);
    CC_IFFREE (zx, double);
    CC_IFFREE (cutval, double);
    return rval;
}

static int comb2dp (int ncount, CCtsp_lpcut_in *comb, CCtsp_lpcut_in **dp)
{
    int rval = 0, i, j, ihandle, *marks = (int *) NULL, ccount;
    int ndomino = 0, hcount, *handle = (int *) NULL, *cliq = (int *) NULL;
    int *Acount = (int *) NULL, **A = (int **) NULL;
    int *Bcount = (int *) NULL, **B = (int **) NULL;

    rval = CCtsp_find_pure_handle (ncount, comb, &ihandle);
    CCcheck_rval (rval, "CCtsp_find_pure_handle failed");
    if (ihandle == -1) {
        fprintf (stderr, "could not find handle in comb\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_clique_to_array (&comb->cliques[ihandle], &handle, &hcount);
    CCcheck_rval (rval, "CCtsp_clique_to array failed");

    CC_MALLOC (marks, ncount, int);
    CCtsp_mark_cut (comb, marks, 0);
    for (i = 0; i < hcount; i++) marks[handle[i]] = 1;

    CC_MALLOC (Acount, comb->cliquecount-1, int);
    CC_MALLOC (Bcount, comb->cliquecount-1, int);
    CC_MALLOC (A, comb->cliquecount-1, int *);
    for (i = 0; i < comb->cliquecount-1; i++) A[i] = (int *) NULL;
    CC_MALLOC (B, comb->cliquecount-1, int *);
    for (i = 0; i < comb->cliquecount-1; i++) B[i] = (int *) NULL;

    for (i = 0; i < comb->cliquecount; i++) {
        if (i != ihandle) {
            rval = CCtsp_clique_to_array (&comb->cliques[i], &cliq, &ccount);
            CCcheck_rval (rval, "CCtsp_clique_to array failed");

            Acount[ndomino] = Bcount[ndomino] = 0;
            for (j = 0; j < ccount; j++) {
                if (marks[cliq[j]] == 0) Acount[ndomino]++;
                else                     Bcount[ndomino]++;
            }
            if (Acount[ndomino] == 0 || Bcount[ndomino] == 0) {
                fprintf (stderr, "bad tooth\n");
                rval = 1; goto CLEANUP;
            }
            CC_MALLOC (A[ndomino], Acount[ndomino], int);
            CC_MALLOC (B[ndomino], Bcount[ndomino], int);
            Acount[ndomino] = Bcount[ndomino] = 0;
            for (j = 0; j < ccount; j++) {
                if (marks[cliq[j]] == 0) A[ndomino][Acount[ndomino]++]=cliq[j];
                else                     B[ndomino][Bcount[ndomino]++]=cliq[j];
            }
            ndomino++;
        }
    }

    rval = CCtsp_build_dp_cut (dp, ndomino, Acount, A, Bcount, B, hcount,
                               handle);
    CCcheck_rval (rval, "CCtsp_build_dp_cut failed");

CLEANUP:
    CC_IFFREE (Acount, int);
    CC_IFFREE (Bcount, int);
    if (A) {
        for (i = 0; i < ndomino; i++) { CC_IFFREE (A[i], int); }
        CC_FREE (A, int *);
    }
    if (B) {
        for (i = 0; i < ndomino; i++) { CC_IFFREE (B[i], int); }
       CC_FREE (B, int *);
    }
    CC_IFFREE (handle, int);
    CC_IFFREE (marks, int);
    return rval;
}

static int grab_sorted_edges (CC_SRKgraph *G, int *oncount, int *oecount,
        int **polist, double **polen, CC_SRKexpinfo *expand)
{
    int rval = 0, i, ncount, ecount;
    int *elist = (int *) NULL, *perm = (int *) NULL, *olist = (int *) NULL;
    double *olen = (double *) NULL, *elen = (double *) NULL;

    rval = CCcut_SRK_grab_edges (G, &ncount, &ecount, &elist, &elen, expand);
    CCcheck_rval (rval, "CCcut_SRK_graph_edges failed");

    CC_MALLOC (perm, ecount, int);
    for (i = 0; i < ecount; i++) {
        perm[i] = i;
        elen[i] = -elen[i];
    }

    CCutil_double_perm_quicksort (perm, elen, ecount);

    CC_MALLOC (olen, ecount, double);
    CC_MALLOC (olist, 2*ecount, int);

    for (i = 0; i < ecount; i++) {
        olen[i] = -elen[perm[i]];
        olist[2*i]   = elist[2*perm[i]];
        olist[2*i+1] = elist[2*perm[i]+1];
    }

    *oncount = ncount;
    *oecount = ecount;
    *polist = olist;
    *polen = olen;

CLEANUP:
    CC_IFFREE (perm, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, double);
    return rval;
}

static int build_dp_cut_expand (CCtsp_lpcut_in **cut, int ndomino, int *Acount,
        int **A, int *Bcount, int **B, int handlecount, int *handle,
        CC_SRKexpinfo *expand, int ncount, int *comb)
{
    int rval = 0, i, valid, large_domino = 0;
    int *new_Acount = (int *) NULL, **new_A = (int **) NULL;
    int *new_Bcount = (int *) NULL, **new_B = (int **) NULL;
    int new_handlecount = 0, *new_handle = (int *) NULL;

    CC_MALLOC (new_Acount, ndomino, int);
    CC_MALLOC (new_Bcount, ndomino, int);
    CC_MALLOC (new_A, ndomino, int *);
    CC_MALLOC (new_B, ndomino, int *);

    for (i = 0; i < ndomino; i++) {
        rval = CCcut_SRK_expand (expand, A[i], Acount[i], &(new_A[i]),
                                 &(new_Acount[i]));
        CCcheck_rval (rval, "CCcut_SRK_expand failed");
        rval = CCcut_SRK_expand (expand, B[i], Bcount[i], &(new_B[i]),
                                 &(new_Bcount[i]));
        CCcheck_rval (rval, "CCcut_SRK_expand failed");
    }

    rval = CCcut_SRK_expand (expand, handle, handlecount, &new_handle,
                             &new_handlecount);
    CCcheck_rval (rval, "CCcut_SRK_expand failed");

    rval = check_raw_domino (ncount, new_handlecount, new_handle, ndomino, 
                             new_Acount, new_A, new_Bcount, new_B, &valid,
                             comb, &large_domino);
    CCcheck_rval (rval, "check_raw_domino failed");
    if (!valid) {
        fprintf (stderr, "Expanded domino is not valid\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCtsp_build_dp_cut (cut, ndomino, new_Acount, new_A, new_Bcount,
              new_B, new_handlecount, new_handle);
    CCcheck_rval (rval, "CCtsp_build_dp_cut failed");
   
CLEANUP:
    CC_IFFREE (new_Acount, int);
    CC_IFFREE (new_Bcount, int);
    if (new_A) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (new_A[i], int);
        }
        CC_IFFREE (new_A, int *);
    }
    if (new_B) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (new_B[i], int);
        }
        CC_IFFREE (new_B, int *);
    }
    CC_IFFREE (new_handle, int);

    return rval;
}

static void free_raw_dominos (int nIneq, int *nHandle, int **Handle, 
        int *nDominoes, int **nAset, int ***Aset, int **nBset, int ***Bset)
{
    int i, j;

    if (Aset) {
        for (i = 0; i < nIneq; i++) {
            for (j = 0; j < nDominoes[i]; j++) {
                CC_IFFREE (Aset[i][j], int);
            }
            CC_IFFREE (Aset[i], int *);
        }
        CC_IFFREE (Aset, int **);
    }
    if (Bset) {
        for (i = 0; i < nIneq; i++) {
            for (j = 0; j < nDominoes[i]; j++) {
                CC_IFFREE (Bset[i][j], int);
            }
            CC_IFFREE (Bset[i], int *);
        }
        CC_IFFREE (Bset, int **);
    }
    if (nAset) {
        for (i = 0; i < nIneq; i++) {
            CC_IFFREE (nAset[i], int);
        }
        CC_IFFREE (nAset, int *);
    }
    if (nBset) {
        for (i = 0; i < nIneq; i++) {
            CC_IFFREE (nBset[i], int);
        }
        CC_IFFREE (nBset, int *);
    }
    if (Handle) {
        for (i = 0; i < nIneq; i++) {
            CC_IFFREE (Handle[i], int); 
        }
        CC_IFFREE (Handle, int *);
    }
    CC_IFFREE (nHandle, int);
    CC_IFFREE (nDominoes, int);
}

int CCtsp_test_dp_cut (int ncount, CCtsp_lpcut_in *c, int *yes_no,
        int *large_domino)
{
    int rval = 0, i, handlecount, ndomino, chkvalid, chkcomb;
    int *Acount = (int *) NULL, **A = (int **) NULL;
    int *Bcount = (int *) NULL, **B = (int **) NULL;
    int *handle = (int *) NULL;

    if (yes_no) *yes_no = 0;
    if (large_domino) *large_domino = 0;

    if (c->dominocount <= 0 || c->cliquecount != 1) {
        printf ("DP-cut has no dominos or no handle\n"); fflush (stdout);
        if (yes_no) *yes_no = 0;
        goto CLEANUP;
    }

    if (c->rhs != CCtsp_DOMINORHS(c)) {
        printf ("DP-cut has wrong RHS\n"); fflush (stdout);
        if (yes_no) *yes_no = 0;
        goto CLEANUP;
    }

    rval = CCtsp_extract_dp_cut (c, &ndomino, &Acount, &A, &Bcount, &B,
                                 &handlecount, &handle);
    CCcheck_rval (rval, "CCtsp_extract_dp_cut failed");

    rval = check_raw_domino (ncount, handlecount, handle, ndomino, 
                  Acount, A, Bcount, B, &chkvalid, &chkcomb, large_domino);
    CCcheck_rval (rval, "check_raw_domino");

    if (chkvalid == 0) {
        printf ("Raw DP-cut is not valid\n"); fflush (stdout);
        if (yes_no) *yes_no = 0;
    } else {
        if (yes_no) *yes_no = 1;
    }

CLEANUP:
    CC_IFFREE (Acount, int);
    CC_IFFREE (Bcount, int);
    if (A) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (A[i], int);
        }
        CC_FREE (A, int *);
    }
    if (B) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (B[i], int);
        }
       CC_FREE (B, int *);
    }
    CC_IFFREE (handle, int);

    return rval;
}

static int check_raw_domino (int ncount, int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B, int *valid, int *comb,
        int *large_domino)
{
    int rval = 0, i, j, marker = 0, ain, aout, bin, bout, *marks = (int *) NULL;

    *valid = 0;
    *comb = 0;
    if (large_domino) *large_domino = 0;

    CC_MALLOC (marks, ncount, int);
    for (i = 0; i < ncount; i++) marks[i] = 0;

    if (dcount % 2 == 0) {
        fprintf (stderr, "Even number of dominos in a DP cut\n");
        goto CLEANUP;
    }

    if (hcount <= 0 || hcount >= ncount) {
        fprintf (stderr, "Handle has %d nodes in a DP cut\n", hcount);
        goto CLEANUP;
    }
    marker++;
    for (i = 0; i < hcount; i++) {
        if (hand[i] < 0 || hand[i] >= ncount) {
            fprintf (stderr, "Bad node in handle of DP cut\n");
            goto CLEANUP;
        }
        if (marks[hand[i]] == marker) {
            fprintf (stderr, "Repeated node in handle of DP cut\n");
            goto CLEANUP;
        }
        marks[hand[i]] = marker;
    }
    for (i = 0; i < dcount; i++) {
        if (Acount[i] + Bcount[i] >= ncount) {
            fprintf (stderr, "Domino covering all nodes in DP cut\n");
            goto CLEANUP;
        }
        if (Acount[i] <= 0 || Bcount[i] <= 0) {
            fprintf (stderr, "Domino with empty side in DP cut\n");
            goto CLEANUP;
        }
        if (Acount[i] + Bcount[i] >= (0.75 * ncount)) {
            printf ("Note: domino with %d nodes\n", Acount[i] + Bcount[i]);
            fflush (stdout);
            if (large_domino) *large_domino = 1;
        }
        
        marker++;
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) {
                fprintf (stderr, "Domino with repeated node in DP cut\n");
                goto CLEANUP;
            }
            marks[A[i][j]] = marker;
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) {
                fprintf (stderr, "Domino with repeated node in DP cut\n");
                goto CLEANUP;
            }
            marks[B[i][j]] = marker;
        }
    }
    *valid = 1;

    /* Test for comb */

    /* First check if dominos intersect handle correctly */

    marker++;
    for (i = 0; i < hcount; i++) marks[hand[i]] = marker;
    for (i = 0; i < dcount; i++) {
        ain = aout = 0;
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) ain++;
            else                          aout++;
        }
        if (ain && aout) goto CLEANUP;

        bin = bout = 0;
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) bin++;
            else                          bout++;
        }
        if (bin && bout) goto CLEANUP;
        if ((aout && bout) || (ain && bin)) goto CLEANUP;
    }

    /* Now check that teeth are disjoint */

    marker++;
    for (i = 0; i < dcount; i++) {
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) goto CLEANUP;
            marks[A[i][j]] = marker;
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) goto CLEANUP;
            marks[B[i][j]] = marker;
        }
    }

    *comb = 1;

CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

static void  print_raw_domino (int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B)
{
    int i, j;

    printf ("Raw Domino\n");
    printf ("    Handle:");
    for (i = 0; i < hcount; i++)  {
        printf (" %d", hand[i]);
        if (i % 15 == 14) printf ("\n           ");
    }
    printf ("\n");
    for (i = 0; i < dcount; i++) {
        printf ("    Domino %d:", i);
        for (j = 0; j < Acount[i]; j++) printf (" %d", A[i][j]); 
        printf ("  |  ");
        for (j = 0; j < Bcount[i]; j++) printf (" %d", B[i][j]); 
        printf ("\n");
    }
    fflush (stdout);
}

#endif
