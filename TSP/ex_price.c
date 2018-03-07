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
/*            ROUTINES TO PRICE EDGES USING BIGGUY ARITHMETIC               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: January 21, 1997; February 12, 2014 (bico)                        */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_exact_price (CCtsp_lp *lp, CCbigguy *bound,                   */
/*      int complete_price, int phase1, int eliminate, int have_exact,      */
/*      int silent)                                                         */
/*    RETURNS a bound that is valid for the entire edge set; the values     */
/*         of the dual variables will be stored in lp->exact_dual unless    */
/*         the existing lp->exact_dual's cutcount agrees with the           */
/*         cutcount for lp                                                  */
/*      -lp is a pointer to the tsp lp                                      */
/*      -bound returns the LP bound                                         */
/*      -if complete_price is nonzero, then price over the complete         */
/*       graph, even if a valid full edge set is present; if eliminating    */
/*       and complete_price is zero and lp->full_edges_valid, then          */
/*       the elimination is based on the lp->fulladj list.  Otherwise       */
/*       the elimination is based on the lp->dat.                           */
/*      -phase1 is either 0 or 1, with 1 indicating that the pricing        */
/*       should be to determine a Farkas' lemma bound to prove that the     */
/*       LP is infeasbile                                                   */
/*      -eliminate is 0 or 1; if 0 only price and if 1 elim edges;          */
/*       uses the bound information to elimination edges and set edges      */
/*       to 1; the remaining edges are placed into lp->fulladj (the old adj */
/*       is freed) and the fixed edges are placed on the list               */
/*       lp->fixededges; the dual values are taken from lp->exact_dual      */
/*      -have_exact is 0 or 1; set to 1 if eliminating and we already have  */
/*       a valid bound in lp->exact_lowerbound                              */
/*      -silent should be nonzero to supress output                         */
/*                                                                          */
/*  int CCtsp_exact_dual (CCtsp_lp *lp)                                     */
/*    RETURNS the dual values as bigguys (used to store the info used       */
/*        to establish the exact lower bound); the values will be           */
/*        stored in lp->exact_dual (and the existing values freed)          */
/*     -lp is the CCtsp_lp                                                  */
/*                                                                          */
/*  int CCtsp_verify_infeasible_lp (CCtsp_lp *lp, int *yesno, int silent)   */
/*    VERIFIES that the lp is infeasible using exact pricing.               */
/*     -yesno is set to 1 if the lp is infeasible and 0 otherwise.          */
/*                                                                          */
/*  int CCtsp_verify_lp_prune (CCtsp_lp *lp, int *yesno, int silent)        */
/*    VERIFIES that the lp bound is less than the lp upperbound - 1.        */
/*     -yesno is set to 1 if bound < upperbound - 1 and 0 otherwise.        */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "bigguy.h"
#include "tsp.h"

#define BIG_PRICE_GEN 5000000

typedef struct bigpredge {
    int ends[2];
    int len;
    CCbigguy rc;
} bigpredge;

typedef struct predge {
    int ends[2];
    int len;
} predge;

static int
    big_pricing_duals (CCtsp_lp *lp, CCbigguy *node_pi, CCbigguy *node_piest,
        CCbigguy *cut_pi, CCbigguy *clique_pi, CCbigguy *rhs_sum),
    big_price_list (CCtsp_lp *lp, int ecount, bigpredge *elist,
        CCbigguy *node_pi, CCbigguy *clique_pi, CCbigguy *cut_pi),
    big_generate_edges (CCtsp_lp *lp, int use_full_edges, CCbigguy *node_piest,
        int nwant, int *gencount, bigpredge *genlist, int *n1, int *n2,
        int *finished, CCbigguy cutoff, int phase1),
    add_to_inlist (CCtsp_lp *lp, int use_full_edges, bigpredge *inlist,
        int *count, int end0, int end1, int phase1),
    add_to_adj (CCtsp_lp *lp, CCtsp_genadj *adj, int end0, int end1,
        int *count),
    test_edge (int end1, int end2, int len, CCbigguy *node_piest,
        CCbigguy cutoff);


#define DIFFTOL 0.1
#define GUESS_MULT 50

int CCtsp_exact_price (CCtsp_lp *lp, CCbigguy *bound, int complete_price,
        int phase1, int eliminate, int have_exact, int silent)
{
    int rval = 0, incount, n1, n2, oldn1, i, j, ek, finished, nbranch = 0;
    int use_full_edges, problem_solved = 0, nfixed = 0, nremain = 0;
    int *newfix = (int *) NULL, newspace = 0;
    int end0, end1, oldnfixed = lp->nfixededges;
    bigpredge *inlist = (bigpredge *) NULL;
    predge *newlist = (predge *) NULL;
    CCbigguy penalty, rhs_sum, biglower, bigupper, cutoff, negcutoff;
    CCbigguy *node_pi = (CCbigguy *) NULL, *node_piest = (CCbigguy *) NULL;
    CCbigguy *clique_pi = (CCbigguy *) NULL, *cut_pi = (CCbigguy *) NULL;
    CCtsp_genadj *adj = (CCtsp_genadj *) NULL;
    CCtsp_genadjobj *adjspace = (CCtsp_genadjobj *) NULL, *pa;
    double szeit;

/*  NEGATIVE RC
    int negative_count = 0;
    CCbigguy negtest;
    FILE *out_rc = (FILE *) NULL;

    negtest = CCbigguy_dtobigguy(-0.1);
    out_rc = fopen ("out.rc", "w");
    if (!out_rc) {
        printf ("Could not open out.rc\n"); exit (1);
    }
*/

    szeit = CCutil_zeit ();
    *bound = CCbigguy_ZERO;

    if (!lp->dat && !lp->full_edges_valid) {
        fprintf (stderr, "must have dat file or full edge set\n");
        rval = 1;  goto CLEANUP;
    }

    if (!lp->dat && complete_price) {
        fprintf (stderr, "must have dat file for complete price\n");
        rval = 1;  goto CLEANUP;
    }

    if (phase1 && eliminate) {
        fprintf (stderr, "cannot eliminate variables in phase1\n");
        rval = 1;  goto CLEANUP;
    }

    if (phase1 && silent<2) { printf ("phase 1 pricing\n"); fflush (stdout); }

    if (complete_price) {
        printf ("Pricing COMPLETE GRAPH\n"); fflush (stdout);
        use_full_edges = 0;
    } else {
        use_full_edges = lp->full_edges_valid;
    }

    /* US50K: With branching elim, this can get called twice, so it can */
    /* end up with the cut counts the same, but with out-of date exact  */
    /* dual.  So always call CCtsp_exact_dual LP (with 1 test)          */

    if (1 || !lp->exact_dual || lp->exact_dual->cutcount != lp->cuts.cutcount) {
        rval = CCtsp_exact_dual (lp);
        CCcheck_rval (rval, "CCtsp_exact_dual failed");
    }

    if (eliminate) {
        if (lp->upperbound == CCtsp_LP_MAXDOUBLE) {
            printf ("Need upper bound for elimination\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }

        if (have_exact &&
            CCbigguy_cmp (lp->exact_lowerbound, CCbigguy_MINBIGGUY) == 0) {
            printf ("Exact lower bound not available\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }

        if (have_exact) biglower = lp->exact_lowerbound;
        else            biglower = CCbigguy_dtobigguy(lp->lowerbound-DIFFTOL);

        cutoff = CCbigguy_dtobigguy (lp->upperbound);
        CCbigguy_sub (&cutoff, biglower);
        CCbigguy_sub (&cutoff, CCbigguy_ONE);
        negcutoff = CCbigguy_ZERO;
        CCbigguy_sub (&negcutoff, cutoff);

        if (!silent) {
            printf ("Edge Elim Cutoff: %f\n", CCbigguy_bigguytod (cutoff));
            fflush (stdout);
        }
        if (CCbigguy_cmp (cutoff, CCbigguy_ZERO) < 0) {
            if (silent < 2) {
                printf ("Cutoff less than ZERO, do not no elim\n");
                fflush (stdout);
            }
            problem_solved = 1;
        } else {
            CC_MALLOC (newlist, GUESS_MULT*lp->graph.ncount, predge);
            newspace = GUESS_MULT*lp->graph.ncount;
            CC_MALLOC (newfix, 2*lp->graph.ncount, int);
            CC_MALLOC (adj, lp->graph.ncount, CCtsp_genadj);
            for (i = 0; i < lp->graph.ncount; i++) { adj[i].deg = 0; }
        }
    } else {
        cutoff = CCbigguy_ZERO;
    }

    for (i = 0; i < lp->branchdepth; i++) {
        if (lp->branchhistory[i].ends[0] != -1) { nbranch++; }
    }

    incount = 0;
    if (BIG_PRICE_GEN >= lp->nfixededges + nbranch) {
        CC_MALLOC (inlist, BIG_PRICE_GEN, bigpredge);
    } else {
        CC_MALLOC (inlist, lp->nfixededges + nbranch, bigpredge);
    }

    CC_MALLOC (node_pi, lp->graph.ncount, CCbigguy);
    CC_MALLOC (node_piest, lp->graph.ncount, CCbigguy);
  
    if (lp->cuts.cliqueend) {
        CC_MALLOC (clique_pi, lp->cuts.cliqueend, CCbigguy);
    }
    if (lp->cuts.cutcount) {
        CC_MALLOC (cut_pi, lp->cuts.cutcount, CCbigguy);
    }

    rval = big_pricing_duals (lp, node_pi, node_piest, cut_pi, clique_pi,
                              &rhs_sum);
    CCcheck_rval (rval, "big_pricing_duals failed");

    finished = 0;
    n1 = 0;
    n2 = (use_full_edges ? 0 : 1);
    penalty = CCbigguy_ZERO;
    oldn1 = n1;

    while (!finished) {
        rval = big_generate_edges (lp, use_full_edges, node_piest,
                 BIG_PRICE_GEN, &incount, inlist, &n1, &n2, &finished, cutoff,
                 phase1);
        CCcheck_rval (rval, "big_generate_edges failed");

        if (!silent && n1 > oldn1 + 10) {
            printf ("n1 = %d\n", n1); fflush (stdout); oldn1 = n1;
        }

        rval = big_price_list (lp, incount, inlist, node_pi, clique_pi, cut_pi);
        CCcheck_rval (rval, "big_price_list failed");

        for (i = 0; i < incount; i++) {
            if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) < 0) {
/* NEGATIVE RC - comment out the following line */
                CCbigguy_add (&penalty, inlist[i].rc);
            }
/* NEGATIVE RC
            if (CCbigguy_cmp (inlist[i].rc, negtest) < 0) {
                fprintf (out_rc, "%d %d %d\n", inlist[i].ends[0],
                                               inlist[i].ends[1],
                                               inlist[i].len);
                negative_count++;
            }
*/
        }

        if (eliminate && problem_solved == 0) {
            for (i = 0; i < incount; i++) {
                if (CCbigguy_cmp (inlist[i].rc, cutoff) <= 0) {
                    adj[inlist[i].ends[0]].deg++;
                    if (nremain >= newspace) {
                        void *tmp_ptr = (void *) newlist;
                        rval = CCutil_reallocrus_count (&tmp_ptr, 
                                             2*newspace, sizeof (predge));
                        CCcheck_rval (rval, "out of memory for newlist");
                        newlist = (predge *) tmp_ptr;
                        newspace = 2*newspace;
                    } 
                    newlist[nremain].ends[0] = inlist[i].ends[0];
                    newlist[nremain].ends[1] = inlist[i].ends[1];
                    newlist[nremain].len = inlist[i].len;
                    nremain++;
                    if (CCbigguy_cmp (inlist[i].rc, negcutoff) < 0) {
                        ek = CCtsp_find_edge (&(lp->graph),inlist[i].ends[0],
                                                           inlist[i].ends[1]);
                        if (ek != -1 && (lp->graph.edges[ek].fixed  == 0 &&
                                         lp->graph.edges[ek].branch == 0)) {
                            newfix[2*nfixed]   = inlist[i].ends[0];
                            newfix[2*nfixed+1] = inlist[i].ends[1];
                            nfixed++;
                        }
                    }
                }
            }
        }
    }

    if (lp->nfixededges + nbranch) {
        /* Adjust bound for fixed/branch edges                            */
        /* If a fixed or a branch=1 edge has positive rc, it can be added */
        /* If a branch=0 edge has negative rc, it should be subtracted    */
        /* from the penalty since it was earlier added.                   */

        incount = 0;
        for (i = 0; i < lp->nfixededges; i++) {
            end0 = lp->fixededges[2*i];
            end1 = lp->fixededges[2*i+1];
            rval = add_to_inlist (lp, use_full_edges, inlist, &incount,
                                  end0, end1, phase1);
            CCcheck_rval (rval, "add_to_inlist failed");
        }
        for (i = 0; i < lp->branchdepth; i++) {
            if (lp->branchhistory[i].ends[0] != -1) {
                end0 = lp->branchhistory[i].ends[0];
                end1 = lp->branchhistory[i].ends[1];
                rval = add_to_inlist (lp, use_full_edges, inlist, &incount,
                                      end0, end1, phase1);
                CCcheck_rval (rval, "add_to_inlist failed");
            }
        }

        rval = big_price_list (lp, incount, inlist, node_pi, clique_pi, cut_pi);
        CCcheck_rval (rval, "big_price_list failed");

        for (i = 0; i < lp->nfixededges; i++) {
            if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) > 0) {
                CCbigguy_add (&penalty, inlist[i].rc);
            }
        }
        for (i = lp->nfixededges, j = 0; i < lp->nfixededges+nbranch; i++) {
            while (lp->branchhistory[j].ends[0] == -1) j++;
            if (lp->branchhistory[j].rhs == 0) {
                if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) < 0) {
                    CCbigguy_sub (&penalty, inlist[i].rc);
                }
            } else {
                if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) > 0) {
                    CCbigguy_add (&penalty, inlist[i].rc);
                }
            }
            j++;
        }
    }

    *bound = rhs_sum;
    CCbigguy_add (bound, penalty);
    if (have_exact && CCbigguy_cmp (*bound, lp->exact_lowerbound) != 0) {
        printf ("Interesting: new bound does not match stored bound\n");
        fflush (stdout);
    }

    if (eliminate && !have_exact) {  
        lp->exact_lowerbound = *bound;
        if (CCbigguy_cmp (*bound, biglower) < 0) {
            if (problem_solved) {
                bigupper = CCbigguy_dtobigguy (lp->upperbound);
                CCbigguy_sub (&bigupper, CCbigguy_ONE);
                if (CCbigguy_cmp (*bound, bigupper) > 0) {
                    printf ("Exact lower bound allows pruning\n");
                    fflush (stdout);
                } else {
                    problem_solved = 0;
                }
            }
            if (!problem_solved) {
                printf ("Re-run Elimination with correct bound\n");
                fflush (stdout);
                CC_IFFREE (newlist, predge);
                CC_IFFREE (newfix, int);
                rval = CCtsp_exact_price (lp, bound, complete_price, phase1,
                                          eliminate, 1, silent);
                CCcheck_rval (rval, "CCtsp_exact_price failed");
                goto CLEANUP;
            }
        }
    }

    if (eliminate && problem_solved == 0) {
        /* leave enough room for the old fixed edges and the branch edges */
        for (i = 0; i < lp->branchdepth; i++) {
            if (lp->branchhistory[i].ends[0] != -1) {
                end0 = lp->branchhistory[i].ends[0];
                end1 = lp->branchhistory[i].ends[1];
                if (end0 < end1) adj[end0].deg++;
                else             adj[end1].deg++;
            }
        }
        for (i = 0; i < lp->nfixededges; i++) {
            end0 = lp->fixededges[2*i];
            end1 = lp->fixededges[2*i+1];
            if (end0 < end1) adj[end0].deg++;
            else             adj[end1].deg++;
        }

        if (nremain + lp->nfixededges + nbranch) {
            CC_MALLOC (adjspace, nremain+lp->nfixededges+nbranch,
                       CCtsp_genadjobj);
        }
        if (nfixed) {
            void *tmp_ptr = (void *) lp->fixededges;
            rval = CCutil_reallocrus_count (&tmp_ptr, 
                      2 * (lp->nfixededges + nfixed), sizeof (int));
            CCcheck_rval (rval, "out of memory for fixed edges");
            lp->fixededges = (int *) tmp_ptr;
        }

        pa = adjspace;
        for (i = 0; i < lp->graph.ncount; i++) {
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

        for (i = 0; i < nfixed; i++) {
            lp->fixededges[2*lp->nfixededges]   = newfix[2*i];
            lp->fixededges[2*lp->nfixededges+1] = newfix[2*i+1];
            lp->nfixededges++;
        }

        /* some old fixed and branched edges may not have gotten into adj */

        for (i = 0; i < lp->branchdepth; i++) {
            if (lp->branchhistory[i].ends[0] != -1) {
                end0 = lp->branchhistory[i].ends[0];
                end1 = lp->branchhistory[i].ends[1];
                rval = add_to_adj (lp, adj, end0, end1, &nremain);
                CCcheck_rval (rval, "add_to_adj failed");
            }
        }
        for (i = 0; i < oldnfixed; i++) {
            rval = add_to_adj (lp, adj, lp->fixededges[2*i],
                                        lp->fixededges[2*i+1], &nremain);
            CCcheck_rval (rval, "add_to_adj failed");
        }

        CC_IFFREE (lp->fulladjspace, CCtsp_genadjobj);
        CC_IFFREE (lp->fulladj, CCtsp_genadj);
        lp->fullcount = nremain;
        lp->fulladjspace = adjspace;
        lp->fulladj = adj;
        lp->full_edges_valid = 1;

        if (!silent) {
            printf ("Remaining Edges: %d (with %d new fixed)\n",nremain,nfixed);
            fflush (stdout);
        }
        rval = CCtsp_eliminate_rebuild (lp, silent);
        CCcheck_rval (rval, "CCtsp_eliminate_rebuild failed");
    }

/* NEGATIVE RC
    printf ("Negative Edges: %d\n", negative_count); fflush (stdout);
    if (out_rc) fclose (out_rc);
*/

    if (!silent) {
        printf ("Exact Price Time: %.2f seconds\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }
    rval = 0;

CLEANUP:
    CC_IFFREE (cut_pi, CCbigguy);
    CC_IFFREE (clique_pi, CCbigguy);
    CC_IFFREE (node_pi, CCbigguy);
    CC_IFFREE (node_piest, CCbigguy);
    CC_IFFREE (inlist, bigpredge);
    CC_IFFREE (newlist, predge);
    CC_IFFREE (newfix, int);
    if (rval) CC_IFFREE (adj, CCtsp_genadj);
    if (rval) CC_IFFREE (adjspace, CCtsp_genadjobj);
    return rval;
}

static int add_to_inlist (CCtsp_lp *lp, int use_full_edges, bigpredge *inlist,
        int *count, int end0, int end1, int phase1)
{
    int rval = 0, j, len = 0;
    CCtsp_genadj *adj = lp->fulladj;

    if (end0 > end1) { CC_SWAP (end0, end1, j); }
    if (!phase1) {
        if (use_full_edges) {
            for (j = 0; j < adj[end0].deg; j++) {
                if (adj[end0].list[j].end == end1) {
                    len = adj[end0].list[j].len; break;
                }
            }
            if (j == adj[end0].deg) {
                fprintf (stderr, "ERROR: edge not in fulladj\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            len = CCutil_dat_edgelen (end0, end1, lp->dat);
        }
    }
    inlist[*count].ends[0] = end0;
    inlist[*count].ends[1] = end1;
    inlist[*count].len     = len;
    (*count)++;

CLEANUP:
    return rval;
}

static int add_to_adj (CCtsp_lp *lp, CCtsp_genadj *adj, int end0, int end1,
        int *count)
{
    int rval = 0, len  = 0, j, k;

    if (end0 > end1) { CC_SWAP (end0, end1, k); }
    for (k = 0; k < adj[end0].deg; k++) {
        if (adj[end0].list[k].end == end1) break;
    }
    if (k == adj[end0].deg) {
        if (lp->full_edges_valid) {
            for (j = 0; j < lp->fulladj[end0].deg; j++) {
                if (lp->fulladj[end0].list[j].end == end1) {
                    len = lp->fulladj[end0].list[j].len; break;
                }
            }
            if (j == lp->fulladj[end0].deg) {
                fprintf (stderr, "ERROR: edge not in fulladj\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            len = CCutil_dat_edgelen (end0, end1, lp->dat);
        }
        adj[end0].list[adj[end0].deg].end = end1;
        adj[end0].list[adj[end0].deg].len = len;
        adj[end0].deg++;
        (*count)++;
    }

CLEANUP:
    return rval;
}

static int big_pricing_duals (CCtsp_lp *lp, CCbigguy *node_pi,
        CCbigguy *node_piest, CCbigguy *cut_pi, CCbigguy *clique_pi,
        CCbigguy *rhs_sum)
{
    int rval = 0, i, j, q, h, s, m, tmp;
    CCbigguy x;
    CCtsp_lpcut *c;

    *rhs_sum = CCbigguy_ZERO;

    if (!lp->exact_dual || lp->exact_dual->cutcount != lp->cuts.cutcount) {
        fprintf (stderr, "no exact_dual in big_pricing_duals\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < lp->graph.ncount; i++) {
        node_pi[i] = lp->exact_dual->node_pi[i];
    }
    for (i = 0; i < lp->cuts.cutcount; i++) {
        cut_pi[i] = lp->exact_dual->cut_pi[i];
    }

    /* Build RHS, adjusting for sparsification via degree equations */

    for (i = 0; i < lp->graph.ncount; i++) {
        CCbigguy_addmult (rhs_sum, node_pi[i], 2);
    }
    for (i = 0; i < lp->cuts.cutcount; i++) {
        x = cut_pi[i];
        CCbigguy_addmult (rhs_sum, x, lp->cuts.cuts[i].rhs);
        for (j = 0; j < lp->cuts.cuts[i].modcount; j++) {
            CCbigguy_addmult (rhs_sum, x,
                          2 * (lp->cuts.cuts[i].mods[j].mult - 128));
        }
    }

    /* Adjust node pi values for the sparsification mods */

    for (i = 0; i < lp->cuts.cutcount; i++) {
        x = cut_pi[i];
        for (j = 0; j < lp->cuts.cuts[i].modcount; j++) {
            CCbigguy_addmult (&(node_pi[lp->cuts.cuts[i].mods[j].node]), x,
                    lp->cuts.cuts[i].mods[j].mult - 128);
        }
    }


    /* Collect the clique pi values */

    for (i = 0; i < lp->cuts.cliqueend; i++) { clique_pi[i] = CCbigguy_ZERO; }
    for (i = 0; i < lp->cuts.cutcount; i++) {
        x = cut_pi[i];
        if (CCtsp_hypergraph_type (&lp->cuts.cuts[i])) {
            for (j = 0; j < lp->cuts.cuts[i].cliquecount; j++) {
                CCbigguy_add (&(clique_pi[lp->cuts.cuts[i].cliques[j]]), x);
            }
        }
    }

    /* Add to node pi values by running through cliques; piest will be used */
    /* to under-estimate the edge reduced costs                             */

    for (i = 0; i < lp->graph.ncount; i++) { node_piest[i] = node_pi[i]; }
    for (i = 0; i < lp->cuts.cliqueend; i++) {
        x = clique_pi[i];
        if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
            CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[i], tmp) {
                CCbigguy_add (&(node_pi[j]), x);
                CCbigguy_add (&(node_piest[j]), x);
            }
        } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
            CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[i], tmp) {
                CCbigguy_add (&(node_pi[j]), x);
            }
        }
    }

    /* Now run through remaining cuts to build node_piest.  The positive  */
    /* cliques get added for both sides of dominos, and one side of the   */
    /* positive semicuts.                                                 */

    /* For triomino cuts, let x = cut_pi[i].  Add to x to piest for each  */
    /* node in each handle, for both sides of dominos, for tsets, and     */
    /* both sides of each semi0 and semi1 (from triominos)                */

    for (i = 0; i < lp->cuts.cutcount; i++) {
        c = &lp->cuts.cuts[i];
        if (CCtsp_hypergraph_type (c) || 
            CCbigguy_cmp (cut_pi[i], CCbigguy_ZERO) == 0) continue; 

        if (CCbigguy_cmp (cut_pi[i], CCbigguy_ZERO) < 0) {
            fprintf (stderr, "YIPES: negative non-hypergraph cut\n");
            rval = 1;  goto CLEANUP;
        }

        m = 1;
        for (q = 0; q < c->cliquecount; q++) {
            if (c->cliquemult) m = c->cliquemult[q];
            x = CCbigguy_ZERO;
            CCbigguy_addmult (&x, cut_pi[i], m);
            if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                h = c->cliques[q];
                CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[h], tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
        for (q = 0; q < c->dominocount; q++) {
            x = cut_pi[i];
            h = c->dominos[q];
            for (s = 0; s < 2; s++) {
                CC_FOREACH_NODE_IN_CLIQUE (j,lp->cuts.dominos[h].sets[s],tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
        m = 1;
        for (q = 0; q < c->semicount; q++) {
            if (c->semimult) m = c->semimult[q];
            x = CCbigguy_ZERO;
            CCbigguy_addmult (&x, cut_pi[i], m);
            if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                h = c->semicuts[q];
                CC_FOREACH_NODE_IN_CLIQUE (j,lp->cuts.dominos[h].sets[0],tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
        if (c->coefcount) {
            printf ("Exact pricing not set up for general LHS coeficients\n");
            rval = 1; goto CLEANUP;
            /* For each node, compute half of the max coef from  */
            /* incidence edges and add to node_piest. Probably   */
            /* should avoid calls here with coefs and valid adj. */
        }

        if (c->TP_handles) {
            x = cut_pi[i];
            for (q = 0; q < 2; q++) {
                h = c->TP_handles[q];
                CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[h], tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
        for (q = 0; q < c->TP_domcount0; q++) {
            x = cut_pi[i];
            h = c->TP_dominos0[q];
            for (s = 0; s < 2; s++) {
                CC_FOREACH_NODE_IN_CLIQUE (j,lp->cuts.dominos[h].sets[s],tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
        for (q = 0; q < c->TP_domcount1; q++) {
            x = cut_pi[i];
            h = c->TP_dominos1[q];
            for (s = 0; s < 2; s++) {
                CC_FOREACH_NODE_IN_CLIQUE (j,lp->cuts.dominos[h].sets[s],tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
        for (q = 0; q < c->TP_tricount; q++) {
            x = cut_pi[i];
            h = c->TP_tsets[q];
            CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[h], tmp) {
                CCbigguy_add (&(node_piest[j]), x);
            }
            h = c->TP_semicuts0[q];
            for (s = 0; s < 2; s++) {
                CC_FOREACH_NODE_IN_CLIQUE (j,lp->cuts.dominos[h].sets[s],tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
            h = c->TP_semicuts1[q];
            for (s = 0; s < 2; s++) {
                CC_FOREACH_NODE_IN_CLIQUE (j,lp->cuts.dominos[h].sets[s],tmp) {
                    CCbigguy_add (&(node_piest[j]), x);
                }
            }
        }
    }

CLEANUP:
    return rval;
}

static int big_price_list (CCtsp_lp *lp, int ecount, bigpredge *elist,
    CCbigguy *node_pi, CCbigguy *clique_pi, CCbigguy *cut_pi)
{
    int rval = 0, i, j, tmp, l, nzlist, nznext, marker = 0;
    int ncount = lp->graph.ncount, ccount = lp->cuts.cliqueend;
    int *temp_elist = (int *) NULL;
    CCtsp_lpadj *adjspace = (CCtsp_lpadj *) NULL, *a;
    CCtsp_lpnode *n = (CCtsp_lpnode *) NULL;
    CCtsp_lpclique *c = lp->cuts.cliques;
    CCtsp_lpcut *cut;
    CCtsp_lpgraph g;
    CCbigguy x;

    CCtsp_init_lpgraph_struct (&g);
    if (ecount == 0) goto CLEANUP;

    CC_MALLOC (n, ncount, CCtsp_lpnode);
    for (i = 0; i < ncount; i++) { n[i].deg = 0; n[i].mark = 0; }
    for (i = 0; i < ecount; i++) {
        elist[i].rc = CCbigguy_itobigguy (elist[i].len);
        CCbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[0]]);
        CCbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[1]]);
        n[elist[i].ends[0]].deg++;
        n[elist[i].ends[1]].deg++;
    }

    CC_MALLOC (adjspace, 2*ecount, CCtsp_lpadj);
    a = adjspace;
    for (i = 0; i < ncount; i++) {
        n[i].adj = a;
        a += n[i].deg;
        n[i].deg = 0;
    }
    for (i = 0; i < ecount; i++) {
        j = elist[i].ends[0];
        n[j].adj[n[j].deg].to = elist[i].ends[1];
        n[j].adj[n[j].deg].edge = i;
        n[j].deg++;
        j = elist[i].ends[1];
        n[j].adj[n[j].deg].to = elist[i].ends[0];
        n[j].adj[n[j].deg].edge = i;
        n[j].deg++;
    }

    for (i = 0; i < ccount; i++) {
        if (CCbigguy_cmp (clique_pi[i], CCbigguy_ZERO)) {
            x = clique_pi[i];
            CCbigguy_add (&x, clique_pi[i]);
            marker++;
            CC_FOREACH_NODE_IN_CLIQUE (j, c[i], tmp) {
                a = n[j].adj;
                for (l = 0; l < n[j].deg; l++) {
                    if (n[a[l].to].mark == marker) {
                        CCbigguy_add (&(elist[a[l].edge].rc), x);
                    }
                }
                n[j].mark = marker;
            }
        }
    }

    /*  Hypergraph cuts, degree equations, and mods have been handled in  */
    /*  the rc field. Now run through remaining cuts, compute their LHS   */
    /*  using the general nzlist code and subtract the pi values from rc. */

    temp_elist = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (temp_elist, "out of memory in price_list");
    for (i = 0; i < ecount; i++) {
        temp_elist[2*i]   = elist[i].ends[0];
        temp_elist[2*i+1] = elist[i].ends[1];
    }
    rval = CCtsp_build_lpgraph (&g,ncount,ecount,temp_elist,(int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
    rval = CCtsp_build_lpadj (&g, 0, ecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    for (i = 0; i < lp->cuts.cutcount; i++) {
        cut = &lp->cuts.cuts[i];
        if (!CCtsp_hypergraph_type (cut)) {
            x = cut_pi[i];
            if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                rval = CCtsp_lpcut_nzlist (&g, cut, lp->cuts.cliques,
                                 lp->cuts.dominos, 0, &nzlist);
                while (nzlist != -1) {
                    nznext = g.edges[nzlist].coefnext;
                    g.edges[nzlist].coefnext = -2;
                    if (g.edges[nzlist].coef) {
                        CCbigguy_addmult (&(elist[nzlist].rc), x,
                                          -g.edges[nzlist].coef);
                        g.edges[nzlist].coef = 0;
                    }
                    nzlist = nznext;
                }
            } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
                fprintf (stderr, "YIPES: negative domino\n");
                rval = 1;  goto CLEANUP;
            }
        }
    }

CLEANUP:
    CC_IFFREE (n, CCtsp_lpnode);
    CC_IFFREE (adjspace, CCtsp_lpadj);
    CC_IFFREE (temp_elist, int);
    CCtsp_free_lpgraph (&g);
    return rval; 
}

static int big_generate_edges (CCtsp_lp *lp, int use_full_edges,
        CCbigguy *node_piest, int nwant, int *gencount, bigpredge *genlist,
        int *n1, int *n2, int *finished, CCbigguy cutoff, int phase1)
{
    int i = *n1, j = *n2, cnt = 0, len = 0, ncount = lp->graph.ncount;
    int end, stop, first = 1;
    CCtsp_genadj *adj = lp->fulladj;
    CCdatagroup *dat = lp->dat;

    *gencount = 0;
    *finished = 0;

    if (!lp->dat && !use_full_edges) {
        fprintf (stderr, "no source of edges in big_generate_edges\n");
        return 1;
    }
    if (i >= ncount) { *finished = 1; return 0; }

    for (; i < ncount; i++) {
        if (use_full_edges) {
            stop = adj[i].deg;
            if (first == 0) j = 0;
        } else {
            stop = ncount;
            if (first == 0) j = i+1;
        }
        first = 0;
        for (; j < stop; j++) {
            if (use_full_edges) {
                end = adj[i].list[j].end;
                if (!phase1) len = adj[i].list[j].len;
            } else {
                end = j;
                if (!phase1) len = CCutil_dat_edgelen (i, j, dat);
            }
            if (test_edge (i, end, len, node_piest, cutoff)) {
                genlist[cnt].ends[0] = i;
                genlist[cnt].ends[1] = end;
                genlist[cnt].len = len;
                cnt++;
                if (cnt == nwant) {
                    *gencount = cnt;
                    *n1 = i; *n2 = j + 1;
                    return 0;
                }
            }
        }
    }

    *finished = 1;
    *gencount = cnt;
    *n1 = ncount; *n2 = ncount;
    return 0;
}

static int test_edge (int end1, int end2, int len, CCbigguy *node_piest,
        CCbigguy cutoff)
{
    CCbigguy rc = CCbigguy_itobigguy (len);

    CCbigguy_sub (&rc, node_piest[end1]);
    CCbigguy_sub (&rc, node_piest[end2]);
    return (CCbigguy_cmp (rc, cutoff) <= 0);
}

int CCtsp_exact_dual (CCtsp_lp *lp)
{
    int rval = 0, i, ncount = lp->graph.ncount, cutcount = lp->cuts.cutcount;
    double *d_node_pi = (double *) NULL, *d_cut_pi = (double *) NULL;
    CCtsp_bigdual *d = (CCtsp_bigdual *) NULL;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
                    (int *) NULL, (int **) NULL, (double **) NULL,
                    (double **) NULL, &d_node_pi, &d_cut_pi);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

    CC_MALLOC (d, 1, CCtsp_bigdual);
    d->cutcount = cutcount;
    d->node_pi = (CCbigguy *) NULL;
    d->cut_pi = (CCbigguy *) NULL;

    CC_MALLOC (d->node_pi, ncount, CCbigguy);
    for (i = 0; i < ncount; i++) {
        d->node_pi[i] = CCbigguy_dtobigguy (d_node_pi[i]);
    }

    if (cutcount) {
        CC_MALLOC (d->cut_pi, cutcount, CCbigguy);
        for (i = 0; i < lp->cuts.cutcount; i++) {
            d->cut_pi[i] = CCbigguy_dtobigguy (d_cut_pi[i]);

            /* Bico October 26, 2004 Correct negative duals */
            if (lp->cuts.cuts[i].sense == 'G') {
                if (CCbigguy_cmp (d->cut_pi[i], CCbigguy_ZERO) < 0) {
                    printf ("CORRECTION: Setting G-dual %f value to 0\n",
                                          CCbigguy_bigguytod (d->cut_pi[i]));
                    fflush (stdout);
                    d->cut_pi[i] = CCbigguy_ZERO;
                }
            }

            if (lp->cuts.cuts[i].sense == 'L') {
                if (CCbigguy_cmp (d->cut_pi[i], CCbigguy_ZERO) > 0) {
                    printf ("CORRECTION: Setting L-dual %f value to 0\n",
                                          CCbigguy_bigguytod (d->cut_pi[i]));
                    fflush (stdout);
                    d->cut_pi[i] = CCbigguy_ZERO;
                }
            }
            /* END Bico October 26, 2004 */
            
        }
    }

    if (lp->exact_dual) { CCtsp_free_bigdual (&lp->exact_dual); }
    lp->exact_dual = d;

CLEANUP:
    CC_IFFREE (d_node_pi, double);
    CC_IFFREE (d_cut_pi, double);
    if (rval) { CCtsp_free_bigdual (&d); }
    return rval;
}

void CCtsp_free_bigdual (CCtsp_bigdual **d)
{
    if (d && *d) {
        CC_IFFREE ((*d)->node_pi, CCbigguy);
        CC_IFFREE ((*d)->cut_pi,  CCbigguy);
        CC_IFFREE ((*d), CCtsp_bigdual);
    }
}

int CCtsp_verify_infeasible_lp (CCtsp_lp *lp, int *yesno, int silent)
{
    int rval = 0;
    CCbigguy exactbound;

    *yesno = 0;
    rval = CCtsp_exact_price (lp, &exactbound, 0, 1, 0, 0, silent);
    CCcheck_rval (rval, "CCtsp_exact_price failed");

    if (!silent) {
        printf ("Exactbound: %f\n", CCbigguy_bigguytod (exactbound));
        fflush (stdout);
    }

    if (CCbigguy_cmp (exactbound, CCbigguy_ZERO) > 0) {
        if (silent < 2) {
            printf ("Problem is shown to be infeasible\n"); fflush (stdout);
        }
        *yesno = 1;
        lp->infeasible = 1;
        lp->lowerbound = CCtsp_LP_MAXDOUBLE;
        lp->exact_lowerbound = CCbigguy_MAXBIGGUY;
    } else {
        printf ("Did not verify an infeasible LP\n"); fflush (stdout);
        *yesno = 0;
    }

CLEANUP:
    return rval;
}

int CCtsp_verify_lp_prune (CCtsp_lp *lp, int *yesno, int silent)
{
    int rval = 0;
    CCbigguy exactbound, bnd;

    *yesno = 0;

    rval = CCtsp_exact_price (lp, &exactbound, 0, 0, 0, 0, silent);
    CCcheck_rval (rval, "CCtsp_exact_price failed");
  
    if (!silent) {
        printf ("Exact LP bound: %f\n", CCbigguy_bigguytod (exactbound));
        fflush (stdout);
    }

    bnd = CCbigguy_dtobigguy (lp->upperbound);
    CCbigguy_sub (&bnd, CCbigguy_ONE);

    if (CCbigguy_cmp (exactbound, bnd) > 0) {
        if (!silent) { printf ("Can prune lp.\n"); fflush (stdout); }
        *yesno = 1;
        lp->exact_lowerbound = exactbound;
    } else {
        if (!silent) { printf ("Cannot prune lp.\n"); fflush (stdout); }
        *yesno = 0;
    }

CLEANUP:
    return rval;
}
