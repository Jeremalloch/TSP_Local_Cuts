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
/*            PROOF ROUTINES TO PRICE EDGES USING BIGGUY ARITHMETIC         */
/*                  MODIFIED VERSION on TSP/ex_price.c                      */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 19, 2007                                                    */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int CCbbproof_elim ()                                                   */
/*  int CCbbproof_elim ()                                                   */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "bigguy.h"
#include "tsp.h"

#define BIG_PRICE_GEN 1000000

typedef struct bigpredge {
    int ends[2];
    int len;
    CCbigguy rc;
} bigpredge;

static int
     big_pricing_duals (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCbigguy **pnode_pi, CCbigguy **pnode_piest, CCbigguy **pnode_domino,
        CCbigguy **pcut_pi, CCbigguy **pclique_pi, CCbigguy *rhs_sum,
        int ncount),
    big_price_list (CCtsp_lpcuts *cuts, int ecount, bigpredge *elist,
        CCbigguy *node_pi, CCbigguy *clique_pi, CCbigguy *cut_pi, int ncount),
    test_edge (int end1, int end2, int len, CCbigguy *node_pi,
        CCbigguy *node_domino, CCbigguy cutoff),
    build_lpgraph (CCtsp_lpgraph *g, int ncount, int ecount, int *elist,
        int *elen),
    build_lpadj (CCtsp_lpgraph *g, int estart, int eend),
    domino_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut *c,
        CCtsp_lpclique *cliques, CCtsp_lpdomino *dominos, int *pnzlist);

static void
    big_generate_edges (CCbigguy *node_piest, CCbigguy *node_domino, int nwant,
        int *gencount, bigpredge *genlist, int *n1, int *n2, int *finished,
        CCbigguy cutoff, CCdatagroup *dat, int ncount),
    lpcut_nonzero_semicut (CCtsp_lpgraph *g, CCtsp_lpdomino *c, int *pnzlist),
    lpcut_nonzero_parity_handle (CCtsp_lpgraph *g, CCtsp_lpclique *c,
        int *pnzlist),
    lpcut_nonzero_domino (CCtsp_lpgraph *g, CCtsp_lpdomino *c, int *pnzlist),
    init_lpgraph_struct (CCtsp_lpgraph *g),
    free_lpgraph (CCtsp_lpgraph *g);

int CCbbproof_price (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCdatagroup *dat, int ncount, int *elist, int ecount, int *fixlist,
        int fcount, CCbigguy *bound, int silent);
int CCbbproof_elim (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual, 
        CCdatagroup *dat, double upperbound, CCbigguy exact_lowerbound,
        int ncount, int **elist, int *ecount, int **fixlist, int *fcount,
        int silent);

int CCbbproof_price (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCdatagroup *dat, int ncount, int *elist, int ecount, int *fixlist,
        int fcount, CCbigguy *bound, int silent)
{
    CCbigguy rhs_sum;
    CCbigguy *node_pi = (CCbigguy *) NULL;
    CCbigguy *node_piest = (CCbigguy *) NULL;
    CCbigguy *node_domino = (CCbigguy *) NULL;
    CCbigguy *clique_pi = (CCbigguy *) NULL;
    CCbigguy *cut_pi = (CCbigguy *) NULL;
    bigpredge *inlist = (bigpredge *) NULL;
    int i, n0, n1, len, incount;
    int rval = 0;

    *bound = CCbigguy_ZERO;

    if (!exact_dual || exact_dual->cutcount != cuts->cutcount) {
        fprintf (stderr, "no exact_dual in CCbbproof_price\n");
        rval = 1; goto CLEANUP;
    }

    rval = big_pricing_duals (cuts, exact_dual, &node_pi, &node_piest,
                &node_domino, &cut_pi, &clique_pi, &rhs_sum, ncount);
    CCcheck_rval (rval, "big_pricing_duals failed");

    incount = 0;
    inlist = CC_SAFE_MALLOC (ecount, bigpredge);
    CCcheck_NULL (inlist, "out of memory for inlist");

    for (i = 0; i < ecount; i++) {
        n0 = elist[2*i];  n1 = elist[2*i+1];
        len = CCutil_dat_edgelen (n0, n1, dat);
        if (test_edge (n0, n1, len, node_piest, node_domino, CCbigguy_ZERO)) {
            inlist[incount].ends[0] = n0;
            inlist[incount].ends[1] = n1;
            inlist[incount].len = len;
            incount++;
        }
    }

    rval = big_price_list (cuts, incount, inlist, node_pi, clique_pi, cut_pi,
                           ncount);
    CCcheck_rval (rval, "big_price_list failed");

    for (i = 0; i < incount; i++) {
        if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) < 0) {
            CCbigguy_add (&rhs_sum, inlist[i].rc);
        }
    }

    /* If a fixed edge has positive rc, it can be added to the bound */

    incount = 0;
    for (i = 0; i < fcount; i++) {
        n0 = fixlist[2*i];  n1 = fixlist[2*i+1];
        len = CCutil_dat_edgelen (n0, n1, dat);
        inlist[incount].ends[0] = n0;
        inlist[incount].ends[1] = n1;
        inlist[incount].len = len;
        incount++;
    }
    rval = big_price_list (cuts, incount, inlist, node_pi, clique_pi, cut_pi,
                           ncount);
    CCcheck_rval (rval, "big_price_list failed");

    for (i = 0; i < incount; i++) {
        if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) > 0) {
            CCbigguy_add (&rhs_sum, inlist[i].rc);
        }
    }

    *bound = rhs_sum;

    if (!silent) {
        printf ("Bound from bb_price: %f\n", CCbigguy_bigguytod (rhs_sum));
        fflush (stdout);
    }

CLEANUP:
    CC_IFFREE (cut_pi, CCbigguy);
    CC_IFFREE (clique_pi, CCbigguy);
    CC_IFFREE (node_pi, CCbigguy);
    CC_IFFREE (node_piest, CCbigguy);
    CC_IFFREE (node_domino, CCbigguy);
    CC_IFFREE (inlist, bigpredge);
    return rval;
}

int CCbbproof_elim (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual, 
        CCdatagroup *dat, double upperbound, CCbigguy exact_lowerbound,
        int ncount, int **elist, int *ecount, int **fixlist, int *fcount,
        int silent)
{
    int incount;
    bigpredge *inlist = (bigpredge *) NULL;
    CCbigguy rhs_sum;
    CCbigguy *node_pi = (CCbigguy *) NULL;
    CCbigguy *node_piest = (CCbigguy *) NULL;
    CCbigguy *node_domino = (CCbigguy *) NULL;
    CCbigguy *clique_pi = (CCbigguy *) NULL;
    CCbigguy *cut_pi = (CCbigguy *) NULL;
    int i, n1, n2, finished, nremain, nfixed;
    CCbigguy cutoff, negcutoff;
    double szeit;
    int *remain = (int *) NULL;
    int *fixed = (int *) NULL;
    int remainsupply = 0;
    int rval = 0;

    /* Do not allow cuts with modcount > 0 or twodomcount > 0 */

    if (CCbigguy_cmp (exact_lowerbound, CCbigguy_MINBIGGUY) == 0) {
        fprintf (stderr, "need an exact lowerbound to run elimination\n");
        rval = 1; goto CLEANUP;
    }
    if (!exact_dual || exact_dual->cutcount != cuts->cutcount) {
        fprintf (stderr, "no exact_dual in CCbbproof_elim\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Eliminating based on lower bound: %.12f\n",
                             CCbigguy_bigguytod (exact_lowerbound));
    fflush (stdout);

    szeit = CCutil_zeit ();

    cutoff = CCbigguy_dtobigguy (upperbound);
    CCbigguy_sub (&cutoff, exact_lowerbound);
    CCbigguy_sub (&cutoff, CCbigguy_ONE);
    negcutoff = CCbigguy_ZERO;
    CCbigguy_sub (&negcutoff, cutoff);
    if (!silent) {
        printf ("Edge Elimination Cutoff: %f\n", CCbigguy_bigguytod (cutoff));
        fflush (stdout);
    }
    if (CCbigguy_cmp (cutoff, CCbigguy_ZERO) < 0) {
        printf ("Cutoff is less than ZERO, do not eliminate\n");
        fflush (stdout);
        return 1;
    }

    remain = CC_SAFE_MALLOC (4*ncount, int); 
    CCcheck_NULL (remain, "out of memory for remain");
    remainsupply = 2*ncount;
    fixed = CC_SAFE_MALLOC (2*ncount, int); 
    CCcheck_NULL (fixed, "out of memory for fixed");

    incount = 0;
    inlist = CC_SAFE_MALLOC (BIG_PRICE_GEN, bigpredge);
    CCcheck_NULL (inlist, "out of memory for inlist");

    rval = big_pricing_duals (cuts, exact_dual, &node_pi, &node_piest,
                       &node_domino, &cut_pi, &clique_pi, &rhs_sum, ncount);
    CCcheck_rval (rval, "big_pricing_duals failed");

    finished = 0;
    nremain = 0;
    nfixed = 0;
    n1 = 0; n2 = 1;

    while (!finished) {
        big_generate_edges (node_piest, node_domino, BIG_PRICE_GEN,
                   &incount, inlist, &n1, &n2, &finished, cutoff, dat, ncount);
        rval = big_price_list (cuts, incount, inlist, node_pi, clique_pi,
                               cut_pi, ncount);
        CCcheck_rval (rval, "big_price_list failed");

        for (i = 0; i < incount; i++) {
            if (CCbigguy_cmp (inlist[i].rc, cutoff) <= 0) {
                if (nremain >= remainsupply) {
                    remain = realloc (remain, (4*remainsupply*sizeof (int)));
                    CCcheck_NULL (remain, "out of memory for remain realloc");
                    remainsupply = 2*remainsupply;
                }
                remain[2*nremain]   = inlist[i].ends[0];
                remain[2*nremain+1] = inlist[i].ends[1];
                nremain++;
                if (CCbigguy_cmp (inlist[i].rc, negcutoff) < 0) {
                    fixed[2*nfixed]   = inlist[i].ends[0];
                    fixed[2*nfixed+1] = inlist[i].ends[1];
                    nfixed++;
                }
            }
            if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) < 0) {
                CCbigguy_add (&rhs_sum, inlist[i].rc);
            }
        }
    }

    if (!silent) {
        printf ("Remaining Edges: %d (with %d fixed)\n", nremain, nfixed);
        printf ("Edge Elimination Time: %.2f seconds\n",
                CCutil_zeit () - szeit);
        fflush (stdout);
    }

    *elist   = remain;
    *ecount  = nremain;
    *fixlist = fixed;
    *fcount  = nfixed;

    printf ("Exactbound from bb_elim: %.12f\n", CCbigguy_bigguytod (rhs_sum));
    fflush (stdout);

/*
    {
        CC_SFILE *f = (CC_SFILE *) NULL;
        f = CCutil_sopen ("joke.bnd", "w");
        rval = CCbigguy_swrite (f, rhs_sum);
        CCcheck_rval (rval, "CCbigguy_swrite failed")
        CCutil_sclose (f);
    }
*/

    if (CCbigguy_cmp (rhs_sum, exact_lowerbound) < 0) {
        fprintf (stderr, "ERROR: Computed lower bound is less than estimate\n");
        fprintf (stderr, "DIFF: %f\n", CCbigguy_bigguytod (rhs_sum) -
                                       CCbigguy_bigguytod (exact_lowerbound));
        rval = 1;  goto CLEANUP;
    }


CLEANUP:

    CC_IFFREE (cut_pi, CCbigguy);
    CC_IFFREE (clique_pi, CCbigguy);
    CC_IFFREE (node_pi, CCbigguy);
    CC_IFFREE (node_piest, CCbigguy);
    CC_IFFREE (node_domino, CCbigguy);
    CC_IFFREE (inlist, bigpredge);
    if (rval) CC_IFFREE (remain, int);
    if (rval) CC_IFFREE (fixed, int);
    return rval;
}

static int big_pricing_duals (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCbigguy **pnode_pi, CCbigguy **pnode_piest, CCbigguy **pnode_domino,
        CCbigguy **pcut_pi, CCbigguy **pclique_pi, CCbigguy *rhs_sum,
        int ncount)
{
    CCbigguy x;
    int i, j, k, s, tmp;
    int rval = 0;
    CCtsp_lpcut *c;
    CCbigguy *node_pi = (CCbigguy *) NULL;
    CCbigguy *node_piest = (CCbigguy *) NULL;
    CCbigguy *node_domino = (CCbigguy *) NULL;
    CCbigguy *clique_pi = (CCbigguy *) NULL;
    CCbigguy *cut_pi = (CCbigguy *) NULL;

    *rhs_sum = CCbigguy_ZERO;

    node_pi    = CC_SAFE_MALLOC (ncount, CCbigguy);
    CCcheck_NULL (node_pi, "out of memory for node_pi");
    node_piest = CC_SAFE_MALLOC (ncount, CCbigguy);
    CCcheck_NULL (node_piest, "out of memory for node_piest");

    node_domino = CC_SAFE_MALLOC (ncount, CCbigguy);
    CCcheck_NULL (node_domino,  "out of memory for node_domino");

    if (cuts->cliqueend) {
        clique_pi = CC_SAFE_MALLOC (cuts->cliqueend, CCbigguy);
        CCcheck_NULL (clique_pi,  "out of memory for clique_pi");
    }
    if (cuts->cutcount) {
        cut_pi = CC_SAFE_MALLOC (cuts->cutcount, CCbigguy);
        CCcheck_NULL (cut_pi,  "out of memory for cut_pi");
    }

    for (i = 0; i < ncount; i++) {
        node_pi[i] = exact_dual->node_pi[i];
        CCbigguy_addmult (rhs_sum, node_pi[i], 2);
    }
    for (i = 0; i < cuts->cutcount; i++) {
        cut_pi[i] = exact_dual->cut_pi[i];
        CCbigguy_addmult (rhs_sum, cut_pi[i], cuts->cuts[i].rhs);
    }

    for (i = 0; i < cuts->cliqueend; i++) {
        clique_pi[i] = CCbigguy_ZERO;
    }
    for (i = 0; i < cuts->cutcount; i++) {
        x = cut_pi[i];
        if (cuts->cuts[i].dominocount <= 0) {
            for (j = 0; j < cuts->cuts[i].cliquecount; j++) {
                CCbigguy_add (&(clique_pi[cuts->cuts[i].cliques[j]]), x);
            }
        }
    }

    for (i = 0; i < ncount; i++) {
        node_piest[i] = node_pi[i];
    }

    for (i = 0; i < cuts->cliqueend; i++) {
        x = clique_pi[i];
        if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
            CC_FOREACH_NODE_IN_CLIQUE (j, cuts->cliques[i], tmp) {
                CCbigguy_add (&(node_pi[j]), x);
                CCbigguy_add (&(node_piest[j]), x);
            }
        } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
            CC_FOREACH_NODE_IN_CLIQUE (j, cuts->cliques[i], tmp) {
                CCbigguy_add (&(node_pi[j]), x);
            }
        }
    }

    for (i = 0; i < ncount; i++) {
        node_domino[i] = CCbigguy_ZERO;
    }
    for (i = 0; i < cuts->cutcount; i++) {
        c = &cuts->cuts[i];
        if (c->dominocount > 0) {
            x = cut_pi[i];
            if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                CC_FOREACH_NODE_IN_CLIQUE (j,
                     cuts->cliques[c->cliques[0]], tmp) {
                    CCbigguy_add (&(node_domino[j]), x);
                }
                for (k = 0; k < c->dominocount; k++) {
                    for (s = 0; s < 2; s++) {
                        CC_FOREACH_NODE_IN_CLIQUE (j,
                             cuts->dominos[c->dominos[k]].sets[s], tmp) {
                            CCbigguy_addmult (&(node_domino[j]), x, 2);
                        }
                    }
                }
            } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
                fprintf (stderr, "YIPES: negative domino\n");
                rval = 1;  goto CLEANUP;
            }
        }
    }

    *pnode_pi = node_pi;
    *pnode_piest = node_piest;
    *pnode_domino = node_domino;
    *pclique_pi = clique_pi;
    *pcut_pi = cut_pi;

CLEANUP:

    return rval;
}

static int big_price_list (CCtsp_lpcuts *cuts, int ecount, bigpredge *elist,
    CCbigguy *node_pi, CCbigguy *clique_pi, CCbigguy *cut_pi, int ncount)
{
    CCtsp_lpadj *adjspace = (CCtsp_lpadj *) NULL;
    CCtsp_lpnode *n = (CCtsp_lpnode *) NULL;
    int *temp_elist = (int *) NULL;
    int i, j, tmp, l, nzlist, nznext;
    CCtsp_lpadj *a;
    int marker = 0;
    CCbigguy x;
    int ccount = cuts->cliqueend;
    CCtsp_lpclique *c = cuts->cliques;
    CCtsp_lpcut *cut;
    CCtsp_lpgraph g;
    int rval = 0;

    init_lpgraph_struct (&g);

    if (ecount == 0) goto CLEANUP;

    n = CC_SAFE_MALLOC (ncount, CCtsp_lpnode);
    CCcheck_NULL (n, "out of memory in big_price_list");
    adjspace = CC_SAFE_MALLOC (2*ecount, CCtsp_lpadj);
    CCcheck_NULL (adjspace, "out of memory in big_price_list");

    for (i = 0; i < ncount; i++) {
        n[i].deg = 0;
        n[i].mark = 0;
    }
    for (i = 0; i < ecount; i++) {
        elist[i].rc = CCbigguy_itobigguy (elist[i].len);
        CCbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[0]]);
        CCbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[1]]);
        n[elist[i].ends[0]].deg++;
        n[elist[i].ends[1]].deg++;
    }
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

    /* Price dominoes, using nzlist */

    if (cuts->dominoend > 0) {
        temp_elist = CC_SAFE_MALLOC (2*ecount, int);
        CCcheck_NULL (temp_elist, "out of memory for temp_elist");
        for (i = 0; i < ecount; i++) {
            temp_elist[2*i]   = elist[i].ends[0];
            temp_elist[2*i+1] = elist[i].ends[1];
        }
        rval = build_lpgraph (&g, ncount, ecount, temp_elist, (int *) NULL);
        CCcheck_rval (rval, "build_lpgraph failed");
        rval = build_lpadj (&g, 0, ecount);
        CCcheck_rval (rval, "build_lpadj failed");
        CC_FREE (temp_elist, int);

        for (i = 0; i < cuts->cutcount; i++) {
            cut = &cuts->cuts[i];
            if (cut->dominocount > 0) {
                x = cut_pi[i];
                if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                    rval = domino_nzlist (&g, cut, cuts->cliques,
                                     cuts->dominos, &nzlist);
                    CCcheck_rval (rval, "domino_nzlist failed");
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
    }

CLEANUP:

    CC_IFFREE (n, CCtsp_lpnode);
    CC_IFFREE (adjspace, CCtsp_lpadj);
    CC_IFFREE (temp_elist, int);
    free_lpgraph (&g);

    return rval; 
}

static void big_generate_edges (CCbigguy *node_piest, CCbigguy *node_domino,
        int nwant, int *gencount, bigpredge *genlist, int *n1, int *n2,
        int *finished, CCbigguy cutoff, CCdatagroup *dat, int ncount)
{
    int i = *n1, j = *n2;
    int len, cnt = 0;

    *gencount = 0;
    *finished = 0;

    if (i >= ncount) { *finished = 1; goto CLEANUP; }

    for (; j < ncount; j++) {
        len = CCutil_dat_edgelen (i, j, dat);
        if (test_edge (i, j, len, node_piest, node_domino, cutoff)) {
            genlist[cnt].ends[0] = i;
            genlist[cnt].ends[1] = j;
            genlist[cnt].len = len;
            cnt++;
            if (cnt == nwant) {
                *finished = 0;
                *gencount = cnt;
                *n1 = i; *n2 = j + 1;
                goto CLEANUP;
            }
        }
    }
    for (i++; i < ncount; i++) {
        for (j = i + 1; j < ncount; j++) {
            len = CCutil_dat_edgelen (i, j, dat);
            if (test_edge (i, j, len, node_piest, node_domino,
                          cutoff)) {
                genlist[cnt].ends[0] = i;
                genlist[cnt].ends[1] = j;
                genlist[cnt].len = len;
                cnt++;
                if (cnt == nwant) {
                    *finished = 0;
                    *gencount = cnt;
                    *n1 = i; *n2 = j + 1;
                    goto CLEANUP;
                }
            }
        }
    }

    *n1 = ncount;
    *n2 = ncount;
    *gencount = cnt;
    *finished = 1;

CLEANUP:

    return;
}

static int test_edge (int end1, int end2, int len, CCbigguy *node_pi,
        CCbigguy *node_domino, CCbigguy cutoff)
{
    CCbigguy rc;

    rc = CCbigguy_itobigguy (len);
    CCbigguy_sub (&rc, node_pi[end1]);
    CCbigguy_sub (&rc, node_pi[end2]);
    if (node_domino) {
        CCbigguy_sub (&rc, node_domino[end1]);
        CCbigguy_sub (&rc, node_domino[end2]);
    }
    if (CCbigguy_cmp (rc, cutoff) <= 0) {
        return 1;
    } else {
        return 0;
    }
}


/*********************  from TSP/tsp_lp.c  **********************************/

static int build_lpgraph (CCtsp_lpgraph *g, int ncount, int ecount,
        int *elist, int *elen)
{
    int i;
    CCtsp_lpnode *n;
    CCtsp_lpedge *e;

    g->ncount = ncount;
    g->ecount = ecount;
    g->nodes = CC_SAFE_MALLOC (ncount, CCtsp_lpnode);
    if (!g->nodes) {
        return 1;
    }
    g->edges = CC_SAFE_MALLOC (ecount, CCtsp_lpedge);
    if (!g->edges) {
        CC_FREE (g->nodes, CCtsp_lpnode);
        return 1;
    }
    g->espace = ecount;
    n = g->nodes;
    e = g->edges;

    for (i = 0; i < ncount; i++) {
        n[i].mark = 0;
    }
    for (i=0; i<ecount; i++) {
        if (elist[2*i] < elist[2*i+1]) {
            e[i].ends[0] = elist[2*i];
            e[i].ends[1] = elist[2*i+1];
        } else {
            e[i].ends[0] = elist[2*i+1];
            e[i].ends[1] = elist[2*i];
        }
        e[i].fixed = 0;
        e[i].branch = 0;
        e[i].age = 0;
        if (elen) {
            e[i].len = elen[i];
        } else {
            e[i].len = 0;
        }
        e[i].coefnext = -2;
        e[i].coef = 0;
    }
    return 0;
}

static int build_lpadj (CCtsp_lpgraph *g, int estart, int eend)
{
    CCtsp_lpadj *a;
    CCtsp_lpnode *n = g->nodes;
    CCtsp_lpedge *e = g->edges;
    int i, j;

    if (g->adjspace) {
        if (g->adjstart == estart && g->adjend == eend) {
            return 0;
        } else {
            CC_FREE (g->adjspace, CCtsp_lpadj);
        }
    }

    if (estart >= eend) {
        g->adjstart = estart;
        g->adjend = eend;
        for (i=0; i<g->ncount; i++) {
            n[i].deg = 0;
            n[i].adj = (CCtsp_lpadj *) NULL;
        }
        return 0;
    }

    g->adjspace = CC_SAFE_MALLOC ((eend - estart)*2, CCtsp_lpadj);
    if (!g->adjspace) {
        return 1;
    }
    a = g->adjspace;
    for (i=0; i<g->ncount; i++) {
        n[i].deg = 0;
    }
    for (i=estart; i<eend; i++) {
        n[e[i].ends[0]].deg++;
        n[e[i].ends[1]].deg++;
    }
    for (i=0; i<g->ncount; i++) {
        n[i].adj = a;
        a += n[i].deg;
        n[i].deg = 0;
    }
    for (i=estart; i<eend; i++) {
        j = e[i].ends[0];
        a = &n[j].adj[n[j].deg];
        a->to = e[i].ends[1];
        a->edge = i;
        n[j].deg++;
        j = e[i].ends[1];
        a = &n[j].adj[n[j].deg];
        a->to = e[i].ends[0];
        a->edge = i;
        n[j].deg++;
    }
    g->adjstart = estart;
    g->adjend = eend;

    return 0;
}

static int domino_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut *c,
        CCtsp_lpclique *cliques, CCtsp_lpdomino *dominos, int *pnzlist)
{
    int i, nzlist = -1;
    int rval = 0;

    *pnzlist = -1;

    if (c->cliquecount != 1 || !c->dominos) {
        printf ("Cut is not a domino-parity inequality\n"); fflush (stdout);
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < c->dominocount; i++) {
        lpcut_nonzero_semicut (g, &dominos[c->dominos[i]], &nzlist);
    }
    lpcut_nonzero_parity_handle (g, &cliques[c->cliques[0]], &nzlist);
    for (i = 0; i < c->dominocount; i++) {
        lpcut_nonzero_domino (g, &dominos[c->dominos[i]], &nzlist);
    }
 
    *pnzlist = nzlist;

CLEANUP:

    return rval;
}

static void lpcut_nonzero_semicut (CCtsp_lpgraph *g, CCtsp_lpdomino *d,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges in the semi-cut (A:B) */

    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[0], tmp) {
        g->nodes[k].mark = nodemarker;
    }

    CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[1], tmp) {
        a = g->nodes[k].adj;
        for (l=0; l<g->nodes[k].deg; l++) {
            if (g->nodes[a[l].to].mark == nodemarker) {
                e = a[l].edge;
                if (g->edges[e].coefnext == -2) {
                    g->edges[e].coefnext = nzlist;
                    nzlist = e;
                }
                g->edges[e].coef++;
            }
        }
    }
    *pnzlist = nzlist;
}

static void lpcut_nonzero_parity_handle (CCtsp_lpgraph *g, CCtsp_lpclique *H,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges F, defined as the edges in an */
    /* odd number of semi-cuts (so an odd value in nzlist) and not in  */
    /* delta(H), or in an even number of semi-cuts and in delta(H).    */

    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    CC_FOREACH_NODE_IN_CLIQUE (k, *H, tmp) {
        g->nodes[k].mark = nodemarker;
    }

    for (e = nzlist; e != -1; e = g->edges[e].coefnext) {
        if (g->edges[e].coef % 2 == 1) {
            if ((g->nodes[g->edges[e].ends[0]].mark != nodemarker &&
                 g->nodes[g->edges[e].ends[1]].mark != nodemarker) || 
                (g->nodes[g->edges[e].ends[0]].mark == nodemarker &&
                 g->nodes[g->edges[e].ends[1]].mark == nodemarker)) {
                    g->edges[e].coef++;
            }
        }
    }

    CC_FOREACH_NODE_IN_CLIQUE (k, *H, tmp) {
        a = g->nodes[k].adj;
        for (l=0; l<g->nodes[k].deg; l++) {
            if (g->nodes[a[l].to].mark != nodemarker) {
                e = a[l].edge;
                if (g->edges[e].coefnext == -2) {
                    g->edges[e].coefnext = nzlist;
                    nzlist = e;
                    g->edges[e].coef++;
                } else if (g->edges[e].coef % 2 == 0) {
                    g->edges[e].coef++;
                }
            }
        }
    }
    *pnzlist = nzlist;
}

static void lpcut_nonzero_domino (CCtsp_lpgraph *g, CCtsp_lpdomino *d,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges in the cut delta(A union B) */

    int nzlist = *pnzlist;
    int nodemarker;
    CCtsp_lpadj *a;
    int e;
    int tmp, i, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    for (i = 0; i < 2; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[i], tmp) {
            g->nodes[k].mark = nodemarker;
        }
    }

    for (i = 0; i < 2; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (k, d->sets[i], tmp) {
            a = g->nodes[k].adj;
            for (l=0; l<g->nodes[k].deg; l++) {
                if (g->nodes[a[l].to].mark != nodemarker) {
                    e = a[l].edge;
                    if (g->edges[e].coefnext == -2) {
                        g->edges[e].coefnext = nzlist;
                        nzlist = e;
                    }
                    g->edges[e].coef++;
                }
            }
        }
    }
    *pnzlist = nzlist;
}

static void init_lpgraph_struct (CCtsp_lpgraph *g)
{
    g->ncount = 0;
    g->ecount = 0;
    g->nodes = (CCtsp_lpnode *) NULL;
    g->edges = (CCtsp_lpedge *) NULL;
    g->adjspace = (CCtsp_lpadj *) NULL;
    g->adjstart = 0;
    g->adjend = 0;
    g->nodemarker = 0;
    g->espace = 0;
}

static void free_lpgraph (CCtsp_lpgraph *g)
{
    CC_IFFREE (g->nodes, CCtsp_lpnode);
    CC_IFFREE (g->edges, CCtsp_lpedge);
    CC_IFFREE (g->adjspace, CCtsp_lpadj);
    g->espace = 0;
}

