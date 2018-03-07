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
/*  Price edges using bigguy arithmetic (modifed version of ex_price.c)     */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, Cook                            */
/*  Date: April 19, 2007                                                    */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int BBprice_price ()                                                    */
/*  int BBprice_elim ()                                                     */
/*                                                                          */
/****************************************************************************/

#include "bbproof.h"

#define BIG_PRICE_GEN 1000000

typedef struct bigpredge {
    int ends[2];
    int len;
    BBbigguy rc;
} bigpredge;

static int
    big_pricing_duals (BBtsp_lpcuts *cuts, BBtsp_bigdual *exact_dual,
        BBbigguy **pnode_pi, BBbigguy **pnode_piest, BBbigguy **pnode_domino,
        BBbigguy **pcut_pi, BBbigguy **pclique_pi, BBbigguy *rhs_sum,
        int ncount),
    big_price_list (BBtsp_lpcuts *cuts, int ecount, bigpredge *elist,
        BBbigguy *node_pi, BBbigguy *clique_pi, BBbigguy *cut_pi, int ncount),
    test_edge (int end1, int end2, int len, BBbigguy *node_pi,
        BBbigguy *node_domino, BBbigguy cutoff),
    build_lpgraph (BBtsp_lpgraph *g, int ncount, int ecount, int *elist,
        int *elen),
    build_lpadj (BBtsp_lpgraph *g, int estart, int eend),
    domino_nzlist (BBtsp_lpgraph *g, BBtsp_lpcut *c,
        BBtsp_lpclique *cliques, BBtsp_lpdomino *dominos, int *pnzlist);

static void
    big_generate_edges (BBbigguy *node_piest, BBbigguy *node_domino, int nwant,
        int *gencount, bigpredge *genlist, int *n1, int *n2, int *finished,
        BBbigguy cutoff, BBdatagroup *dat, int ncount),
    lpcut_nonzero_semicut (BBtsp_lpgraph *g, BBtsp_lpdomino *c, int *pnzlist),
    lpcut_nonzero_parity_handle (BBtsp_lpgraph *g, BBtsp_lpclique *c,
        int *pnzlist),
    lpcut_nonzero_domino (BBtsp_lpgraph *g, BBtsp_lpdomino *c, int *pnzlist),
    init_lpgraph_struct (BBtsp_lpgraph *g),
    free_lpgraph (BBtsp_lpgraph *g);

int BBprice_price (BBtsp_lpcuts *cuts, BBtsp_bigdual *exact_dual,
        BBdatagroup *dat, int ncount, int *elist, int ecount, int *fixlist,
        int fcount, BBbigguy *bound, int silent)
{
    int i, n0, n1, len, incount, rval = 0;
    BBbigguy rhs_sum;
    BBbigguy *node_pi = (BBbigguy *) NULL;
    BBbigguy *node_piest = (BBbigguy *) NULL;
    BBbigguy *node_domino = (BBbigguy *) NULL;
    BBbigguy *clique_pi = (BBbigguy *) NULL;
    BBbigguy *cut_pi = (BBbigguy *) NULL;
    bigpredge *inlist = (bigpredge *) NULL;

    *bound = BBbigguy_ZERO;

    if (!exact_dual || exact_dual->cutcount != cuts->cutcount) {
        fprintf (stderr, "no exact_dual in BBprice_price\n");
        rval = 1; goto CLEANUP;
    }

    rval = big_pricing_duals (cuts, exact_dual, &node_pi, &node_piest,
                &node_domino, &cut_pi, &clique_pi, &rhs_sum, ncount);
    BBcheck_rval (rval, "big_pricing_duals failed");

    incount = 0;
    inlist = BB_SAFE_MALLOC (ecount, bigpredge);
    BBcheck_NULL (inlist, "out of memory for inlist");

    for (i = 0; i < ecount; i++) {
        n0 = elist[2*i];  n1 = elist[2*i+1];
        len = BButil_dat_edgelen (n0, n1, dat);
        if (test_edge (n0, n1, len, node_piest, node_domino, BBbigguy_ZERO)) {
            inlist[incount].ends[0] = n0;
            inlist[incount].ends[1] = n1;
            inlist[incount].len = len;
            incount++;
        }
    }

    rval = big_price_list (cuts, incount, inlist, node_pi, clique_pi, cut_pi,
                           ncount);
    BBcheck_rval (rval, "big_price_list failed");

    for (i = 0; i < incount; i++) {
        if (BBbigguy_cmp (inlist[i].rc, BBbigguy_ZERO) < 0) {
            BBbigguy_add (&rhs_sum, inlist[i].rc);
        }
    }

    /* If a fixed edge has positive rc, it can be added to the bound */

    incount = 0;
    for (i = 0; i < fcount; i++) {
        n0 = fixlist[2*i];  n1 = fixlist[2*i+1];
        len = BButil_dat_edgelen (n0, n1, dat);
        inlist[incount].ends[0] = n0;
        inlist[incount].ends[1] = n1;
        inlist[incount].len = len;
        incount++;
    }
    rval = big_price_list (cuts, incount, inlist, node_pi, clique_pi, cut_pi,
                           ncount);
    BBcheck_rval (rval, "big_price_list failed");

    for (i = 0; i < incount; i++) {
        if (BBbigguy_cmp (inlist[i].rc, BBbigguy_ZERO) > 0) {
            BBbigguy_add (&rhs_sum, inlist[i].rc);
        }
    }
    *bound = rhs_sum;

    if (!silent) {
        printf ("Bound from bb_price: %f\n", BBbigguy_bigguytod (rhs_sum));
        fflush (stdout);
    }

CLEANUP:
    BB_IFFREE (cut_pi, BBbigguy);
    BB_IFFREE (clique_pi, BBbigguy);
    BB_IFFREE (node_pi, BBbigguy);
    BB_IFFREE (node_piest, BBbigguy);
    BB_IFFREE (node_domino, BBbigguy);
    BB_IFFREE (inlist, bigpredge);
    return rval;
}

int BBprice_elim (BBtsp_lpcuts *cuts, BBtsp_bigdual *exact_dual, 
        BBdatagroup *dat, double upperbound, BBbigguy exact_lowerbound,
        int ncount, int **elist, int *ecount, int **fixlist, int *fcount,
        int silent)
{
    int remainsupply = 0, rval = 0;
    int i, n1, n2, finished, nremain, nfixed, incount;
    int *remain = (int *) NULL, *fixed = (int *) NULL;
    double szeit;
    bigpredge *inlist = (bigpredge *) NULL;
    BBbigguy rhs_sum;
    BBbigguy *node_pi = (BBbigguy *) NULL;
    BBbigguy *node_piest = (BBbigguy *) NULL;
    BBbigguy *node_domino = (BBbigguy *) NULL;
    BBbigguy *clique_pi = (BBbigguy *) NULL;
    BBbigguy *cut_pi = (BBbigguy *) NULL;
    BBbigguy cutoff, negcutoff;

    if (BBbigguy_cmp (exact_lowerbound, BBbigguy_MINBIGGUY) == 0) {
        fprintf (stderr, "need an exact lowerbound to run elimination\n");
        rval = 1; goto CLEANUP;
    }
    if (!exact_dual || exact_dual->cutcount != cuts->cutcount) {
        fprintf (stderr, "no exact_dual in BBprice_elim\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Eliminating based on lower bound: %.12f\n",
                             BBbigguy_bigguytod (exact_lowerbound));
    fflush (stdout);

    szeit = BButil_zeit ();

    cutoff = BBbigguy_dtobigguy (upperbound);
    BBbigguy_sub (&cutoff, exact_lowerbound);
    BBbigguy_sub (&cutoff, BBbigguy_ONE);
    negcutoff = BBbigguy_ZERO;
    BBbigguy_sub (&negcutoff, cutoff);
    if (!silent) {
        printf ("Edge Elimination Cutoff: %f\n", BBbigguy_bigguytod (cutoff));
        fflush (stdout);
    }
    if (BBbigguy_cmp (cutoff, BBbigguy_ZERO) < 0) {
        printf ("Cutoff is less than ZERO, do not eliminate\n");
        fflush (stdout);
        return 1;
    }

    remain = BB_SAFE_MALLOC (4*ncount, int); 
    BBcheck_NULL (remain, "out of memory for remain");
    remainsupply = 2*ncount;
    fixed = BB_SAFE_MALLOC (2*ncount, int); 
    BBcheck_NULL (fixed, "out of memory for fixed");

    incount = 0;
    inlist = BB_SAFE_MALLOC (BIG_PRICE_GEN, bigpredge);
    BBcheck_NULL (inlist, "out of memory for inlist");

    rval = big_pricing_duals (cuts, exact_dual, &node_pi, &node_piest,
                       &node_domino, &cut_pi, &clique_pi, &rhs_sum, ncount);
    BBcheck_rval (rval, "big_pricing_duals failed");

    finished = 0;
    nremain = 0;
    nfixed = 0;
    n1 = 0; n2 = 1;

    while (!finished) {
        big_generate_edges (node_piest, node_domino, BIG_PRICE_GEN,
                   &incount, inlist, &n1, &n2, &finished, cutoff, dat, ncount);
        rval = big_price_list (cuts, incount, inlist, node_pi, clique_pi,
                               cut_pi, ncount);
        BBcheck_rval (rval, "big_price_list failed");

        for (i = 0; i < incount; i++) {
            if (BBbigguy_cmp (inlist[i].rc, cutoff) <= 0) {
                if (nremain >= remainsupply) {
                    remain = realloc (remain, (4*remainsupply*sizeof (int)));
                    BBcheck_NULL (remain, "out of memory for remain realloc");
                    remainsupply = 2*remainsupply;
                }
                remain[2*nremain]   = inlist[i].ends[0];
                remain[2*nremain+1] = inlist[i].ends[1];
                nremain++;
                if (BBbigguy_cmp (inlist[i].rc, negcutoff) < 0) {
                    fixed[2*nfixed]   = inlist[i].ends[0];
                    fixed[2*nfixed+1] = inlist[i].ends[1];
                    nfixed++;
                }
            }
            if (BBbigguy_cmp (inlist[i].rc, BBbigguy_ZERO) < 0) {
                BBbigguy_add (&rhs_sum, inlist[i].rc);
            }
        }
    }

    if (!silent) {
        printf ("Remaining Edges: %d (with %d fixed)\n", nremain, nfixed);
        printf ("Edge Elimination Time: %.2f seconds\n",
                BButil_zeit () - szeit);
        fflush (stdout);
    }

    *elist   = remain;
    *ecount  = nremain;
    *fixlist = fixed;
    *fcount  = nfixed;

    printf ("Exactbound from bb_elim: %.12f\n", BBbigguy_bigguytod (rhs_sum));
    fflush (stdout);

    if (BBbigguy_cmp (rhs_sum, exact_lowerbound) < 0) {
        fprintf (stderr, "ERROR: Computed lower bound is less than estimate\n");
        fprintf (stderr, "DIFF: %f\n", BBbigguy_bigguytod (rhs_sum) -
                                       BBbigguy_bigguytod (exact_lowerbound));
        rval = 1;  goto CLEANUP;
    }

CLEANUP:
    BB_IFFREE (cut_pi, BBbigguy);
    BB_IFFREE (clique_pi, BBbigguy);
    BB_IFFREE (node_pi, BBbigguy);
    BB_IFFREE (node_piest, BBbigguy);
    BB_IFFREE (node_domino, BBbigguy);
    BB_IFFREE (inlist, bigpredge);
    if (rval) BB_IFFREE (remain, int);
    if (rval) BB_IFFREE (fixed, int);
    return rval;
}

static int big_pricing_duals (BBtsp_lpcuts *cuts, BBtsp_bigdual *exact_dual,
        BBbigguy **pnode_pi, BBbigguy **pnode_piest, BBbigguy **pnode_domino,
        BBbigguy **pcut_pi, BBbigguy **pclique_pi, BBbigguy *rhs_sum,
        int ncount)
{
    int i, j, k, s, tmp, rval = 0;
    BBbigguy x;
    BBtsp_lpcut *c;
    BBbigguy *node_pi = (BBbigguy *) NULL;
    BBbigguy *node_piest = (BBbigguy *) NULL;
    BBbigguy *node_domino = (BBbigguy *) NULL;
    BBbigguy *clique_pi = (BBbigguy *) NULL;
    BBbigguy *cut_pi = (BBbigguy *) NULL;

    *rhs_sum = BBbigguy_ZERO;
    node_pi    = BB_SAFE_MALLOC (ncount, BBbigguy);
    BBcheck_NULL (node_pi, "out of memory for node_pi");
    node_piest = BB_SAFE_MALLOC (ncount, BBbigguy);
    BBcheck_NULL (node_piest, "out of memory for node_piest");
    node_domino = BB_SAFE_MALLOC (ncount, BBbigguy);
    BBcheck_NULL (node_domino,  "out of memory for node_domino");

    if (cuts->cliqueend) {
        clique_pi = BB_SAFE_MALLOC (cuts->cliqueend, BBbigguy);
        BBcheck_NULL (clique_pi,  "out of memory for clique_pi");
    }
    if (cuts->cutcount) {
        cut_pi = BB_SAFE_MALLOC (cuts->cutcount, BBbigguy);
        BBcheck_NULL (cut_pi,  "out of memory for cut_pi");
    }

    for (i = 0; i < ncount; i++) {
        node_pi[i] = exact_dual->node_pi[i];
        BBbigguy_addmult (rhs_sum, node_pi[i], 2);
    }
    for (i = 0; i < cuts->cutcount; i++) {
        cut_pi[i] = exact_dual->cut_pi[i];
        BBbigguy_addmult (rhs_sum, cut_pi[i], cuts->cuts[i].rhs);
    }

    for (i = 0; i < cuts->cliqueend; i++) {
        clique_pi[i] = BBbigguy_ZERO;
    }
    for (i = 0; i < cuts->cutcount; i++) {
        x = cut_pi[i];
        if (cuts->cuts[i].dominocount <= 0) {
            for (j = 0; j < cuts->cuts[i].cliquecount; j++) {
                BBbigguy_add (&(clique_pi[cuts->cuts[i].cliques[j]]), x);
            }
        }
    }

    for (i = 0; i < ncount; i++) {
        node_piest[i] = node_pi[i];
    }

    for (i = 0; i < cuts->cliqueend; i++) {
        x = clique_pi[i];
        if (BBbigguy_cmp (x, BBbigguy_ZERO) > 0) {
            BB_FOREACH_NODE_IN_CLIQUE (j, cuts->cliques[i], tmp) {
                BBbigguy_add (&(node_pi[j]), x);
                BBbigguy_add (&(node_piest[j]), x);
            }
        } else if (BBbigguy_cmp (x, BBbigguy_ZERO) < 0) {
            BB_FOREACH_NODE_IN_CLIQUE (j, cuts->cliques[i], tmp) {
                BBbigguy_add (&(node_pi[j]), x);
            }
        }
    }

    for (i = 0; i < ncount; i++) {
        node_domino[i] = BBbigguy_ZERO;
    }
    for (i = 0; i < cuts->cutcount; i++) {
        c = &cuts->cuts[i];
        if (c->dominocount > 0) {
            x = cut_pi[i];
            if (BBbigguy_cmp (x, BBbigguy_ZERO) > 0) {
                BB_FOREACH_NODE_IN_CLIQUE (j,
                     cuts->cliques[c->cliques[0]], tmp) {
                    BBbigguy_add (&(node_domino[j]), x);
                }
                for (k = 0; k < c->dominocount; k++) {
                    for (s = 0; s < 2; s++) {
                        BB_FOREACH_NODE_IN_CLIQUE (j,
                             cuts->dominos[c->dominos[k]].sets[s], tmp) {
                            BBbigguy_addmult (&(node_domino[j]), x, 1);  /*2*/
                        }
                    }
                }
            } else if (BBbigguy_cmp (x, BBbigguy_ZERO) < 0) {
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

static int big_price_list (BBtsp_lpcuts *cuts, int ecount, bigpredge *elist,
    BBbigguy *node_pi, BBbigguy *clique_pi, BBbigguy *cut_pi, int ncount)
{
    int i, j, tmp, l, nzlist, nznext, marker = 0, rval = 0;
    int ccount = cuts->cliqueend;
    int *temp_elist = (int *) NULL;
    BBbigguy x;
    BBtsp_lpadj *adjspace = (BBtsp_lpadj *) NULL;
    BBtsp_lpnode *n = (BBtsp_lpnode *) NULL;
    BBtsp_lpadj *a;
    BBtsp_lpclique *c = cuts->cliques;
    BBtsp_lpcut *cut;
    BBtsp_lpgraph g;

    init_lpgraph_struct (&g);
    if (ecount == 0) goto CLEANUP;

    n = BB_SAFE_MALLOC (ncount, BBtsp_lpnode);
    BBcheck_NULL (n, "out of memory in big_price_list");
    adjspace = BB_SAFE_MALLOC (2*ecount, BBtsp_lpadj);
    BBcheck_NULL (adjspace, "out of memory in big_price_list");

    for (i = 0; i < ncount; i++) {
        n[i].deg = 0;
        n[i].mark = 0;
    }
    for (i = 0; i < ecount; i++) {
        elist[i].rc = BBbigguy_itobigguy (elist[i].len);
        BBbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[0]]);
        BBbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[1]]);
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
        if (BBbigguy_cmp (clique_pi[i], BBbigguy_ZERO)) {
            x = clique_pi[i];
            BBbigguy_add (&x, clique_pi[i]);
            marker++;
            BB_FOREACH_NODE_IN_CLIQUE (j, c[i], tmp) {
                a = n[j].adj;
                for (l = 0; l < n[j].deg; l++) {
                    if (n[a[l].to].mark == marker) {
                        BBbigguy_add (&(elist[a[l].edge].rc), x);
                    }
                }
                n[j].mark = marker;
            }
        }
    }

    /* Price dominoes, using nzlist */

    if (cuts->dominoend > 0) {
        temp_elist = BB_SAFE_MALLOC (2*ecount, int);
        BBcheck_NULL (temp_elist, "out of memory for temp_elist");
        for (i = 0; i < ecount; i++) {
            temp_elist[2*i]   = elist[i].ends[0];
            temp_elist[2*i+1] = elist[i].ends[1];
        }
        rval = build_lpgraph (&g, ncount, ecount, temp_elist, (int *) NULL);
        BBcheck_rval (rval, "build_lpgraph failed");
        rval = build_lpadj (&g, 0, ecount);
        BBcheck_rval (rval, "build_lpadj failed");
        BB_FREE (temp_elist, int);

        for (i = 0; i < cuts->cutcount; i++) {
            cut = &cuts->cuts[i];
            if (cut->dominocount > 0) {
                x = cut_pi[i];
                if (BBbigguy_cmp (x, BBbigguy_ZERO) > 0) {
                    rval = domino_nzlist (&g, cut, cuts->cliques,
                                     cuts->dominos, &nzlist);
                    BBcheck_rval (rval, "domino_nzlist failed");
                    while (nzlist != -1) {
                        nznext = g.edges[nzlist].coefnext;
                        g.edges[nzlist].coefnext = -2;
                        if (g.edges[nzlist].coef) {
                            BBbigguy_addmult (&(elist[nzlist].rc), x,
                                              -g.edges[nzlist].coef);
                            g.edges[nzlist].coef = 0;
                        }
                        nzlist = nznext;
                    }
                } else if (BBbigguy_cmp (x, BBbigguy_ZERO) < 0) {
                    fprintf (stderr, "YIPES: negative domino\n");
                    rval = 1;  goto CLEANUP;
                }
            }
        }
    }

CLEANUP:
    BB_IFFREE (n, BBtsp_lpnode);
    BB_IFFREE (adjspace, BBtsp_lpadj);
    BB_IFFREE (temp_elist, int);
    free_lpgraph (&g);

    return rval; 
}

static void big_generate_edges (BBbigguy *node_piest, BBbigguy *node_domino,
        int nwant, int *gencount, bigpredge *genlist, int *n1, int *n2,
        int *finished, BBbigguy cutoff, BBdatagroup *dat, int ncount)
{
    int i = *n1, j = *n2;
    int len, cnt = 0;

    *gencount = 0;
    *finished = 0;

    if (i >= ncount) { *finished = 1; goto CLEANUP; }

    for (; j < ncount; j++) {
        len = BButil_dat_edgelen (i, j, dat);
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
            len = BButil_dat_edgelen (i, j, dat);
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

static int test_edge (int end1, int end2, int len, BBbigguy *node_pi,
        BBbigguy *node_domino, BBbigguy cutoff)
{
    BBbigguy rc;

    rc = BBbigguy_itobigguy (len);
    BBbigguy_sub (&rc, node_pi[end1]);
    BBbigguy_sub (&rc, node_pi[end2]);
    if (node_domino) {
        BBbigguy_sub (&rc, node_domino[end1]);
        BBbigguy_sub (&rc, node_domino[end2]);
    }
    if (BBbigguy_cmp (rc, cutoff) <= 0) {
        return 1;
    } else {
        return 0;
    }
}


/*********************  from TSP/tsp_lp.c  **********************************/

static int build_lpgraph (BBtsp_lpgraph *g, int ncount, int ecount,
        int *elist, int *elen)
{
    int i;
    BBtsp_lpnode *n;
    BBtsp_lpedge *e;

    g->ncount = ncount;
    g->ecount = ecount;
    g->nodes = BB_SAFE_MALLOC (ncount, BBtsp_lpnode);
    if (!g->nodes) {
        return 1;
    }
    g->edges = BB_SAFE_MALLOC (ecount, BBtsp_lpedge);
    if (!g->edges) {
        BB_FREE (g->nodes, BBtsp_lpnode);
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

static int build_lpadj (BBtsp_lpgraph *g, int estart, int eend)
{
    BBtsp_lpadj *a;
    BBtsp_lpnode *n = g->nodes;
    BBtsp_lpedge *e = g->edges;
    int i, j;

    if (g->adjspace) {
        if (g->adjstart == estart && g->adjend == eend) {
            return 0;
        } else {
            BB_FREE (g->adjspace, BBtsp_lpadj);
        }
    }

    if (estart >= eend) {
        g->adjstart = estart;
        g->adjend = eend;
        for (i=0; i<g->ncount; i++) {
            n[i].deg = 0;
            n[i].adj = (BBtsp_lpadj *) NULL;
        }
        return 0;
    }

    g->adjspace = BB_SAFE_MALLOC ((eend - estart)*2, BBtsp_lpadj);
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

static int domino_nzlist (BBtsp_lpgraph *g, BBtsp_lpcut *c,
        BBtsp_lpclique *cliques, BBtsp_lpdomino *dominos, int *pnzlist)
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

static void lpcut_nonzero_semicut (BBtsp_lpgraph *g, BBtsp_lpdomino *d,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges in the semi-cut (A:B) */

    int nzlist = *pnzlist;
    int nodemarker;
    BBtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    BB_FOREACH_NODE_IN_CLIQUE (k, d->sets[0], tmp) {
        g->nodes[k].mark = nodemarker;
    }

    BB_FOREACH_NODE_IN_CLIQUE (k, d->sets[1], tmp) {
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

static void lpcut_nonzero_parity_handle (BBtsp_lpgraph *g, BBtsp_lpclique *H,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges F, defined as the edges in an */
    /* odd number of semi-cuts (so an odd value in nzlist) and not in  */
    /* delta(H), or in an even number of semi-cuts and in delta(H).    */

    int nzlist = *pnzlist;
    int nodemarker;
    BBtsp_lpadj *a;
    int e;
    int tmp, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    BB_FOREACH_NODE_IN_CLIQUE (k, *H, tmp) {
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

    BB_FOREACH_NODE_IN_CLIQUE (k, *H, tmp) {
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

static void lpcut_nonzero_domino (BBtsp_lpgraph *g, BBtsp_lpdomino *d,
        int *pnzlist)
{
    /* Add in the nonzeros for the edges in the cut delta(A union B) */

    int nzlist = *pnzlist;
    int nodemarker;
    BBtsp_lpadj *a;
    int e;
    int tmp, i, k, l;

    g->nodemarker++;
    nodemarker = g->nodemarker;

    for (i = 0; i < 2; i++) {
        BB_FOREACH_NODE_IN_CLIQUE (k, d->sets[i], tmp) {
            g->nodes[k].mark = nodemarker;
        }
    }

    for (i = 0; i < 2; i++) {
        BB_FOREACH_NODE_IN_CLIQUE (k, d->sets[i], tmp) {
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

static void init_lpgraph_struct (BBtsp_lpgraph *g)
{
    g->ncount = 0;
    g->ecount = 0;
    g->nodes = (BBtsp_lpnode *) NULL;
    g->edges = (BBtsp_lpedge *) NULL;
    g->adjspace = (BBtsp_lpadj *) NULL;
    g->adjstart = 0;
    g->adjend = 0;
    g->nodemarker = 0;
    g->espace = 0;
}

static void free_lpgraph (BBtsp_lpgraph *g)
{
    BB_IFFREE (g->nodes, BBtsp_lpnode);
    BB_IFFREE (g->edges, BBtsp_lpedge);
    BB_IFFREE (g->adjspace, BBtsp_lpadj);
    g->espace = 0;
}
