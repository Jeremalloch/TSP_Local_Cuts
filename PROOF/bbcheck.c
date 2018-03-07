/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2007 by David Applegate, Robert Bixby,              */
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
/*           A PROGRAM TO VERIFY TSP PROOF INFORMATION                      */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 2, 2007                                                       */
/*                                                                          */
/****************************************************************************/

#include "bbproof.h"

static char *prooffname = (char *) NULL;
static char *tsplibfname = (char *) NULL;
static int verify_hk = 0;

static int
    checkcuts (int ncount, int ncuts, BBtsp_lpcuts *cuts, BBcuttype *cutlist,
        int doheld),
    check_cutproof (BBtsp_lpcut_in *c, BBcutproof *p, int ncount, int ind),
    cut_compare (BBtsp_lpcut_in *cut1, BBtsp_lpcut_in *cut2),
    checktree (BBtsp_lpcuts *cuts, BBproofnode *lproof, int ncount,
        BBdatagroup *dat, int ecount, int *elist, int fixcount, int *fixlist,
        double upbound, int n, int side, int depth, int silent),
    add_branch (BBtsp_lpcuts *cuts, BBtsp_branchobj *b, int side, int depth,
        int silent),
    remove_branch (BBtsp_lpcuts *cuts, BBtsp_branchobj *b, int side),
    parseargs (int ac, char **av);

static void
    print_branch (BBtsp_lpclique *c, char sense, int rhs, int depth),
    usage (char *fname);

int main (int ac, char **av)
{
    int i, ecount, fixcount, ncount = 0, rval = 0;
    char *probname = (char *) NULL;
    BBdatagroup dat;
    double upbound;
    int *elist = (int *) NULL;
    int *fixlist = (int *) NULL;
    BBtsp_lpcuts cuts;
    double szeit, tzeit;
    BBcuttype *cutlist = (BBcuttype *) NULL;
    BBcutproof *p = (BBcutproof *) NULL;
    int ncuts = 0;

    BButil_init_datagroup (&dat);
    BButil_init_tsp_lpcuts_struct (&cuts);

    rval = BButil_print_command (ac, av);
    BBcheck_rval (rval, "BButil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (tsplibfname == (char *) NULL) {
        fprintf (stderr, "Need TSPLIB file\n");
        rval = 1; goto CLEANUP;
    }

    szeit = BButil_zeit ();

    rval = BButil_gettsplib (tsplibfname, &i, &dat); 
    BBcheck_rval (rval, "BButil_gettsplib failed");

    rval = BBio_read_short_proof (&p, &cutlist, prooffname);
    BBcheck_rval (rval, "read_short_proof failed");

    if (i != p->ncount) {
        fprintf (stderr, "TSPLIB file not match proof file\n");
        rval = 1;  goto CLEANUP;
    }
    ncount = p->ncount;
    ncuts = p->cuts->cutcount;

    upbound = 0.0;
    for (i = 0; i < ncount-1; i++) {
       upbound += (double) BButil_dat_edgelen (p->tour[i],p->tour[i+1],&dat);
    }
    upbound += (double) BButil_dat_edgelen (p->tour[ncount-1],p->tour[0],&dat);
    printf ("Tour Value: %.0f\n", upbound); fflush (stdout);

    rval = BButil_datagroup_perm (ncount, &dat, p->tour);
    BBcheck_rval (rval, "BButil_datagroup_perm failed")

    tzeit =  BButil_zeit ();
    rval = checkcuts (ncount, ncuts, p->cuts, cutlist, verify_hk);
    BBcheck_rval (rval, "checkcuts failed");
    printf ("Cut Check Time: %.2f seconds\n", BButil_zeit () - tzeit);

    rval = BBprice_elim (p->cuts, p->lproof[0].exact_dual, &dat, upbound,
              p->lpbound, ncount, &elist, &ecount, &fixlist, &fixcount, 0);
    BBcheck_rval (rval, "BBprice_elim failed");

    tzeit =  BButil_zeit ();
    rval = checktree (p->cuts, p->lproof, ncount, &dat, ecount, elist,
                      fixcount, fixlist, upbound, 0, -1, 0, 0);
    BBcheck_rval (rval, "checktree failed");
    printf ("Tree Check Time: %.2f seconds\n", BButil_zeit () - tzeit);

    printf ("Verified optimal value: %.0f\n", upbound); 
    printf ("Total Running Time: %.2f seconds\n", BButil_zeit () - szeit);
    fflush (stdout);

CLEANUP:
    if (cutlist) {
        for (i = 0; i < ncuts; i++) {
            if (cutlist[i].proof) {
                BBio_free_cutproof (cutlist[i].proof);
                BB_FREE (cutlist[i].proof, BBcutproof);
            }
        }
        BB_FREE (cutlist, BBcuttype);
    }
    BB_IFFREE (probname, char);
    BButil_freedatagroup (&dat);
    if (p) BBio_free_cutproof (p);
    BB_IFFREE (p, BBcutproof);
    BB_IFFREE (elist, int);
    BB_IFFREE (fixlist, int);

    return rval;
}

static int checkcuts (int ncount, int ncuts, BBtsp_lpcuts *cuts, 
        BBcuttype *cutlist, int doheld)
{
    int i, m, type, no_outside, niso = 0, nhk = 0, nver = 0, rval = 0;
    int nsubtour = 0, ncomb = 0, nstar = 0, nbip = 0, ndomino = 0;
    BBtsp_lpcut_in cut, icut, bro, ibro;

    if (!cutlist) {
        for (i = 0; i < ncuts; i++) {
            rval = BButil_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &cut);
            BBcheck_rval (rval, "BButil_lpcut_to_lpcut_in failed");
            m = (BB_TYPE_SUBTOUR | BB_TYPE_COMB |
                 BB_TYPE_STAR | BB_TYPE_BIPARTITION);
            rval = BBcuts_verify (&cut, ncount, m, &type, (BBdatagroup *) NULL);
            BBcheck_rval (rval, "BBcuts_verify failed");
            switch (type) {
                case BB_TYPE_SUBTOUR:     nsubtour++;  break;
                case BB_TYPE_COMB:        ncomb++;  break;
                case BB_TYPE_STAR:        nstar++;  break;
                case BB_TYPE_BIPARTITION: nbip++;  break;
                default:  
                    fprintf (stderr, "unknown class %d\n", type);
                    rval = 1;  goto CLEANUP;
            }
            BButil_free_lpcut_in (&cut);
        }
        printf ("  LP has %d subtours, %d combs, %d stars, %d bipartition\n",
                   nsubtour, ncomb, nstar, nbip);
        fflush (stdout);
    } else {
        for (i = 0; i < ncuts; i++) {
            rval = BButil_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &cut);
            BBcheck_rval (rval, "BButil_lpcut_to_lpcut_in failed");

            if (cutlist[i].class != BB_NOCLASS) {
                rval = BBcuts_verify (&cut, ncount, cutlist[i].class, &type, 
                         (BBdatagroup *) NULL);
                BBcheck_rval (rval, "BBcuts_verify failed");
                if (type != cutlist[i].class) {
                    fprintf (stderr, "verified wrong class!\n");
                    rval = 1;  goto CLEANUP;
                }
                switch (type) {
                    case BB_TYPE_SUBTOUR:     nsubtour++;  break;
                    case BB_TYPE_COMB:        ncomb++;  break;
                    case BB_TYPE_STAR:        nstar++;  break;
                    case BB_TYPE_BIPARTITION: nbip++;  break;
                    case BB_TYPE_DOMINO:      ndomino++;  break;
                }
            } else if (cutlist[i].isomorph >= 0) {
                rval = BButil_lpcut_to_lpcut_in (cuts,
                             &(cuts->cuts[cutlist[i].isomorph]), &icut);
                BBcheck_rval (rval, "BButil_lpcut_to_lpcut_in failed");
                rval = BBcuts_backbone (&cut, &bro, &no_outside);
                BBcheck_rval (rval, "BBcuts_backbone failed");
                rval = BBcuts_backbone (&icut, &ibro, &no_outside);
                BBcheck_rval (rval, "BBcuts_backbode failed");

                if (cut_compare (&bro, &ibro)) {
                    fprintf (stderr, "cuts %d and %d are not isomorphs\n",
                                                   i, cutlist[i].isomorph);
                    rval = 1;  goto CLEANUP;
                }
              
                BButil_free_lpcut_in (&icut);
                BButil_free_lpcut_in (&bro);
                BButil_free_lpcut_in (&ibro);
                niso++;
            } else if (cutlist[i].hk == 1) {
                if (doheld) {
                    printf ("Call Held-Karp for cut %d ...\n", i);
                    fflush (stdout);
                    rval = BBcuts_verify (&cut, ncount, BB_TYPE_OTHER, &type,
                                                      (BBdatagroup *) NULL);
                    BBcheck_rval (rval, "BBcuts_verify failed for Held-Karp");
                }
                nhk++;
            } else {
                if (cutlist[i].proof) {
                    printf ("Verify proof for cut %d ...\n", i);
                    fflush (stdout);
                    rval = check_cutproof (&cut, cutlist[i].proof, ncount, i);
                    BBcheck_rval (rval, "check_cutproof failed");
                    nver++;
                } else {
                    printf ("Need a proof for cut %d\n", i); fflush (stdout);
                }
            }
            BButil_free_lpcut_in (&cut);
        }
        if (nhk > 0 && !doheld) {
            printf ("Need to verify %d cuts with Held-Karp\n", nhk);
            fflush (stdout);
        }
        printf ("Distribution of cuts\n");
        printf ("    %4d subtours\n", nsubtour);
        printf ("    %4d combs\n", ncomb);
        printf ("    %4d stars\n", nstar);
        printf ("    %4d bipartition inequalities\n", nbip);
        printf ("    %4d domino-parity inequalities\n", ndomino);
        printf ("    %4d other cuts (%d isomorphs, %d Held-Karp, %d proofs)\n",
                          niso+nhk+nver, niso, nhk, nver);
        fflush (stdout);
    }

CLEANUP:
    return rval;
}

static int check_cutproof (BBtsp_lpcut_in *c, BBcutproof *p, int ncount,
        int ind)
{
    int i, ncuts, ecount, fixcount, rval = 0;
    BBdatagroup vdat;
    double upbound;
    int *elist = (int *) NULL, *fixlist = (int *) NULL;


    BButil_init_datagroup (&vdat);

    if (p->ncount != c->skel.atomcount) {
        fprintf (stderr, "cut proof does not match atomcount");
        rval = 1;  goto CLEANUP;
    }
    ncuts = p->cuts->cutcount;

    rval = BBcuts_verify (c, ncount, 0, &i, &vdat);
    BBcheck_rval (rval, "BBcuts_verify failed for getdat");

    upbound = 0.0;
    for (i = 0; i < p->ncount-1; i++) {
       upbound += (double) BButil_dat_edgelen (p->tour[i], p->tour[i+1],
                                               &vdat);
    }
    upbound += (double) BButil_dat_edgelen (p->tour[p->ncount-1], p->tour[0],
                                           &vdat);

    if (c->rhs != (int) upbound) {
        fprintf (stderr, "Cut RHS %d, but proof gives %f\n", c->rhs, upbound);
        rval = 1; goto CLEANUP;
    }

    rval = BButil_datagroup_perm (p->ncount, &vdat, p->tour);
    BBcheck_rval (rval, "BButil_datagroup_perm failed")

    rval = checkcuts (ncount, ncuts, p->cuts, (BBcuttype *) NULL, 0);
    BBcheck_rval (rval, "checkcuts failed");

    rval = BBprice_elim (p->cuts, p->lproof[0].exact_dual, &vdat, upbound,
            p->lpbound, p->ncount, &elist, &ecount, &fixlist, &fixcount, 1);
    BBcheck_rval (rval, "BBprice_elim failed");

    rval = checktree (p->cuts, p->lproof, p->ncount, &vdat, ecount, elist,
                      fixcount, fixlist, upbound, 0, -1, 0, 1);
    BBcheck_rval (rval, "checktree failed");

    printf ("Verified RHS %d for cut %d\n", c->rhs, ind);  fflush (stdout); 

CLEANUP:
    BButil_freedatagroup (&vdat);
    BB_IFFREE (elist, int);
    BB_IFFREE (fixlist, int);
    return rval;
}

static int cut_compare (BBtsp_lpcut_in *cut1, BBtsp_lpcut_in *cut2)
{
    int i, j;
    BBtsp_lpclique *a, *b;

    if (cut1->cliquecount != cut2->cliquecount) return 1;
    if (cut1->rhs != cut2->rhs) return 1;
    if (cut1->sense != cut2->sense) return 1;
    for (i = 0; i < cut1->cliquecount; i++) {
        a = &cut1->cliques[i];
        b = &cut2->cliques[i];
        if (a->segcount != b->segcount) return 1;
        for (j = 0; j < a->segcount; j++) {
            if (a->nodes[j].lo != b->nodes[j].lo) return 1;
            if (a->nodes[j].hi != b->nodes[j].hi) return 1;
        }
    }
    return 0;
}

static int checktree (BBtsp_lpcuts *cuts, BBproofnode *lproof,
        int ncount, BBdatagroup *dat, int ecount, int *elist, int fixcount,
        int *fixlist, double upbound, int n, int side, int depth, int silent)
{
    int rval = 0;
    int parent;
    BBbigguy bound, ubound;

    if (!silent) {
        printf ("Node %d: Side=%d Child0=%d Child1=%d Parent=%d\n",
                n, lproof[n].side, lproof[n].child0, lproof[n].child1,
                lproof[n].parent);
    }

    if (side != lproof[n].side) {
        printf ("Error: sides do not agree (%d)\n", side);
        rval = 1;  goto CLEANUP;
    }

    if ((lproof[n].child0 == -1 && lproof[n].child1 != -1) &&
        (lproof[n].child1 == -1 && lproof[n].child0 != -1)) {
        printf ("Error: node with only one child\n");
        rval = 1;  goto CLEANUP;
    }

    parent = lproof[n].parent;

    if (depth > 0) {
        rval = add_branch (cuts, lproof[parent].branch, lproof[n].side,
                           depth, silent);
        BBcheck_rval (rval, "add_branch failed");
    }

    if (depth == 0 || lproof[n].child0 == -1) {
        if (!lproof[n].exact_dual) {
            printf ("Error: Missing exact_dual data\n"); fflush (stdout);
            rval = 1;  goto CLEANUP;
        }
        if (lproof[n].exact_dual->cutcount != cuts->cutcount) {
            printf ("Error: Mismatch in cutcount\n"); fflush (stdout);
            rval = 1;  goto CLEANUP;
        }
        rval = BBprice_price (cuts, lproof[n].exact_dual, dat,
                     ncount, elist, ecount, fixlist, fixcount, &bound, silent);
        BBcheck_rval (rval, "BBprice_price failed");

        if (!silent) {
            printf ("Exactbound: %f", BBbigguy_bigguytod (bound));
            fflush (stdout);
        }
        if (depth > 0) {
            ubound = BBbigguy_dtobigguy (upbound);
            BBbigguy_sub (&ubound, BBbigguy_ONE);
            if (BBbigguy_cmp (bound, ubound) > 0) {
                if (!silent) { printf ("  PRUNED\n"); fflush (stdout); }
            } else {
                printf ("\nERROR: Cannot verify pruned node\n");
                fflush (stdout);
                rval = 1; goto CLEANUP;
            }
        } else {
           if (!silent) { printf ("  Root\n"); fflush (stdout); }
        }
    }

    if (lproof[n].child0 != -1) {
        rval = checktree (cuts, lproof, ncount, dat, ecount, elist, fixcount,
                  fixlist, upbound, lproof[n].child0, 0, depth+1, silent);
        BBcheck_rval (rval, "checktree failed");
    }
    if (lproof[n].child1 != -1) {
        rval = checktree (cuts, lproof, ncount, dat, ecount, elist, fixcount,
                  fixlist, upbound, lproof[n].child1, 1, depth+1, silent);
        BBcheck_rval (rval, "checktree failed");
    }

    if (depth > 0) {
        rval = remove_branch (cuts, lproof[parent].branch, lproof[n].side);
        BBcheck_rval (rval, "remove_branch failed");
    }

CLEANUP:
    return rval;
}

static int add_branch (BBtsp_lpcuts *cuts, BBtsp_branchobj *b, int side,
        int depth, int silent)
{
    int ar[2], rval = 0;
    BBtsp_lpcut *c;

    if (b == NULL) {
        fprintf (stderr, "add_branch called without a branchobj\n");
        rval = 1; goto CLEANUP;
    }
    if (!silent) { printf ("Add Side=%d ", side); fflush (stdout); }

    if (cuts->cutcount >= cuts->cutspace ||
                          cuts->cliqueend >= cuts->cliquespace) {
        fprintf (stderr, "out of room for branching cut or clique");
        rval = 1;  goto CLEANUP;
    }
    c = &cuts->cuts[cuts->cutcount++];
    BButil_init_lpcut (c);
    c->cliques = BB_SAFE_MALLOC (1, int);
    BBcheck_NULL (c->cliques, "out of memory for c->cliques");
    c->cliquecount = 1;

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        rval = BButil_array_to_lpclique (ar,2,&cuts->cliques[cuts->cliqueend]);
        BBcheck_rval (rval, "BButil_array_to_lpclique failed");
        if (side == 1) {
            c->rhs = 2.0;
            c->sense = 'L';
        } else {
            c->rhs = 4.0;
            c->sense = 'G';
        }
    } else {
        if (!b->clique) {
            fprintf (stderr, "BBtsp_branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        rval = BButil_copy_lpclique (b->clique,
                                 &cuts->cliques[cuts->cliqueend]);
        BBcheck_rval (rval, "BBtsp_copy_lpclique failed");
        if (side == 0) {
            c->rhs = 2.0;
            c->sense = 'L';
        } else {
            c->rhs = 4.0;
            c->sense = 'G';
        }
    }
    c->cliques[0] = cuts->cliqueend++;
    c->branch = 1;

    if (!silent) {
        print_branch (&cuts->cliques[c->cliques[0]], c->sense, c->rhs, depth);
    }

CLEANUP:
    return rval;
}

static int remove_branch (BBtsp_lpcuts *cuts, BBtsp_branchobj *b, int side)
{
    int rval = 0;
    int ar[2];
    int k;
    char sense;
    int rhs, num;
    BBtsp_lpcut *cu;
    BBtsp_lpclique *c = (BBtsp_lpclique *) NULL;

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        c = BB_SAFE_MALLOC (1, BBtsp_lpclique);
        BBcheck_NULL (c, "out of memory for c");

        rval = BButil_array_to_lpclique (ar, 2, c);
        BBcheck_rval (rval, "BButil_array_to_lpclique failed");

        if (side == 1) {
            sense = 'L';
            rhs = 2;
        } else {
            sense = 'G';
            rhs = 4;
        }
    } else {
        c = b->clique;
        if (side == 0) {
            sense = 'L';
            rhs = 2;
        } else {
            sense = 'G';
            rhs = 4;
        }
    }

    num = cuts->cutcount - 1;

    /* Check that last cut in the LP is really the branching cut */

    cu = &(cuts->cuts[num]);
    if (cu->cliquecount != 1 || cu->sense != sense || cu->rhs != rhs) {
        printf ("Error: Last LP row does not match branch\n");
        printf ("Count = %d, Sense = %c, RHS = %d\n", cu->cliquecount,
                 cu->sense, cu->rhs);
        rval = 1;  goto CLEANUP;
    }

    BButil_lpclique_compare (&(cuts->cliques[cu->cliques[0]]),c, &k);
    if (k != 0) {
        printf ("Error: Last LP row clique does not match branch\n");
        rval = 1;  goto CLEANUP;
    }

    BB_IFFREE (cu->cliques, int);
    cuts->cutcount--;

    if (b->ends[0] != -1) {
        BButil_free_lpclique (c);
        BB_FREE (c, BBtsp_lpclique);
    }

CLEANUP:
    return rval;
}

static void print_branch (BBtsp_lpclique *c, char sense, int rhs, int depth)
{
    int i;

    printf ("Depth=%d ", depth);
    printf ("Clique ");
    for (i = 0; i < c->segcount; i++) {
        printf ("%d->%d ", c->nodes[i].lo, c->nodes[i].hi);
    }
    if (sense == 'L') {
        printf ("<= %d\n", rhs);
    } else {
        printf (">= %d\n",rhs);
    }
    fflush (stdout);
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = BButil_bix_getopt (ac, av, "h?", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'h':
            verify_hk= 1;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind != ac-2) {
        usage (av[0]);
        return 1;
    }

    tsplibfname = av[boptind++];
    prooffname = av[boptind++];

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-h] TSPLIB_file proof_file\n", fname);
    fprintf (stderr, "   -h    verify non-classified cuts with Held-Karp\n");
}

