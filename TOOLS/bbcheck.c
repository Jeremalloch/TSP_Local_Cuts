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
/*               A PROGRAM TO VERIFY PROOF INFORMATION                      */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 2, 2007                                                       */
/*                                                                          */
/*  SEE short describtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "lp.h"

typedef struct proofnode {
    int number; 
    int parent;
    int side;
    int child0;
    int child1;
    CCtsp_branchobj *branch;
    CCtsp_bigdual *exact_dual;
} proofnode;

typedef struct cutproof  {
    int ncount;
    int bbcount;
    int *tour;
    char *probname;
    proofnode *lproof;
    CCtsp_lpcuts *cuts;
    CCbigguy lpbound;
} cutproof;

typedef struct cuttype {
    int class;
    int isomorph;
    int hk;
    cutproof *proof;
} cuttype;


static char *prooffname = (char *) NULL;
static char *tsplibfname = (char *) NULL;
static int run_silent = 1;
static int verify_hk = 0;

int BBverify_cut (CCtsp_lpcut_in *cut, int ncount, int check_types, int *type,
        int use_tsp, int *did_use_tsp, char *cname, CCdatagroup *getdat);

static int
    checkcuts (int ncount, int ncuts, CCtsp_lpcuts *cuts, cuttype *cutlist,
        int doheld),
    check_cutproof (CCtsp_lpcut_in *c, cutproof *p, int ncount, int ind),
    checktree (CCtsp_lpcuts *cuts, proofnode *lproof, int ncount,
        CCdatagroup *dat, int ecount, int *elist, int fixcount, int *fixlist,
        double upbound, int n, int side, int depth, int silent),
    add_branch (CCtsp_lpcuts *cuts, CCtsp_branchobj *b, int side, int depth,
        int silent),
    remove_branch (CCtsp_lpcuts *cuts, CCtsp_branchobj *b, int side),
    read_short_proof (cutproof **pp, cuttype **pcutlist, char *fname),
    read_short_proof_work (CC_SFILE *f, cutproof *p, int silent),
    read_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *cuts),
    read_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount),
    read_lpdomino (CC_SFILE *f, CCtsp_lpdomino *d, int ncount),
    parseargs (int ac, char **av);

static void
    init_proofnode (proofnode *p),
    free_proofnode (proofnode *p),
    free_lpcuts (CCtsp_lpcuts *cuts),
    init_cutproof (cutproof *p),
    free_cutproof (cutproof *p),
    print_branch (CCtsp_lpclique *c, char sense, int rhs, int depth),
    usage (char *fname);

int CCbbproof_price (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCdatagroup *dat, int ncount, int *elist, int ecount, int *fixlist,
        int fcount, CCbigguy *bound, int silent);
int CCbbproof_elim (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCdatagroup *dat, double upperbound, CCbigguy exact_lowerbound,
        int ncount, int **elist, int *ecount, int **fixlist, int *fcount,
        int silent);

#define BB_NOCLASS 99999
#define BB_TYPE_SUBTOUR 1
#define BB_TYPE_COMB 2
#define BB_TYPE_STAR 4
#define BB_TYPE_BIPARTITION 8
#define BB_TYPE_OTHER 16
#define BB_TYPE_DOMINO 64



int main (int ac, char **av)
{
    int i, ecount, fixcount, ncount = 0, rval = 0;
    char *probname = (char *) NULL;
    int *ptour = (int *) NULL;
    CCrandstate rstate;
    CCdatagroup dat;
    int seed = 99;
    double upbound;
    int *elist = (int *) NULL;
    int *fixlist = (int *) NULL;
    CCtsp_lpcuts cuts;
    double szeit;
    cuttype *cutlist = (cuttype *) NULL;
    cutproof *p = (cutproof *) NULL;
    int ncuts = 0;

    CCutil_init_datagroup (&dat);
    CCtsp_init_tsp_lpcuts_struct (&cuts);

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (tsplibfname == (char *) NULL) {
        fprintf (stderr, "Need TSPLIB file\n");
        rval = 1; goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    rval = CCutil_gettsplib (tsplibfname, &i, &dat); 
    CCcheck_rval (rval, "CCutil_gettsplib failed");

    rval = read_short_proof (&p, &cutlist, prooffname);
    CCcheck_rval (rval, "read_short_proof failed");

    if (i != p->ncount) {
        fprintf (stderr, "TSPLIB file not match proof file\n");
        rval = 1;  goto CLEANUP;
    }
    ncount = p->ncount;
    ncuts = p->cuts->cutcount;

    upbound = 0.0;
    for (i = 0; i < ncount-1; i++) {
       upbound += (double) CCutil_dat_edgelen (p->tour[i],p->tour[i+1],&dat);
    }
    upbound += (double) CCutil_dat_edgelen (p->tour[ncount-1],p->tour[0],&dat);
    printf ("Tour Value: %.0f\n", upbound); fflush (stdout);

    rval = CCutil_datagroup_perm (ncount, &dat, p->tour);
    CCcheck_rval (rval, "CCutil_datagroup_perm failed")

    rval = checkcuts (ncount, ncuts, p->cuts, cutlist, verify_hk);
    CCcheck_rval (rval, "checkcuts failed");

    rval = CCbbproof_elim (p->cuts, p->lproof[0].exact_dual, &dat, upbound,
              p->lpbound, ncount, &elist, &ecount, &fixlist, &fixcount, 0);
    CCcheck_rval (rval, "CCbbproof_elim failed");

    rval = checktree (p->cuts, p->lproof, ncount, &dat, ecount, elist,
                      fixcount, fixlist, upbound, 0, -1, 0, 0);
    CCcheck_rval (rval, "checktree failed");

    printf ("Verified optimal value: %.0f\n", upbound); 
    printf ("Total Running Time: %.2f seconds\n", CCutil_zeit () - szeit);
    fflush (stdout);


CLEANUP:

    if (cutlist) {
        for (i = 0; i < ncuts; i++) {
            if (cutlist[i].proof) free_cutproof (cutlist[i].proof);
        }
        CC_FREE (cutlist, cuttype);
    }
    CC_IFFREE (ptour, int);
    CC_IFFREE (probname, char);
    CCutil_freedatagroup (&dat);
    if (p) free_cutproof (p);
    CC_IFFREE (p, cutproof);
    CC_IFFREE (elist, int);
    CC_IFFREE (fixlist, int);

    return rval;
}

static int checkcuts (int ncount, int ncuts, CCtsp_lpcuts *cuts, 
        cuttype *cutlist, int doheld)
{
    int i, type, niso = 0, nhk = 0, nver = 0, rval = 0;
    int nsubtour = 0, ncomb = 0, nstar = 0, nbip = 0, ndomino = 0;
    CCtsp_lpcut_in cut;

    if (!cutlist) {
        fprintf (stderr, "Need a cutlist to verify the cuts\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncuts; i++) {
        rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &cut);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

        if (cutlist[i].class != BB_NOCLASS) {
            rval = BBverify_cut (&cut, ncount, cutlist[i].class, &type, 0,
                     (int *) NULL, (char *) NULL, (CCdatagroup *) NULL);
            CCcheck_rval (rval, "BBverify_cut failed");
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
            niso++;
        } else if (cutlist[i].hk == 1) {
            if (doheld) {
                printf ("Call Held-Karp for cut %d ...\n", i);
                fflush (stdout);
                rval = BBverify_cut (&cut, ncount, BB_TYPE_OTHER, &type, 0,
                     (int *) NULL, (char *) NULL, (CCdatagroup *) NULL);
                CCcheck_rval (rval, "BBverify_cut failed for Held-Karp");
            }
            nhk++;
        } else {
            if (cutlist[i].proof) {
                printf ("Verify proof for cut %d ...\n", i);
                fflush (stdout);
                rval = check_cutproof (&cut, cutlist[i].proof, ncount, i);
                CCcheck_rval (rval, "check_cutproof failed");
                nver++;
            } else {
                printf ("Need a proof for cut %d\n", i); fflush (stdout);
            }
        }
        CCtsp_free_lpcut_in (&cut);
    }
    if (!doheld) {
        printf ("Need to verify %d cuts with Held-Karp\n", nhk);
        fflush (stdout);
    }


    printf ("Distribution of cuts\n");
    printf ("    %4d subtours\n", nsubtour);
    printf ("    %4d combs\n", ncomb);
    printf ("    %4d stars\n", nstar);
    printf ("    %4d bipartition inequalities\n", nbip);
    printf ("    %4d other cuts (%d isomorphs, %d Held-Karp, %d proofs)\n",
                      niso+nhk+nver, niso, nhk, nver);
    fflush (stdout);

CLEANUP:
    return rval;
}


static int check_cutproof (CCtsp_lpcut_in *c, cutproof *p, int ncount, int ind)
{
    int i, ecount, fixcount, rval = 0;
    CCdatagroup vdat;
    double upbound;
    int *elist = (int *) NULL, *fixlist = (int *) NULL;

    CCutil_init_datagroup (&vdat);

    if (p->ncount != c->skel.atomcount) {
        fprintf (stderr, "cut proof does not match atomcount");
        rval = 1;  goto CLEANUP;
    }

    rval = BBverify_cut (c, ncount, 0, &i, 0, (int *) NULL,
                        (char *) NULL, &vdat);
    CCcheck_rval (rval, "BBverify_cut failed for getdat");

    upbound = 0.0;
    for (i = 0; i < p->ncount-1; i++) {
       upbound += (double) CCutil_dat_edgelen (p->tour[i], p->tour[i+1],
                                               &vdat);
    }
    upbound += (double) CCutil_dat_edgelen (p->tour[p->ncount-1], p->tour[0],
                                           &vdat);

    if (c->rhs != (int) upbound) {
        fprintf (stderr, "Cut RHS %d, but proof gives %f\n", c->rhs, upbound);
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_datagroup_perm (p->ncount, &vdat, p->tour);
    CCcheck_rval (rval, "CCutil_datagroup_perm failed")

    rval = CCbbproof_elim (p->cuts, p->lproof[0].exact_dual, &vdat, upbound,
            p->lpbound, p->ncount, &elist, &ecount, &fixlist, &fixcount, 1);
    CCcheck_rval (rval, "CCbbproof_elim failed");

    rval = checktree (p->cuts, p->lproof, p->ncount, &vdat, ecount, elist,
                      fixcount, fixlist, upbound, 0, -1, 0, 1);
    CCcheck_rval (rval, "checktree failed");

    printf ("Verified RHS %d for cut %d\n", c->rhs, ind);  fflush (stdout); 

CLEANUP:
    CCutil_freedatagroup (&vdat);
    return rval;
}

static int checktree (CCtsp_lpcuts *cuts, proofnode *lproof,
        int ncount, CCdatagroup *dat, int ecount, int *elist, int fixcount,
        int *fixlist, double upbound, int n, int side, int depth, int silent)
{
    int rval = 0;
    int parent;
    CCbigguy bound, ubound;

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
        CCcheck_rval (rval, "add_branch failed");
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
        rval = CCbbproof_price (cuts, lproof[n].exact_dual, dat,
                     ncount, elist, ecount, fixlist, fixcount, &bound, silent);
        CCcheck_rval (rval, "CCbbproof_price failed");

        if (!silent) {
            printf ("Exactbound: %f", CCbigguy_bigguytod (bound));
            fflush (stdout);
        }
        if (depth > 0) {
            ubound = CCbigguy_dtobigguy (upbound);
            CCbigguy_sub (&ubound, CCbigguy_ONE);
            if (CCbigguy_cmp (bound, ubound) > 0) {
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
        CCcheck_rval (rval, "checktree failed");
    }
    if (lproof[n].child1 != -1) {
        rval = checktree (cuts, lproof, ncount, dat, ecount, elist, fixcount,
                  fixlist, upbound, lproof[n].child1, 1, depth+1, silent);
        CCcheck_rval (rval, "checktree failed");
    }

    if (depth > 0) {
        rval = remove_branch (cuts, lproof[parent].branch, lproof[n].side);
        CCcheck_rval (rval, "remove_branch failed");
    }

CLEANUP:

    return rval;
}

static int add_branch (CCtsp_lpcuts *cuts, CCtsp_branchobj *b, int side,
        int depth, int silent)
{
    int ar[2], rval = 0;
    CCtsp_lpcut *c;

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
    CCtsp_init_lpcut (c);
    c->cliques = CC_SAFE_MALLOC (1, int);
    CCcheck_NULL (c->cliques, "out of memory for c->cliques");
    c->cliquecount = 1;

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        rval = CCtsp_array_to_lpclique (ar, 2, &cuts->cliques[cuts->cliqueend]);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        if (side == 1) {
            c->rhs = 2.0;
            c->sense = 'L';
        } else {
            c->rhs = 4.0;
            c->sense = 'G';
        }
    } else {
        if (!b->clique) {
            fprintf (stderr, "CCtsp_branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_copy_lpclique (b->clique, &cuts->cliques[cuts->cliqueend]);
        CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
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

static int remove_branch (CCtsp_lpcuts *cuts, CCtsp_branchobj *b, int side)
{
    int rval = 0;
    int ar[2];
    int k;
    char sense;
    int rhs, num;
    CCtsp_lpcut *cu;
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;

    /* printf ("Remove Side=%d ", side); fflush (stdout); */

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        CCcheck_NULL (c, "out of memory for c");

        rval = CCtsp_array_to_lpclique (ar, 2, c);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

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

    CCtsp_lpclique_compare (&(cuts->cliques[cu->cliques[0]]),c, &k);
    if (k != 0) {
        printf ("Error: Last LP row clique does not match branch\n");
        rval = 1;  goto CLEANUP;
    }

    CC_IFFREE (cu->cliques, int);
    cuts->cutcount--;

    if (b->ends[0] != -1) {
        CCtsp_free_lpclique (c);
        CC_FREE (c, CCtsp_lpclique);
    }

CLEANUP:
    return rval;
}

static void init_proofnode (proofnode *p) 
{
    if (p) {
        p->number = -1;
        p->parent = -1;
        p->child0 = -1;
        p->child1 = -1;
        p->side   = -1;
        p->branch = (CCtsp_branchobj *) NULL;
        p->exact_dual = (CCtsp_bigdual *) NULL;
    }
}

static void free_proofnode (proofnode *p) 
{
    if (p) {
        if (p->branch) {
            CCtsp_free_branchobj (p->branch);
            CC_FREE (p->branch, CCtsp_branchobj);
        }
        if (p->exact_dual) {
            CCtsp_free_bigdual (&p->exact_dual);
        }
        init_proofnode (p);
    }
}

static void print_branch (CCtsp_lpclique *c, char sense, int rhs, int depth)
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

static int read_short_proof (cutproof **pp, cuttype **pcutlist, char *fname)
{
    int rval = 0;
    CC_SFILE *f = (CC_SFILE *) NULL;
    int i, ind, cutcount, proofcount;
    cutproof *p;
    cuttype *clist = (cuttype *) NULL;

    *pcutlist = (cuttype *) NULL;
    *pp = (cutproof *) NULL;

    p = CC_SAFE_MALLOC (1, cutproof);
    CCcheck_NULL (p, "out of memory for p");
    init_cutproof (p);

    f = CCutil_sopen (fname, "r");
    if (!f) {
        fprintf (stderr, "could not safe open %s for writing\n", fname);
        rval =1; goto CLEANUP;
    }  

    rval = read_short_proof_work (f, p, 0);
    CCcheck_rval (rval, "read_short_proof_work failed");

    /* Get cutlist if present */

    rval = CCutil_sread_int (f, &cutcount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    if (cutcount) {
        if (cutcount != p->cuts->cutcount) {
            fprintf (stderr, "cutlist does not match the proof\n");
            rval = 1; goto CLEANUP;
        }

        clist = CC_SAFE_MALLOC (cutcount, cuttype);
        CCcheck_NULL (clist, "out of memory for clist");
        for (i = 0; i < cutcount; i++) {
            rval = CCutil_sread_int (f, &clist[i].class);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            rval = CCutil_sread_int (f, &clist[i].isomorph);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            rval = CCutil_sread_int (f, &clist[i].hk);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            clist[i].proof = (cutproof *) NULL;
        }
        *pcutlist = clist;

        rval = CCutil_sread_int (f, &proofcount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        printf ("File has %d cut proofs\n", proofcount); fflush (stdout);
        for (i = 0; i < proofcount; i++) {
            rval = CCutil_sread_int (f, &ind);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            clist[ind].proof = CC_SAFE_MALLOC (1, cutproof);
            CCcheck_NULL (clist[ind].proof, "out of memory for proof");
            printf ("Read proof for cut %d\n", ind); fflush (stdout);
            rval = read_short_proof_work (f, clist[ind].proof, 1);
            CCcheck_rval (rval, "read_short_proof_work failed for cut");
        }
    }
    *pp = p;

CLEANUP:
    if (f) CCutil_sclose (f);
    return rval;
}

static int read_short_proof_work (CC_SFILE *f, cutproof *p, int silent)
{
    int rval = 0;
    int i, j, gotbranch, segcount, cutcount;
    CCtsp_bigdual *d = (CCtsp_bigdual *) NULL;
    CCtsp_branchobj *b = (CCtsp_branchobj *) NULL;
    proofnode *t;

    p->probname = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    CCcheck_NULL (p->probname, "out of memory for probname");

    for (i = 0; i < CCtsp_PROB_FILE_NAME_LEN; i++) {
        rval = CCutil_sread_char (f, &p->probname[i]);
        CCcheck_rval (rval, "CCutil_sread_char failed");
    }

    rval = CCutil_sread_int (f, &p->ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    p->cuts = CC_SAFE_MALLOC (1, CCtsp_lpcuts);
    CCcheck_NULL (p->cuts, "out of memory for cuts");
    CCtsp_init_tsp_lpcuts_struct (p->cuts);

    rval = read_cuts (f, p->ncount, p->cuts);
    CCcheck_rval (rval, "read_cuts failed");
 
    if (!silent) {
        printf ("Number of cuts: %d\n", p->cuts->cutcount); fflush (stdout);
    }

    p->tour = CC_SAFE_MALLOC (p->ncount, int);
    CCcheck_NULL (p->tour, "out of memory for tour");

    for (i = 0; i < p->ncount; i++) {
        rval = CCutil_sread_int (f, &p->tour[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    rval = CCbigguy_sread (f, &p->lpbound);
    CCcheck_rval (rval, "CCbigguy_sread failed");

    if (!silent) {
        printf ("LP bound: %12f\n", CCbigguy_bigguytod (p->lpbound));
        fflush (stdout);
    }

    rval = CCutil_sread_int (f, &p->bbcount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    t = CC_SAFE_MALLOC (p->bbcount, proofnode);
    CCcheck_NULL (t, "out of memory for lproof");
    p->lproof = t;

    for (i = 0; i < p->bbcount; i++) {
        rval = CCutil_sread_int (f, &t[i].number);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &t[i].parent);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &t[i].side);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &t[i].child0);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &t[i].child1);
        CCcheck_rval (rval, "CCutil_sread_int failed");

        rval = CCutil_sread_int (f, &gotbranch);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        if (gotbranch) {
            b = CC_SAFE_MALLOC (1, CCtsp_branchobj);
            CCcheck_NULL (b, "out of memory for branchobj");
            CCtsp_init_branchobj (b);
            rval = CCutil_sread_int (f, &b->ends[0]);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            rval = CCutil_sread_int (f, &b->ends[1]);
            CCcheck_rval (rval, "CCutil_sread_int failed");

            rval = CCutil_sread_int (f, &segcount);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            if (segcount) {
                b->clique = CC_SAFE_MALLOC (1, CCtsp_lpclique);
                CCcheck_NULL (b->clique, "out or memory for clique");
                CCtsp_init_lpclique (b->clique);
                b->clique->segcount = segcount;
                b->clique->nodes = CC_SAFE_MALLOC (segcount, CCtsp_segment);
                CCcheck_NULL (b->clique->nodes, "out of memory for nodes");
                for (j = 0; j < segcount; j++) {
                    rval = CCutil_sread_int (f, &b->clique->nodes[j].lo);
                    CCcheck_rval (rval, "CCutil_sread_int failed");
                    rval = CCutil_sread_int (f, &b->clique->nodes[j].hi);
                    CCcheck_rval (rval, "CCutil_sread_int failed");
                }
            } else {
                b->clique = (CCtsp_lpclique *) NULL;
            }
            t[i].branch = b;
        } else {
            t[i].branch = (CCtsp_branchobj *) NULL;
        }
        rval = CCutil_sread_int (f, &cutcount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        if (cutcount) {
            d = CC_SAFE_MALLOC (1, CCtsp_bigdual);
            CCcheck_NULL (d, "out of memory for d");
            d->cutcount = cutcount;
            d->node_pi = CC_SAFE_MALLOC (p->ncount, CCbigguy);
            CCcheck_NULL (d->node_pi, "out of memory for node_pi");
            d->cut_pi = CC_SAFE_MALLOC (cutcount, CCbigguy);
            CCcheck_NULL (d->cut_pi, "out of memory for cut_pi");
            for (j = 0; j < p->ncount; j++) {
                rval = CCbigguy_sread (f, &d->node_pi[j]);
                CCcheck_rval (rval, "CCbigguy_sread failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = CCbigguy_sread (f, &d->cut_pi[j]);
                CCcheck_rval (rval, "CCbigguy_sread failed");
            }

            t[i].exact_dual = d;
        } else {
            t[i].exact_dual = (CCtsp_bigdual *) NULL;
        }
    }
    if (!silent) {
        printf ("BBcount: %d\n", p->bbcount); fflush (stdout);
    }

CLEANUP:
    return rval;
}

#define BB_EXTRA_CLIQUES 1000

static int read_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *c)
{
    int rval = 0;
    int i, j, nbits, cbits, dbits, cliqcount, dominocount, cutcount;
    CCtsp_lpcut *u;
    char version;

    rval = CCutil_sread_char (f, &version);
    CCcheck_rval (rval, "CCutil_sread_char failed");
    if (version != 3) {
        fprintf (stderr, "Possibly old cut version %d\n", (unsigned) version);
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_sread_int (f, &i);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    if (i != ncount) {
        fprintf (stderr, "Cuts do not match the proof\n");
        rval = 1;  goto CLEANUP;
    }
    nbits = CCutil_sbits (ncount);

    rval = CCutil_sread_int (f, &cliqcount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    c->cliques = CC_SAFE_MALLOC (cliqcount+BB_EXTRA_CLIQUES, CCtsp_lpclique);
    CCcheck_NULL (c->cliques, "out of memory for c->cliques");
    c->cliquespace = cliqcount+BB_EXTRA_CLIQUES;
    for (i = 0; i < cliqcount; i++) {
        rval = read_lpclique (f, &c->cliques[i], ncount);
        CCcheck_rval (rval, "read_lpclique failed");
    }
    c->cliqueend = cliqcount;

    rval = CCutil_sread_int (f, &dominocount);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    if (dominocount > 0) {
        c->dominos  = CC_SAFE_MALLOC (dominocount, CCtsp_lpdomino);
        for (i = 0; i < dominocount; i++) {
            rval = read_lpdomino (f, &c->dominos[i], ncount);
            CCcheck_rval (rval, "read_lpdomino failed");
        }
    }
    c->dominoend = dominocount;

    cbits = CCutil_sbits (cliqcount);
    dbits = CCutil_sbits (dominocount);
    rval = CCutil_sread_int (f, &cutcount);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    c->cuts = CC_SAFE_MALLOC (cutcount+BB_EXTRA_CLIQUES, CCtsp_lpcut);
    CCcheck_NULL (c->cuts, "out of memory for c->cuts");
    c->cutspace = cutcount+BB_EXTRA_CLIQUES;
    for (i = 0; i < cutcount; i++) {
        u = &c->cuts[i];
        rval = CCutil_sread_int (f, &u->cliquecount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &u->twodom_cliquecount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &u->dominocount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &u->twodomcount[0]);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &u->twodomcount[1]);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        if (u->twodom_cliquecount > 0 || u->twodomcount[0] > 0 ||
                                         u->twodomcount[1] > 0 ) {
            fprintf (stderr, "Not set up for 2-dom cuts\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCutil_sread_int (f, &u->rhs);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_char (f, &u->sense);
        u->branch = 0;
        u->cliques = CC_SAFE_MALLOC (u->cliquecount, int);
        CCcheck_NULL (u->cliques, "out of memory for u->cliques");
        for (j = 0; j < u->cliquecount; j++) {
            rval = CCutil_sread_bits (f, &u->cliques[j], cbits);
            CCcheck_rval (rval, "CCutil_sread_bits failed");
        }
        if (u->dominocount > 0) {
            u->dominos = CC_SAFE_MALLOC (u->dominocount, int);
            CCcheck_NULL (u->dominos, "out of memory for u->dominos");
            for (j = 0; j < u->dominocount; j++) {
                rval = CCutil_sread_bits (f, &u->dominos[j], dbits);
                CCcheck_rval (rval, "CCutil_sread_bits failed");
            }
        } else {
            u->dominos = (int *) NULL;
        }
        rval = CCtsp_read_skeleton (f, &u->skel, ncount);
        CCcheck_rval (rval, "CCtsp_read_skeleton failed");
    }
    c->cutcount = cutcount;

CLEANUP:
    return rval;
}

static int read_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount)
{
    int nbits = CCutil_sbits (ncount);
    int i, lo, hi, rval = 0;

    rval = CCutil_sread_bits (f, &c->segcount, nbits);
    CCcheck_rval (rval, "CCutil_sread_bits failed");
    c->nodes = CC_SAFE_MALLOC (c->segcount, CCtsp_segment);
    CCcheck_NULL (c->nodes, "out of memory for c->nodes");

    for (i = 0; i < c->segcount; i++) {
        rval = CCutil_sread_bits (f, &lo, nbits);
        CCcheck_rval (rval, "CCutil_sread_bits failed");
        rval = CCutil_sread_bits (f, &hi, nbits);
        CCcheck_rval (rval, "CCutil_sread_bits failed");
        c->nodes[i].lo = lo;
        c->nodes[i].hi = hi;
    }
    c->hashnext = c->refcount = 0;

CLEANUP:
    return rval;
}

static int read_lpdomino (CC_SFILE *f, CCtsp_lpdomino *d, int ncount)
{
    int k, rval = 0;

    for (k = 0; k < 2; k++) {
        rval = read_lpclique (f, &(d->sets[k]), ncount);
        CCcheck_rval (rval, "read_lpclique failed");
    }

CLEANUP:
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "D:hv?", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'D':
            tsplibfname = boptarg;
            break;
        case 'h':
            verify_hk= 1;
            break;
        case 'v':
            run_silent = 0;
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
    fprintf (stderr, "Usage: %s [-pv] TSPLIB_file proof_file\n", fname);
    fprintf (stderr, "   -h    verify non-classified cuts with Held-Karp\n");
    fprintf (stderr, "   -v    verbose\n");
}

static void free_lpcuts (CCtsp_lpcuts *cuts)
{
    int i;

    if (cuts->cuts) {
        for (i=0; i < cuts->cutcount; i++) {
            CC_IFFREE (cuts->cuts[i].cliques, int);
            CC_IFFREE (cuts->cuts[i].dominos, int);
            CCtsp_free_skeleton (&cuts->cuts[i].skel);
        }
        CC_FREE (cuts->cuts, CCtsp_lpcut);
    }

    if (cuts->cliques) {
        for (i=0; i < cuts->cliqueend; i++) {
            CC_IFFREE (cuts->cliques[i].nodes, CCtsp_segment);
        }
        CC_FREE (cuts->cliques, CCtsp_lpclique);
    }
    if (cuts->dominos) {
        for (i=0; i < cuts->dominoend; i++) {
            CCtsp_free_lpdomino (&cuts->dominos[i]);
        }
        CC_FREE (cuts->dominos, CCtsp_lpdomino);
    }
    CC_IFFREE (cuts->cliquehash, int);
    CC_IFFREE (cuts->dominohash, int);
}

static void init_cutproof (cutproof *p)
{
    if (p) {
        p->ncount = 0;
        p->bbcount = 0;
        p->tour = (int *) NULL;
        p->probname = (char *) NULL;
        p->lproof = (proofnode *) NULL;
        p->cuts = (CCtsp_lpcuts *) NULL;
        p->lpbound = CCbigguy_ZERO ;
    }
}

static void free_cutproof (cutproof *p)
{
    int i;

    if (p) {
        CC_IFFREE (p->tour, int);
        CC_IFFREE (p->probname, char);
        if (p->lproof) {
            for (i = 0; i < p->bbcount; i++) {
                free_proofnode (&p->lproof[i]);
            }
            CC_FREE (p->lproof, proofnode);
        }
        if (p->cuts) free_lpcuts (p->cuts);
        CC_IFFREE (p->cuts, CCtsp_lpcuts);
        init_cutproof (p);
    }
}


