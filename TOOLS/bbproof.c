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
/*               A PROGRAM TO VERIFY PROOF INFORMATION                      */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 1, 2007                                                     */
/*                                                                          */
/*  SEE short describtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "lp.h"
#include "verify.h"

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
static char *masterfname = (char *) NULL;
static char *ntourfname = (char *) NULL;
static char *opttourfname = (char *) NULL;
static char *rootfname = (char *) NULL;
static char *cutfname = (char *) NULL;
static char *tsplibfname = (char *) NULL;
static double initial_upbound = -1.0;
static double initial_lowerbound = -1.0;
static int run_silent = 1;
static int tracenode = -1;
static int writeshort = -1;
static int readshort = -1;
static int verifycuts = 1;

static int
    checkcuts (CCtsp_lpcuts *pcuts, int ncount, char *cutfile,
        cuttype **cutlist),
    check_cutproof (CCtsp_lpcuts *pcuts, int ind, int ncount,
        cutproof **cproof),
    cut_compare (CCtsp_lpcut_in *cut1, CCtsp_lpcut_in *cut2),
    checktree (CCtsp_lpcuts *cuts, proofnode *lproof, int ncount,
        CCdatagroup *dat, int ecount, int *elist, int fixcount, int *fixlist,
        double upbound, int n, int side, int depth),
    add_branch (CCtsp_lpcuts *cuts, CCtsp_branchobj *b, int side, int depth),
    add_cut (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *d),
    remove_branch (CCtsp_lpcuts *cuts, CCtsp_branchobj *b, int side),
    read_proof (char *fname, int *pncount, int *pbbcount, proofnode **plproof),
    read_lp (CCtsp_lp **lp, char *probname, char *fname, int ncount,
        CCdatagroup *dat, int *ptour, int silent),
    read_probfile (CCtsp_lp *lp, char *fname, char *probloc, int *ncount,
        int silent),
    write_perm_probfile (CCtsp_lp *lp, int ncount, int *perm, char *fname),
    write_perm_proof (CCtsp_lp *lp, char *probname, int ncount, int bbcount,
        proofnode *lproof, int *ptour, char *fname),
    write_short_proof (CCtsp_lpcuts *cuts, char *probname, int ncount,
        int bbcount, proofnode *lproof, int *tour, CCbigguy lpbound,
        cuttype *cutlist),
    write_short_proof_work (CC_SFILE *f, CCtsp_lpcuts *cuts,
        char *probname, int ncount, int bbcount, proofnode *lproof, int *tour,
        CCbigguy lpbound),
    read_short_proof (cutproof **cproof, char *fname),
    parseargs (int ac, char **av);

static void
    init_proofnode (proofnode *p),
    free_proofnode (proofnode *p),
    free_lpcuts (CCtsp_lpcuts *cuts),
    print_branch (CCtsp_lpclique *c, char sense, int rhs, int depth),
    init_cuttype (cuttype *c), 
    init_cutproof (cutproof *p),
    free_cutproof (cutproof *p),
    usage (char *fname);

int CCbbproof_price (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCdatagroup *dat, int ncount, int *elist, int ecount, int *fixlist,
        int fcount, CCbigguy *bound, int silent);
int CCbbproof_elim (CCtsp_lpcuts *cuts, CCtsp_bigdual *exact_dual,
        CCdatagroup *dat, double upperbound, CCbigguy exact_lowerbound,
        int ncount, int **elist, int *ecount, int **fixlist, int *fcount,
        int silent);
int BBverify_cut (CCtsp_lpcut_in *cut, int ncount, int check_types, int *type,
        int use_tsp, int *did_use_tsp, char *cname, CCdatagroup *getdat);
int BBverify_backbone_cut (CCtsp_lpcut_in *cut, CCtsp_lpcut_in *new,
        int *no_outside, int ncount);
int BBadd_to_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcut_in *cut, int *hit);


int main (int ac, char **av)
{
    int rval = 0;
    char *probname = (char *) NULL;
    int *ptour = (int *) NULL;
    int ncount = 0;
    int i, bbcount = 0;
    CCrandstate rstate;
    CCdatagroup dat;
    int seed = 99;
    double upbound;
    proofnode *lproof = (proofnode *) NULL;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    int *elist = (int *) NULL;
    int *fixlist = (int *) NULL;
    int ecount, fixcount;
    CCtsp_lpcuts cuts;
    CCtsp_lpcuts *pcuts = (CCtsp_lpcuts *) NULL;
    CCbigguy lpbound;
    cuttype *cutlist = (cuttype *) NULL;
    cutproof *cproof = (cutproof *) NULL;

    CCutil_init_datagroup (&dat);
    CCtsp_init_tsp_lpcuts_struct (&cuts);

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    CCutil_sprand (seed, &rstate);

    if (readshort == -1) {
        if (tsplibfname) {
            if (opttourfname == (char *) NULL) {
                fprintf (stderr, "Need optimal tour file\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            if (masterfname == (char *) NULL) {
                fprintf (stderr, "Need TSPLIB or master file\n");
                rval = 1; goto CLEANUP;
            }
        }

        if (rootfname == (char *) NULL) {
            fprintf (stderr, "Need root LP to start the proof\n");
            rval = 1; goto CLEANUP;
        }

        if (initial_upbound == -1 && opttourfname == (char *) NULL) {
            fprintf (stderr, "Need optimal tour or tour length\n");
            rval = 1; goto CLEANUP;
        } 

        if (initial_lowerbound == -1) {
            fprintf (stderr, "Need a lower bound\n");
            rval = 1; goto CLEANUP;
        }

        if (tsplibfname) {
            probname = CCtsp_problabel (tsplibfname);
        } else {
            probname = CCtsp_problabel (masterfname);
        }

        rval = read_proof (prooffname, &ncount, &bbcount, &lproof);
        CCcheck_rval (rval, "read_proof failed");

        if (tsplibfname) {
            rval = CCutil_gettsplib (tsplibfname, &i, &dat); 
            CCcheck_rval (rval, "CCutil_gettsplib failed");
            if (i != ncount) {
                fprintf (stderr, "TSPLIB file not match proof file\n");
                rval = 1;  goto CLEANUP;
            }
            ptour = CC_SAFE_MALLOC (ncount, int);
            CCcheck_NULL (ptour, "out of memory for ptour");
            rval = CCutil_getcycle_tsplib (ncount, opttourfname, ptour); 
            CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");

            rval = CCutil_datagroup_perm (ncount, &dat, ptour);
            CCcheck_rval (rval, "CCutil_datagroup_perm failed")
        } else {
            rval = CCutil_getmaster (masterfname, &i, &dat, &ptour);
            CCcheck_rval (rval, "CCutil_getmaster failed");
            if (i != ncount) {
                fprintf (stderr, "master file not match proof file\n");
                rval = 1;  goto CLEANUP;
            }
        }

        if (opttourfname) {
            upbound = 0.0;
            for (i = 0; i < ncount-1; i++) {
               upbound += (double) CCutil_dat_edgelen (i, i+1, &dat);
            }
            upbound += (double) CCutil_dat_edgelen (ncount-1, 0, &dat);
        } else {
            upbound = initial_upbound;
        }

        rval = read_lp (&rootlp, probname, rootfname, ncount, &dat, ptour,
                       run_silent);
        CCcheck_rval (rval, "read_lp failed");

        lpbound = CCbigguy_dtobigguy (initial_lowerbound);
/*
        {
            CC_SFILE *f = (CC_SFILE *) NULL;
            f = CCutil_sopen ("bnd.pla8", "r");
            rval = CCbigguy_sread (f, &lpbound);
            CCcheck_rval (rval, "CCbigguy_sread failed")
            CCutil_sclose (f);
        }
*/

        pcuts = &rootlp->cuts;
    } else {
        if (tsplibfname == (char *) NULL) {
            fprintf (stderr, "Need TSPLIB file\n");
            rval = 1; goto CLEANUP;
        }

        rval = read_short_proof (&cproof, prooffname);
        CCcheck_rval (rval, "read_short_proof failed");

        pcuts = cproof->cuts;
        probname = cproof->probname;
        ncount = cproof->ncount;
        bbcount = cproof->bbcount;
        ptour = cproof->tour;
        lproof = cproof->lproof;
        lpbound = cproof->lpbound;

        rval = CCutil_gettsplib (tsplibfname, &i, &dat); 
        CCcheck_rval (rval, "CCutil_gettsplib failed");
        if (i != ncount) {
            fprintf (stderr, "TSPLIB file not match proof file\n");
            rval = 1;  goto CLEANUP;
        }

        upbound = 0.0;
        for (i = 0; i < ncount-1; i++) {
           upbound += (double) CCutil_dat_edgelen (ptour[i], ptour[i+1], &dat);
        }
        upbound += (double) CCutil_dat_edgelen (ptour[ncount-1], ptour[0],
                                               &dat);
        printf ("Tour Value: %.0f\n", upbound);

        rval = CCutil_datagroup_perm (ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_datagroup_perm failed")
    }

    if (ntourfname) {
        rval = write_perm_probfile (rootlp, ncount, ptour, ntourfname);
        CCcheck_rval (rval, "write_perm_propfile failed");
        rval = write_perm_proof (rootlp, probname, ncount, bbcount, lproof,
                                 ptour, ntourfname);
        CCcheck_rval (rval, "write_perm_proof failed");
        printf ("Wrote permuted proof\n"); fflush (stdout);
        goto CLEANUP;
    }

    if (verifycuts == 1) {
        rval = checkcuts (pcuts, ncount, cutfname, &cutlist);
        CCcheck_rval (rval, "checkcuts failed");
    }

    if (writeshort == 1) {
        rval = write_short_proof (pcuts, probname, ncount, bbcount,
                                  lproof, ptour, lpbound, cutlist);
        CCcheck_rval (rval, "write_short_proof failed");
        printf ("Wrote short proof\n"); fflush (stdout);
        goto CLEANUP;
    }


    {
        double szeit;
        CCbigguy b;

        szeit = CCutil_zeit ();

        rval = CCbbproof_elim (pcuts, lproof[0].exact_dual, &dat, 
                    upbound, lpbound, ncount, &elist, &ecount,
                    &fixlist, &fixcount, 0);
        CCcheck_rval (rval, "CCbbproof_elim failed");

        rval = CCbbproof_price (pcuts, lproof[0].exact_dual, &dat,
                     ncount, elist, ecount, fixlist, fixcount, &b, 0);
        CCcheck_rval (rval, "CCbbproof_price failed");

        printf ("Sparse exact bound: %f\n", CCbigguy_bigguytod (b));
        fflush (stdout);
    }
    

    rval = checktree (pcuts, lproof, ncount, &dat, ecount, elist,
                      fixcount, fixlist, upbound, 0, -1, 0);
    CCcheck_rval (rval, "checktree failed");

CLEANUP:

    CC_IFFREE (ptour, int);
    CC_IFFREE (probname, char);
    CC_IFFREE (cutlist, cuttype);  /* Should free small proofs */
    CCutil_freedatagroup (&dat);
    if (lproof) {
        for (i = 0; i < bbcount; i++) {
            free_proofnode (&lproof[i]);
        }
        CC_FREE (lproof, proofnode);
    }
    if (!rootlp) free_lpcuts (&cuts);
    return rval;
}

#define BB_NOCLASS 99999
#define BB_TYPE_PROOF (CC_TYPE_SUBTOUR     | CC_TYPE_COMB  | CC_TYPE_STAR | \
                       CC_TYPE_BIPARTITION)

static int checkcuts (CCtsp_lpcuts *pcuts, int ncount, char *cutfile, 
        cuttype **cutlist)
{
    CCtsp_lpcut_in cut, back;
    int type, nsubtour = 0, ncomb = 0, nstar = 0;
    int nbipartition = 0, ndomino = 0, n2p = 0, nother = 0;
    int i, ncuts = pcuts->cutcount, hkcnt, done, missed;
    int rval = 0;
    CCtsp_lpcuts *opool = (CCtsp_lpcuts *) NULL;
    int no_outside, hit, cnt = 0; 
    int *other = (int *) NULL;
    cuttype *list = (cuttype *) NULL;
    FILE *in = (FILE *) NULL;
/*
    int hardcuts[5];
*/
    cutproof *cproof = (cutproof *) NULL;
    
    printf ("Verify Cuts\n"); fflush (stdout);

    if (cutlist) *cutlist = (cuttype *) NULL;

    if (cutfile) {
        int id;
        in = fopen (cutfile, "r");
        if (!in) {
            fprintf (stderr, "could not open %s for reading\n", cutfile);
            rval = 1; goto CLEANUP;
        }

        list = CC_SAFE_MALLOC (ncuts, cuttype);
        CCcheck_NULL (list, "out of memory for list");

        fscanf (in, "%d", &i);
        if (i != ncuts) {
            fprintf (stderr, "cut info file has wrong number of cuts\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncuts; i++) {
            fscanf (in, "%d %d %d %d", &id, &list[i].class,
                                       &list[i].isomorph, &list[i].hk);
            if (id != i) {
                fprintf (stderr, "problem in cut info file\n");
                rval = 1; goto CLEANUP;
            }
        }
        fclose (in); in = (FILE *) NULL;

        hkcnt = 0;
        for (i = 0; i < ncuts; i++) {
            if (list[i].class == BB_NOCLASS && list[i].isomorph == -1 &&
                                                     list[i].hk == 0) {
                printf ("%d\n", i); fflush (stdout);
                hkcnt++;
                rval = check_cutproof (pcuts, i, ncount, &cproof);
                CCcheck_rval (rval, "check_cutproof failed");
                list[i].proof = cproof;
            }
        }
        printf ("Need to verify %d cuts with HK\n", hkcnt); fflush (stdout);
/*
        done = missed = 0;
        for (i = 0; i < ncuts; i++) {
            if (list[i].class == BB_NOCLASS && list[i].isomorph == -1 &&
                list[i].hk == 0) {
                rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[i]), &cut);
                CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
                rval = BBverify_cut (&cut, ncount, CC_TYPE_OTHER, &type, 0,
                                    (int *) NULL, (char *) NULL),
                                    (CCdatagroup *) NULL);
                if (rval)  {
                    list[i].hk = 0;
                    printf ("Failed cut %d\n", i); fflush (stdout);
                    missed++;
                } else {
                    list[i].hk = 1;
                }
                done++;
                printf ("%d done, %d failures\n", done, missed);
                fflush (stdout);
                CCtsp_free_lpcut_in (&cut);
            }
        }
*/
        if (cutlist) *cutlist = list;
        goto CLEANUP;
    }

    list = CC_SAFE_MALLOC (ncuts, cuttype);
    CCcheck_NULL (list, "out of memory for list");
    for (i = 0; i < ncuts; i++) {
        init_cuttype (&list[i]);
    }

    for (i = 0; i < ncuts; i++) {
        rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[i]), &cut);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

        rval = BBverify_cut (&cut, ncount, BB_TYPE_PROOF, &type, 1,
                     (int *) NULL, (char *) NULL, (CCdatagroup *) NULL);
        if (rval) {
            list[i].class = BB_NOCLASS;
            nother++;
        } else {
            switch (type) {
            case CC_TYPE_SUBTOUR:
                list[i].class = CC_TYPE_SUBTOUR;
                nsubtour++;  break;
            case CC_TYPE_COMB:
                list[i].class = CC_TYPE_COMB;
                ncomb++;  break;
            case CC_TYPE_STAR:
                list[i].class = CC_TYPE_STAR;
                nstar++;  break;
            case CC_TYPE_BIPARTITION:
                list[i].class = CC_TYPE_BIPARTITION;
                nbipartition++;  break;
            case CC_TYPE_DOMINO:
                list[i].class = CC_TYPE_DOMINO;
                ndomino++;  break;
            case CC_TYPE_2P:
                list[i].class = CC_TYPE_2P;
                n2p++;  break;
            default:
                printf ("UNKOWN CUT TYPE: %d\n", type); 
                fflush (stdout);
            }
        }
        CCtsp_free_lpcut_in (&cut);
    }
    printf ("Distribution of cuts\n");
    printf ("    %4d subtours\n", nsubtour);
    printf ("    %4d combs\n", ncomb);
    printf ("    %4d stars\n", nstar);
    printf ("    %4d bipartition inequalities\n", nbipartition);
    printf ("    %4d domino parity inequalities\n", ndomino);
    printf ("    %4d 2-domino parity inequalities\n", n2p);
    printf ("    %4d other cuts\n", nother);
    fflush (stdout);

    if (nother == 0) goto DONE;

    other = CC_SAFE_MALLOC (nother, int);
    CCcheck_NULL (other, "out of memory for other");

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &opool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    opool->cuts = CC_SAFE_MALLOC (nother+1, CCtsp_lpcut);
    CCcheck_NULL (opool->cuts, "out of memory for opool.cuts");
    opool->cutspace = nother+1;

    for (i = 0; i < ncuts; i++) {
        if (list[i].class == BB_NOCLASS) {
            rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[i]), &cut);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            rval = BBverify_backbone_cut (&cut, &back, &no_outside, ncount);
            CCcheck_rval (rval, "CCverify_backbone_cut failed");

            if (no_outside) {
                fprintf (stderr, "backbone cut has no outside atom\n");
                rval = 1; goto CLEANUP;
            } else {
                rval = BBadd_to_cutpool (opool, &back, &hit);
                CCcheck_rval (rval, "BBadd_to_cutpool failed");
                if (hit == -1) {
                    other[cnt++] = i;
                } else {
                    list[i].isomorph = other[hit];  
                    {
                        CCtsp_lpcut_in bro, sis, bbro, bsis;
                        rval = CCtsp_lpcut_to_lpcut_in (pcuts,
                                       &(pcuts->cuts[other[hit]]), &bro);
                        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed\n");
                        rval = CCtsp_lpcut_to_lpcut_in (pcuts,
                                   &(pcuts->cuts[i]), &sis);
                        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed\n");
                        rval = BBverify_backbone_cut (&bro, &bbro,
                                               &no_outside, ncount);
                        CCcheck_rval (rval,
                                 "BBverify_backbone_cut failed");
                        rval = BBverify_backbone_cut (&sis, &bsis,
                                               &no_outside, ncount);
                        CCcheck_rval (rval, "BBverify_backbone_cut failed");
                        if  (cut_compare (&bbro, &bsis)) {
                            fprintf (stderr, "cuts not isomorphs\n");
                            rval = 1;  goto CLEANUP;
                        }
                        CCtsp_free_lpcut_in (&bro);
                        CCtsp_free_lpcut_in (&sis);
                        CCtsp_free_lpcut_in (&bbro);
                        CCtsp_free_lpcut_in (&bsis);
                    }
                }
            }
            CCtsp_free_lpcut_in (&back);
            CCtsp_free_lpcut_in (&cut);
        }
    }

    printf ("Number of non-isomorphic other cuts: %d\n", opool->cutcount);
    fflush (stdout);
    if (cnt != opool->cutcount) {
        fprintf (stderr, "lost a cut in pool\n");
        rval = 1;  goto CLEANUP;
    }

/*
    hardcuts[0] = other[47];
    hardcuts[1] = other[67];
    hardcuts[2] = other[68];
    hardcuts[3] = other[88];
    hardcuts[4] = other[518];

    for (i = 0; i < 5; i++) printf ("%d ", hardcuts[i]);
    printf ("\n"); fflush (stdout);
*/

/*
    {
        CCtsp_lpcuts *jpool = (CCtsp_lpcuts *) NULL;
        char buf[1028];

        hkcnt = 0;
        for (i = 0; i < ncuts; i++) {
            if (list[i].class == BB_NOCLASS && list[i].isomorph == -1) {
                hkcnt++;
            }
        }
        printf ("Save  %d HK cuts\n", hkcnt); fflush (stdout);


        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &jpool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed");
        jpool->cuts = CC_SAFE_MALLOC (hkcnt+1, CCtsp_lpcut);
        CCcheck_NULL (jpool->cuts, "out of memory for jpool.cuts");
        jpool->cutspace = hkcnt+1;

        for (i = 0; i < ncuts; i++) {
            if (list[i].class == BB_NOCLASS && list[i].isomorph == -1) {
                rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[i]),
                                                &cut);
                CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
                rval = BBadd_to_cutpool (jpool, &cut, &hit);
                CCcheck_rval (rval, "BBadd_to_cutpool failed");
                CCtsp_free_lpcut_in (&cut);
           }
       }

        sprintf (buf, "allhk.cuts");
        rval = CCtsp_write_cutpool (ncount, buf, jpool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");

        CCtsp_free_cutpool (&jpool);
    }
*/

/*
    {
        CCtsp_lpcuts *jpool = (CCtsp_lpcuts *) NULL;
        char buf[1028];

        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &jpool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed");
        jpool->cuts = CC_SAFE_MALLOC (6, CCtsp_lpcut);
        CCcheck_NULL (jpool->cuts, "out of memory for jpool.cuts");
        jpool->cutspace = 6;

        for (i = 0; i < 5; i++) {
            rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[hardcuts[i]]),
                                            &cut);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            rval = BBadd_to_cutpool (jpool, &cut, &hit);
            CCcheck_rval (rval, "BBadd_to_cutpool failed");
            CCtsp_free_lpcut_in (&cut);
       }

        sprintf (buf, "hard.cuts");
        rval = CCtsp_write_cutpool (ncount, buf, jpool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");

        CCtsp_free_cutpool (&jpool);
    }
*/

/*
    {
        char buf[1028];
        sprintf (buf, "junk.remain");
        printf ("Saving the %d cuts to %s\n", opool->cutcount, buf);
        fflush (stdout);

        rval = CCtsp_write_cutpool (ncount, buf, opool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }
*/


    hkcnt = 0;
    for (i = 0; i < ncuts; i++) {
        if (list[i].class == BB_NOCLASS && list[i].isomorph == -1) {
            hkcnt++;
        }
    }
    printf ("Need to verify %d cuts with Held-Karp\n", hkcnt); fflush (stdout);

    done = missed = 0;
    for (i = 0; i < ncuts; i++) {
        if (list[i].class == BB_NOCLASS && list[i].isomorph == -1) {
            rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[i]), &cut);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            rval = BBverify_cut (&cut, ncount, CC_TYPE_OTHER, &type, 0,
                  (int *) NULL, (char *) NULL, (CCdatagroup *) NULL);
            if (rval)  {
                list[i].hk = 0;
                printf ("Failed cut %d\n", i); fflush (stdout);
                missed++;
            } else {
                list[i].hk = 1;
            }
            done++;
            printf ("%d done, %d failures\n", done, missed);
            fflush (stdout);
            CCtsp_free_lpcut_in (&cut);
        }
    }

/*
    {
        for (i = 0; i < ncuts; i++) {
            if (list[i].class == BB_NOCLASS && list[i].isomorph == -1) {
                list[i].hk = 1;
            }
        }

        for (i = 0; i < 5; i++) list[hardcuts[i]].hk = 0;
    }
*/


DONE:
    {
        FILE *out = (FILE *) NULL;

        out = fopen ("cut.list", "w");
        if (!out) {
            fprintf (stderr, "could not open file for writing\n");
            rval = 1; goto CLEANUP;
        }

        fprintf (out, "%d\n", ncuts);
        for (i = 0; i < ncuts; i++) {
            fprintf (out, "%d %d %d %d\n", i, list[i].class, list[i].isomorph,
                                              list[i].hk);
        }
        fclose (out);
    }

    if (cutlist) *cutlist = list;

CLEANUP:
    CC_IFFREE (other, int);
    if (rval) {
        CC_IFFREE (list, cuttype);
    }
    if (opool) CCtsp_free_cutpool (&opool);
    if (in) fclose (in);
    return rval;
}

static int check_cutproof (CCtsp_lpcuts *pcuts, int ind, int ncount,
        cutproof **cproof)
{
    int i, type, ecount, fixcount, rval = 0;
    CCtsp_lpcut_in cut;
    CCdatagroup vdat;
    char buf[1028];
    CCtsp_lpcuts vcuts;
    CCbigguy b;
    char *vprobname = (char *) NULL;
    int *vtour = (int *) NULL, *elist = (int *) NULL, *fixlist = (int *) NULL;
    double upbound;
    cutproof *p = (cutproof *) NULL;

    CCtsp_init_tsp_lpcuts_struct (&vcuts);
    CCutil_init_datagroup (&vdat);
    CCtsp_init_lpcut_in (&cut);

    rval = CCtsp_lpcut_to_lpcut_in (pcuts, &(pcuts->cuts[ind]), &cut);
    CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
    rval = BBverify_cut (&cut, ncount, 0, &type, 0, (int *) NULL,
                        (char *) NULL, &vdat);
    CCcheck_rval (rval, "BBverify_cut failed for getdat");

    sprintf (buf, "v%d.sproof", ind);

    rval = read_short_proof (&p, buf);
    CCcheck_rval (rval, "read_short_proof failed");

    if (p->ncount != cut.skel.atomcount) {
        fprintf (stderr, "cut proof does not match atomcount");
        rval = 1;  goto CLEANUP;
    }

    upbound = 0.0;
    for (i = 0; i < p->ncount-1; i++) {
       upbound += (double) CCutil_dat_edgelen (p->tour[i], p->tour[i+1],
                                               &vdat);
    }
    upbound += (double) CCutil_dat_edgelen (p->tour[p->ncount-1], p->tour[0],
                                           &vdat);
    printf ("Cut tour Value: %.0f\n", upbound);

    rval = CCutil_datagroup_perm (p->ncount, &vdat, p->tour);
    CCcheck_rval (rval, "CCutil_datagroup_perm failed")

    rval = CCbbproof_elim (p->cuts, p->lproof[0].exact_dual, &vdat, upbound,
            p->lpbound, p->ncount, &elist, &ecount, &fixlist, &fixcount, 0);
    CCcheck_rval (rval, "CCbbproof_elim failed");

    rval = CCbbproof_price (p->cuts, p->lproof[0].exact_dual, &vdat, p->ncount,
                           elist, ecount, fixlist, fixcount, &b, 0);
    CCcheck_rval (rval, "CCbbproof_price failed");

    printf ("Sparse exact bound: %f\n", CCbigguy_bigguytod (b));
    fflush (stdout);

    rval = checktree (p->cuts, p->lproof, p->ncount, &vdat, ecount, elist,
                      fixcount, fixlist, upbound, 0, -1, 0);
    CCcheck_rval (rval, "checktree failed");

    if (cut.rhs == (int) upbound) {
        printf ("Verified RHS value %d for cut %d\n", cut.rhs, ind);
        fflush (stdout);
    } else {
        fprintf (stderr, "Cut RHS %d, but proof gives %f\n", cut.rhs, upbound);
        rval = 1; goto CLEANUP;
    }

    if (cproof == (cutproof **) NULL) {
        fprintf (stderr, "need to pass in cproof\n");
        rval = 1; goto CLEANUP;
    }
    *cproof = p;

CLEANUP:
    if (rval) free_cutproof (p);
    CCutil_freedatagroup (&vdat);
    CCtsp_free_lpcut_in (&cut);
    CC_IFFREE (vprobname, char);
    CC_IFFREE (vtour, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (fixlist, int);
    return rval;
}

static int cut_compare (CCtsp_lpcut_in *cut1, CCtsp_lpcut_in *cut2)
{
    int i, j;
    CCtsp_lpclique *a, *b;

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

static int checktree (CCtsp_lpcuts *cuts, proofnode *lproof,
        int ncount, CCdatagroup *dat, int ecount, int *elist, int fixcount,
        int *fixlist, double upbound, int n, int side, int depth)
{
    int rval = 0;
    int parent;
    CCbigguy bound, ubound;

    printf ("Node %d: Side=%d Child0=%d Child1=%d Parent=%d\n",
            n, lproof[n].side, lproof[n].child0, lproof[n].child1,
            lproof[n].parent);

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
                           depth);
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
                     ncount, elist, ecount, fixlist, fixcount, &bound, 0);
        CCcheck_rval (rval, "CCbbproof_price failed");

        printf ("Exactbound: %f", CCbigguy_bigguytod (bound));
        fflush (stdout);
        if (depth > 0) {
            ubound = CCbigguy_dtobigguy (upbound);
            CCbigguy_sub (&ubound, CCbigguy_ONE);
            if (CCbigguy_cmp (bound, ubound) > 0) {
                printf ("  PRUNED\n");
                fflush (stdout);
            } else {
                printf ("\nERROR: Cannot verify pruned node\n");
                fflush (stdout);
                rval = 1; goto CLEANUP;
            }
        } else {
           printf ("  Root\n"); fflush (stdout);
        }
    }

    if (lproof[n].child0 != -1) {
        rval = checktree (cuts, lproof, ncount, dat, ecount, elist, fixcount,
                  fixlist, upbound, lproof[n].child0, 0, depth+1);
        CCcheck_rval (rval, "checktree failed");
    }
    if (lproof[n].child1 != -1) {
        rval = checktree (cuts, lproof, ncount, dat, ecount, elist, fixcount,
                  fixlist, upbound, lproof[n].child1, 1, depth+1);
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
        int depth)
{
    int rval = 0;
    int ar[2];
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;
    CCtsp_lpcut_in d;

    CCtsp_init_lpcut_in (&d);

    if (b == NULL) {
        fprintf (stderr, "add_branch called without a branchobj\n");
        rval = 1; goto CLEANUP;
    }
    printf ("Add Side=%d ", side); fflush (stdout);

    c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    CCcheck_NULL (c, "out of memory for c");

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        rval = CCtsp_array_to_lpclique (ar, 2, c);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        if (side == 1) {
            d.rhs = 2.0;
            d.sense = 'L';
        } else {
            d.rhs = 4.0;
            d.sense = 'G';
        }
    } else {
        if (!b->clique) {
            fprintf (stderr, "CCtsp_branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_copy_lpclique (b->clique, c);
        CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
        if (side == 0) {
            d.rhs = 2.0;
            d.sense = 'L';
        } else {
            d.rhs = 4.0;
            d.sense = 'G';
        }
    }

    d.cliquecount = 1;
    d.branch = 1;
    d.cliques = c;

    print_branch (c, d.sense, d.rhs, depth);

    rval = add_cut (cuts, &d);
    CCcheck_rval (rval, "add_cut failed");

    CCtsp_free_lpcut_in (&d);

CLEANUP:

    return rval;
}

static int add_cut (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *d)
{
    int rval = 0;
    CCtsp_lpcut new;

    CCtsp_init_lpcut (&new);

    new.rhs         = d->rhs;
    new.sense       = d->sense;
    new.branch      = d->branch;
    rval = CCtsp_register_cliques (cuts, d, &new);
    CCcheck_rval (rval, "CCtsp_register_cliques failed");
    rval = CCtsp_register_dominos (cuts, d, &new);
    CCcheck_rval (rval, "CCtsp_register_dominos failed");

    /* No need to build a skeleton */

    if (cuts->cutcount >= cuts->cutspace) {
        void *tmp_ptr = (void *) cuts->cuts;
        if (CCutil_reallocrus_scale (&tmp_ptr, &cuts->cutspace,
                cuts->cutcount + 1, 1.3, sizeof (CCtsp_lpcut))) {
            return -1;
        }
        cuts->cuts = (CCtsp_lpcut *) tmp_ptr;
    }
    cuts->cuts[cuts->cutcount] = new;
    cuts->cutcount++;

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
        if (!b->clique) {
            printf ("Error: Un-branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        c = b->clique;

        if (side == 0) {
            sense = 'L';
            rhs = 2;
        } else {
            sense = 'G';
            rhs = 4;
        }
    }

    /* print_branch (c, sense, rhs, depth); */

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

    /* Remove the last cut from the LP */

    CCtsp_unregister_cliques (cuts, &(cuts->cuts[num]));
    CCtsp_unregister_dominos (cuts, &(cuts->cuts[num]));
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


static int read_proof (char *fname, int *pncount, int *pbbcount,
        proofnode **plproof)
{
    int rval = 0;
    proofnode *lproof = (proofnode *) NULL;
    int i, j, gotbranch, segcount, cutcount, ncount = 0, bbcount = 0;
    CCtsp_branchobj *b = (CCtsp_branchobj *) NULL;
    CCtsp_bigdual *d = (CCtsp_bigdual *) NULL;
    CC_SFILE *f = (CC_SFILE *) NULL;

    f = CCutil_sopen (fname, "r");
    if (!f) {
        fprintf (stderr, "could not safe open %s for reading\n", fname);
        rval = 1; goto CLEANUP;
    }  

    rval = CCutil_sread_int (f, &ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    rval = CCutil_sread_int (f, &bbcount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    lproof = CC_SAFE_MALLOC (bbcount, proofnode);
    for (i = 0; i < bbcount; i++) init_proofnode (&lproof[i]);

    for (i = 0; i < bbcount; i++) {
        rval = CCutil_sread_int (f, &lproof[i].number);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &lproof[i].parent);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &lproof[i].side);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &lproof[i].child0);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &lproof[i].child1);
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
            lproof[i].branch = b;
        } else {
            lproof[i].branch = (CCtsp_branchobj *) NULL;
        }
        rval = CCutil_sread_int (f, &cutcount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        if (cutcount) {
            d = CC_SAFE_MALLOC (1, CCtsp_bigdual);
            CCcheck_NULL (d, "out of memory for d");
            d->cutcount = cutcount;
            d->node_pi = CC_SAFE_MALLOC (ncount, CCbigguy);
            CCcheck_NULL (d->node_pi, "out of memory for node_pi");
            d->cut_pi = CC_SAFE_MALLOC (cutcount, CCbigguy);
            CCcheck_NULL (d->cut_pi, "out of memory for cut_pi");
            for (j = 0; j < ncount; j++) {
                rval = CCbigguy_sread (f, &d->node_pi[j]);
                CCcheck_rval (rval, "CCbigguy_sread failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = CCbigguy_sread (f, &d->cut_pi[j]);
                CCcheck_rval (rval, "CCbigguy_sread failed");
            }

            lproof[i].exact_dual = d;
        } else {
            lproof[i].exact_dual = (CCtsp_bigdual *) NULL;
        }

    }
    *pncount = ncount;
    *pbbcount = bbcount;
    *plproof = lproof;

CLEANUP:

    if (f) CCutil_sclose (f);
    if (rval) {
        if (lproof) {
            for (i = 0; i < bbcount; i++) {
                free_proofnode (&lproof[i]);
            }
            CC_FREE (lproof, proofnode);
        }
    }

    return rval;
}

static int write_perm_proof (CCtsp_lp *lp, char *probname, int ncount,
        int bbcount, proofnode *lproof, int *perm, char *fname)
{
    int rval = 0;
    int i, j;
    CC_SFILE *f = (CC_SFILE *) NULL;
    CCtsp_lpclique *c;
    char buf[1024];
    CCtsp_bigdual *d;
    int *nperm = (int *) NULL;
    int *invperm = (int *) NULL;
    int *ninvperm = (int *) NULL;

    nperm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for nperm");
    rval = CCutil_getcycle_tsplib (ncount, fname, nperm); 
    CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");

    ninvperm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for ninvperm");
    for (i = 0; i < ncount; i++) ninvperm[nperm[i]] = i;

    invperm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for invperm");
    for (i = 0; i < ncount; i++) invperm[perm[i]] = i;


    sprintf (buf, "%s.newproof", probname);
    f = CCutil_sopen (buf, "w");
    if (!f) {
        fprintf (stderr, "could not safe open %s for writing\n", buf);
        rval =1; goto CLEANUP;
    }  

    rval = CCutil_swrite_int (f, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (f, bbcount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < bbcount; i++) {
        rval = CCutil_swrite_int (f, lproof[i].number);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].parent);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].side);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].child0);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].child1);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        if (lproof[i].branch) {
            rval = CCutil_swrite_int (f, 1);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            if (lproof[i].branch->ends[0] != -1) {
                rval = CCutil_swrite_int (f,
                   ninvperm[perm[lproof[i].branch->ends[0]]]);
            } else {
                rval = CCutil_swrite_int (f, lproof[i].branch->ends[0]);
            }
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            if (lproof[i].branch->ends[1] != -1) {
                rval = CCutil_swrite_int (f,
                   ninvperm[perm[lproof[i].branch->ends[1]]]);
            } else {
                rval = CCutil_swrite_int (f, lproof[i].branch->ends[1]);
            }
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            if (lproof[i].branch->clique) {
                CCtsp_lpclique con;
                int *arr = (int *) NULL;
                int acount;

                CCtsp_init_lpclique (&con);

                c = lproof[i].branch->clique;
                rval = CCtsp_clique_to_array (c, &arr, &acount);
                CCcheck_rval (rval, "CCtsp_clique_to_array failed");
    
                for (j = 0; j < acount; j++) {
                    arr[j] = ninvperm[perm[arr[j]]];
                }       
        
                rval = CCtsp_array_to_lpclique (arr, acount, &con);
                CCcheck_rval (rval, "CCtsp_array_lpclique failed")


                rval = CCutil_swrite_int (f, con.segcount);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
                for (j = 0; j < con.segcount; j++) {
                    rval = CCutil_swrite_int (f, con.nodes[j].lo);
                    CCcheck_rval (rval, "CCutil_swrite_int failed");
                    rval = CCutil_swrite_int (f, con.nodes[j].hi);
                    CCcheck_rval (rval, "CCutil_swrite_int failed");
                }
                CC_IFFREE (arr, int);
                CCtsp_free_lpclique (&con);
            } else {
                rval = CCutil_swrite_int (f, 0);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
            }
        } else {
            rval = CCutil_swrite_int (f, 0);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
        if (lproof[i].exact_dual) {
            d = lproof[i].exact_dual;
            rval = CCutil_swrite_int (f, d->cutcount);
            CCcheck_rval (rval, "CCutil_swrite_int failed");

            /* Adjust node_pi for deleting sparsify mods */
            for (j = 0; j < d->cutcount; j++) {
                int k;
                CCbigguy t, x;

                x = d->cut_pi[j];
                for (k = 0; k < lp->cuts.cuts[j].modcount; k++) {
                    t = d->node_pi[lp->cuts.cuts[j].mods[k].node];
                    CCbigguy_addmult (&t, x,
                        (((int) lp->cuts.cuts[j].mods[k].mult) - 128));
                    d->node_pi[lp->cuts.cuts[j].mods[k].node] = t;
                }
            }

            for (j = 0; j < ncount; j++) {
                rval = CCbigguy_swrite (f, d->node_pi[invperm[nperm[j]]]);
                CCcheck_rval (rval, "CCbigguy_swrite failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = CCbigguy_swrite (f, d->cut_pi[j]);
                CCcheck_rval (rval, "CCbigguy_swrite failed");
            }
        } else {
            rval = CCutil_swrite_int (f, 0);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
    }

CLEANUP:

    if (f) CCutil_sclose (f);
    CC_IFFREE (nperm, int);
    CC_IFFREE (invperm, int);
    CC_IFFREE (ninvperm, int);

    return rval;
}

static int write_short_proof (CCtsp_lpcuts *cuts, char *probname, int ncount,
        int bbcount, proofnode *lproof, int *tour, CCbigguy lpbound,
        cuttype *cutlist)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    char buf[1024];
    int i, hkcnt, rval = 0;
    cutproof *p;

    sprintf (buf, "%s.shortproof", probname);
    f = CCutil_sopen (buf, "w");
    if (!f) {
        fprintf (stderr, "could not safe open %s for writing\n", buf);
        rval = 1; goto CLEANUP;
    }  

    rval = write_short_proof_work (f, cuts, probname, ncount, bbcount, lproof,
                                   tour, lpbound);
    CCcheck_rval (rval, "write_short_proof work failed");

    if (cutlist) {
        rval = CCutil_swrite_int (f, cuts->cutcount);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        hkcnt = 0;
        for (i = 0; i < cuts->cutcount; i++) {
            if (cutlist[i].class == BB_NOCLASS && cutlist[i].isomorph == -1 &&
                                                  cutlist[i].hk == 0) {
                hkcnt++;
            }
            rval = CCutil_swrite_int (f, cutlist[i].class);
            CCcheck_rval (rval, "CCutil_swrite_int failed (class)");
            rval = CCutil_swrite_int (f, cutlist[i].isomorph);
            CCcheck_rval (rval, "CCutil_swrite_int failed (isomorph)");
            rval = CCutil_swrite_int (f, cutlist[i].hk);
            CCcheck_rval (rval, "CCutil_swrite_int failed (hk))");
        }

        printf ("Have %d cut proofs to write\n", hkcnt); fflush (stdout);
        rval = CCutil_swrite_int (f, hkcnt);
        CCcheck_rval (rval, "CCutil_swrite_int failed");

        for (i = 0; i < cuts->cutcount; i++) {
            if (cutlist[i].class == BB_NOCLASS && cutlist[i].isomorph == -1 &&
                                                  cutlist[i].hk == 0) {
                rval = CCutil_swrite_int (f, i);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
                if (cutlist[i].proof == (cutproof *) NULL) {
                    fprintf (stderr, "no proof for cut %d\n", i);
                    rval = 1; goto CLEANUP;
                }
                p = cutlist[i].proof;
                rval = write_short_proof_work (f, p->cuts, p->probname,
                        p->ncount, p->bbcount, p->lproof, p->tour, p->lpbound);
                CCcheck_rval (rval, "write_short_proof failed");
                printf ("%d (bbcount = %d, ncount = %d)\n", i,
                     cutlist[i].proof->bbcount, cutlist[i].proof->ncount);
                fflush (stdout);
            }
        }
    } else {
        rval = CCutil_swrite_int (f, 0);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

CLEANUP:
    if (f) CCutil_sclose (f);
    return rval;
}

static int write_short_proof_work (CC_SFILE *f, CCtsp_lpcuts *cuts,
        char *probname, int ncount, int bbcount, proofnode *lproof, int *tour,
        CCbigguy lpbound)
{
    int i, j,  rval = 0;
    CCtsp_lpclique *c;
    CCtsp_bigdual *d;

    /* Probname and ncount */

    for (i = 0; probname[i] != '\0' && i < CCtsp_PROB_FILE_NAME_LEN-1; i++) {
        rval = CCutil_swrite_char (f, probname[i]);
        CCcheck_rval (rval, "CCutil_swrite_char failed");
    }
    for (; i < CCtsp_PROB_FILE_NAME_LEN; i++) {
        rval = CCutil_swrite_char (f, '\0');
        CCcheck_rval (rval, "CCutil_swrite_char failed");
    }

    rval = CCutil_swrite_int (f, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    /* Cuts for full LP */

    rval = CCtsp_write_cuts (f, ncount, cuts, 0);
    CCcheck_rval (rval, "CCtsp_write_cuts failed");


    /* Optimal tour */

    for (i = 0; i < ncount; i++) {
        rval = CCutil_swrite_int (f, tour[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    /* LP bound for elimination */

    rval = CCbigguy_swrite (f, lpbound);
    CCcheck_rval (rval, "CCbigguy_swrite failed");

    /* BB tree */

    rval = CCutil_swrite_int (f, bbcount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < bbcount; i++) {
        rval = CCutil_swrite_int (f, lproof[i].number);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].parent);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].side);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].child0);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].child1);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        if (lproof[i].branch) {
            rval = CCutil_swrite_int (f, 1);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            rval = CCutil_swrite_int (f, lproof[i].branch->ends[0]);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            rval = CCutil_swrite_int (f, lproof[i].branch->ends[1]);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            if (lproof[i].branch->clique) {
                c = lproof[i].branch->clique;
                rval = CCutil_swrite_int (f, c->segcount);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
                for (j = 0; j < c->segcount; j++) {
                    rval = CCutil_swrite_int (f, c->nodes[j].lo);
                    CCcheck_rval (rval, "CCutil_swrite_int failed");
                    rval = CCutil_swrite_int (f, c->nodes[j].hi);
                    CCcheck_rval (rval, "CCutil_swrite_int failed");
                }
            } else {
                rval = CCutil_swrite_int (f, 0);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
            }
        } else {
            rval = CCutil_swrite_int (f, 0);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
        if (lproof[i].exact_dual) {
            d = lproof[i].exact_dual;
            rval = CCutil_swrite_int (f, d->cutcount);
            CCcheck_rval (rval, "CCutil_swrite_int failed");

            for (j = 0; j < ncount; j++) {
                rval = CCbigguy_swrite (f, d->node_pi[j]);
                CCcheck_rval (rval, "CCbigguy_swrite failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = CCbigguy_swrite (f, d->cut_pi[j]);
                CCcheck_rval (rval, "CCbigguy_swrite failed");
            }
        } else {
            rval = CCutil_swrite_int (f, 0);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
    }

CLEANUP:
    return rval;
}

static int read_short_proof (cutproof **cproof, char *fname)
{
    int rval = 0;
    int i, j, gotbranch, segcount, cutcount;
    CC_SFILE *f = (CC_SFILE *) NULL;
    CCtsp_bigdual *d = (CCtsp_bigdual *) NULL;
    CCtsp_branchobj *b = (CCtsp_branchobj *) NULL;
    proofnode *p;
    CCtsp_lpcuts *cuts;
    char *probname;
    int ncount, bbcount;
    proofnode *lproof;
    int *tour;
    CCbigguy lpbound;

    printf ("read_short_proof ...\n"); fflush (stdout);

    f = CCutil_sopen (fname, "r");
    if (!f) {
        fprintf (stderr, "could not safe open %s for reading\n", fname);
        rval =1; goto CLEANUP;
    }  

    /* Probname and ncount */

    probname = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    CCcheck_NULL (probname, "out of memory for probname");

    for (i = 0; i < CCtsp_PROB_FILE_NAME_LEN; i++) {
        rval = CCutil_sread_char (f, &probname[i]);
        CCcheck_rval (rval, "CCutil_sread_char failed");
    }

    rval = CCutil_sread_int (f, &ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    printf ("Prob Name: %s\n", probname);
    printf ("Number of Nodes: %d\n", ncount);
    fflush (stdout);

    printf ("CCtsp_read_cuts ..."); fflush (stdout);

    cuts = CC_SAFE_MALLOC (1, CCtsp_lpcuts);
    CCcheck_NULL (cuts, "out of memory for cuts");
    CCtsp_init_tsp_lpcuts_struct (cuts);

    rval = CCtsp_init_cliquehash (cuts, ncount);
    CCcheck_rval (rval, "CCtsp_init_cliquehash failed");

    rval = CCtsp_init_dominohash (cuts, ncount);
    CCcheck_rval (rval, "CCtsp_init_dominohash failed");

    rval = CCtsp_read_cuts (f, &j, cuts, 0, 0);
    CCcheck_rval (rval, "CCtsp_write_cuts failed");

    printf ("Read %d cuts\n", cuts->cutcount); fflush (stdout);

    tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (tour, "out of memory for tour");

    for (i = 0; i < ncount; i++) {
        rval = CCutil_sread_int (f, &tour[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    printf ("Read the tour\n"); fflush (stdout);

    rval = CCbigguy_sread (f, &lpbound);
    CCcheck_rval (rval, "CCbigguy_sread failed");

    printf ("LP bound: %12f\n", CCbigguy_bigguytod (lpbound));
    fflush (stdout);

    rval = CCutil_sread_int (f, &bbcount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    lproof = CC_SAFE_MALLOC (bbcount, proofnode);
    CCcheck_NULL (lproof, "out of memory for lproof");

    p = lproof;

    for (i = 0; i < bbcount; i++) {
        rval = CCutil_sread_int (f, &p[i].number);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &p[i].parent);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &p[i].side);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &p[i].child0);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        rval = CCutil_sread_int (f, &p[i].child1);
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
            p[i].branch = b;
        } else {
            p[i].branch = (CCtsp_branchobj *) NULL;
        }
        rval = CCutil_sread_int (f, &cutcount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
        if (cutcount) {
            d = CC_SAFE_MALLOC (1, CCtsp_bigdual);
            CCcheck_NULL (d, "out of memory for d");
            d->cutcount = cutcount;
            d->node_pi = CC_SAFE_MALLOC (ncount, CCbigguy);
            CCcheck_NULL (d->node_pi, "out of memory for node_pi");
            d->cut_pi = CC_SAFE_MALLOC (cutcount, CCbigguy);
            CCcheck_NULL (d->cut_pi, "out of memory for cut_pi");
            for (j = 0; j < ncount; j++) {
                rval = CCbigguy_sread (f, &d->node_pi[j]);
                CCcheck_rval (rval, "CCbigguy_sread failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = CCbigguy_sread (f, &d->cut_pi[j]);
                CCcheck_rval (rval, "CCbigguy_sread failed");
            }

            p[i].exact_dual = d;
        } else {
            p[i].exact_dual = (CCtsp_bigdual *) NULL;
        }
    }

    printf ("BBcount: %d\n", bbcount); fflush (stdout);

    *cproof = CC_SAFE_MALLOC (1, cutproof);
    CCcheck_NULL ((*cproof), "out of memory for proof");
    (*cproof)->lproof = lproof;
    (*cproof)->ncount = ncount;
    (*cproof)->bbcount = bbcount;
    (*cproof)->tour = tour;
    (*cproof)->probname = probname;
    (*cproof)->lproof = lproof;
    (*cproof)->cuts = cuts;
    (*cproof)->lpbound = lpbound;

CLEANUP:

    if (f) CCutil_sclose (f);
    return rval;
}


static int read_lp (CCtsp_lp **lp, char *probname, char *fname, int ncount,
        CCdatagroup *dat, int *ptour, int silent)
{
    int rval = 0;

    *lp = CC_SAFE_MALLOC (1, CCtsp_lp);
    if ( !(*lp) ) {
       rval = 1; goto CLEANUP;
    }

    CCtsp_init_tsp_lp_struct (*lp);
    (*lp)->perm = ptour;
    (*lp)->dat = dat;

    rval = read_probfile (*lp, fname, probname, &ncount, silent);
    CCcheck_rval (rval, "read_probfile failed");

    CCtsp_free_bigdual (&((*lp)->exact_dual));

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "C:D:I:l:M:R:T:t:u:vwx?", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'C':
            ntourfname = boptarg;
            break;
        case 'D':
            tsplibfname = boptarg;
            break;
        case 'I':
            cutfname = boptarg;
            break;
        case 'l':
            initial_lowerbound = atof (boptarg);
            break;
        case 'M':
            masterfname = boptarg;
            break;
        case 'R':
            rootfname = boptarg;
            break;
        case 'T':
            opttourfname = boptarg;
            break;
        case 't':
            tracenode = atoi (boptarg);
            break;
        case 'u':
            initial_upbound = atof (boptarg);
            break;
        case 'v':
            run_silent = 0;
            break;
        case 'w':
            writeshort = 1;
            break;
        case 'x':
            readshort = 1;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }

    prooffname = av[boptind++];

    if (boptind < ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] proof_fname\n", fname);
    fprintf (stderr, "   -D f  TSPLIB data file\n");
    fprintf (stderr, "   -I f  cut info file\n");
    fprintf (stderr, "   -M f  master file (if no tour and TSPLIB files)\n");
    fprintf (stderr, "   -R f  root LP for proof\n");
    fprintf (stderr, "   -T f  tour file\n");
    fprintf (stderr, "   -u #  optimal tour length (if no tour file)\n");
    fprintf (stderr, "   -l #  lower bound (required for now)\n");
    fprintf (stderr, "   -v    verbose\n");
    fprintf (stderr, "   -w    write a short combined proof\n");
    fprintf (stderr, "   -x    short combined proof file\n");
    fprintf (stderr, "   -C f  TSPLIB tour to permute master file (build new proof)\n");
}

static int read_probfile (CCtsp_lp *lp, char *fname, char *probloc, int *ncount,
        int silent)
{
    int tncount;
    int rval = 0;
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;

    p = CCtsp_prob_read_name (fname);
    if (!p) {
        fprintf (stderr, "could not open %s for reading\n", fname);
        rval = 1;  goto CLEANUP;
    }

    lp->problabel = CCtsp_problabel (fname);
    if (lp->problabel == (char *) NULL) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        CCtsp_prob_rclose (p);
        rval = 1;  goto CLEANUP;
    }

    if (probloc == (char *) NULL) {
        lp->probloc = CCutil_strdup (lp->problabel);
    } else {
        lp->probloc = CCutil_strdup (probloc);
    }

    rval = CCtsp_prob_getnnodes (p, &tncount);
    if (rval == -1) goto CLEANUP;
    if (rval == 0) {
        if (ncount != (int *) NULL && *ncount != 0 && tncount != *ncount) {
            fprintf (stderr, "node counts differ in probfile and input\n");
            rval = 1; goto CLEANUP;
        }
        if (ncount != (int *) NULL && *ncount == 0) {
            *ncount = tncount;
        }
    } else {
        if (ncount == (int *) NULL || *ncount == 0) {
            fprintf (stderr, "node count not present in probfile or input\n");
            rval = 1; goto CLEANUP;
        } else {
            tncount = *ncount;
        }
    }

    rval = CCtsp_init_cliquehash (&lp->cuts, 2*tncount);
    if (rval) return rval;

    rval = CCtsp_init_dominohash (&lp->cuts, 2*tncount);
    if (rval) return rval;

    lp->problabel = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    if (!lp->problabel) {
        fprintf (stderr, "out of memory in read_probfile\n");
        rval = 1;
        goto CLEANUP;
    }
    rval = CCtsp_prob_getname (p, lp->problabel);
    if (rval == -1) goto CLEANUP;
    if (!silent) {
        printf ("Prob Name: %s\n", lp->problabel); fflush (stdout);
    }

    rval = CCtsp_prob_getcuts (p, &tncount, &(lp->cuts), silent);
    if (rval == -1) goto CLEANUP;

    lp->graph.ncount = tncount;

    rval = CCtsp_prob_getfixed (p, tncount, &(lp->nfixededges),
            &(lp->fixededges), silent);
    if (rval == -1) goto CLEANUP;

    printf ("Fixed edges in LP file: %d\n", lp->nfixededges);

    rval = CCtsp_prob_getfulladj (p, tncount, &(lp->fullcount),
                             &(lp->fulladj), &(lp->fulladjspace), silent);
    if (rval == -1) {
        fprintf (stderr, "CCtsp_prob_getfulladj failed\n");
        goto CLEANUP;
    }
    if (!rval) {
        if (!silent) {
            printf ("Read LP full adj\n"); fflush (stdout);
        }
        if (lp->fullcount) {
            lp->full_edges_valid = 1;
        }
    }
    printf ("Remaining edges in LP file: %d\n", lp->fullcount);
    fflush (stdout);

    rval = 0;

CLEANUP:

    if (CCtsp_prob_rclose (p)) {
        fprintf (stderr, "CCtsp_prob_rclose failed\n");
        return 1;
    }

    if (!silent) {
        printf ("Done with read_probfile\n"); fflush (stdout);
    }
    return rval;
}


static int write_perm_probfile (CCtsp_lp *lp, int ncount, int *perm,
        char *fname)
{
    CCtsp_PROB_FILE *q = (CCtsp_PROB_FILE *) NULL;
    char nambuf[1024];
    int *nperm = (int *) NULL;
    int *ninvperm = (int *) NULL;
    int *invperm = (int *) NULL;
    int *globalnames = (int *) NULL;
    int i;
    int rval = 0;
    CCtsp_lp *newlp = (CCtsp_lp *) NULL;
    CCtsp_lpcut_in cin, con;

    nperm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for nperm");
    rval = CCutil_getcycle_tsplib (ncount, fname, nperm); 
    CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");

    ninvperm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for ninvperm");
    for (i = 0; i < ncount; i++) ninvperm[nperm[i]] = i;

    invperm = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for invperm");
    for (i = 0; i < ncount; i++) invperm[perm[i]] = i;

    globalnames = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (nperm, "out of memory for globalnames");
    for (i = 0; i < ncount; i++) globalnames[i] = i;

    newlp = CC_SAFE_MALLOC (1, CCtsp_lp);
    CCcheck_NULL (newlp, "out of memory for newlp");
    CCtsp_init_tsp_lp_struct (newlp);
    rval = CCtsp_init_cliquehash (&newlp->cuts, 2*ncount);
    CCcheck_rval (rval, "CCtsp_init_cliquehash failed");
    rval = CCtsp_init_dominohash (&newlp->cuts, 2*ncount);
    CCcheck_rval (rval, "CCtsp_init_dominohash failed");
  
    for (i = 0; i < lp->cuts.cutcount; i++) {   
        rval = CCtsp_lpcut_to_lpcut_in (&lp->cuts, &lp->cuts.cuts[i], &cin);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_int failed");
        rval = CCtsp_copy_lpcut_in_global (&cin, &con, globalnames, perm,
                                          ninvperm, ncount);
        CCcheck_rval (rval, "CCtsp_copy_lpcut_in_global failed");
        rval = CCtsp_add_cut (newlp, &con, (CCtsp_lprow *) NULL, 0);
        CCcheck_rval (rval, "CCtsp_add_cut failed");
        CCtsp_free_lpcut_in (&cin);
        CCtsp_free_lpcut_in (&con);
    }

    sprintf (nambuf, "short_%s", lp->problabel);

    q = CCtsp_prob_write (nambuf, 100);
    if (!q) {
        fprintf (stderr, "could not open %s for writing\n", nambuf);
        rval = 1;  goto CLEANUP;
    }

    if (CCtsp_prob_putname (q, lp->problabel)) {
        fprintf (stderr, "CCtsp_prob_putname failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putnnodes (q, ncount)) {
        fprintf (stderr, "CCtsp_prob_putnodes failed\n"); goto CLEANUP;
    }
    if (CCtsp_prob_putcuts (q, ncount, &(newlp->cuts))) {
        fprintf (stderr, "CCtsp_prob_putcuts failed\n"); goto CLEANUP;
    }
    if (lp->nfixededges > 0) { 
        int nbits = CCutil_sbits (ncount);
        int k;

        q->offsets.fix = CCutil_stell (q->f);
        if (CCutil_swrite_char (q->f, 1)) return 1;
        if (CCutil_swrite_int (q->f, ncount)) return 1;
        if (CCutil_swrite_int (q->f, lp->nfixededges)) return 1;
        for (i = 0; i < 2*lp->nfixededges; i++) {
            k = ninvperm[perm[lp->fixededges[i]]];
            if (CCutil_swrite_bits (q->f, k, nbits)) return 1;
        }
    }

    if (lp->fullcount > 0 && lp->full_edges_valid) {  
        int j, k, m;
        CCtsp_genadj *adj = lp->fulladj;
        CCtsp_genadjobj *adjspace = (CCtsp_genadjobj *) NULL;
        CCtsp_genadjobj *pa;
        CCtsp_genadj *nadj = (CCtsp_genadj *) NULL;


        nadj = CC_SAFE_MALLOC (ncount, CCtsp_genadj);
        CCcheck_NULL (nadj, "out of memory for nadj");
        for (i = 0; i < ncount; i++) nadj[i].deg = 0;

        for (i = 0; i < ncount; i++) {
            k = invperm[nperm[i]];
            for (j = 0; j < adj[k].deg; j++) {
                m = ninvperm[perm[adj[k].list[j].end]];
                if (m < i) nadj[m].deg++;
                else       nadj[i].deg++;
            }
        }

        k = 0;
        for (i = 0; i < ncount; i++) {
            k += nadj[i].deg;
        }
        if (k != lp->fullcount) {
            printf ("Lost some edges\n"); fflush (stdout); exit (1);
        }

        adjspace = CC_SAFE_MALLOC (lp->fullcount, CCtsp_genadjobj);
        CCcheck_NULL (adjspace, "out of memory for adjspace");

        pa = adjspace;
        for (i = 0; i < ncount; i++) {
            nadj[i].list = pa;
            pa += nadj[i].deg;
            nadj[i].deg = 0;
        }

        for (i = 0; i < ncount; i++)  {
            k = invperm[nperm[i]];
            for (j = 0; j < adj[k].deg; j++) {
                m = ninvperm[perm[adj[k].list[j].end]];
                if (m < i) {
                    nadj[m].list[nadj[m].deg].end = i;
                    nadj[m].list[nadj[m].deg].len = adj[k].list[j].len;
                    nadj[m].deg++;
                } else { 
                    nadj[i].list[nadj[i].deg].end = m;
                    nadj[i].list[nadj[i].deg].len = adj[k].list[j].len;
                    nadj[i].deg++;
                }
            }
        }

        rval = CCtsp_prob_putfulladj (q, ncount, lp->fullcount, nadj);
        CCcheck_rval (rval, "CCtsp_prob_putfulladj failed");
    }

CLEANUP:

    if (q) CCtsp_prob_wclose (q);
    CC_IFFREE (nperm, int);
    CC_IFFREE (globalnames, int);
    CC_IFFREE (ninvperm, int);
    CC_IFFREE (invperm, int);
    if (newlp) CCtsp_free_tsp_lp_struct (&newlp);

    return rval;
}

static void free_lpcuts (CCtsp_lpcuts *cuts)
{
    int i;

    if (cuts->cuts) {
        for (i=0; i < cuts->cutcount; i++) {
            CC_IFFREE (cuts->cuts[i].cliques, int);
            CC_IFFREE (cuts->cuts[i].twodom_cliques, int);
            CC_IFFREE (cuts->cuts[i].dominos, int);
            CC_IFFREE (cuts->cuts[i].twodom[0], int);
            CC_IFFREE (cuts->cuts[i].twodom[1], int);
            CC_IFFREE (cuts->cuts[i].mods, CCtsp_sparser);
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

static void init_cuttype (cuttype *c) 
{
    if (c) {
        c->class = BB_NOCLASS;
        c->isomorph = -1;
        c->hk = 0;
        c->proof = (cutproof *) NULL;
    }
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
        init_cutproof (p);
    }
}

