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


#include "machdefs.h"
#include "tsp.h"
#include "util.h"
#include "verify.h"

static int gather_cuts = 0;
static int edge_clone = 0;
static int backbone_cuts = 0;
static int textin = 0;
static int probin = 0;
static int poolin = 0;
static char *cutname = (char *) NULL;
static int in_ncount = 0;
static int run_silently = 0;


static int
    get_text_pool (char *poolname, int nodecount, CCtsp_lpcuts **pool),
    get_problem_lp (char *probname, CCtsp_lp **lp, int *ncount,
         int silent),
    verify_pool (CCtsp_lpcuts *pool, int ncount, int just_gather,
        char *probname),
    edgeclone_pool (CCtsp_lpcuts *pool, int ncount, char *probname),
    backbone_pool (CCtsp_lpcuts *pool, int ncount, char *probname),
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int rval;
    double szeit;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    int ncount = 0;
    char *probname = (char *) NULL;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;


    rval = parseargs (ac, av);
    if (rval) return -1;

    szeit = CCutil_zeit ();

    probname = CCtsp_problabel (cutname);
    CCcheck_NULL (probname, "CCtsp_problabel failed");

    if (textin) {
        rval = get_text_pool (cutname, in_ncount, &pool);
        CCcheck_rval (rval, "get_text_pool failed");
        ncount = in_ncount;
    } else if (probin) {
        rval = get_problem_lp (cutname, &lp, &ncount, run_silently);
        CCcheck_rval (rval, "get_problem_lp failed");
        pool = &lp->cuts;
    } else if (poolin) {
        rval = CCtsp_init_cutpool (&ncount, cutname, &pool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    }

    if (edge_clone) {
        rval = edgeclone_pool (pool, ncount, probname);
        CCcheck_rval (rval, "edgeclone_pool failed");
        printf ("Completed in %.2f seconds\n", CCutil_zeit () - szeit);
        fflush (stdout);
    } else if (backbone_cuts) {
        rval = backbone_pool (pool, ncount, probname);
        CCcheck_rval (rval, "backbone_pool failed");
        printf ("Completed in %.2f seconds\n", CCutil_zeit () - szeit);
        fflush (stdout);
    } else {
        rval = verify_pool (pool, ncount, gather_cuts, probname);

        printf ("Verification completed (%s) in %.2f seconds\n",
                rval ? "unsuccessful" : "successful",
                CCutil_zeit () - szeit);
        fflush (stdout);
    }

CLEANUP:

    CC_IFFREE (probname, char);
    if (pool && !lp) CCtsp_free_cutpool (&pool);
    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int get_text_pool (char *poolname, int nodecount, CCtsp_lpcuts **pool)
{
    int rval = 0;
    int *tour = (int *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *cnext;
    int cutcount;
    int i;

    rval = CCtsp_init_cutpool (&nodecount, (char *) NULL, pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    tour = CC_SAFE_MALLOC (nodecount, int);
    CCcheck_NULL (tour, "out of memory for tour");
    for (i=0; i<nodecount; i++) tour[i] = i;
    
    rval = CCtsp_file_cuts (poolname, &cuts, &cutcount, nodecount, tour,
                           (double *) NULL);
    CCcheck_rval (rval, "CCtsp_file_cuts failed");

    while (cuts) {
        cnext = cuts->next;
        rval = CCtsp_add_to_cutpool_lpcut_in (*pool, cuts);
        CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
        CCtsp_free_lpcut_in (cuts);
        CC_FREE (cuts, CCtsp_lpcut_in);
        cuts = cnext;
    }

CLEANUP:

    CC_IFFREE (tour, int);
    while (cuts) {
        cnext = cuts->next;
        CCtsp_free_lpcut_in (cuts);
        CC_FREE (cuts, CCtsp_lpcut_in);
        cuts = cnext;
    }
    return rval;
}

static int get_problem_lp (char *probname, CCtsp_lp **lp, int *ncount,
         int silent)
{
    int rval = 0;

    *lp = CC_SAFE_MALLOC (1, CCtsp_lp);
    CCcheck_NULL (*lp, "out of memory for lp");

    CCtsp_init_tsp_lp_struct (*lp);
    rval = CCtsp_read_probfile (*lp, probname, (char *) NULL, ncount, silent);
    CCcheck_rval (rval, "CCtsp_read_probfile failed");

    rval = CCtsp_build_lpadj (&(*lp)->graph, 0, (*lp)->graph.ecount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");
    
    rval = CCtsp_add_branchhistory_to_lp (*lp);
    CCcheck_rval (rval, "CCtsp_add_branchhistory_to_lp failed");

CLEANUP:

    return rval;
}

static int verify_pool (CCtsp_lpcuts *pool, int ncount, int just_gather,
        char *probname)
{
    int rval = 0;
    int i, type;
    CCtsp_lpcut_in cut;
    int nfail = 0;
    int nsuccess = 0;
    int nsubtour = 0;
    int ncomb = 0;
    int dirtycomb = 0;
    int nstar = 0;
    int nbipartition = 0;
    int ndomino = 0;
    int n2p = 0;
    int nother = 0;
    CCtsp_lpcuts *opool = (CCtsp_lpcuts *) NULL;

    if (just_gather) {
        rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &opool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    }

    for (i=0; i<pool->cutcount; i++) {
        if (pool->cuts[i].branch == 0) {
            rval = CCtsp_lpcut_to_lpcut_in (pool, &(pool->cuts[i]), &cut);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            if (just_gather) {
                rval = CCverify_cut (&cut, ncount, CC_TYPE_SIMPLE, &type, 1,
                                    (int *) NULL, (char *) NULL);
            } else {
                rval = CCverify_cut (&cut, ncount, CC_TYPE_ALL, &type, 0,
                                    (int *) NULL, (char *) NULL);
            }
            if (rval) {
                if (just_gather) {
                    nother++;
                    rval = CCtsp_add_to_cutpool_lpcut_in (opool, &cut);
                    CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
                } else {
                    fprintf (stderr, "CCverify_cut failed\n");
                    nfail++;
                }
            } else {
                nsuccess++;
                switch (type) {
                case CC_TYPE_SUBTOUR:
                    nsubtour++;  break;
                case CC_TYPE_COMB:
                    ncomb++;  break;
                case CC_TYPE_DIRTY_COMB:
                    dirtycomb++;  break;
                case CC_TYPE_STAR:
                    nstar++;  break;
                case CC_TYPE_BIPARTITION:
                    nbipartition++;  break;
                case CC_TYPE_DOMINO:
                    ndomino++;  break;
                case CC_TYPE_2P:
                    n2p++;  break;
                case CC_TYPE_OTHER:
                    if (just_gather) {
                        fprintf (stderr, "A non-simple cut in gather\n");
                        rval = 1;  goto CLEANUP;
                    } else {
                        nother++;  break;
                    }
                default:
                    printf ("UNKOWN CUT TYPE: %d\n", type);
                }
            }
            CCtsp_free_lpcut_in (&cut);
        }
    }

    printf ("%d of %d cuts failed verification\n", nfail, nfail + nsuccess);
    printf ("Distribution of cuts\n");
    printf ("    %4d subtours\n", nsubtour);
    printf ("    %4d combs\n", ncomb);
    printf ("    %4d dirty combs\n", dirtycomb);
    printf ("    %4d stars\n", nstar);
    printf ("    %4d bipartition inequalities\n", nbipartition);
    printf ("    %4d domino parity inequalities\n", ndomino);
    printf ("    %4d 2-domino parity inequalities\n", n2p);
    printf ("    %4d other cuts\n", nother);
    fflush (stdout);

    if (just_gather) {
        if (nother) {
            char buf[1024];
            sprintf (buf, "%s.other", probname);
            printf ("Writing non-simple cuts to %s\n", buf); fflush (stdout);
            rval = CCtsp_write_cutpool (ncount, buf, opool);
            CCcheck_rval (rval, "CCtsp_write_cutpool failed");
        } else {
            printf ("No non-simple cuts\n"); fflush (stdout);
        }
    }

    if (nfail == 0) rval = 0;
    else rval = -1;

CLEANUP:

    if (opool) CCtsp_free_cutpool (&opool);
    return rval;
}

static int edgeclone_pool (CCtsp_lpcuts *pool, int ncount, char *probname)
{
    int rval = 0;
    CCtsp_lpcut_in cut, cloned;
    CCtsp_lpcuts *opool = (CCtsp_lpcuts *) NULL;
    int i, yesno, clone_cnt = 0;
    char buf[1024];

    printf ("Number of cuts: %d\n", pool->cutcount);
    fflush (stdout);

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &opool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    for (i=0; i<pool->cutcount; i++) {
        if (pool->cuts[i].branch == 0) {
            rval = CCtsp_lpcut_to_lpcut_in (pool, &(pool->cuts[i]), &cut);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
 
            rval = CCverify_edgeclone (&cut, ncount, &yesno, &cloned);
            CCcheck_rval (rval, "CCverify_edgeclone failed");

            if (yesno == 0) {
                rval = CCtsp_add_to_cutpool_lpcut_in (opool, &cut);
                CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
            } else {
                rval = CCtsp_add_to_cutpool_lpcut_in (opool, &cloned);
                CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
                clone_cnt++;
            }

            CCtsp_free_lpcut_in (&cut);
            if (yesno) CCtsp_free_lpcut_in (&cloned);
        }
    }

    printf ("Number of edge-clone reductions: %d\n", clone_cnt);
    fflush (stdout);
    printf ("Number of non-isomporhic cuts: %d\n", opool->cutcount);
    fflush (stdout);

    sprintf (buf, "%s.cloned", probname);
    printf ("Writing clone-reduced cuts to %s\n", buf); fflush (stdout);
    rval = CCtsp_write_cutpool (ncount, buf, opool);
    CCcheck_rval (rval, "CCtsp_write_cutpool failed");

CLEANUP:

    if (opool) CCtsp_free_cutpool (&opool);
    return rval;
}

static int backbone_pool (CCtsp_lpcuts *pool, int ncount, char *probname)
{
    int rval = 0;
    CCtsp_lpcut_in back, cut;
    CCtsp_lpcuts *opool = (CCtsp_lpcuts *) NULL;
    int i, no_outside, no_cnt = 0;
    char buf[1024];

    printf ("Number of cuts: %d\n", pool->cutcount);
    fflush (stdout);

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &opool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    
    for (i=0; i<pool->cutcount; i++) {
        if (pool->cuts[i].branch == 0) {
            rval = CCtsp_lpcut_to_lpcut_in (pool, &(pool->cuts[i]), &cut);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
 
            rval = CCverify_backbone_cut (&cut, &back, &no_outside, ncount);
            CCcheck_rval (rval, "CCverify_backbone_cut failed");

            if (no_outside) {
                no_cnt++;
                rval = CCtsp_add_to_cutpool_lpcut_in (opool, &cut);
                CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
            } else {
                rval = CCtsp_add_to_cutpool_lpcut_in (opool, &back);
                CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
            }

            CCtsp_free_lpcut_in (&back);
            CCtsp_free_lpcut_in (&cut);
        }
    }
    printf ("Number of cuts with no outside node: %d\n", no_cnt);
    fflush (stdout);
    printf ("Number of non-isomporhic cuts: %d\n", opool->cutcount);
    fflush (stdout);

    sprintf (buf, "%s.iso", probname);
    printf ("Writing non-isomorphic cuts to %s\n", buf); fflush (stdout);
    rval = CCtsp_write_cutpool (ncount, buf, opool);
    CCcheck_rval (rval, "CCtsp_write_cutpool failed");

CLEANUP:

    if (opool) CCtsp_free_cutpool (&opool);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "egiNPt:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'e': edge_clone = 1; break;
        case 'g': gather_cuts = 1; break;
        case 'i': backbone_cuts = 1; break;
        case 'N': probin = 1; break;
        case 'P': poolin = 1; break;
        case 't': textin = 1; in_ncount = atoi(boptarg); break;
        default: usage (av[0]); return 1;
        }
    }
    if (textin == 0 && probin == 0 && poolin == 0) {
        usage (av[0]);
        return 1;
    }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    cutname = av[boptind++];
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [NPt:] cutfile\n", fname);
    fprintf (stderr, "   -e        try edge-clone reductions\n");
    fprintf (stderr, "   -g        gather unusual cuts into a pool file\n");
    fprintf (stderr, "   -i        just run cut isomorphism\n");
    fprintf (stderr, "   -N        cutfile is problem file\n");
    fprintf (stderr, "   -P        cutfile is a cut pool\n");
    fprintf (stderr, "   -t n      cutfile is text file, n is the nodecount\n");
    fprintf (stderr, "   One of -N, -P, or -t must be specified\n");
}
