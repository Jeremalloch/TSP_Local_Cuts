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
/*      BUILD AN LP FILE FROM A SAV FILE (WITHOUT STAR MODIFICATIONS)       */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: November 15, 2005                                                 */
/*                                                                          */
/*  CREATES an LP file from a .sav file.  The LP is used to build GMI and   */
/*  other general IP cuts.  To make it easier to verify the cuts the LP is  */
/*  created without using multiples of the degree equations to sparsify the */
/*  constraint matrix.                                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static char *fname = (char *) NULL;
static char *masterfname = (char *) NULL;
static int seed = 0;
static int run_silently = 1;
static int usefull = 0;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *f);


int main (int ac, char **av)
{
    int rval = 0;
    int i, j, k, ncount;
    CCrandstate rstate;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_PROB_FILE *p = (CCtsp_PROB_FILE *) NULL;
    char lpname[256];
    CCtsp_predge *plist = (CCtsp_predge *) NULL;
    CCdatagroup dat;
    int *perm = (int *) NULL;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!fname) {
        usage (av[0]);
        goto CLEANUP;
    }

    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    p = CCtsp_prob_read_name (fname);
    if (!p) {
        fprintf (stderr, "could not open %s for reading\n", fname);
        rval = 1;  goto CLEANUP;
    }

    rval = CCtsp_prob_getnnodes (p, &ncount);
    CCcheck_rval (rval, "CCtsp_prob_getnnodes failed");

    rval = CCtsp_prob_rclose (p);
    CCcheck_rval (rval, "CCtsp_prob_rclose failed");

    lp = CC_SAFE_MALLOC (1, CCtsp_lp);
    CCcheck_NULL (lp, "out of memory for lp");
    CCtsp_init_tsp_lp_struct (lp);

    rval = CCtsp_read_probfile (lp, fname, (char *) NULL, &ncount,
                                run_silently);
    CCcheck_rval (rval, "CCtsp_read_probfile failed");
    printf ("Prob Name: %s\n", lp->problabel); fflush (stdout);

    if (usefull) {
        int ntry = 0;
        if (lp->full_edges_valid == 0) {
            fprintf (stderr, "Error: full edge set is not labeled as valid\n");
            rval = 1;  goto CLEANUP;
        }
        if (!masterfname) {
            fprintf (stderr, "need to specify a master file for full edges\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCutil_getmaster (masterfname, &ntry, &dat, &perm);
        CCcheck_rval (rval, "CCutil_getmaster failed");
        if (ntry != ncount) {
            fprintf (stderr, "master does not match the LP\n");
            rval = 1;  goto CLEANUP;
        }
    }

    for (i = 0; i < lp->cuts.cutcount; i++) {
        lp->cuts.cuts[i].modcount = 0;
    }

    if (lp->nfixededges) {
        printf ("Restoring bounds to %d fixed edges\n", lp->nfixededges);
        fflush (stdout);

        for (i = 0; i < lp->nfixededges; i++) {
            k = CCtsp_find_edge (&(lp->graph), lp->fixededges[2*i],
                                               lp->fixededges[2*i+1]);
            if (k != -1) {
                lp->graph.edges[k].fixed = 0;
            } else {
                printf ("WARNING: File want's to fix a non-lp edge\n");
                fflush (stdout);
            }
        }
        lp->nfixededges = 0;
        CC_IFFREE (lp->fixededges, int);
    }

    rval = CClp_init (&(lp->lp));
    CCcheck_rval (rval, "CClp_init failed");

    if (lp->warmstart) {
        printf ("Free warmstart information\n");  fflush (stdout);
        CClp_free_warmstart (&(lp->warmstart));
    }

    rval = CCtsp_load_lp (lp, &rstate, run_silently);
    CCcheck_rval (rval, "CCtsp_load_lp failed");

    if (usefull) {
        int cnt = 0;
        int incnt = 0;
        CCtsp_genadj *adj = lp->fulladj;

        for (i = 0; i < ncount; i++) {
            cnt += adj[i].deg;
        }

        plist = CC_SAFE_MALLOC (cnt, CCtsp_predge);
        CCcheck_NULL (plist, "out of memory for plist");

        cnt = 0;
        for (i = 0; i < ncount; i++) {
            for (j = 0; j < adj[i].deg; j++) {
                k = CCtsp_find_edge (&(lp->graph), i, adj[i].list[j].end);
                if (k == -1) {
                    plist[cnt].ends[0] = i;
                    plist[cnt].ends[1] = adj[i].list[j].end;
                    plist[cnt].len = CCutil_dat_edgelen (i, adj[i].list[j].end,
                                                         &dat);
                    cnt++;
                } else {
                    if (lp->graph.edges[k].len !=
                        CCutil_dat_edgelen (i, adj[i].list[j].end, &dat)) {
                        fprintf (stderr, "edge length error\n");
                        rval = 1;  goto CLEANUP;
                    }
                    incnt++;
                }
            }
        }
        printf ("%d edges in LP, %d extra edges\n", incnt, cnt);
        fflush (stdout);

        if (cnt) {
            rval = CCtsp_add_vars_to_lp (lp, plist, cnt);
            CCcheck_rval (rval, "CCtsp_add_vars_to_lp failed");
        }
    }


    sprintf (lpname, "%s.lp", lp->problabel);
    printf ("Writing LP to %s\n", lpname); fflush (stdout);

    rval = CClp_dump_lp (lp->lp, (const char *) lpname);
    CCcheck_rval (rval, "CClp_dump_lp failed");

CLEANUP:

    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    CC_IFFREE (plist, CCtsp_predge);
    CC_IFFREE (perm, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "fM:s:v", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'f':
            usefull = 1;
            break;
        case 'M':
            masterfname = boptarg;
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'v':
            run_silently = 0;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        fname = av[boptind++];
    }
    if (boptind < ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-flags below] [filename]\n", f);
    fprintf (stderr, "   -f    add full edge set to LP\n");
    fprintf (stderr, "   -M f  master file (needed if -f is used)\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -v    verbose (turn on lots of messages)\n");
}
