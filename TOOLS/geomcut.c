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
/*                         CALL THE SNOWCUT GENERATOR                       */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: July 8, 2016                                                      */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

static char *masterfname = (char *) NULL;
static char *rootfname = (char *) NULL;
static char *tourfname = (char *) NULL;
static int seed = 0;

int main (int ac, char **av);
static int lp_value (CCtsp_lp *lp, double *val);
static int lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x);
static int parseargs (int ac, char **av);
static void usage (char *fname);

int main (int ac, char **av)
{
    int  i, ncount, rval = 0, infeasible = 0, xcount = 0, cutcount = 0;
    int *ptour = (int *) NULL, *xlist = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dominopool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    int *snowtour = (int *) NULL, *invperm = (int *) NULL;
    double val, z, *x = (double *) NULL;
    
    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname) {
        fprintf (stderr, "Must specify a master file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    CCutil_printlabel ();
    CCutil_signal_init ();
    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    if (tourfname) {
        CC_MALLOC (snowtour, ncount, int);
        CC_MALLOC (invperm, ncount, int);
        for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;

        rval = CCutil_getcycle (ncount, tourfname, snowtour, 0);
        CCcheck_rval (rval, "CCutil_getcycle failed");
        for (i = 0; i < ncount; i++) {
            snowtour[i] = invperm[snowtour[i]];
        }
    }

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &dominopool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");

    rval = CCtsp_init_lp (&rootlp, (char *) NULL, -1, rootfname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               pool, dominopool, 0, &rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    rval = lp_value (rootlp, &val);
    CCcheck_rval (rval, "lp_value failed");
    printf ("LP Value: %f\n", val); fflush (stdout);

    rval = lp_x (rootlp, &xcount, &xlist, &x);
    CCcheck_rval (rval, "lp_x failed");

    rval = CCtsp_geom_cuts (&cuts, &cutcount, rootlp, xcount, xlist, x,
                            snowtour, &z, &rstate);
    CCcheck_rval (rval, "CCtsp_geom_cuts failed");

    printf ("Found %d Snow Cuts\n", cutcount);
    fflush (stdout);

    {
        int cut_added = 0, tighten = 1;

        CCtsp_add_cuts_to_queue (rootlp, &cuts);
        rval = CCtsp_process_cuts (rootlp, &cut_added, tighten, 0, &rstate,
                                  (double *) NULL, &infeasible);
        CCcheck_rval (rval, "CCtsp_process_cuts failed");
        if (infeasible) {
            printf ("Infeasible LP after processing cuts\n");  fflush (stdout);
            goto CLEANUP;
        }
    }

CLEANUP:
    if (pool) { CCtsp_free_cutpool (&pool); }
    if (dominopool) { CCtsp_free_cutpool (&dominopool); }
    CC_IFFREE (ptour, int);
    CCutil_freedatagroup (&dat);
    if (rootlp) CCtsp_free_tsp_lp_struct (&rootlp);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    CC_IFFREE (snowtour, int);
    CC_IFFREE (invperm, int);

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


static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "M:s:t:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'M': masterfname = boptarg; break;
        case 's': seed = atoi(boptarg); break;
        case 't': tourfname =  boptarg;  break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        rootfname = av[boptind++];
    } else{
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] savfile_root\n", fname);
    fprintf (stderr, "   -M f  specify a master file (required)\n");
    fprintf (stderr, "   -s #  seed\n");
}
