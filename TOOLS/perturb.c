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
/*                TOOL TO CALL CUTTING ROUTINES FOR AN X-VECTOR             */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: May 15, 2013                                                      */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

static char *edgefname = (char *) NULL;
static char *tourfname = (char *) NULL;
static int seed = 0;

int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int rval = 0, i, ncount, ecount;
    int *ptour = (int *) NULL, *elist = (int *) NULL, *elen = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    double szeit;
    char *probname = (char *) NULL, buf[256]; 
    
    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!tourfname) {
        fprintf (stderr, "Must specify a permutation file\n");
        usage (av[0]); goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    probname = CCtsp_problabel (edgefname);
    printf ("probname = %s\n", probname); fflush (stdout);

    rval = CCutil_getedgelist_n (&ncount, edgefname, &ecount, &elist, &elen, 0);
    CCcheck_rval (rval, "CCutil_getedgelist_n failed");

    ptour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (ptour, "out of memory for ptour");
    rval = CCutil_getcycle (ncount, tourfname, ptour, 0);
    CCcheck_rval (rval, "CCutil_getcycle failed");

    printf ("ncount = %d, ecount = %d\n", ncount, ecount);
    fflush (stdout);

#if 0
    /* For E100k with long edge lengths */
    for (i = 0; i < ecount; i++) {
        elen[i] += ((CCutil_lprand (&rstate) % 201) - 100);
    }
#endif

    for (i = 0; i < ecount; i++) {
        if (CCutil_lprand (&rstate) % 2 == 1) {
            (elen[i])++;
        } else {
            (elen[i])--;
        }
    }

    rval = CCutil_dat_setnorm (&dat, CC_SPARSE);
    CCcheck_rval (rval, "CCutil_dat_setnorm failed");
    rval = CCutil_build_sparse_dat (ncount, ecount, elist, elen, &dat,
                                    1000000);
    CCcheck_rval (rval, "CCutil_build_sparse_dat failed");

    rval = CCutil_datagroup_perm (ncount, &dat, ptour);
    CCcheck_rval (rval, "CCutil_datagroup_perm");
    
    sprintf (buf, "%s_perturb_%d.edg", probname, seed);
    rval = CCutil_writeedges_int (ncount, buf, ecount, elist, elen, 0);
    CCcheck_rval (rval, "CCutil_writeedges_int failed");

    sprintf (buf, "%s_perturb_%d.mas", probname, seed);
    rval = CCutil_putmaster (buf, ncount, &dat, ptour);
    CCcheck_rval (rval, "CCutil_putmaster failed");


CLEANUP:
    CC_IFFREE (ptour, int);
    CCutil_freedatagroup (&dat);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "s:t:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 's':
            seed = atoi(boptarg);
            break;
        case 't':
            tourfname  = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        edgefname = av[boptind++];
    } else{
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] edge_file\n", fname);
    fprintf (stderr, "   -s #  seed\n");
    fprintf (stderr, "   -t f  specify a permutation file (required)\n");
}

