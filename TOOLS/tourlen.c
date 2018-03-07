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
/*          A PROGRAM TO COMPUTE THE LENGTH OF A TOUR IN A FILE             */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 3, 1999                                                       */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *cycfname  = (char *) NULL;
static char *tspfname  = (char *) NULL;
static char *edgefname = (char *) NULL;
static int eformat = 0;
static int tsplib_in = 1;
static int simpletour = 0;
static int norm = CC_EUCLIDEAN;
static int seed = 0;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int ncount = 0, rval = 0, *tour = (int *) NULL;
    double val;
    CCdatagroup dat;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    if ((!edgefname && !tspfname) || (edgefname && tspfname)) {
        usage (av[0]);
        goto CLEANUP;
    }

    CCutil_sprand (seed, &rstate);

    if (tspfname) {
        if (tsplib_in) {
            rval = CCutil_gettsplib (tspfname, &ncount, &dat);
            CCcheck_rval (rval, "CCutil_gettsplib failed");
        } else {
            rval = CCutil_getdata (tspfname, 0, norm, &ncount, &dat, 1, 0,
                                   &rstate);
            CCcheck_rval (rval, "CCutil_getdata failed");
        }
    } else {
        rval = CCutil_getdata (edgefname, 0, CC_SPARSE, &ncount, &dat,
                               0, 0, &rstate);
        CCcheck_rval (rval, "CCutil_getdata failed");
    }
    CC_MALLOC (tour, ncount, int);

    if (eformat == 0) {
        if (simpletour) {
            rval = CCutil_getcycle (ncount, cycfname, tour, 0);
            CCcheck_rval (rval, "CCutil_getcycle failed");
        } else {
            rval = CCutil_getcycle_tsplib (ncount, cycfname, tour);
            CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");
        }
    } else {
        rval = CCutil_getcycle_edgelist (ncount, cycfname, tour, 0);
        CCcheck_rval (rval, "CCutil_getcycle_edgelist failed");
    }

    {
        int *chk = (int *) NULL;
        int i;

        CC_MALLOC (chk, ncount, int);
        for (i = 0; i < ncount; i++) chk[i] = 0;
        for (i = 0; i < ncount; i++) {
            if (chk[tour[i]] == 1) {
                fprintf (stderr, "duplicate node in tour: %d, position %d\n",
                                  tour[i], i);
                rval = 1;  goto CLEANUP;
            }
            chk[tour[i]] = 1;
        }
        CC_IFFREE (chk, int);
    }

    if (edgefname) {
        int istour;
        rval = CCutil_sparse_real_tour (ncount, &dat, tour, &istour);
        CCcheck_rval (rval, "CCutil_sparse_real_tour failed");
        if (istour == 0) {
            printf ("Tour is not contained in the sparse edge set\n");
            fflush (stdout);
            goto CLEANUP;
        }
    }

    CCutil_cycle_len (ncount, &dat, tour, &val);
    printf ("Tour Length: %.0f\n", val); fflush (stdout);

CLEANUP:
    CC_IFFREE (tour, int);
    CCutil_freedatagroup (&dat);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, boptind = 1, inorm = 0;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "eE:N:tT:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'e': eformat = 1; break;
        case 'E': edgefname = boptarg; break;
        case 'N':
            inorm = atoi(boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 3: norm = CC_EUCLIDEAN_3D; break;
            case 4: norm = CC_USER; break;
            case 5: norm = CC_ATT; break;
            case 6: norm = CC_GEOGRAPHIC; break;
            case 7: norm = CC_MATRIXNORM; break;
            case 8: norm = CC_DSJRANDNORM; break;
            case 9: norm = CC_CRYSTAL; break;
            case 10: norm = CC_SPARSE; break;
            case 11: norm = CC_RHMAP1; break;
            case 12: norm = CC_RHMAP2; break;
            case 13: norm = CC_RHMAP3; break;
            case 14: norm = CC_RHMAP4; break;
            case 15: norm = CC_RHMAP5; break;
            case 16: norm = CC_EUCTOROIDAL; break;
            case 17: norm = CC_GEOM; break;
            case 18: norm = CC_EUCLIDEAN_CEIL; break;
            case 20: norm = CC_ROAD; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case 't': simpletour = 1; break;
        case 'T': tspfname = boptarg; break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        cycfname = av[boptind++];
    } else {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] tour_file\n", fname);
    fprintf (stderr, "   -e    tour file is an edge list (default is TSPLIB)\n");
    fprintf (stderr, "   -t    tour file in concorde permutation format (default TSPLIB)\n"); 
    fprintf (stderr, "   -E f  edge file to specify lengths\n");
    fprintf (stderr, "   -T f  TSPLIB (dat) file to specify lengths\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON 20=ROAD\n");
}
