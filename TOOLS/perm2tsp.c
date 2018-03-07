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
/*         Convert a Concorde Permutation File to TSPLIB Tour File          */
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

static char *infname  = (char *) NULL;
static char *outfname  = (char *) NULL;
static char *tspfname  = (char *) NULL;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, ncount, rval = 0;
    int *tour = (int *) NULL;
    double val;
    CCdatagroup dat;
    FILE *out = (FILE *) NULL;

    CCutil_init_datagroup (&dat);

    rval = parseargs (ac, av);
    if (rval) return 1;

    if ((!tspfname)) {
        usage (av[0]);
        return 1;
    }

    rval = CCutil_gettsplib (tspfname, &ncount, &dat);
    if (rval) {
        fprintf (stderr, "CCutil_gettsplib failed\n"); goto CLEANUP;
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        fprintf (stderr, "out of memory in main\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_getcycle (ncount, infname, tour, 0);
    if (rval) {
        fprintf (stderr, "CCutil_getcycle failed\n"); goto CLEANUP;
    }

    CCutil_cycle_len (ncount, &dat, tour, &val);
    printf ("Tour Length: %.0f\n", val); fflush (stdout);

    out = fopen (outfname, "w");
    if (!out) {
        fprintf (stderr, "could not open %s for writing\n", outfname);
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "NAME : %s\n", outfname);
    fprintf (out, "COMMENT : Length %.0f\n", val);
    fprintf (out, "TYPE : TOUR\n");
    fprintf (out, "DIMENSION : %d\n", ncount);
    fprintf (out, "TOUR_SECTION\n");
    for (i = 0; i < ncount; i++) {
        fprintf (out, "%d\n", tour[i]+1); 
    }
    fprintf (out, "-1\n");
    fprintf (out, "EOF\n");

    fclose (out);

CLEANUP:

    CC_IFFREE (tour, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "T:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'T':
            tspfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac+1) {
        infname = av[boptind++];
        outfname = av[boptind++];
    } else {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s -T tspfile in_perm_file out_tour_file\n",
                      fname);
}

