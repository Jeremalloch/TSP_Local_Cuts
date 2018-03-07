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
/*                  A PROGRAM TO INTERSECT EDGE FILES                       */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 17, 2005                                                  */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static int nfiles = 0;
static int sdifference = 0;
static int ssubtract = 0;
static char **filelist = (char **) NULL;
static char  *outfname = (char *) NULL;

int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int ncount, ecount, rval = 0;
    int *elist = (int *) NULL, *elen  = (int *) NULL;

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    if (nfiles != 2) {
        fprintf (stderr, "must specify two edge files\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_getedgelist_n (&ncount, filelist[0], &ecount, &elist,
                                 &elen, 0);
    CCcheck_rval (rval, "CCutil_getedgelist_n failed");

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    if (sdifference == 0 && ssubtract == 0) {
        rval = CCutil_edge_file_intersect (ncount, filelist[0], filelist[1],
                                &ecount, &elist, &elen);
        CCcheck_rval (rval, "CCutil_edge_file_intersect failed");
        printf ("Intersected Edge List: %d edges\n", ecount); fflush (stdout);
    } else if (ssubtract == 0) {
        rval = CCutil_edge_file_difference (ncount, filelist[0], filelist[1],
                                &ecount, &elist, &elen, 0);
        CCcheck_rval (rval, "CCutil_edge_file_difference failed");
        printf ("Symmetric Difference: %d edges\n", ecount); fflush (stdout);
    } else {
        rval = CCutil_edge_file_difference (ncount, filelist[0], filelist[1],
                                &ecount, &elist, &elen, 1);
        CCcheck_rval (rval, "CCutil_edge_file_difference failed");
        printf ("Set Subtraction: %d edges\n", ecount); fflush (stdout);
    }
      
    if (outfname) {
        rval = CCutil_writeedges_int (ncount, outfname, ecount, elist,
                                      elen, 0);
        CCcheck_rval (rval, "CCutil_writeedges_int failed");
    }

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "o:ST", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'o':
            outfname = boptarg;
            break;
        case 'S':
            sdifference = 1;
            break;
        case 'T':
            ssubtract = 1;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (ac <= boptind) {
        usage (av[0]);
        return 1;
    }

    nfiles = ac - boptind;
    filelist = &(av[boptind]);
    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] edge_file1 edge_file2\n", fname);
    fprintf (stderr, "   -o f  write intersection edge file\n");
    fprintf (stderr, "   -S    symmetric difference\n");
    fprintf (stderr, "   -T    only A - B (first set minus second set)\n");
}
