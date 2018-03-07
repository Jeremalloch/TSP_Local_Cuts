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
/*          Check the Euclidean distances for proximity to 0.5              */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: July 20, 2007                                                     */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *tspfname  = (char *) NULL;
static int silent = 1;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, j, ncount, norm, rval = 0;
    double val, t, ti, f, d, t1, t2, xi, yi;
    double *x, *y;
    CCdatagroup dat;
    double maxx = 0.0, minx = CCutil_MAXDOUBLE;

    CCutil_init_datagroup (&dat);

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!tspfname) {
        usage (av[0]);
        return 1;
    }

    rval = CCutil_gettsplib (tspfname, &ncount, &dat);
    CCcheck_rval (rval, "CCutil_gettsplib failed");

    CCutil_dat_getnorm (&dat, &norm);
    if (norm != CC_EUCLIDEAN_CEIL) {
        fprintf (stderr, "TSPLIB instance is not Ceil-Euclidean\n");
        rval = 1; goto CLEANUP;
    }

    x = dat.x;
    y = dat.y;
    val = 1.0;
    for (i = 0; i < ncount; i++) {
        xi = x[i];
        yi = y[i];
        if (xi > maxx) maxx = xi;
        if (yi > maxx) maxx = yi;
        if (xi < minx) minx = xi;
        if (yi < minx) minx = yi;
        for (j = i+1; j < ncount; j++) {
            if (x[j] > maxx) maxx = x[j];
            if (y[j] > maxx) maxx = x[j];
            if (x[j] < minx) minx = x[j];
            if (y[j] < minx) minx = x[j];
            t1 = xi - x[j];
            t2 = yi - y[j];          
            t  = sqrt (t1 * t1 + t2 * t2);
            ti = (int) (t);
            f = (double) ti;
            if (f != t) {
                f = t - ti;
                if (f > 0.5) d = 1.0 -f;
                else         d = f;
                if (d < val)  {
                    val = d;
                    printf ("%d %d: %.16f\n", i, j, val); 
                    fflush (stdout);
                }
            } else {
                f = f * f;
                d = t1 * t1 + t2 * t2;
                if (f != d) {
                    printf ("bad zero: %d %d: %.0f %.0f\n", i, j, t1, t2); 
                    fflush (stdout);
                }
            } 
        }
    }

    printf ("Min coord: %f\n", minx); fflush (stdout); 
    printf ("Max coord: %f\n", maxx); fflush (stdout); 


CLEANUP:

    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "v", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'v':
            silent = 0;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        tspfname = av[boptind++];
    } else {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] TSPLIB_file\n", fname);
    fprintf (stderr, "   -v    turn on more output\n");
}
