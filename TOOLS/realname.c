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
/*          Convert the permuted name of vertices to actual name            */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 26, 2006                                                    */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *masterfname  = (char *) NULL;
static char *vertexfname = (char *) NULL;

int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, k, n, ncount, rval = 0;
    int *perm = (int *) NULL;
    CCdatagroup dat;
    FILE *in = (FILE *) NULL;

    CCutil_init_datagroup (&dat);

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname) {
        fprintf (stderr, "Must specify a master file\n");
        usage (av[0]);
        return 1;
    }

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &perm);       
    CCcheck_rval (rval, "CCutil_getmaster failed");

    in = fopen (vertexfname, "r");
    if (!in) {
        fprintf (stderr, "could not open vertex file %s\n", vertexfname);
        rval = 1; goto CLEANUP;
    }

    if ((fscanf (in, "%d", &k) != 1 || k <= 0)) {
        fprintf (stderr, "vertex file in wrong format\n"); 
    }
    printf ("Cities: ");
    for (i = 0; i < k; i++) {
        if ((fscanf (in, "%d", &n) != 1 || n < 0 || n >= ncount)) {
            fprintf (stderr, "vertex file in wrong format\n"); 
        }
        printf ("%d ", perm[n]);
    }
    printf ("\n");


CLEANUP:

    if (in) fclose (in);
    CC_IFFREE (perm, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "M:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'M':
            masterfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        vertexfname = av[boptind++];
    } else {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] vertex_file\n", fname);
    fprintf (stderr, "   -M f  master file\n");
}
