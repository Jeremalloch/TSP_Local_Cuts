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
/*                       Demo of Shrink Routines                            */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 4, 2004                                                     */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"


static int seed = 0;
char *fname = (char *) NULL;

int
    main (int ac, char **av);

static void
    usage (char *f);

static int
    parseargs (int ac, char **av),
    test_shrink (int ncount, int ecount, int *elist, double *elen);

int main (int ac, char **av)
{
    int rval = 0;
    double szeit;
    int ncount, ecount;
    int *elist = (int *) NULL;
    double *elen = (double *) NULL;
    CCrandstate rstate;

    seed = (int) CCutil_real_zeit ();
    if (parseargs (ac, av)) goto CLEANUP;
    CCutil_sprand (seed, &rstate);

    szeit = CCutil_zeit ();
    
    rval = CCutil_getedges_double (&ncount, fname, &ecount, &elist, &elen, 0);
    CCcheck_rval (rval, "CCutil_getedges_double failed");

    rval = test_shrink (ncount, ecount, elist, elen);
    CCcheck_rval (rval, "test_shrink failed");

    printf ("Total time: %0.2f seconds\n", CCutil_zeit () - szeit);
    fflush (stdout);


CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, double);
    return 0;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "s:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 's':
            seed = atoi (boptarg);
            break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }

    fname = av[boptind++];
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- below -] x_file\n", f);
    fprintf (stderr, "    s #: random seed\n");
}


static int test_shrink (int ncount, int ecount, int *elist, double *elen)
{
    int rval = 0;
    CC_SRKgraph G;
    int k;
    int s_ncount, s_ecount;
    int *s_elist = (int *) NULL;
    double *s_elen = (double *) NULL;
    
    CCcut_SRK_init_graph (&G);

    printf ("Initial Ncount = %d, Ecount = %d\n", ncount, ecount);
    fflush (stdout);

    rval = CCcut_SRK_buildgraph (&G, ncount, ecount, elist, elen);
    CCcheck_rval (rval, "CCcut_SRK_buildgraph failed");

    rval = CCcut_SRK_defluff (&G);
    CCcheck_rval (rval, "CCcut_SRK_defluff failed");

    CCcut_SRK_identify_paths_to_edges (&G, &k, 0);
    printf ("Identify paths to edges: %d (nodes remaining)\n", k);
    fflush (stdout);
    rval = CCcut_SRK_grab_edges (&G, &s_ncount, &s_ecount, &s_elist, &s_elen,
                                 (CC_SRKexpinfo *) NULL);
    CCcheck_rval (rval, "CCcut_SRK_grab_edges failed");
    printf ("Ncount = %d, Ecount = %d\n", s_ncount, s_ecount);
    fflush (stdout);
    CC_IFFREE (s_elist, int);
    CC_IFFREE (s_elen, double);

    CCcut_SRK_identify_one_triangles (&G, &k, (CC_SRKnode *) NULL, 0.001, 2.0,
                                      0);
    printf ("Identify one triangles: %d (edges shrunk)\n", k); fflush (stdout);
    rval = CCcut_SRK_grab_edges (&G, &s_ncount, &s_ecount, &s_elist, &s_elen,
                                 (CC_SRKexpinfo *) NULL);
    CCcheck_rval (rval, "CCcut_SRK_grab_edges failed");
    printf ("Ncount = %d, Ecount = %d\n", s_ncount, s_ecount);
    fflush (stdout);
    CC_IFFREE (s_elist, int);
    CC_IFFREE (s_elen, double);


    CCcut_SRK_identify_triangle_square (&G, &k, 0.001, 0);
    printf ("Identify triangles in squares: %d (triangles shrunk)\n", k);
    fflush (stdout);

    rval = CCcut_SRK_grab_edges (&G, &s_ncount, &s_ecount, &s_elist, &s_elen,
                                 (CC_SRKexpinfo *) NULL);
    CCcheck_rval (rval, "CCcut_SRK_grab_edges failed");
    printf ("Ncount = %d, Ecount = %d\n", s_ncount, s_ecount);
    fflush (stdout);
    CC_IFFREE (s_elist, int);
    CC_IFFREE (s_elen, double);


    CCcut_SRK_identify_one_square (&G, &k, 0.001, 1.0, 0);
    printf ("Identify ones in squares: %d (pairs of ones)\n", k);

    rval = CCcut_SRK_grab_edges (&G, &s_ncount, &s_ecount, &s_elist, &s_elen,
                                 (CC_SRKexpinfo *) NULL);
    CCcheck_rval (rval, "CCcut_SRK_grab_edges failed");

    printf ("Final Ncount = %d, Ecount = %d\n", s_ncount, s_ecount);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (s_elist, int);
    CC_IFFREE (s_elen, double);
    CCcut_SRK_free_graph (&G);
    return rval;
}

