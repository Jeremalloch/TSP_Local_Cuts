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
/*           Solve small (say 50 city) ATSP by reduction to TSP             */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 20, 2015                                                    */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

static char *matrixfname = (char *) NULL;
static char *permfname = (char *) NULL;
static int seed = 0;

#define NMAX 500

int main (int ac, char **av);
static int solve_atsp (int ncount, int **dist, double *pval,
    CCrandstate *rstate, int *tour);
static int read_matrix (char *fname, int *pncount, int ***pdist);
static int test_perm (char *fname, int ncount, int **dist);
static int parseargs (int ac, char **av);
static void usage (char *fname);

int main (int ac, char **av)
{
    int  rval = 0, i, ncount, **dist = (int **) NULL, *tour = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    double val, tval, szeit;
    
    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!matrixfname) {
        fprintf (stderr, "Must specify a distance-matrix file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    rval = read_matrix (matrixfname, &ncount, &dist);
    CCcheck_rval (rval, "read_matrix failed");

    if (ncount > NMAX) {
        printf ("only set up for instances with at most %d cities\n", NMAX);
        rval = 1; goto CLEANUP;
    }

    if (permfname) {
        rval = test_perm (permfname, ncount, dist);
        CCcheck_rval (rval, "test_perm failed");
    }

    CC_MALLOC (tour, ncount, int);

    rval = solve_atsp (ncount, dist, &val, &rstate, tour);
    CCcheck_rval (rval, "solve_atsp failed");

    tval = 0.0;
    for (i = 0; i < ncount; i++) {
        tval += (int) (dist[tour[i]][tour[(i+1)%ncount]]);
    }
    printf ("Optimal Tour Length: %.0f\n", tval); fflush (stdout);
    printf ("Optimal Tour Order:\n");
    for (i = 0; i < ncount; i++) {
        printf ("%2d ", tour[i]);
        if (i%10 == 9) printf ("\n");
    }
    if (i%10) printf ("\n");


CLEANUP:
    CCutil_freedatagroup (&dat);
    if (dist) {
        for (i = 0; i < ncount; i++) {
            CC_IFFREE (dist[i], int);
        }
        CC_IFFREE (dist, int *);
    }
    CC_IFFREE (tour, int);

    return rval;
}

static int solve_atsp (int ncount, int **dist, double *pval,
        CCrandstate *rstate, int *tour)
{
    int rval = 0, success = 0, foundtour = 0, a, b, i, k, reverse = 0;
    int ecount = 0, *elist = (int *) NULL, *elen = (int *) NULL;
    int ancount = 3*ncount, *atour = (int *) NULL;

    *pval = 0.0;

    /* Split each node v into v_in - v_middle - v_opt */

    ecount = ncount * (ncount + 1);
    CC_MALLOC (elist, 2*ecount, int);
    CC_MALLOC (elen, ecount, int);

    for (a = 0, k = 0; a < ncount; a++) {
        for (b = 0; b < ncount; b++) {
            if (a == b) continue;
            elist[2*k]   = 3*a+2;
            elist[2*k+1] = 3*b;
            elen[k] = dist[a][b];
            k++;
        }
        elist[2*k]   = 3*a;
        elist[2*k+1] = 3*a+1;
        elen[k] = 0;
        k++;
        elist[2*k]   = 3*a+1;
        elist[2*k+1] = 3*a+2;
        elen[k] = 0;
        k++;
    }

    CC_MALLOC (atour, ancount, int);

    rval = CCtsp_solve_sparse (ancount, ecount, elist, elen, (int *) NULL,
               atour, (double *) NULL, pval, &success, &foundtour,
               (char *) NULL, (double *) NULL, (int *) NULL, 
               (CCtsp_cutselect *) NULL, 0, 1, rstate);
    CCcheck_rval (rval, "CCtsp_solve_sparse failed");

    if (success == 0) {
        printf ("tsp search hit a limit and was terminated\n");
        rval = 1;  goto CLEANUP;
    }

    /* find direction of tour by looking at node 0 */

    for (i = 0; i < ancount; i++) {
        if (atour[i] == 0) {
            if (     atour[(i+1)%ancount] == 1) reverse = 0;
            else if (atour[(i+ancount-1)%ancount] == 1) reverse = 1;
            else {
                printf ("returned tour does not match reduction (reverse)\n");
                rval = 1; goto CLEANUP;
            }
            break;
        }
    }

    for (i = 0; i < ancount; i++) {
        if (atour[i] % 3 == 0) {
            if (reverse == 0) {
                if (atour[(i+1)%ancount] != atour[i]+1 ||
                    atour[(i+2)%ancount] != atour[i]+2) {
                    printf ("returned tour does not match reduction process\n");
                    printf ("%d %d %d\n", atour[i], atour[(i+1)%ancount],
                                                 atour[(i+2)%ancount]);
                   rval = 1; goto CLEANUP;
                }
            } else {
                if (atour[(i+ancount-1)%ancount] != atour[i]+1 ||
                    atour[(i+ancount-2)%ancount] != atour[i]+2) {
                    printf ("returned tour does not match reduction process\n");
                    printf ("%d %d %d\n", atour[i],
                                          atour[(i+ancount-1)%ancount],
                                          atour[(i+ancount-2)%ancount]);
                    rval = 1; goto CLEANUP;
      
                }
            }
        }
    }

    k = 0;
    if (reverse == 0) {
        for (i = 0; i < ancount; i++) {
            if (atour[i] % 3 == 0) {
                if (tour) tour[k++] = atour[i]/3;
            }
        }
    } else {
        for (i = ancount-1; i  >= 0; i--) {
            if (atour[i] % 3 == 0) {
                if (tour) tour[k++] = atour[i]/3;
            }
        }
    }

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    return rval;
}

static int read_matrix (char *fname, int *pncount, int ***pdist)
{
    int rval = 0, ncount = 0, i, j, len, **dist = (int **) NULL;
    FILE *f = (FILE *) NULL;

    f = fopen (fname, "r");
    if (f == (FILE *) NULL) {
        fprintf (stderr, "Unable to open %s for input\n", fname);
        rval = 1; goto CLEANUP;
    }

    fscanf (f, "%d", &ncount);

    CC_MALLOC (dist, ncount, int *);
    for (i = 0; i < ncount; i++) dist[i] = (int *) NULL;
    for (i = 0; i < ncount; i++) {
        CC_MALLOC (dist[i], ncount, int);
    }

    for (i = 0; i < ncount; i++) {
        for (j = 0; j < ncount; j++) {
            fscanf (f, "%d", &len);
            dist[i][j] = len;
        }
    }

    fclose (f);
        
    *pncount = ncount;
    *pdist = dist;

CLEANUP:
    return rval;
}

static int test_perm (char *fname, int ncount, int **dist)
{
    int rval = 0, i, *perm = NULL;
    int t = 0;
    double val = 0.0;

    CC_MALLOC (perm, ncount, int);
    rval = CCutil_getcycle (ncount, fname, perm, 0);
    CCcheck_rval (rval, "CCutil_getcycle failed");

    for (i = 0; i < ncount; i++) {
        val += (double) (dist[perm[i]][perm[(i+1)%ncount]]);
        printf ("%2d %2d %9d  %.0f\n", perm[i], perm[(i+1)%ncount],
                           dist[perm[i]][perm[(i+1)%ncount]], val);
    }
    printf ("Test tour length: %.0f\n", val);

CLEANUP:
    CC_IFFREE (perm, int);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "P:s:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'P':
            permfname = boptarg;
            break;
        case 's':
            seed = atoi(boptarg);
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        matrixfname = av[boptind++];
    } else{
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] distance_matrix_file\n", fname);
    fprintf (stderr, "   -s #  seed\n");
    fprintf (stderr, "   -P f  permuation file for starting tour\n");
}


