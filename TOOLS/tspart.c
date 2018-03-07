/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2012 by David Applegate, Robert Bixby,              */
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
/*                  RUN TSP ART Routines (ala Kaplan-Bosch)                 */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 14, 2012                                                    */
/*  Changes:                                                                */
/*                                                                          */
/*  CREDIT: The swarm method adopts the algorithm and code of               */
/*  Jim Bumgardner's Stipple Cam project.  Jim's cool project is described  */
/*  on the page http://joyofprocessing.com/blog/2011/11/stipple-cam/        */
/*                                                                          */
/****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "tsp.h"

static char *matrixfname = (char *) NULL;
static char *outfname = (char *) NULL;
static char *outtour = (char *) NULL;
static int usegrid = 0;
static int useselect = 0;
static int uselloyd = 0;
static int usebmatch = 0;
static int probsize = 1;
static int seed = 0;

static int read_matrix (char *filename, int *width, int *height,
    int ***gmatrix);
static int write_xy (char *filename, int ncount, double *x, double *y);
static int write_tour (char *filename, int ncount, int *tour, CCdatagroup *dat);
static int write_edges (char *filename, int ncount, int ecount, int *elist,
    int *elen);
static int parseargs (int ac, char **av);
static void usage (char *f);
int main (int ac, char **av);

int CCutil_stipple_swarm (int probsize, int width, int height, int **gmatrix,
        int *ncount, double **xout, double **yout, CCrandstate *rstate);
int CCutil_stipple_grid (int probsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate);
int CCutil_stipple_select (int probsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate);
int CCutil_stipple_lloyd (int probsize, int width, int height, int **gmatrix,
    int *nount, double **xout, double **yout, CCrandstate *rstate);
int CCutil_stipple_bmatch (int probsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, int **degreeout,
    CCrandstate *rstate);

static int read_matrix (char *filename, int *width, int *height, int ***gmatrix)
{
    int rval = 0;
    FILE *in = (FILE *) NULL;
    int i, j, r, g, b;
    int **m;

    if ((in = fopen (filename, "r")) == NULL) {
        fprintf (stderr, "could not open %s for input\n", filename);
        rval = 1;  goto CLEANUP;
    }

    fscanf (in, "%d %d\n", width, height);
    *gmatrix = (int **) malloc (*height * sizeof (int *));
    if (!(*gmatrix)) {
        fprintf (stderr, "out of memory for gmatrix\n");
        rval = 1;  goto CLEANUP;
    }
    m = *gmatrix;

    for (i = 0; i < *height; i++) m[i] = (int *) NULL;
    for (i = 0; i < *height; i++) {
        m[i] = malloc (*width * sizeof (int));
        if (!m[i]) {
             fprintf (stderr, "out of memory for row of gmatrix\n");
             rval = 1;  goto CLEANUP;
        }
    }

    for (i = 0; i < *height; i++) {
        for (j = 0; j < *width; j++) {
            fscanf (in, "%d %d %d", &r, &g, &b);
            m[i][j] = (30*r + 59*g + 11*b) / 100;
        }
    }

    printf ("height = %d,  width = %d\n", *height, *width);
    fflush (stdout);

CLEANUP:

    if (in) fclose (in);
    return rval;
}

static int write_xy (char *filename, int ncount, double *x, double *y)
{
    int i, rval = 0;
    FILE *out = (FILE *) NULL;

    if ((out = fopen (filename, "w")) == NULL) {
        fprintf (stderr, "could not open %s for output\n", filename);
        rval = 1;  goto CLEANUP;
    }

    fprintf (out, "%d\n", ncount);
    for (i = 0; i < ncount; i++) {
        fprintf (out, "%f %f\n", x[i], y[i]);
    }

CLEANUP:
    if (out) fclose (out);
    return rval;
}

static int write_tour (char *filename, int ncount, int *tour, CCdatagroup *dat)
{
    int i, rval = 0;
    FILE *out = (FILE *) NULL;

    if ((out = fopen (filename, "w")) == NULL) {
        fprintf (stderr, "could not open %s for output\n", filename);
        rval = 1;  goto CLEANUP;
    }

    fprintf (out, "%d %d\n", ncount, ncount);
    for (i = 1; i < ncount; i++) {
        fprintf (out, "%d %d %d\n", tour[i-1], tour[i],
           (int) CCutil_dat_edgelen(tour[i-1], tour[i], dat));
    }
    fprintf (out, "%d %d %d\n", tour[ncount-1], tour[0],
           (int) CCutil_dat_edgelen(tour[ncount-1], tour[0], dat));

CLEANUP:
    if (out) fclose (out);
    return rval;
}

static int write_edges (char *filename, int ncount, int ecount, int *elist,
        int *elen)
{
    int i, rval = 0;
    FILE *out = (FILE *) NULL;

    if ((out = fopen (filename, "w")) == NULL) {
        fprintf (stderr, "could not open %s for output\n", filename);
        rval = 1;  goto CLEANUP;
    }
    fprintf (out, "%d %d\n", ncount, ecount);
    for (i = 0; i < ecount; i++) {
        fprintf (out, "%d %d %d\n", elist[2*i], elist[2*i+1], elen[i]);
    }

CLEANUP:
    if (out) fclose (out);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "abghlmo:p:s:", &boptind,
                                                          &boptarg )) != EOF)
	switch (c) {
	case 'a':
            probsize = 0;
	    break;
	case 'b':
            probsize = 2;
	    break;
	case 'g':
            usegrid = 1;
	    break;
	case 'h':
            useselect = 1;
	    break;
        case 'l':
            uselloyd = 1;
	    break;
        case 'm':
            usebmatch = 1;
	    break;
        case 'o':
            outfname = boptarg;
            break;
        case 'p':
            outtour = boptarg;
            break;
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
    matrixfname = av[boptind++];
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [flags below] matrix_file\n", f);
    fprintf (stderr, "    -a    use fewer points\n");
    fprintf (stderr, "    -b    use more points\n");
    fprintf (stderr, "    -g    use grid-based rendering\n");
    fprintf (stderr, "    -h    use random-selection rendering\n");
    fprintf (stderr, "    -l    use weighted Lloyd rendering\n");
    fprintf (stderr, "    -m    create b-matching data\n");
    fprintf (stderr, "    -o f  specify an output xy-file f\n");
    fprintf (stderr, "    -p f  specify an output edge-file f\n");
    fprintf (stderr, "    -s #  specify an integer for random seed\n");
}

int main (int ac, char **av)
{
    int rval = 0;
    int i, ncount = 0, width = 0, height = 0;
    int *tour = (int *) NULL, *degree = (int *) NULL;
    int mcount, *melist = (int *) NULL, *melen = (int *) NULL;
    int **gmatrix = (int **) NULL;
    double *x = (double *) NULL, *y = (double *) NULL;
    double tval;
    CCrandstate rstate;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;
    CCutil_sprand (seed, &rstate);

    rval = read_matrix (matrixfname, &width, &height, &gmatrix);
    if (rval) {
        fprintf (stderr, "read_matrix failed\n"); goto CLEANUP;
    }

    if (usebmatch) {
        rval = CCutil_stipple_bmatch (probsize, width, height, gmatrix, &ncount,
                                      &x, &y, &degree, &rstate);
        CCcheck_rval (rval, "CCutil_stipple_bmatch failed");
    } else if (usegrid) {
        rval = CCutil_stipple_grid (probsize, width, height, gmatrix, &ncount,
                             &x, &y, &rstate);
        CCcheck_rval (rval, "CCutil_stipple_grid failed");
    } else if (useselect) {
        rval = CCutil_stipple_select (probsize, width, height, gmatrix, &ncount,
                             &x, &y, &rstate);
        CCcheck_rval (rval, "CCutil_stipple_select failed");
    } else if (uselloyd) {
        rval = CCutil_stipple_lloyd (probsize, width, height, gmatrix, &ncount,
                             &x, &y, &rstate);
        CCcheck_rval (rval, "CCutil_stipple_lloyd failed");
    } else {
        rval = CCutil_stipple_swarm (probsize, width, height, gmatrix, &ncount,
                             &x, &y, &rstate);
        CCcheck_rval (rval, "CCutil_stipple_swarm failed");
    }

    if (outfname) {
        rval = write_xy (outfname, ncount, x, y);
        if (rval) {
            fprintf (stderr, "write_xy failed\n"); goto CLEANUP;
        }
    }

    if (usebmatch) {
        FILE *out = fopen ("bmatch.degree", "w");

        fprintf (out, "%d\n", ncount);
        for (i = 0; i < ncount; i++) {
            fprintf (out, "%d\n", degree[i]);
        }
        fclose (out);
    }

    if (outtour) {
        CCutil_dat_setnorm (&dat, CC_EUCLIDEAN);
        dat.x = CC_SAFE_MALLOC (ncount, double);
        CCcheck_NULL (dat.x, "out of memory for dat.x");
        dat.y = CC_SAFE_MALLOC (ncount, double);
        CCcheck_NULL (dat.y, "out of memory for dat.y");
        for (i = 0; i < ncount; i++) {
            dat.x[i] = x[i];
            dat.y[i] = y[i];
        }

        if (usebmatch) {
            rval = CCbmat_bmatch_geom (ncount, &dat, degree, &tval, &mcount,
                          &melist, &melen, &rstate);
            CCcheck_rval (rval, "CCbmat_bmatch_geom failed");

            rval = write_edges (outtour, ncount, mcount, melist, melen);
            CCcheck_rval (rval, "write_tour failed");
        } else {
            tour = CC_SAFE_MALLOC (ncount, int);
            CCcheck_NULL (tour, "out of memory for tour");
 
            rval = CCtsp_call_linkern (ncount, &dat, tour, &tval, 0.0, 1,
                                       &rstate);
            CCcheck_rval (rval, "CCtsp_call_linkern failed");

            rval = write_tour (outtour, ncount, tour, &dat);
            CCcheck_rval (rval, "write_tour failed");
        } 
    }

CLEANUP:
    if (gmatrix) {
        for (i = 0; i < height; i++) {
            if (gmatrix[i]) free (gmatrix[i]);
        }
        free (gmatrix);
        gmatrix = (int **) NULL;
    }
    CC_IFFREE (x, double);
    CC_IFFREE (y, double);
    CC_IFFREE (tour, int);
    CC_IFFREE (degree, int);
    CC_IFFREE (melist, int);
    CC_IFFREE (melen, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

