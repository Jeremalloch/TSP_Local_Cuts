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
/*            STIPPLING ROUINES FOR BOSCH-KAPLAN-STYLE TSP ART              */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 14, 2012                                                    */
/*  Changes:                                                                */
/*                                                                          */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int CCutil_stipple_swarm (int probsize, int width, int height,          */
/*        int **gmatrix, int *ncount, double **xout, double **yout,         */
/*        CCrandstate *rstate);                                             */
/*    -Uses Jim Bumgardner's swarm algorithm for point placement            */
/*    -probsize should be 0, 1, or 2, specifying a desired smaller set of   */
/*     point (0), a default set of points (1), or a larger set of points    */
/*     (2)                                                                  */
/*    -width and height give the dimenstions of gmatrix                     */
/*    -gmatrix gives the grayscale for each pixel, from 0 to 255            */
/*    -ncount returns the number of points                                  */
/*    -xout and yout return the coordinates of the points (the arrays       */
/*     are allocated by the routine and should be freed by the calling      */
/*     function                                                             */
/*                                                                          */
/*  int CCutil_stipple_grid (int probsize, int width, int height,           */
/*        int **gmatrix, int *ncount, double **xout, double **yout,         */
/*        CCrandstate *rstate);                                             */
/*    -Uses the grid-based method for point placement similar to the        */
/*     techniques described by Bosch-Herman and Kaplan-Bosch                */
/*                                                                          */
/*  int CCutil_stipple_select (int probsize, int width, int height,         */
/*        int **gmatrix, int *ncount, double **xout, double **yout,         */
/*        CCrandstate *rstate);                                             */
/*    -Uses the random point placement based on gray scale (usually gives   */
/*     poor renderings                                                      */
/*     techniques described by Bosch-Herman and Kaplan-Bosch                */
/*                                                                          */
/*  int CCutil_stipple_lloyd (int probsize, int width, int height,          */
/*        int **gmatrix, int *ncount, double **xout, double **yout,         */
/*        CCrandstate *rstate);                                             */
/*    -Uses a method similar to the Secord's weighted version of Lloyd's    */
/*     algorithm.  Secord's algorithm was proposed for TSP Art by Kaplan    */
/*     and Bosch.                                                           */
/*                                                                          */
/*  CREDIT: The swarm method adopts the algorithm and code of               */
/*  Jim Bumgardner's Stipple Cam project.  Jim's cool project is described  */
/*  on the page http://joyofprocessing.com/blog/2011/11/stipple-cam/        */
/*                                                                          */
/****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "kdtree.h"
#include "edgegen.h"

#define DENSITY 0.50       /* default 0.5, higher = more points */
#define DAMPING 0.7        
#define KSPEED  3.0        
#define MINDISTFACTOR 2.5  /* default 2.5, smaller = faster rendering */

#define EPSDIST 0.001      /* epsilon for positive values */
#define SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

typedef struct particle {
    double  x, y, vx, vy, rad;
    double  fx, fy, wt;
} particle;

char *matrixfname = (char *) NULL;
char *outfname = (char *) NULL;
char *outtour = (char *) NULL;
int usegrid = 0;
int useselect = 0;
int uselloyd = 0;
int bigprob=0;
int smallprob=0;
int seed = 0;

static int stipple_swarm (int gridsize, int width, int height, int **gmatrix, 
    int *ncount, double **xout, double **yout, CCrandstate *rstate);
static int stipple_grid (int gridsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate);
static int stipple_select (int target, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate);
static int stipple_lloyd (int gridsize, int width, int height, int **gmatrix,
    int *outcount, double **xout, double **yout, CCrandstate *rstate);

static void init_particle (particle *p);
static int setup_particles (int width, int height, int cellsize, int boxwidth,
    int boxheight, double *minRadius, double *maxRadius, int *npart,
    particle **particles, CCrandstate *rstate);
static int move_particles (int width, int height, int **gmatrix, int cellsize,
    int boxwidth, int boxheight, int **gbox, double minRadius, 
    double maxRadius, int npart, particle *particles, CCrandstate *rstate);
static int buildbox (int gridsize, int width, int height, int **gmatrix,
    int *cellsize, int *boxwidth, int *boxheight, int ***gbox);
static int sort_corners (int ncount, int *triangles, int *ccount,
    int **corners, int *edge);
static void triangle_weight (double *points, int width, int height,
    int **gmatrix, double *tcenterx, double *tcentery, double *tweight);
static void build_equation (double ax, double ay, double bx, double by,
        double cx, double cy, double *coef, double *rhs);
static int remove_duplicates (int *ncount, double *x, double *y);
static double realrand (CCrandstate *rstate);
#if 0
static void perm_quicksort (int *perm, double *len, int n);
#endif
static void point_quicksort (double *x, double *y, int n);

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

int CCutil_stipple_swarm (int probsize, int width, int height, int **gmatrix,
        int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int gridsize, rval = 0;

    switch (probsize) {
    case 0:  gridsize = 100; break;
    case 1:  gridsize = 150; break;
    case 2:  gridsize = 200; break;
    default: gridsize = 150; break;
    }

    rval = stipple_swarm (gridsize, width, height, gmatrix, ncount, xout,
                          yout, rstate);
    CCcheck_rval (rval, "stipple_swarm failed");

CLEANUP:
    return rval;
}

int CCutil_stipple_grid (int probsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int gridsize, rval = 0;

    printf ("Grid-based Rendering\n"); fflush (stdout);

    switch (probsize) {
    case 0:  gridsize =  75; break;
    case 1:  gridsize = 100; break;
    case 2:  gridsize = 125; break;
    default: gridsize = 100; break;
    }

    rval = stipple_grid (gridsize, width, height, gmatrix, ncount, xout,
                         yout, rstate);
    CCcheck_rval (rval, "stipple_grid failed");

CLEANUP:
    return rval;
}

int CCutil_stipple_select (int probsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int target, rval = 0;

    printf ("Random Selection Rendering\n"); fflush (stdout);

    switch (probsize) {
    case 0:  target = 10000; break;
    case 1:  target = 15000; break;
    case 2:  target = 20000; break;
    default: target = 15000; break;
    }

    rval = stipple_select (target, width, height, gmatrix, ncount, xout,
                           yout, rstate);
    CCcheck_rval (rval, "stipple_select failed");

CLEANUP:
    return rval;
}

int CCutil_stipple_lloyd (int probsize, int width, int height, int **gmatrix,
    int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int gridsize, rval = 0;

    printf ("Weighted Lloyd's Algorithm Rendering\n"); fflush (stdout);

    switch (probsize) {
    case 0:  gridsize =  75; break;
    case 1:  gridsize = 100; break;
    case 2:  gridsize = 125; break;
    default: gridsize = 100; break;
    }

    rval = stipple_lloyd (gridsize, width, height, gmatrix, ncount, xout,
                           yout, rstate);
    CCcheck_rval (rval, "stipple_select failed");

CLEANUP:
    return rval;
}

#define MAXBMATCH 9 

int CCutil_stipple_bmatch (int probsize, int width, int height, int **gmatrix,
        int *ncount, double **xout, double **yout, int **degreeout,
        CCrandstate *rstate)
{
    int gridsize, rval = 0;
    int i, j, num;
    int count = 0, interval = 255 / (MAXBMATCH + 1);
    double *x = (double *) NULL, *y = (double *) NULL;
    int *degree = (int *) NULL;
    int **gbox = (int **) NULL;
    int cellsize, boxwidth, boxheight, total;

    *ncount = 0;
    *xout = (double *) NULL;
    *yout = (double *) NULL;
    *degreeout = (int *) NULL;

    printf ("Create Stippling for b-Matching \n"); fflush (stdout);

    switch (probsize) {
    case 0:  gridsize =  75; break;
    case 1:  gridsize = 100; break;
    case 2:  gridsize = 125; break;
    default: gridsize = 100; break;
    }

    rval = buildbox (gridsize, width, height, gmatrix, &cellsize, &boxwidth,
                     &boxheight, &gbox);
    CCcheck_rval (rval, "buildbox failed");

    count = 0;
    for (i = 0; i < boxheight; i++) {
        for (j = 0; j < boxwidth; j++) {
            num = (255 - gbox[i][j]) / interval;
            if (num) count++;
        }
    } 
    printf ("Number of points: %d\n", count);

    x = CC_SAFE_MALLOC (count, double);
    CCcheck_NULL (x, "out of memory for x");
    y = CC_SAFE_MALLOC (count, double);
    CCcheck_NULL (y, "out of memory for y");
    degree = CC_SAFE_MALLOC (count, int);
    CCcheck_NULL (degree, "out of memory for degree");
    
    count = 0, total = 0;
    for (i = 0; i < boxheight; i++) {
        for (j = 0; j < boxwidth; j++) {
            num = (255 - gbox[i][j]) / interval;
            /* num = (num * num) / 3; */
            if (num) {
                degree[count] = num;
                x[count] = j*cellsize + (CCutil_lprand(rstate) % cellsize);
                y[count] = i*cellsize + (CCutil_lprand(rstate) % cellsize); 
                count++;
                total += num;
            } 
        }
    }
    if (total % 2) degree[count-1]++;

    *ncount = count;
    *xout = x;
    *yout = y;
    *degreeout = degree; 

CLEANUP:
    if (gbox) {
        for (i = 0; i < boxheight; i++) {
            CC_IFFREE (gbox[i], int);
        }
        CC_FREE (gbox, int *);
    }

    if (rval) {
        CC_IFFREE (x, double);
        CC_IFFREE (y, double);
        CC_IFFREE (degree, int);
    }
    return rval;
}

static int stipple_swarm (int gridsize, int width, int height, int **gmatrix, 
        int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int i, rval = 0, npart = 0;
    double *x = (double *) NULL, *y = (double *) NULL;
    double minRadius = 0.0, maxRadius = 0.0;
    particle *particles = (particle *) NULL;
    int **gbox = (int **) NULL;
    int cellsize, boxwidth, boxheight;
    int rounds = 50;  /* default=50*/
    
    *ncount = 0;
    *xout = (double *) NULL;
    *yout = (double *) NULL;

    rval = buildbox (gridsize, width, height, gmatrix, &cellsize, &boxwidth,
                     &boxheight, &gbox);
    CCcheck_rval (rval, "buildbox failed");

    rval = setup_particles (width, height, cellsize, boxwidth, boxheight,
                       &minRadius, &maxRadius, &npart, &particles, rstate);
    CCcheck_rval (rval, "setup_particles failed");

    for (i = 0; i < rounds; i++) {
        rval = move_particles (width, height, gmatrix, cellsize, boxwidth,
                     boxheight, gbox, minRadius, maxRadius, npart, particles,
                     rstate);
        CCcheck_rval (rval, "move_particles failed");
        printf ("."); fflush (stdout);
        if (i % 25 == 24) printf ("\n");
    }
    if (i % 25) printf ("\n");

    x = CC_SAFE_MALLOC (npart, double);
    CCcheck_NULL (x, "out of memory for x");
    y = CC_SAFE_MALLOC (npart, double);
    CCcheck_NULL (y, "out of memory for y");
    for (i = 0; i < npart; i++) {
        x[i] = particles[i].x;
        y[i] = particles[i].y;
    }

    *ncount = npart;
    *xout = x;
    *yout = y;

CLEANUP:
    CC_IFFREE (particles, particle);
    if (gbox) {
        for (i = 0; i < boxheight; i++) {
            CC_IFFREE (gbox[i], int);
        }
        CC_FREE (gbox, int *);
    }
    if (rval) {
        CC_IFFREE (x, double);
        CC_IFFREE (y, double);
    }
    return rval;
}

static int setup_particles (int width, int height, int cellsize, int boxwidth,
        int boxheight, double *minRadius, double *maxRadius, int *npart,
        particle **particles, CCrandstate *rstate)
{
    int rval = 0;
    int dx, dy;
    double medArea, medRadius;
    int maxParticles = boxwidth*boxheight;
    particle *p;

    *particles = (particle *) malloc (maxParticles * sizeof (particle));
    if (!(*particles)) {
        fprintf (stderr, "out of memory for particles\n");
        rval = 1; goto CLEANUP;
    }
    p = *particles;

    *npart = 0;
    for (dy = 0; dy < boxheight; dy++) {
        for (dx = 0; dx < boxwidth; dx++) {
            if (realrand(rstate) < DENSITY) {   /* could work with intensity */
                init_particle (&p[*npart]);
                p[*npart].x = dx*cellsize;
                p[*npart].y = dy*cellsize;
                (*npart)++;
            }
        }
    }
    printf ("Number of particles: %d\n", *npart);

    medArea = (double) (width*height) / (double) (*npart + 1);
    medRadius = sqrt(medArea/M_PI);
    *minRadius = medRadius; /* using medRadius > 1 improves black areas */
    *maxRadius = medRadius*medRadius;
    printf ("minRadius = %f maxRadius = %f\n", *minRadius, *maxRadius);

CLEANUP:
    return rval;
}

static int try_the_move (int i, int j, void *pass_param)
{
    particle *particles = (particle *) pass_param;
    particle *p  = &particles[i];
    particle *pj = &particles[j];
    double dx, dy, distance, maxDist, diff;

    dx = p->x - pj->x;
    dy = p->y - pj->y;
    distance = sqrt(dx*dx+dy*dy);
      
    maxDist = (p->rad + pj->rad);
    diff = maxDist - distance;
    if (diff > EPSDIST && distance > EPSDIST && maxDist > EPSDIST) {
        double scle = diff/maxDist;
        scle = scle*scle;
        p->wt += scle;
        pj->wt += scle;
        scle = scle*KSPEED/distance;
        p->fx += dx*scle;
        p->fy += dy*scle;
        pj->fx -= dx*scle;
        pj->fy -= dy*scle;
    }
    return 0;
}

static int move_particles (int width, int height, int **gmatrix, int cellsize,
        int boxwidth, int boxheight, int **gbox, double minRadius, 
        double maxRadius, int npart, particle *particles, CCrandstate *rstate)
{
    int i, px, py, rval = 0;
    double v;
    int *perm = (int *) NULL;
    double *xarray = (double *) NULL;
    int usesnap = 0;

    if (usesnap) {
        for (i = 0; i < npart; i++) {
            px = (int) ((particles[i].x / cellsize) + 0.5);
            py = (int) ((particles[i].y / cellsize) + 0.5);
            if (px >= 0 && px < boxwidth && py >= 0 && py < boxheight) {
                v = (double) gbox[py][px];
                particles[i].rad = minRadius +
                     ((v/255.0) * (maxRadius-minRadius));
            }
        }
    } else {
        for (i = 0; i < npart; i++) {
            px = (int) (particles[i].x + 0.5);
            py = (int) (particles[i].y + 0.5);
            if (px < 0) px = 0;
            if (py < 0) py = 0;
            if (px > width-1)  px = width-1;
            if (py > height-1) py = height-1;
            v = (double) gmatrix[py][px];
            particles[i].rad = minRadius + ((v/255.0) * (maxRadius-minRadius));
        }
    }
  
    for (i = 0; i < npart; ++i) {
        particle *p = &particles[i];
        p->fx = p->fy = p->wt = 0;
        p->vx *= DAMPING;
        p->vy *= DAMPING;
    }

    /* Processing by increasing values of x can cut off the search early.  */
    /* So sort the x values, with order placed in the perm array.          */
    /* Note: could likely save more time with a kd-tree-based ball search. */

#if 0
    perm = (int *) malloc (npart * sizeof (int));
    CCcheck_NULL (perm, "out of memory for perm");
    xarray = (double *) malloc (npart * sizeof (double));
    CCcheck_NULL (xarray, "out of memory for xarray");
    for (i = 0; i < npart; i++) {
        perm[i] = i;
        xarray[i] = particles[i].x;
    }
    perm_quicksort (perm, xarray, npart);
#endif

{
    CCdatagroup swarmdat;
    CCkdtree K;

    CCutil_init_datagroup (&swarmdat);
    rval = CCutil_dat_setnorm (&swarmdat, CC_EUCLIDEAN);
    CCcheck_rval (rval, "CCutil_setnorm failed");
    swarmdat.x = CC_SAFE_MALLOC (npart, double);
    CCcheck_NULL (swarmdat.x, "out of memory for swarmdat.x");
    swarmdat.y = CC_SAFE_MALLOC (npart, double);
    CCcheck_NULL (swarmdat.y, "out of memory for swarmdat.y");
    for (i = 0; i < npart; i++) {
        swarmdat.x[i] = particles[i].x;
        swarmdat.y[i] = particles[i].y;
    }
    rval = CCkdtree_build (&K, npart, &swarmdat, (double *) NULL, rstate);
    CCcheck_rval (rval, "CCkdtree_build failed");

    for (i = 0; i < npart; i++) {
        rval = CCkdtree_fixed_radius_nearest (&K, &swarmdat, (double *) NULL,
                  i, particles[i].rad*MINDISTFACTOR, try_the_move,
                  (void *) particles);
        CCcheck_rval (rval, "CCkdtree_fixed_radius_nearest failed");
        CCkdtree_delete (&K, i);
    }

    CCutil_freedatagroup (&swarmdat);
    CCkdtree_free (&K);
}

#if 0
    /* particle interactions */
    for (i = 0; i < npart-1; ++i) {
        particle *p = &particles[perm[i]];
        int j;
        for (j = i+1; j < npart; ++j) {
            double dx, dy, distance, maxDist, diff;
            particle *pj = &particles[perm[j]];
            if (   (pj->x - p->x) > p->rad*MINDISTFACTOR) break;
            if (abs(pj->y - p->y) > p->rad*MINDISTFACTOR) continue; 

            dx = p->x - pj->x;
            dy = p->y - pj->y;
            distance = sqrt(dx*dx+dy*dy);
      
            maxDist = (p->rad + pj->rad);
            diff = maxDist - distance;
            if (diff > EPSDIST && distance > EPSDIST && maxDist > EPSDIST) {
                double scle = diff/maxDist;
                scle = scle*scle;
                p->wt += scle;
                pj->wt += scle;
                scle = scle*KSPEED/distance;
                p->fx += dx*scle;
                p->fy += dy*scle;
                pj->fx -= dx*scle;
                pj->fy -= dy*scle;
            }
        }
    }
#endif
  
    for (i = 0; i < npart; ++i) {
        particle *p = &particles[i];
        double dx, dy, distance, scle, diff;
        double maxDist = p->rad;

        /* keep within edges */
        /* left edge */
        distance = dx = p->x;
        diff = maxDist - distance;
        if (diff > EPSDIST && distance > EPSDIST && maxDist > EPSDIST) {
	    scle = diff/maxDist;
	    scle = scle*scle;
	    p->wt += scle;
	    scle = scle*KSPEED/distance;
            p->fx += 2*dx*scle;
        }
        /* right edge */
        dx = p->x - width; 
        distance = -dx;
        diff = maxDist - distance;
        if (diff > EPSDIST && distance > EPSDIST && maxDist > EPSDIST) {
	    scle = diff/maxDist;
	    scle = scle*scle;
	    p->wt += scle;
	    scle = scle*KSPEED/distance;
            p->fx += 2*dx*scle;
        }
        /* top edge */
        distance = dy = p->y; 
        diff = maxDist - distance;
        if (diff > EPSDIST && distance > EPSDIST && maxDist > EPSDIST) {
	    scle = diff/maxDist;
    	    scle = scle*scle;
	    p->wt += scle;
	    scle = scle*KSPEED/distance;
            p->fy += 2*dy*scle;
        }
        /* bot edge */
        dy = p->y - height;
        distance = -dy;
        diff = maxDist - distance;
        if (diff > EPSDIST && distance > EPSDIST && maxDist > EPSDIST) {
    	    scle = diff/maxDist;
	    scle = scle*scle;
	    p->wt += scle;
	    scle = scle*KSPEED/distance;
            p->fy += 2*dy*scle;
        }
        if (p->wt > EPSDIST) {
            p->vx += p->fx/p->wt;
            p->vy += p->fy/p->wt;
        }
        p->x += p->vx;
        p->y += p->vy;
        if (p->x < 0.0) p->x = 0.0;
        if (p->x > width) p->x = width;
        if (p->y < 0.0) p->y = 0.0;
        if (p->y > height) p->y = height;
    }

CLEANUP:
    if (perm) free (perm);
    if (xarray) free (xarray);
    return rval;
}

static void init_particle (particle *p)
{
    if (p) {
        p->vx = 0.0;
        p->vy = 0.0;
        p->rad = 1.0;
    }
}

static int buildbox (int gridsize, int width, int height, int **gmatrix,
        int *cellsize, int *boxwidth, int *boxheight, int ***gbox)
{
    int rval = 0;
    int i, j;
    int **b;

    if (width < height) {
        if (width >  gridsize) *cellsize = width  / gridsize;
    } else {
        if (height > gridsize) *cellsize = height / gridsize;
    }
    /* UPDATE: use a minimum cellsize, say cellsize = 3 */
    if (*cellsize < 3) *cellsize = 3;

    *boxwidth = width/(*cellsize);    /* leave off last fractional col */
    *boxheight = height/(*cellsize);  /* leave off last fractional row */

    *gbox = (int **) malloc (*boxheight * sizeof (int *));
    if (!(*gbox)) {
        fprintf (stderr, "out of memory for gbox\n");
        rval = 1; goto CLEANUP;
    }
    b = *gbox;

    for (i = 0; i < *boxheight; i++) b[i] = (int *) NULL;
    for (i = 0; i < *boxheight; i++) {
        b[i] = (int *) malloc (*boxwidth * sizeof (int));
        if (!b[i]) {
            fprintf (stderr, "out of memory for row of gbox\n"); 
            rval = 1; goto CLEANUP;
        }
    }
    /* Put average grayscale in each cell in gbox */
    for (i = 0; i < *boxheight; i++) {
        int row = i*(*cellsize);
        for (j = 0; j < *boxwidth; j++) {
            int col = j*(*cellsize);
            int gtotal = 0;
            int ri, rj;
            for (ri = 0; ri < (*cellsize); ri++) {
                for (rj = 0; rj < (*cellsize); rj++) {
                    gtotal += gmatrix[row+ri][col+rj];
                }
            }
            b[i][j] = (gtotal / ((*cellsize) * (*cellsize)));
        }
    }

CLEANUP:
    return rval;
}

#define MAXINCELL 3

static int stipple_grid (int gridsize, int width, int height, int **gmatrix,
        int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int i, j, k, num, rval = 0;
    int count = 0, interval = 255 / (MAXINCELL + 1);
    double *x = (double *) NULL, *y = (double *) NULL;
    int **gbox = (int **) NULL;
    int cellsize, boxwidth, boxheight;

    *ncount = 0;
    *xout = (double *) NULL;
    *yout = (double *) NULL;

    rval = buildbox (gridsize, width, height, gmatrix, &cellsize, &boxwidth,
                     &boxheight, &gbox);
    CCcheck_rval (rval, "buildbox failed");

    for (i = 0; i < boxheight; i++) {
        for (j = 0; j < boxwidth; j++) {
            num = (255 - gbox[i][j]) / interval;
            count += ((num * num) / 3);
        }
    } 
    printf ("Number of points: %d\n", count);

    x = CC_SAFE_MALLOC (count, double);
    CCcheck_NULL (x, "out of memory for x");
    y = CC_SAFE_MALLOC (count, double);
    CCcheck_NULL (y, "out of memory for y");
    
    count = 0;
    for (i = 0; i < boxheight; i++) {
        for (j = 0; j < boxwidth; j++) {
            num = (255 - gbox[i][j]) / interval;
            num = (num * num) / 3;
            for (k = 0; k < num; k++) {
                x[count] = j*cellsize + (CCutil_lprand(rstate) % cellsize);
                y[count] = i*cellsize + (CCutil_lprand(rstate) % cellsize); 
                count++;
            }
        }
    }

    *ncount = count;
    *xout = x;
    *yout = y;

CLEANUP:
    if (gbox) {
        for (i = 0; i < boxheight; i++) {
            CC_IFFREE (gbox[i], int);
        }
        CC_FREE (gbox, int *);
    }

    if (rval) {
        CC_IFFREE (x, double);
        CC_IFFREE (y, double);
    }
    return rval;
}

static int stipple_select (int target, int width, int height, int **gmatrix,
        int *ncount, double **xout, double **yout, CCrandstate *rstate)
{
    int rval = 0;
    int count = 0;
    double *x = (double *) NULL, *y = (double *) NULL;
    double v, dw, total = 0.0;
    int i, j;

    *ncount = 0;
    *xout = (double *) NULL;
    *yout = (double *) NULL;

    x = CC_SAFE_MALLOC (2*target, double);
    CCcheck_NULL (x, "out of memory for x");
    y = CC_SAFE_MALLOC (2*target, double);
    CCcheck_NULL (y, "out of memory for y");

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            total += (double) (255 - gmatrix[i][j]);
        }
    }
    if (total < 1.0) {
        fprintf (stderr, "image is white\n");
        goto CLEANUP;
    }

    dw = (double) target;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            v = (double) (255 - gmatrix[i][j]);
            if (realrand(rstate) <= (dw * v) / total) {
                x[count] = (double) j;
                y[count] = (double) i;
                count++;
                if (count == 2*target) break;
            }
        }
    }
    printf ("Selected %d points\n", count);

    *ncount = count;
    *xout = x;
    *yout = y;

CLEANUP:
    return rval;
}

static int stipple_lloyd (int gridsize, int width, int height, int **gmatrix,
        int *outcount, double **xout, double **yout, CCrandstate *rstate)
{
    int rval = 0;
    int i, j, z, tricount = 0, ncount = 0, newcount;
    int *triangles = (int *) NULL;
    double *x = (double *) NULL, *y = (double *) NULL;
    double points[6];
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);
    *outcount = 0;
    *xout = (double *) NULL;
    *yout = (double *) NULL;

    rval = stipple_grid (gridsize, width, height, gmatrix, &ncount, &x, &y,
                         rstate);
    CCcheck_rval (rval, "stipple_grid failed");

    rval = remove_duplicates (&ncount, x, y);
    CCcheck_rval (rval, "remove_duplicates failed");

    for (z = 0; z < 10; z++) {  /* 10 rounds is okay */
        int **corners = (int **) NULL;
        int *ccount = (int *) NULL;
        int *edge = (int *) NULL;
        double *xcenters = (double *) NULL;
        double *ycenters = (double *) NULL;

        CCutil_init_datagroup (&dat);
        printf ("Starting with %d points\n", ncount);

        CCutil_dat_setnorm (&dat, CC_EUCLIDEAN);
        dat.x = CC_SAFE_MALLOC (ncount, double);
        CCcheck_NULL (dat.x, "out of memory for dat.x");
        dat.y = CC_SAFE_MALLOC (ncount, double);
        CCcheck_NULL (dat.y, "out of memory for dat.y");

        for (i = 0; i < ncount; i++) {
            dat.x[i] = x[i];
            dat.y[i] = y[i];
        }
  
        rval = CCedgegen_delaunay (ncount, &dat, 0, (int *) NULL,
                   (int **) NULL, &tricount, &triangles);
        CCcheck_rval (rval, "CCedgegen_delaunay failed");
        printf ("Found %d triangles\n", tricount);

        /* For each point create a list of its triangles */

        ccount = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (ccount, "out of memory for ccount");
        corners = CC_SAFE_MALLOC (ncount, int *);
        CCcheck_NULL (ccount, "out of memory for corners");
        for (i = 0; i < ncount; i++) {
            ccount[i] = 0;
            corners[i] = (int *) NULL;
        }
        for (i = 0; i < tricount; i++) {
            ccount[triangles[3*i]]++;
            ccount[triangles[3*i+1]]++;
            ccount[triangles[3*i+2]]++;
        }
    
        for (i = 0; i < ncount; i++)  {
            if (ccount[i] == 0) {
                fprintf (stderr, "point in no triangle\n");
                rval = 1; goto CLEANUP;
            }
            corners[i] = CC_SAFE_MALLOC (ccount[i], int);
            CCcheck_NULL (corners[i], "out of memory for corners[i]");
            ccount[i] = 0;
        }
        for (i = 0; i < tricount; i++) {
            corners[triangles[3*i  ]][ccount[triangles[3*i  ]]++] = i;
            corners[triangles[3*i+1]][ccount[triangles[3*i+1]]++] = i;
            corners[triangles[3*i+2]][ccount[triangles[3*i+2]]++] = i;
        }

        /* Put each points list of triangles in circular order */
    
        edge = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (edge, "out of memory for edge");
        rval = sort_corners (ncount, triangles, ccount, corners, edge);
        CCcheck_rval (rval, "sort_corners failed");

        /* Compute circumcircle of each triangle */

        xcenters = CC_SAFE_MALLOC(tricount, double);
        CCcheck_NULL(xcenters, "out of memory for xcenters");
        ycenters = CC_SAFE_MALLOC(tricount, double);
        CCcheck_NULL(ycenters, "out of memory for ycenters");

        for (i = 0; i < tricount; i++) {
            double ax, ay, bx, by, cx, cy, d;
            int p, q;
            ax = x[triangles[3*i]]; 
            ay = y[triangles[3*i]]; 
            bx = x[triangles[3*i+1]]; 
            by = y[triangles[3*i+1]]; 
            cx = x[triangles[3*i+2]]; 
            cy = y[triangles[3*i+2]]; 
            d = 2.0*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by));
            if (d == 0.0) {
                printf ("WHAT %d?\n", i);
                printf ("a: %f %f  b: %f %f  c: %f %f\n", ax,ay,bx,by,cx,cy);
                printf ("%d %d %d\n", triangles[3*i], triangles[3*i+1],
                                      triangles[3*i+2]);
                exit (1);
            }
            xcenters[i] = ((ax*ax + ay*ay)*(by-cy) + (bx*bx + by*by)*(cy-ay) +
                           (cx*cx + cy*cy)*(ay-by))/d;
            ycenters[i] = ((ax*ax + ay*ay)*(cx-bx) + (bx*bx + by*by)*(ax-cx) +
                           (cx*cx + cy*cy)*(bx-ax))/d;

            p = (int) xcenters[i];
            q = (int) ycenters[i];
        }

        /* Move each point to its weighted center */

        newcount = 0;
        for (i = 0; i < ncount; i++) {
            double tcenterx=0.0, tcentery=0.0, tweight=0.0;
            double sum, xtotal, ytotal;
            double *cweights = (double *) NULL;
            double *centerx = (double *) NULL;
            double *centery = (double *) NULL;
            int newx, newy;

            cweights = CC_SAFE_MALLOC (ccount[i], double);
            CCcheck_NULL (cweights, "out of memory for cweights");
            centerx = CC_SAFE_MALLOC (ccount[i], double);
            CCcheck_NULL (centerx, "out of memory for centerx");
            centery = CC_SAFE_MALLOC (ccount[i], double);
            CCcheck_NULL (centery, "out of memory for centery");
            points[0] = x[i];
            points[1] = y[i];
            for (j = 1; j < ccount[i]; j++) {
                points[2] = xcenters[corners[i][j-1]];
                points[3] = ycenters[corners[i][j-1]];
                points[4] = xcenters[corners[i][j]];
                points[5] = ycenters[corners[i][j]];
                triangle_weight (points, width, height, gmatrix, &tcenterx,
                                 &tcentery, &tweight);
                cweights[j-1] = tweight; 
                centerx[j-1] = tcenterx;
                centery[j-1] = tcentery;
            }

            if (edge[i]) {
                sum = 0.0; xtotal = 0.0;  ytotal = 0.0;
            } else {
                points[2] = xcenters[corners[i][j-1]];
                points[3] = ycenters[corners[i][j-1]];
                points[4] = xcenters[corners[i][0]];
                points[5] = ycenters[corners[i][0]];
                triangle_weight (points, width, height, gmatrix, &tcenterx,
                                 &tcentery, &tweight);
                sum = tweight;
                xtotal = tweight*tcenterx;
                ytotal = tweight*tcentery;
            }

            for (j = 0; j < ccount[i]-1; j++) {
                sum += cweights[j];
                xtotal += cweights[j] * centerx[j];
                ytotal += cweights[j] * centery[j];
            }
            newx = (int) (xtotal / sum + 0.5); 
            newy = (int) (ytotal / sum + 0.5); 
            if (newx >= 0 && newx <= (double) width - 1 &&
                newy >= 0 && newy <= (double) height - 1) {
                x[newcount] = (double) newx;
                y[newcount] = (double) newy;
                newcount++;
            }
            
            CC_IFFREE (cweights, double);
            CC_IFFREE (centerx, double);
            CC_IFFREE (centery, double);
        }
        ncount = newcount;

        CC_IFFREE (ccount, int);
        if (corners) {
            for (i = 0; i < ncount; i++) {
                CC_IFFREE (corners[i], int);
            }
            CC_FREE (corners, int *);
        }
        CC_IFFREE (xcenters, double);
        CC_IFFREE (ycenters, double);
        CC_IFFREE (edge, int);
        CC_IFFREE (triangles, int);
        CCutil_freedatagroup (&dat);

        rval = remove_duplicates (&ncount, x, y);
        CCcheck_rval (rval, "remove_duplicates failed");
    }

    *outcount = ncount;
    *xout = x;
    *yout = y;

CLEANUP:
    CCutil_freedatagroup (&dat);
    CC_IFFREE (triangles, int);
    return rval;
}

typedef struct triangle
{
    int a;
    int b;
    int c;
} triangle;

static int sort_corners (int ncount, int *triangles, int *ccount,
        int **corners, int *edge)
{
    int rval = 0;
    int i, j, n, t, b, c, temp;
    int *hit = (int *) NULL, *order = (int *) NULL, *done = (int *) NULL;
    triangle *u = (triangle *) NULL;

    hit = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (hit, "out of memory for hit");
    for (n = 0; n < ncount; n++) hit[n] = 0;

    for (i = 0; i < ncount; i++) edge[i] = 0;

    for (n = 0; n < ncount; n++) {
        u = CC_SAFE_MALLOC (ccount[n], triangle);
        CCcheck_NULL (u, "out of memory for u");
        order = CC_SAFE_MALLOC (ccount[n], int);
        CCcheck_NULL (order, "out of memory for order");
        done = CC_SAFE_MALLOC (ccount[n], int);
        CCcheck_NULL (done, "out of memory for done");
        for (i = 0; i < ccount[n]; i++) {
            done[i] = 0;
            t = corners[n][i];
            u[i].a = triangles[3*t];
            u[i].b = triangles[3*t+1];
            u[i].c = triangles[3*t+2];
            if (u[i].b == n) {
                SWAP(u[i].a, u[i].b, temp);
            } else if (u[i].c == n) {
                SWAP(u[i].a, u[i].c, temp);
            }
            if (u[i].a != n) {
                fprintf (stderr, "triangle without its corner\n");
                rval = 1; goto CLEANUP;
            }
            hit[u[i].b]++;
            hit[u[i].c]++;
        }
        for (i = 0; i < ccount[n]; i++) {
            if (hit[u[i].b] == 1 || hit[u[i].c] == 1) break;
        }
        if (i < ccount[n]) {
            order[0] = i;
            done[i] = 1;
            edge[n] = 1;
        } else {
            order[0] = 0;
            done[0] = 1;
        }

        for (i = 1; i < ccount[n]; i++) {
            b = u[order[i-1]].b;
            c = u[order[i-1]].c;
            for (j = 0; j < ccount[n]; j++) {
                if (done[j]) continue;
                if (u[j].b == b || u[j].b == c || u[j].c == b || u[j].c == c)
                    break;
            }
            if (j == ccount[n]) {
                fprintf (stderr, "triangles not in circular order\n");
                rval = 1; goto CLEANUP;
            }
            done[j] = 1;
            order[i] = j;
        }
        for (i = 0; i < ccount[n]; i++) {
            done[i] = corners[n][order[i]];
        }
        for (i = 0; i < ccount[n]; i++) {
            corners[n][i] = done[i];
        }
        for (i = 0; i < ccount[n]; i++) {
            hit[u[i].b] = 0;
            hit[u[i].c] = 0;
        }
        CC_IFFREE (u, triangle);
        CC_IFFREE (done, int);
        CC_IFFREE (order, int);
    }
  
CLEANUP:
    CC_IFFREE (hit, int);
    CC_IFFREE (u, triangle);
    CC_IFFREE (done, int);
    CC_IFFREE (order, int);
    return rval;
}

static void triangle_weight (double *points, int width, int height,
        int **gmatrix, double *tcenterx, double *tcentery, double *tweight)
{
    int ipoints[6];
    double ax = points[0], ay = points[1];
    double bx = points[2], by = points[3];
    double cx = points[4], cy = points[5];
    double xcenter, ycenter, w;
    double M[3][2], b[3];
    int i, x, y, xmax, xmin, ymax, ymin;
    double xsum, ysum;

    for (i = 0; i < 6; i++) {
        ipoints[i] = (int) (points[i] + 0.5);
    }

    build_equation (ax, ay, bx, by, cx, cy, M[0], &b[0]);
    build_equation (ax, ay, cx, cy, bx, by, M[1], &b[1]);
    build_equation (bx, by, cx, cy, ax, ay, M[2], &b[2]);

    xmax = ipoints[0];  xmin = ipoints[0];
    ymax = ipoints[1];  ymin = ipoints[1];
    for (i = 1; i <= 2; i++) {
        if (ipoints[2*i]   > xmax) xmax = ipoints[2*i];
        if (ipoints[2*i]   < xmin) xmin = ipoints[2*i];
        if (ipoints[2*i+1] > ymax) ymax = ipoints[2*i+1];
        if (ipoints[2*i+1] < ymin) ymin = ipoints[2*i+1];
    }
    if (xmin < 0) xmin = 0;
    if (xmax >= width)  xmax = width-1;
    if (ymin < 0) ymin = 0;
    if (ymax >= height) ymax = height-1;

    w = 1.0;
    xsum = 1.0;
    ysum = 1.0;
    for (x = xmin; x <= xmax; x++) {
        for (y = ymin; y <= ymax; y++) {
            for (i = 0; i < 3; i++) {
                if (M[i][0]*(double)x + M[i][1]*(double)y > b[i]) break;
            }
            if (i == 3) {
                w += (double) (255 - gmatrix[y][x]);
                xsum += (x * (double) (255 - gmatrix[y][x]));
                ysum += (y * (double) (255 - gmatrix[y][x]));
            }
        }
    }

    xcenter = xsum / w;
    ycenter = ysum / w;

    xcenter = (ax + bx + cx) / 3.0;
    ycenter = (ay + by + cy) / 3.0;

    for (i = 0; i < 3; i++) {
        if (M[i][0]*xcenter + M[i][1]*ycenter > b[i] + 0.001) break;
    }
    if (i < 3) printf ("+");

    x = (int) (xcenter + 0.5);
    y = (int) (ycenter + 0.5);
    if (x < 0) x = 0;
    if (x >= width) x = width-1;
    if (y < 0) y = 0;
    if (y >= height) y = height-1;

    *tcenterx = (double) x;
    *tcentery = (double) y;
    *tweight = w;
}

static void build_equation (double ax, double ay, double bx, double by,
        double cx, double cy, double *coef, double *rhs)
{
    if (ax == bx) {
        coef[0] = 1.0;
        coef[1] = 0.0;
        rhs[0]  = ax;
    } else {
        coef[0] = by-ay;
        coef[1] = ax-bx;
        rhs[0]  = ax*(by-ay) - ay*(bx-ax);
    }

    if (coef[0]*bx + coef[1]*by > rhs[0] + 0.001 ||
        coef[0]*bx + coef[1]*by < rhs[0] - 0.001) {
        printf ("-"); fflush (stdout);
    }

    if (coef[0]*cx + coef[1]*cy > rhs[0]) {
        coef[0] = -coef[0];
        coef[1] = -coef[1];
        rhs[0]  = -rhs[0];
    }
}

static int remove_duplicates (int *ncount, double *x, double *y)
{
    int i, rval = 0, dupcount = 0, n = *ncount;
    double *xdup = (double *) NULL, *ydup = (double *) NULL;
    
    xdup = CC_SAFE_MALLOC (n, double);
    CCcheck_NULL (xdup, "out of memory for xdup");
    ydup = CC_SAFE_MALLOC (n, double);
    CCcheck_NULL (ydup, "out of memory for ydup");
    for (i = 0; i < n; i++) {
        xdup[i] = x[i];
        ydup[i] = y[i];
    }

    point_quicksort (xdup, ydup, n);
    for (i = 0; i < n-1; i++) {
        if (xdup[i] != xdup[i+1] || ydup[i] != ydup[i+1]) {
            x[dupcount] = xdup[i];
            y[dupcount] = ydup[i];
            dupcount++;
        }
    }
    x[dupcount] = xdup[n-1];
    y[dupcount] = ydup[n-1];
    dupcount++;
    if (dupcount < n) {
        printf ("Removed %d duplicate points\n", n-dupcount);
        fflush (stdout);
    }
    *ncount = dupcount;

CLEANUP:
    CC_IFFREE(xdup, double);
    CC_IFFREE(ydup, double);
    return rval;
}


#define THE_MULT 10000000

static double realrand (CCrandstate *rstate)
{
    return (double) (CCutil_lprand (rstate) % THE_MULT) / ((double) THE_MULT);
}

#if 0
static void perm_quicksort (int *perm, double *len, int n)
{
    int i, j, temp;
    double t;

    if (n <= 1)
        return;

    SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        SWAP (perm[i], perm[j], temp);
    }
    SWAP (perm[0], perm[j], temp);
    perm_quicksort (perm, len, j);
    perm_quicksort (perm + i, len, n - i);
}
#endif

static void point_quicksort (double *x, double *y, int n)
{
    int i, j;
    double tx, ty, temp;

    if (n <= 1)
        return;

    CC_SWAP (x[0], x[(n - 1)/2], temp);
    CC_SWAP (y[0], y[(n - 1)/2], temp);

    i = 0;
    j = n;
    tx = x[0];
    ty = y[0];

    while (1) {
        do i++; while (i < n && (x[i] < tx || (x[i] == tx && y[i] < ty)));
        do j--; while (x[j] > tx || (x[j] == tx && y[j] > ty));
        if (j < i) break;
        CC_SWAP (x[i], x[j], temp);
        CC_SWAP (y[i], y[j], temp);
    }
    CC_SWAP (x[0], x[j], temp);
    CC_SWAP (y[0], y[j], temp);

    point_quicksort (x, y, j);
    point_quicksort (x + i, y + i, n - i);

}

