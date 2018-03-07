#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "edgegen.h"
#include "bmatch.h"

static char *fname = (char *) NULL;
static char *outfname = (char *) NULL;
static char *pricefname = (char *) NULL;
static char *degreefname = (char *) NULL;
static char *pricedatname = (char *) NULL;
static int b_degree = 2;
static int seed = 0;
static int nnodes_want = 0;
static int gridsize = 0;
static int tsplib_in = 1;
static int norm = CC_EUCLIDEAN;

static int get_degrees (char *f, int ncount, int **pdegrees);
static int parseargs (int ac, char **av);
static void usage (char *f);

int main (int ac, char **av)
{
    int i, k, rval = 0;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    int *melist = (int *) NULL, *melen = (int *) NULL;
    int *degrees = (int *) NULL;
    int ncount, ecount, mcount, allow_dups, use_gridsize;
    double t;
    long val;
    FILE *out = (FILE *) NULL;
    CCrandstate rstate;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;
    CCutil_sprand (seed, &rstate);

    if ((!nnodes_want && !fname) || (tsplib_in && !fname)) {
        usage (av[0]);
        goto CLEANUP;
    }

    if (tsplib_in) {
        rval = CCutil_gettsplib (fname, &ncount, &dat);
        CCcheck_rval (rval, "CCutil_gettsplib failed");
        CCutil_dat_getnorm (&dat, &norm);
    } else {
        ncount = nnodes_want;
        if (gridsize < 0) {
            use_gridsize = -gridsize;
            allow_dups = 0;
        } else if (gridsize > 0) {
            use_gridsize = gridsize;
            allow_dups = 1;
        } else {
            use_gridsize = nnodes_want;
            allow_dups = 0;
        }
        rval = CCutil_getdata (fname, 0, norm, &ncount, &dat, use_gridsize,
                               allow_dups, &rstate);
        CCcheck_rval (rval, "CCutil_getdata failed");
    }

    if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "code only set up for kd-tree norms\n");
        rval = 1; goto CLEANUP;
    }
    if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "code only set up for 2-dimensional norms\n");
        rval = 1; goto CLEANUP;
    }

    if (degreefname) {
        rval = get_degrees (degreefname, ncount, &degrees);
        CCcheck_rval (rval, "get_degrees failed");
        for (i = 0, k = 0; i < ncount; i++) k += degrees[i];
        printf ("Sum of b-degrees: %d\n", k);
        fflush (stdout);
        if (k % 2) {
            fprintf (stderr, "Sum of b-degrees is odd\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        if ((b_degree * ncount) % 2) {
            fprintf (stderr, "Sum of b-degrees is odd\n");
            rval = 1; goto CLEANUP;
        }
    }

    t = CCutil_zeit ();

    if (degrees) {
        printf ("Computing optimal b-matching ....\n"); fflush (stdout);
        rval = CCbmat_bmatch_geom (ncount, &dat, degrees, &val, &mcount,
                       &melist, &melen, &rstate);
        CCcheck_rval (rval, "CCbmat_bmatch failed");
    } else if (b_degree == 2) {
        printf ("Computing optimal 2-matching ....\n"); fflush (stdout);
        rval = CCbmat_twomatch_geom (ncount, &dat, &val, &mcount, &melist,
                       &melen, &rstate);
        CCcheck_rval (rval, "CCbmat_twomatch failed");
    } else {
        printf ("Computing optimal 1-matching ....\n"); fflush (stdout);
        rval = CCbmat_onematch_geom (ncount, &dat, &val, &mcount, &melist,
                       &melen, &rstate);
        CCcheck_rval (rval, "CCbmat_onematch failed");
    } 
    printf ("Total Running Time: %.2f\n", CCutil_zeit () - t);

    if (outfname) {
        out = fopen (outfname, "w");
        if (!out) {
            fprintf (stderr, "unable to open %s for output\n", outfname);
            rval = 1; goto CLEANUP;
        } 
        fprintf (out, "%d %d\n", ncount, mcount);
        for (i = 0; i < mcount; i++) {
            fprintf (out, "%d %d %d\n", melist[2*i], melist[2*i+1], melen[i]);
        }
    }

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (melist, int);
    CC_IFFREE (melen, int);
    CC_IFFREE (degrees, int);
    CCutil_freedatagroup (&dat); 
    if (out) fclose (out);
    return rval;
}

#if 0
static int get_edges (CCdatagroup *dat, int ncount, int *pecount, int **pelist,
        int **pelen, int bdegree, int *degrees, CCrandstate *rstate)
{
    int i, k, ecount, end1, end2, w, rval = 0;
    CCedgegengroup plan;
    int *elist = (int *) NULL, *elen = (int *) NULL;

    CCedgegen_init_edgegengroup (&plan);

    if (degrees) {
        k = 0;
        for (i = 0; i < ncount; i++) {
            if (degrees[i] > k) k = degrees[i];
        }
    } else {
        k = bdegree;
    }

    plan.f2match_nearest.number = 5*k;
    plan.f2match_nearest.basic = 0;
    plan.f2match_nearest.priced = 1;
    plan.tour.greedy = 1;
 
    printf ("Generating %d fractional 2-matching neighbors\n", 5*k);
    fflush (stdout);

    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &ecount,
                            &elist, 0, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");
    
    elen = CC_SAFE_MALLOC (ecount, int);
    CCcheck_NULL (elen, "out of memory for elen");
    for (i = 0; i < ecount; i++) {
        elen[i] = CCutil_dat_edgelen (elist[2*i], elist[2*i+1], dat);
    }

    *pecount = ecount;
    *pelist = elist;
    *pelen = elen;

CLEANUP:
    if (rval) {
        CC_IFFREE (elist, int);
        CC_IFFREE (elen, int);
    }
    return rval;
}
#endif

static int get_degrees (char *f, int ncount, int **pdegrees)
{
    int i, incount, k, rval = 0;
    int *degrees = (int *) NULL;
    FILE *in = (FILE *) NULL;

    if ((in = fopen (f, "r")) == NULL) {
        fprintf (stderr, "Unable to open %s for input\n", f);
        rval = 1; goto CLEANUP;
    }
    if (fscanf (in, "%d", &incount) != 1) {
        fprintf (stderr,"input file has invalid format\n");
        rval = 1; goto CLEANUP;
    }
    if (incount != ncount) {
        fprintf (stderr, "degree file does not match problem data\n");
        rval = 1; goto CLEANUP;
    }

    degrees = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (degrees, "out of memory for degrees");

    for (i=0; i<ncount; i++) {
        if (fscanf(in,"%d", &k) != 1) {
            fprintf (stderr,"input file has invalid format\n");
            rval = 1; goto CLEANUP;
        }
        degrees[i] = k;
    }

    *pdegrees = degrees;

CLEANUP:
    if (in) fclose (in);
    if (rval) {
        CC_IFFREE (degrees, int);
    }
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt(ac, av, "b:1o:k:G:s:N:", &boptind, &boptarg)) != EOF) switch (c) {
        case 'b':
            degreefname = boptarg;
            break;
        case '1':
            b_degree = 1;
            break;
        case 'o':
            outfname = boptarg;
            break;
        case 'k':
            nnodes_want = atoi (boptarg);
            tsplib_in = 0;
            break;
        case 'G':
            gridsize = atoi(boptarg);
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'N':
            inorm = atoi(boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 18: norm = CC_EUCLIDEAN_CEIL; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage(av[0]); return 1;
    }
    if (boptind < ac) {
        fname = av[boptind++];
    }
    if (boptind != ac) {
        usage (av[0]); return 1;
    }
    return 0;
}

static void usage (char *f)
{
fprintf (stderr, "Usage: %s [- see below -] [tsplib_file or dat_file]\n", f);
    fprintf (stderr, "   -1    compute 1-matching (default 2-matching)\n");
    fprintf (stderr, "   -b f  vector file for b-matching\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -G #  use #x# grid for random points, no duplicates if # < 0\n"
);
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -o f  output file for matching edges\n");
    fprintf (stderr, "   -N #  norm (must specify if not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 18=JOHNSON\n");
}

