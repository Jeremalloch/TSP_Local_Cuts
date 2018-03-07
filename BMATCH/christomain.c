#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "edgegen.h"
#include "bmatch.h"

static char *fname = (char *) NULL;
static char *outfname = (char *) NULL;
static int seed = 0;
static int nnodes_want = 0;
static int gridsize = 0;
static int tsplib_in = 1;
static int norm = CC_EUCLIDEAN;

static int parseargs (int ac, char **av);
static void usage (char *f);

int CCbmat_christofides (int ncount, CCdatagroup *dat, int *val, int *tour,
    CCrandstate *rstate);

int main (int ac, char **av)
{
    int rval = 0;
    int ncount, allow_dups, use_gridsize, val;
    int *tour = (int *) NULL;
    double t;
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

    t = CCutil_zeit ();

    printf ("Total Running Time: %.2f\n", CCutil_zeit () - t);

    tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (tour, "out of memory for tour");

    rval = CCbmat_christofides (ncount, &dat, &val, tour, &rstate);
    CCcheck_rval (rval, "CCbmat_christofides failed");

#if 0
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
#endif

CLEANUP:
    CCutil_freedatagroup (&dat); 
    if (out) fclose (out);
    CC_IFFREE (tour, int);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt(ac, av, "o:k:G:s:N:", &boptind, &boptarg)) != EOF) switch (c) {
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
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -G #  use #x# grid for random points, no duplicates if # < 0\n"
);
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -o f  output file for tour edges\n");
    fprintf (stderr, "   -N #  norm (must specify if not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 18=JOHNSON\n");
}

