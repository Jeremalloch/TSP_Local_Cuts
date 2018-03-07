#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "bmatch.h"

static char *fname = (char *) NULL;
static char *outfname = (char *) NULL;
static char *pricefname = (char *) NULL;
static char *degreefname = (char *) NULL;
static char *pricedatfname = (char *) NULL;
static int b_degree = 2;

static int get_edges (char *f, int *pncount, int *pecount, int **pelist,
    int **pelen);
static int get_degrees (char *f, int ncount, int **pdegrees);
static int parseargs (int ac, char **av);
static void usage (char *f);

int main (int ac, char **av)
{
    int i, k, rval = 0;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    int *melist = (int *) NULL, *melen = (int *) NULL;
    int *degrees = (int *) NULL;
    int ncount, ecount, mcount;
    double t;
    long val;
    FILE *out = (FILE *) NULL;
    CCdatagroup D, *dat = (CCdatagroup *) NULL;

    CCutil_init_datagroup (&D);

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    rval = get_edges (fname, &ncount, &ecount, &elist, &elen);
    CCcheck_rval (rval, "get_edges failed");

    printf ("Nodes: %d, Edges %d\n", ncount, ecount);
    fflush (stdout);

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

    if (pricedatfname) {
        int inorm, incount;
        rval = CCutil_gettsplib (pricedatfname, &incount, &D);
        CCcheck_rval (rval, "CCutil_gettsplib failed");
        if (incount != ncount) {
            fprintf (stderr, "TSPLIB file does not match edge file\n");
            rval = 1; goto CLEANUP;
        }
        CCutil_dat_getnorm (&D, &inorm);
        if ((inorm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
            fprintf (stderr, "TSPLIB file does not support kd-tree norms\n");
            rval = 1; goto CLEANUP;
        }
        if ((inorm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
            fprintf (stderr, "code only set up for 2-dimensional norms\n");
            rval = 1; goto CLEANUP;
        }
        dat = &D;
    }

    t = CCutil_zeit ();
    if (degrees) {
        printf ("Computing optimal b-matching ....\n"); fflush (stdout);
        rval = CCbmat_bmatch (ncount, ecount, elist, elen, degrees, &val,
                        dat, pricefname, &mcount, &melist, &melen);
        CCcheck_rval (rval, "CCbmat_bmatch failed");
    } else if (b_degree == 2) {
        printf ("Computing optimal 2-matching ....\n"); fflush (stdout);
        rval = CCbmat_twomatch (ncount, ecount, elist, elen, &val, dat,
                       pricefname, &mcount, &melist, &melen);
        CCcheck_rval (rval, "CCbmat_twomatch failed");
    } else {
        printf ("Computing optimal 1-matching ....\n"); fflush (stdout);
        rval = CCbmat_onematch (ncount, ecount, elist, elen, &val, dat,
                       pricefname, &mcount, &melist, &melen);
        CCcheck_rval (rval, "CCbmat_onematch failed");
    } 
    printf ("Running Time: %.2f\n", CCutil_zeit () - t);

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
    CCutil_freedatagroup (&D); 
    if (out) fclose (out);
    return rval;
}

static int get_edges (char *f, int *pncount, int *pecount, int **pelist,
        int **pelen)
{
    int i, ncount, ecount, end1, end2, w, rval = 0;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    FILE *in = (FILE *) NULL;

    if ((in = fopen (f, "r")) == NULL) {
        fprintf (stderr, "Unable to open %s for input\n", f);
        rval = 1; goto CLEANUP;
    }
    if (fscanf (in, "%d %d", &ncount, &ecount) != 2) {
        fprintf (stderr,"input file has invalid format\n");
        rval = 1; goto CLEANUP;
    }

    elist = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (elist, "out of memory for elist");
    elen = CC_SAFE_MALLOC (ecount, int);
    CCcheck_NULL (elen, "out of memory for elen");

    for (i=0; i<ecount; i++) {
        if (fscanf(in,"%d %d %d",&end1, &end2, &w) != 3) {
            fprintf (stderr,"input file has invalid format\n");
            rval = 1; goto CLEANUP;
        }
        elist[2*i] = end1;
        elist[2*i+1] = end2;
        elen[i] = w;
    }

    *pncount = ncount;
    *pecount = ecount;
    *pelist = elist;
    *pelen = elen;

CLEANUP:
    if (in) fclose (in);
    if (rval) {
        CC_IFFREE (elist, int);
        CC_IFFREE (elen, int);
    }
    return rval;
}

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
    int c;
    extern int optind;
    extern char *optarg;

    while ((c = getopt(ac, av, "b:1p:T:o:")) != EOF) switch (c) {
        case 'b': degreefname = optarg; break;
        case '1': b_degree = 1; break;
        case 'p': pricefname = optarg; break;
        case 'T': pricedatfname = optarg; break;
        case 'o': outfname = optarg; break;
        default: usage(av[0]); return 1;
    }
    if (optind >= ac) {
        usage(av[0]); return 1;
    }
    fname = av[optind++];
    if (optind != ac) {
        usage (av[0]); return 1;
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-b:1T:p:o:] edge_file\n",f);
    fprintf (stderr, "   -1    compute 1-matching (default 2-matching)\n");
    fprintf (stderr, "   -b f  vector file for b-matching\n");
    fprintf (stderr, "   -T f  TSPLIB file (Euclidean or MAX norm) for pricing edges on complete graph\n");
    fprintf (stderr, "   -p f  file of additional edges to price out\n");
    fprintf (stderr, "   -o f  output file for matching edges\n");
}

