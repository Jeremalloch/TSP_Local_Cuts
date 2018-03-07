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
/*                       SLOW (EXACT) TWO-OPT                               */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 9, 2005                                                     */
/*                                                                          */
/*  For a short describtion see usage ()                                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "linkern.h"
#include "util.h"
#include "kdtree.h"
#include "edgegen.h"
#include "macrorus.h"

#define BIGDOUBLE (1e30)

static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int tsplib_in = 1;
static int run_silently = 0;

static char *cycfname = (char *) NULL;
static char *nodefile = (char *) NULL;
static char *edgecycfname = (char *) NULL;
static char *saveit_final = (char *) NULL;


int
   main (int, char **);
static int
   print_command (int ac, char **av),
   twoopt_tour (int ncount, CCdatagroup *dat, int *incycle,
        int *outcycle, double *val, int silent, CCrandstate *rstate),
   swapalot (int ncount, int *tcyc, CCdatagroup *dat),
   parseargs (int, char **);
static void
   randcycle (int ncount, int *cyc, CCrandstate *rstate),
   usage (char *f);
static double
   cycle_length (int ncount, int *cyc, CCdatagroup *dat);


int main (int ac, char **av)
{
    int ncount;
    double val, startzeit;
    int *incycle = (int *) NULL, *outcycle = (int *) NULL;
    CCdatagroup dat;
    int rval = 0;
    CCrandstate rstate;

    CCutil_printlabel ();
    CCutil_init_datagroup (&dat);

    rval = print_command (ac, av);
    CCcheck_rval (rval, "print_command failed");

    seed = (int) CCutil_real_zeit ();
    if (parseargs (ac, av))
        return 1;
    CCutil_sprand (seed, &rstate);

    printf ("Two-Opt with seed %d\n", seed);
    fflush (stdout);

    if (!nodefile) {
        usage (av[0]);
        return 1;
    }

    startzeit = CCutil_zeit ();

    if (tsplib_in) {
        if (CCutil_gettsplib (nodefile, &ncount, &dat)) {
            fprintf (stderr, "could not read the TSPLIB file\n");
            rval = 1;
            goto CLEANUP;
        }
        CCutil_dat_getnorm (&dat, &norm);
    }

    incycle = CC_SAFE_MALLOC (ncount, int);
    if (!incycle) {
        rval = 1;
        goto CLEANUP;
    }
    if (cycfname) {
        rval = CCutil_getcycle (ncount, cycfname, incycle, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle failed\n");
            goto CLEANUP;
        }
    } else if (edgecycfname) {
        rval = CCutil_getcycle_edgelist (ncount, edgecycfname, incycle, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle_edgelist failed\n");
            goto CLEANUP;
        }
    } else {
        printf ("Use random starting cycle\n");
        randcycle (ncount, incycle, &rstate);
    }

    outcycle = CC_SAFE_MALLOC (ncount, int);
    if (!outcycle) {
        rval = 1; goto CLEANUP;
    }

    rval = twoopt_tour (ncount, &dat, incycle, outcycle, &val, run_silently,
                        &rstate);
    if (rval) {
        fprintf (stderr, "CClinkern_tour failed\n");
        goto CLEANUP;
    }

    if (saveit_final) {
        rval = CCutil_writecycle_edgelist (ncount, saveit_final, outcycle,
                                           &dat, 0);
        if (rval) {
            fprintf (stderr, "could not write the cycle\n");
            goto CLEANUP;
        }
    }

    printf ("Final Cycle: %.0f\n", val);
    fflush (stdout);
    printf ("Total Running Time: %.2f\n", CCutil_zeit () - startzeit);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (incycle, int);
    CC_IFFREE (outcycle, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static void randcycle (int ncount, int *cyc, CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++)
        cyc[i] = i;

    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }
}


static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "o:Qs:y:Y:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'o':
            saveit_final = boptarg;
            break;
        case 'Q':
            run_silently++;
            break;
        case 's':
            seed = atoi (boptarg);
            break;
        case 'y':
            cycfname = boptarg;
            break;
        case 'Y':
            edgecycfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind < ac)
        nodefile = av[boptind++];

    if (boptind > ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- see below -] tsplib_file\n", f);
    fprintf (stderr, "   -s #  random number seed\n");
    fprintf (stderr, "   -o f  save final tour\n");
    fprintf (stderr, "   -y f  starting cycle\n");
    fprintf (stderr, "   -Y f  starting cycle (as an edgelist)\n");
    fprintf (stderr, "   -Q    run silently\n");
}

static int print_command (int ac, char **av)
{
    int rval = 0;
    int i, cmdlen = 0;
    char *cmdout = (char *) NULL;

    for (i=0; i<ac; i++) {
        cmdlen += strlen(av[i]) + 1;
    }
    cmdout = CC_SAFE_MALLOC (cmdlen, char);
    CCcheck_NULL (cmdout, "out of memory in print_command");

    cmdlen = 0;
    for (i=0; i<ac; i++) {
        strcpy (cmdout + cmdlen, av[i]);
        cmdlen += strlen(av[i]);
        cmdout[cmdlen] = ' ';
        cmdlen++;
    }
    cmdout[cmdlen-1] = '\0';
    printf ("%s\n", cmdout); fflush (stdout);

CLEANUP:

    CC_IFFREE (cmdout, char);
    return rval;
}


#define FLIP(aprev, a, b, bnext, f, x) {                                   \
    CClinkern_flipper_flip ((x),(a), (b));                                 \
    (f)->stack[(f)->counter].first = (a);                                  \
    (f)->stack[(f)->counter++].last = (b);                                 \
}

#define UNFLIP(aprev, a, b, bnext, f, x) {                                 \
    CClinkern_flipper_flip ((x), (b), (a));                                \
    (f)->counter--;                                                        \
}

static int twoopt_tour (int ncount, CCdatagroup *dat, int *incycle,
        int *outcycle, double *val, int silent, CCrandstate *rstate)
{
    int rval = 0;
    int i;
    int *tcyc = (int *) NULL;
    double startzeit;

    if (silent == 0) {
        printf ("twoopt ...\n"); fflush (stdout);
    }
    startzeit = CCutil_zeit ();

    tcyc = CC_SAFE_MALLOC (ncount, int);
    if (tcyc == (int *) NULL) {
        fprintf (stderr, "out of memory in linkern\n");
        rval = 1; goto CLEANUP;
    }

    if (incycle) {
        for (i = 0; i < ncount; i++) tcyc[i] = incycle[i];
    } else {
        randcycle (ncount, tcyc, rstate);
    }
    *val = cycle_length (ncount, tcyc, dat);
    if (silent == 0) {
        printf ("Starting Cycle: %.0f\n", *val); fflush (stdout);
    }

    rval = swapalot (ncount, tcyc, dat);
    if (rval) {
        fprintf (stderr, "swapalot failed\n");
        goto CLEANUP;
    }

    if (silent == 0) {
        *val = cycle_length (ncount, tcyc, dat);
        printf ("Best cycle length: %.0f\n", *val);
        printf ("Two-Opt Running Time: %.2f\n",
                  CCutil_zeit () - startzeit);
        fflush (stdout);
    }

    if (outcycle) {
        for (i = 0; i < ncount; i++) outcycle[i] = tcyc[i];
    }

CLEANUP:

    CC_IFFREE (tcyc, int);
    return rval;
}

static int swapalot (int ncount, int *tcyc, CCdatagroup *dat)
{
    int rval = 0;
    int a, b, c, d, improve;
    int k = 0;
    double sub1, sub2, add1, add2, total = 0.0;
    CClk_flipper F;

    rval = CClinkern_flipper_init (&F, ncount, tcyc);
    if (rval) {
        fprintf (stderr, "CClinkern_flipper_init failed\n");
        goto CLEANUP;
    }

    do {
       improve = 0;
       a = 0;
       do {
           b = CClinkern_flipper_next (&F, a);
           sub1 = (double) CCutil_dat_edgelen (a, b, dat);
           c = CClinkern_flipper_next (&F, b);
           d = CClinkern_flipper_next (&F, c);
           do {
               sub2 = (double) CCutil_dat_edgelen (c, d, dat);
               add1 = (double) CCutil_dat_edgelen (a, c, dat);
               add2 = (double) CCutil_dat_edgelen (b, d, dat);
               if (add1 + add2 < sub1 + sub2) {
                   printf ("Improve %d: %.0f\n", k, sub1 + sub2 - add1 - add2);
                   fflush (stdout);
                   CClinkern_flipper_flip (&F, b, c);
                   improve = 1;
                   k++;
                   total += (sub1 + sub2 - add1 - add2);
               }
               c = d;
               d = CClinkern_flipper_next (&F, d);
           } while (!improve && d != a);  
           a = b;
        } while (!improve && a != 0);
    } while (improve);

    CClinkern_flipper_cycle (&F, tcyc);
    CClinkern_flipper_finish (&F);
 
    printf ("Total Improvement: %.0f\n", total); fflush (stdout);

CLEANUP:

    return rval;
}

static double cycle_length (int ncount, int *cyc, CCdatagroup *dat)
{
    int i;
    double val = 0.0;

    for (i = 1; i < ncount; i++) {
        val += (double) CCutil_dat_edgelen (cyc[i-1], cyc[i], dat);
    }
    val += (double) CCutil_dat_edgelen (cyc[0], cyc[ncount-1], dat);

    return val;
}

