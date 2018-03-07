#include <stdio.h>
#include "machdefs.h"
#include "util.h"

static int
    receive_tspfile (CC_SFILE *s, char *name, int id, int *out_ncount);

static double
    child_zeit (void);

int main (int ac, char **av)
{
    char *bosshost = (char *) NULL;
    double newcost = -1.0, rtime = 0.0;
    char probname[CCutil_FILE_NAME_LEN];
    int i, ncount, id = -1, rval = 0;
    CC_SFILE *s = (CC_SFILE *) NULL;
    double szeit;
    char buf[1024], buf2[1024], key[1024], tbuf[1024];
    char *tb;
    FILE *p = (FILE *) NULL;
    FILE *out = (FILE *) NULL;

    if (ac < 2) {
        printf ("usage %s: boss\n", *av);
        rval = 1;  goto CLEANUP;
    }

    CCutil_printlabel ();

    bosshost = av[1];

    while (1) {
        s = CCutil_snet_open (bosshost, CC_SUBDIV_PORT);
        if (!s) {
            fprintf (stderr, "CCutil_snet_open failed\n");
            rval = 1;  goto CLEANUP;
        }
       
        rval = CCutil_swrite_int (s, id);
        CCcheck_rval (rval, "CCutil_swrite_int failed (id)");
        rval = CCutil_swrite_double (s, rtime);
        CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");
        rval = CCutil_swrite_double (s, newcost);
        CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");

        rval = CCutil_sread_int (s, &id);
        CCcheck_rval (rval, "CCutil_sread_int failed (id)");

        if (id == -1) {
            CCutil_sclose (s);
            s = (CC_SFILE *) NULL;
            goto DONE;
        }

        rval = CCutil_sread_string (s, probname, CCutil_FILE_NAME_LEN);
        CCcheck_rval (rval, "CCutil_sread_string failed (probname)");

        rval = receive_tspfile (s, probname, id, &ncount);
        CCcheck_rval (rval, "receive_tspfile failed");

        CCutil_sclose (s);

        sprintf (tbuf, "%s_temp.%d", probname, id);
        out = fopen (tbuf, "w");
        if (!out) {
            fprintf (stderr, "could not open %s for output\n", tbuf);
            rval = 1; goto CLEANUP;
        }
        fprintf (out, "NAME : perm with %d cities\n", ncount);
        fprintf (out, "TYPE : TOUR\n");
        fprintf (out, "DIMENSION : %d\n", ncount);
        fprintf (out, "TOUR_SECTION\n");
        for (i = 1; i <= ncount; i++) { fprintf (out, "%d\n", i); }
        fprintf (out, "-1\n");
        fprintf (out, "EOF\n");
        fclose (out);  out = (FILE *) NULL;

        sprintf (buf, "%s_par.%d", probname, id);
        out = fopen (buf, "w");
        if (!out) {
            fprintf (stderr, "could not open %s for output\n", buf);
            rval = 1; goto CLEANUP;
        }

        fprintf (out, "PROBLEM_FILE = %s_lkh_tsp.%d\n", probname, id);
        fprintf (out, "OUTPUT_TOUR_FILE = %s_tour.%d\n", probname, id);
        fprintf (out, "INITIAL_TOUR_FILE = %s\n", tbuf);
        fprintf (out, "RUNS = 1\n");
        fprintf (out, "TRACE_LEVEL = 0\n");

        /* For Million-City Star */

        fprintf (out, "INITIAL_PERIOD = 1000\n");
        fprintf (out, "CANDIDATE_SET_TYPE = DELAUNAY\n");
        fprintf (out, "MAX_CANDIDATES = 4\n");
        fprintf (out, "EXTRA_CANDIDATE_SET_TYPE = QUADRANT\n");
        fprintf (out, "EXTRA_CANDIDATES = 4\n");
        fprintf (out, "MOVE_TYPE = 10\n");
        fprintf (out, "PATCHING_C = 10\n");
        fprintf (out, "PATCHING_A = 9\n");
        fprintf (out, "MAX_TRIALS = 1\n");


/*
        if (ncount > 10000) {
            fprintf (out, "MAX_TRIALS = %d\n", 10000);
        } else {
            fprintf (out, "MAX_TRIALS = %d\n", ncount);
        }
        if (ncount > 10000) {
            fprintf (out, "CANDIDATE_SET_TYPE = DELAUNAY\n");
            fprintf (out, "INITIAL_PERIOD = 1000\n");
        }
*/

        fclose (out);  out = (FILE *) NULL;

        szeit = child_zeit ();

        sprintf (buf, "lkh_diamond %s_par.%d", probname, id);
        p = popen (buf, "r");
        if (!p) {
            perror (buf);
            fprintf (stderr, "popen failed\n");
            rval = 1; goto CLEANUP;
        }

        while ((fgets (buf2, sizeof (buf2), p)) != NULL) {
            buf2[sizeof (buf2) - 1] = '\0';
            fputs (buf2, stdout);
            fflush (stdout);
    
            if (sscanf (buf2, "%s", key) != EOF) {
                if (!strcmp (key, "Cost.min")) {
                    tb = buf2;
                    tb += strlen (key);
                    while (*tb == ' ' || *tb == '=') tb++;
                    if (sscanf (tb, "%lf", &newcost) == EOF) {
                        fprintf (stderr, "Could not read tourlen\n");
                        rval = 1;  goto CLEANUP;
                    }
                }
            }
        }

        pclose (p);

        if (newcost == -1.0) {
            fprintf (stderr, "failed to produce a tourlen\n");
            rval = 1; goto CLEANUP;
        }
 
        rtime = child_zeit () - szeit;
        printf ("New Tour Cost: %0.0f\n", newcost); fflush (stdout);

        rval = unlink (tbuf);
        if (rval) {
            fprintf (stderr, "Unable to delete temp path file %s\n", tbuf);
            rval = 0;
        }
    }

DONE:

    printf ("No work available.  Shutting down.\n"); fflush (stdout);

CLEANUP:

    if (out) fclose (out);
    return rval;
}

static int receive_tspfile (CC_SFILE *s, char *name, int id, int *out_ncount)
{
    int rval = 0;
    int i, norm, ncount;
    double *x = (double *) NULL;
    double *y = (double *) NULL;
    FILE *out = (FILE *) NULL;
    char buf[1024];

    rval = CCutil_sread_int (s, &norm);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    rval = CCutil_sread_int (s, &ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    if (out_ncount) *out_ncount = ncount;

    x = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (x, "out of memory in receive_tspfile");
    y = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (x, "out of memory in receive_tspfile");

    for (i = 0; i < ncount; i++) {
        rval = CCutil_sread_double (s, &x[i]);
        CCcheck_rval (rval, "CCutil_sread_double failed");
        rval = CCutil_sread_double (s, &y[i]);
        CCcheck_rval (rval, "CCutil_sread_double failed");
    }

    sprintf (buf, "%s_lkh_tsp.%d", name, id);
    out = fopen (buf, "w");
    if (!out) {
        fprintf (stderr, "unable to open %s for output\n", buf);
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "NAME: %s\n", buf);
    fprintf (out, "TYPE: TSP\n");
    fprintf (out, "COMMENT: Generated by lkhgrunt\n");
    fprintf (out, "DIMENSION: %d\n", ncount);

    fprintf (out, "EDGE_WEIGHT_TYPE: ");
    switch (norm) {
    case CC_MAXNORM:
        fprintf (out, "MAX_2D\n");
        break;
    case CC_MANNORM:
        fprintf (out, "MAN_2D\n");
        break;
    case CC_EUCLIDEAN_CEIL:
        fprintf (out, "CEIL_2D\n");
        break;
    case CC_EUCLIDEAN:
        fprintf (out, "EUC_2D\n");
        break;
    case CC_GEOM:
        fprintf (out, "GEOM\n");
        break;
    default:
        fprintf (stderr, "illegal NORM\n");
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "NODE_COORD_SECTION\n");
    for (i = 0; i < ncount; i++) {
        fprintf (out, "%d %f %f\n", i+1, x[i], y[i]);
    }

    fprintf (out, "FIXED_EDGES_SECTION\n");
    fprintf (out, "%d %d\n", ncount, 1);
    fprintf (out, "-1\n");
    fclose (out);  out = (FILE *) NULL;


CLEANUP:

    CC_IFFREE (x, double);
    CC_IFFREE (y, double);
    if (out) fclose (out);

    return rval;
}

static double child_zeit (void)
{
    struct rusage ru;

    getrusage (RUSAGE_CHILDREN, &ru);

    return ((double) ru.ru_utime.tv_sec) +
            ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

