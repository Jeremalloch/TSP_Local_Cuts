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
/*         A PROGRAM TO COMPUTE CONVERT A TOUR FILE TO AN EDGE FILE         */
/*                 or CONVERT TO A TSPLIB TOUR FILE                         */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: January 24, 2001                                                  */
/*        May 7, 2002 (bico)                                                */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *cycfname = (char *) NULL;
static char *tspfname = (char *) NULL;
static char *outfname = (char *) NULL;
static int tsplib_in = 1;
static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int simpletour = 0;
static int tsplibtour = 0;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static char
    *get_problabel (const char *probloc);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int rval = 0, i, j, ncount = 0, *tour = (int *) NULL;
    double val;
    CCdatagroup dat;
    CCrandstate rstate;
    char *name = (char *) NULL, *p, buf[256], key[256], field[256];
    FILE *out = (FILE *) NULL, *in = (FILE *) NULL;

    CCutil_init_datagroup (&dat);
    seed = (int) CCutil_real_zeit ();
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    CCutil_sprand (seed, &rstate);

    if (tspfname) {
        if (tsplib_in) {
            rval = CCutil_gettsplib (tspfname, &ncount, &dat);
            CCcheck_rval (rval, "CCutil_gettsplib failed");
        } else {
            rval = CCutil_getdata (tspfname, 0, norm, &ncount, &dat, 1, 0,
                                   &rstate);
            CCcheck_rval (rval, "CCutil_getdata failed");
        }
    } else {
        printf ("Edges will not have length fields\n"); fflush (stdout);

        /* get ncount from top of cycle file */
        in = fopen (cycfname, "r");
        if (!in) {
            fprintf (stderr, "could not open %s for reading\n", cycfname);
            rval = 1; goto CLEANUP;
        }
        if (simpletour) {
            fscanf (in, "%d", &ncount);
        } else {
            while (ncount == 0 && fgets (buf, 254, in) != (char *) NULL) {
                p = buf;
                while (*p != '\0') {
                    if (*p == ':') *p = ' ';
                    p++;
                }
                p = buf;
                if (sscanf (p, "%s", key) != EOF) {
                    p += strlen (key);
                    while (*p == ' ') p++;
                    if (!strcmp (key, "DIMENSION")) {
                        if (sscanf (p, "%s", field) == EOF) {
                            fprintf (stderr, "ERROR in DIMENSION line\n");
                            rval = 1; goto CLEANUP;
                        }
                        ncount = atoi (field);
                    }
                }
            } 
        }
        fclose (in);  in = (FILE *) NULL;
        printf ("ncount = %d\n", ncount);
        if (ncount == 0) {
            fprintf (stderr, "unable to grab ncount from cycle file\n");
            rval = 1; goto CLEANUP;
        }
    }

    CC_MALLOC (tour, ncount, int);
    if (simpletour) {
        rval = CCutil_getcycle (ncount, cycfname, tour, 0);
        CCcheck_rval (rval, "CCutil_getcycle failed");
    } else {
        rval = CCutil_getcycle_tsplib (ncount, cycfname, tour);
        CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");
    }

    if (tspfname) {
        CCutil_cycle_len (ncount, &dat, tour, &val);
        printf ("Tour Length: %.0f\n", val); fflush (stdout);
    }

    if (outfname) {
        out = fopen (outfname, "w");
        if (!out) {
            fprintf (stderr, "could not open %s for writing\n", outfname);
            rval = 1; goto CLEANUP;
        }
        
        if (tsplibtour) {
            name = get_problabel (tspfname);
            fprintf (out, "NAME: %s\n", name);
            fprintf (out, "COMMENT: Tour length %.0f\n", val);
            fprintf (out, "TYPE: TOUR\n");
            fprintf (out, "DIMENSION: %d\n", ncount);
            fprintf (out, "TOUR_SECTION\n");
            for (i = 0; i < ncount; i++) {
                fprintf (out, "%d\n", tour[i] + 1);
            }
            fprintf (out, "-1\n");
            fprintf (out, "EOF\n");
        } else {
            if (tspfname) {
                fprintf (out, "%d %d\n", ncount, ncount);
                for (i = 0; i < ncount; i++) {
                    j = (i+1) % ncount;
                    fprintf (out, "%d %d %d\n", tour[i], tour[j],
                            CCutil_dat_edgelen (tour[i], tour[j], &dat));
                }
            } else {
                for (i = 0; i < ncount; i++) {
                    j = (i+1) % ncount;
                    fprintf (out, "[%d, %d],\n", tour[i], tour[j]);
                }
            }
        }
    }

CLEANUP:
    CC_IFFREE (tour, int);
    CC_IFFREE (name, char);
    CCutil_freedatagroup (&dat);
    if (out) fclose (out);
    if (in) fclose (in);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, boptind = 1, inorm = 0;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "N:o:StT:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'o': outfname = boptarg; break;
        case 'N':
            inorm = atoi(boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 3: norm = CC_EUCLIDEAN_3D; break;
            case 4: norm = CC_USER; break;
            case 5: norm = CC_ATT; break;
            case 6: norm = CC_GEOGRAPHIC; break;
            case 7: norm = CC_MATRIXNORM; break;
            case 8: norm = CC_DSJRANDNORM; break;
            case 9: norm = CC_CRYSTAL; break;
            case 10: norm = CC_SPARSE; break;
            case 11: norm = CC_RHMAP1; break;
            case 12: norm = CC_RHMAP2; break;
            case 13: norm = CC_RHMAP3; break;
            case 14: norm = CC_RHMAP4; break;
            case 15: norm = CC_RHMAP5; break;
            case 16: norm = CC_EUCTOROIDAL; break;
            case 17: norm = CC_GEOM; break;
            case 18: norm = CC_EUCLIDEAN_CEIL; break;
            case 20: norm = CC_ROAD; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case 'S': tsplibtour = 1; break;
        case 't': simpletour = 1; break;
        case 'T': tspfname = boptarg; break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]); return 1;
        }
    }

    if (boptind < ac) {
        cycfname = av[boptind++];
    } else {
        usage (av[0]); return 1;
    }

    return 0;
}

static char *get_problabel (const char *probloc)
{
    const char *p;
    const char *problabel = probloc;
    char *probcopy = (char *) NULL;
    char *p2;

    p = CCutil_strrchr_c (problabel, ':');
    if (p != (const char *) NULL) problabel = p+1;
    p = CCutil_strrchr_c (problabel, '/');
    if (p != (const char *) NULL) problabel = p+1;
    probcopy = CCutil_strdup (problabel);
    if (probcopy == (char *) NULL) return (char *) NULL;
    p2 = CCutil_strchr (probcopy, '.');
    if (p2 != (char *) NULL) *p2 = '\0';
    return probcopy;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] tour_file\n", fname);
    fprintf (stderr, "   -t    tour file in concorde format (default TSPLIB)\n");
    fprintf (stderr, "   -o f  output file (for the edge list)\n");
    fprintf (stderr, "   -S    write a TSPLIB tour file with tour length\n");
    fprintf (stderr, "   -T f  TSPLIB (dat) file to specify lengths\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON 20=ROAD\n");
}
