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
/*     A PROGRAM TO COMPUTE SUBDIVIDE A TSP INSTANCE FOR PARALLEL CODE      */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: September 18, 2017 (starting from tdivide.c)                      */
/*                                                                          */
/*  Splits EUCLID_3D instance into subproblems for cutting planes.          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "kdtree.h"
#include "macrorus.h"
#include "tsp.h"
#include "linkern.h"

static char *masterfname = (char *) NULL;
static char *probname  = (char *) NULL;
static int seed = 0;
static int bound_part = 0;
static int bound_count = 0;

int
    main (int ac, char **av);

static int
    build_partition (int ncount, CCdatagroup *dat, int partsize,
        int *p_scount, CCsubdiv **p_sublist, int ***p_partlist),
    create_subproblems (int ncount, CCdatagroup *dat, int *invtour,
        char *pname, int tcount, CCsubdiv *trac, int **slist),
    build_subproblem (CCdatagroup *dat, int *invtour, char *pname, int id,
        int scount, int *slist, char *hit),
    parseargs (int ac, char **av);

static char
    *get_problabel (const char *probloc);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, ncount, pcount = 0, rval = 0, inorm = -1;
    int *tour = (int *) NULL, *invtour = (int *) NULL;
    int **plist = (int **) NULL;
    char *name = (char *) NULL;
    CCdatagroup dat, permdat;
    CCsubdiv *trac = (CCsubdiv *) NULL;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);
    CCutil_init_datagroup (&permdat);
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname) {
        fprintf (stderr, "master file must be specified\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (!bound_part && !bound_count) {
        fprintf (stderr, "Must use -k or -n to specify # of subproblems\n");
        usage (av[0]);
        goto CLEANUP;
    }

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")
    CCutil_printlabel ();

    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    if (probname) {
        name = get_problabel (probname);
    } else {
        name = get_problabel (masterfname);
    }

    printf ("Name: %s\n", name); fflush (stdout);

    rval = CCutil_getmaster (masterfname, &ncount, &permdat, &tour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    CCutil_dat_getnorm (&permdat, &inorm);
    if (inorm != CC_EUCLIDEAN_3D) {
        fprintf (stderr, "Only set up for CC_EUCLIDEAN_3D norm\n");
        rval = 1;  goto CLEANUP;
    }

    CC_MALLOC (invtour, ncount, int);
    for (i = 0; i < ncount; i++) invtour[tour[i]] = i;

    /* Subdivision uses unpermuted data, so unpermute the dat struct */

    rval = CCutil_copy_datagroup (ncount, &permdat, &dat);
    CCcheck_rval (rval, "CCutil_copy_datagroup failed");

    rval = CCutil_datagroup_perm (ncount, &dat, invtour); 
    CCcheck_rval (rval, "CCutil_datagroup_perm failed");

    if (bound_count > 0) {
        bound_part = (int) (((double) ncount / (double) bound_count) + 0.5);
    }

    rval = build_partition (ncount, &dat, bound_part, &pcount, &trac, &plist);
    CCcheck_rval (rval, "build_partition failed");

    rval = create_subproblems (ncount, &dat, invtour, name, pcount, trac,
                               plist);
    CCcheck_rval (rval, "create_subproblems failed");

    rval = CCutil_write_subdivision_index (name, ncount, pcount, trac);
    CCcheck_rval (rval, "CCutil_write_subdivision_index failed");

CLEANUP:
    CC_IFFREE (tour, int);
    CC_IFFREE (invtour, int);
    CC_IFFREE (name, char);
    CC_IFFREE (trac, CCsubdiv);
    CCutil_freedatagroup (&dat);
    CCutil_freedatagroup (&permdat);
    if (plist) {
        for (i = 0; i < pcount; i++) { CC_IFFREE (plist[i], int); }
        CC_FREE (plist, int *);
    }
    return rval;
}

static int create_subproblems (int ncount, CCdatagroup *dat, int *invtour,
        char *pname, int tcount, CCsubdiv *trac, int **slist)
{
    int i, rval = 0;
    char *hit = (char *) NULL;

    CC_MALLOC (hit, ncount, char);
    for (i = 0; i < ncount; i++) hit[i] = 0;

    for (i = 0; i < tcount; i++) {
        rval = build_subproblem (dat, invtour, pname, trac[i].id, trac[i].cnt,
                                 slist[i], hit); 
        CCcheck_rval (rval, "build_subproblem failed");
    }

    for (i = 0; i < ncount; i++) {
        if (hit[i] == 0) {
            fprintf (stderr, "missed a node in partition\n");
            rval = 1; goto CLEANUP;
        }
    }

CLEANUP:
    CC_IFFREE (hit, char);
    return rval;
}

static int build_subproblem (CCdatagroup *dat, int *invtour, char *pname,
        int id, int scount, int *slist, char *hit)
{
    int rval = 0, i;
    int *tpos = (int *) NULL, *perm = (int *) NULL, *p_slist = (int *) NULL;
    char buf[CCutil_FILE_NAME_LEN];
    CC_SFILE *out = (CC_SFILE *) NULL;
    FILE *fout = (FILE *) NULL;

    sprintf (buf, "%s_%d.mas", pname, id);
    printf ("Create %s, scount = %d\n", buf, scount);
    fflush (stdout);

    if (hit) {
        for (i = 0; i < scount; i++) {
            if (hit[slist[i]]) {
                fprintf (stderr, "duplicate node in partition");
                rval = 1;  goto CLEANUP;
            }
            hit[slist[i]] = 1;
        }
    }

    /* Build the permutation from the exisiting tour */

    CC_MALLOC (tpos, scount, int);
    CC_MALLOC (perm, scount, int);
    for (i = 0; i < scount; i++) {
        tpos[i] = invtour[slist[i]];
        perm[i] = i;
    }
    CCutil_int_perm_quicksort (perm, tpos, scount);
    CC_FREE (tpos, int);

    CC_MALLOC (p_slist, scount, int);
    for (i = 0; i < scount; i++) p_slist[i] = slist[perm[i]];

    /* The header information */

    out = CCutil_sopen (buf, "w");
    if (out == (CC_SFILE *) NULL) {
        fprintf (stderr, "Could not open %s for output\n", buf);
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_swrite_int (out, scount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (out, CC_MASTER_DAT);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    rval = CCutil_swrite_int (out, CC_EUCLIDEAN_3D);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    /* The permuted datagroup */

    for (i = 0; i < scount; i++) {
        rval = CCutil_swrite_double (out, dat->x[p_slist[i]]);
        CCcheck_rval (rval, "CCutil_swrite_double failed");
        rval = CCutil_swrite_double (out, dat->y[p_slist[i]]);
        CCcheck_rval (rval, "CCutil_swrite_double failed");
        rval = CCutil_swrite_double (out, dat->z[p_slist[i]]);
        CCcheck_rval (rval, "CCutil_swrite_double failed");
    }

    /* The permutation */

    for (i = 0; i < scount; i++) {
        rval = CCutil_swrite_int (out, perm[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    /* Original node names (not permuted) */

    /* Write a names file with same info as above */
    sprintf (buf, "%s_%d.nam", pname, id);
    printf ("Create %s\n", buf); fflush (stdout);

    fout = fopen (buf, "w");
    if (fout == (FILE *) NULL) {
        fprintf (stderr, "Could not open %s for output\n", buf);
        rval = 1;  goto CLEANUP;
    }
    fprintf (fout, "%d\n", scount);
    for (i = 0; i < scount; i++) {
        fprintf (fout, "%d\n", slist[i]);
    }
    fclose (fout);
    fout = (FILE *) NULL;

CLEANUP:
    CC_IFFREE (tpos, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (p_slist, int);
    if (out) CCutil_sclose (out);
    if (fout) fclose (fout);
    return rval;
}

static int build_partition (int ncount, CCdatagroup *dat, int partsize,
        int *p_scount, CCsubdiv **p_sublist, int ***p_partlist)
{
    int i, k, rval = 0, totalcnt = 0, cnt, scount, ind;
    int **partlist = (int **) NULL, *perm = (int *) NULL;
    CCsubdiv *sublist = (CCsubdiv *) NULL;

    *p_scount = 0;
    *p_sublist = (CCsubdiv *) NULL;
    *p_partlist = (int **) NULL;

    if (dat->norm != CC_EUCLIDEAN_3D) {
        fprintf (stderr, "Only set up for CC_EUCLIDEAN_3D norm\n");
        rval = 1;  goto CLEANUP;
    }

    if (ncount < 2*partsize) {
        fprintf (stderr, "two few nodes for partition size\n");
        rval = 1;  goto CLEANUP;
    }

    CC_MALLOC (perm, ncount, int);
    for (i = 0; i < ncount; i++) perm[i] = i;
    CCutil_double_perm_quicksort (perm, dat->z, ncount);

    scount = ncount / partsize;

    CC_MALLOC (partlist, scount, int *);
    for (i = 0; i < scount; i++) partlist[i] = (int *) NULL;
    CC_MALLOC (sublist, scount, CCsubdiv);

    for (ind = 0, k = 0; k < scount; k++) {
        if (k < scount-1) cnt = partsize;
        else              cnt = ncount - totalcnt;

        CC_MALLOC (partlist[k], cnt, int);
        for (i = 0; i < cnt; i++) {
            partlist[k][i] = perm[ind++];
        }

        sublist[k].id = k;
        sublist[k].xrange[0] = 0.0;
        sublist[k].xrange[1] = 0.0;
        sublist[k].yrange[0] = 0.0;
        sublist[k].yrange[1] = 0.0;
        sublist[k].cnt = cnt;
        sublist[k].bound = -1.0;
        totalcnt += cnt;
    }

    *p_scount = scount;
    *p_sublist = sublist;
    *p_partlist = partlist;

CLEANUP:
    CC_IFFREE (perm, int);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "k:M:n:P:s:?", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'k': bound_part = atoi(boptarg); break;
        case 'M': masterfname  = boptarg; break;
        case 'n': bound_count = atoi(boptarg); break;
        case 'P': probname = boptarg; break;
        case 's': seed = atoi (boptarg); break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        fprintf (stderr, "extra items\n");
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below]\n", fname);
    fprintf (stderr, "   -k #  create bound subproblems with bucketsize #\n");
    fprintf (stderr, "   -n #  create bound subproblems with bucketsize ncities/#\n");
    fprintf (stderr, "   -M f  specify a master file\n");
    fprintf (stderr, "   -P s  specify a name for the output master files\n");
    fprintf (stderr, "   -s #  random seed\n");
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

