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
/*  Date: March 12, 2002                                                    */
/*                                                                          */
/*  Splits the instance into subproblems for lb, lkh, or linkern.  For      */
/*  lb, code uses m depot nodes (m is set to 4*sqrt(nodes) by default).     */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "kdtree.h"
#include "macrorus.h"
#include "tsp.h"
#include "linkern.h"

static char *tourfname = (char *) NULL;
static char *masterfname = (char *) NULL;
static char *tspfname  = (char *) NULL;
static char *savfname  = (char *) NULL;
static char *probname  = (char *) NULL;
static char *in_subdiv_name = (char *) NULL;
static int simpletour = 0;
static int depotcount = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int seed = 0;
static int bound_part = 0;
static int bound_count = 0;
static int lkh_part = 0;
static int tour_part = 0;
static int first_bucket = 0;
static int set_for_cuts = 0;
static int global_norm = CC_EUCLIDEAN;

#define RCMULT 100.0

typedef struct rcinfo {
    CCtsp_lp *lp;
    int ecount;     
    int *elist;      /* good representation, like quad-nearest 25 */
    int *rlen;       /* rc*RCMULT rounded to an int */
} rcinfo;

int
    main (int ac, char **av);

static int
    create_tour_problems (CCdatagroup *dat, char *pname, int tcount,
        CCsubdiv *trac, int **slist),
    build_tour_subproblem (CCdatagroup *dat, char *pname, int id,
        int scount, int *slist),
    create_subproblems (int ncount, CCdatagroup *dat, int *invtour,
        char *pname, int ndepot, int tcount, CCsubdiv *trac, int **slist,
        CCrandstate *rstate, int inorm, int just_cuts, rcinfo *R),
    build_subproblem (CCdatagroup *dat, CCkdtree *kt, int *invtour,
        char *pname, int ndepot, int id, int scount, int *slist, char *hit,
        int inorm, CCsubdiv *sbox, rcinfo *R),
    lkh_subproblem (int id, CCdatagroup *dat, int scount, int *sub, int start,
        CCsubdiv_lkh *plist, char *pname),
    geom_len (double x00, double y00, double x11, double y11, int *len),
    create_lkh (char *name, int ncount, CCdatagroup *dat, int *tour,
        int bucketsize, int first, double tourlen),
    border_opt (int ncount, CCdatagroup *dat, char *name, int pcount,
        CCsubdiv *intrac, int *invtour, CCrandstate *rstate),
    get_names (char *fname, int *names, int ncount),
    find_tour (int ncount, CCdatagroup *dat, int *perm, CCrandstate *rstate),
    build_rcinfo (rcinfo *R, int ncount, CCdatagroup *permdat,
       int *tour, int *invtour, char *fname, CCrandstate *rstate, int pcount,
       CCsubdiv *trac, int **plist),
    build_rc_edges (int ncount, CCdatagroup *dat, CCtsp_lp *lp,
       CCrandstate *rstate, int *outcount, int **outlist, int **outlen,
       int *cutmarks),
    parseargs (int ac, char **av);

static char
    *get_problabel (const char *probloc);

static void
    no_overlap (int tcount, CCsubdiv *trac, int *yesno),
    in_box (double x, double y, CCsubdiv *b, int *yesno, int i, int j,
        CCsubdiv *trac),
    init_rcinfo (rcinfo *R),
    free_rcinfo (rcinfo *R),
    usage (char *fname);


int main (int ac, char **av)
{
    int i, lap, ncount, pcount, rval = 0;
    int *tour = (int *) NULL, *invtour = (int *) NULL, **plist = (int **) NULL;
    double val;
    char *name = (char *) NULL,buf[1024];
    CCdatagroup dat, permdat;
    CCsubdiv *trac = (CCsubdiv *) NULL;
    CCsubdiv *in_trac = (CCsubdiv *) NULL;
    rcinfo *R = (rcinfo *) NULL;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);
    CCutil_init_datagroup (&permdat);
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!tspfname && !masterfname) {
        fprintf (stderr, "TSPLIB/DAT file or master file must be specified\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (set_for_cuts && (!bound_part && !bound_count)) {
        fprintf (stderr, "Must specify bound subproblems to set up for cuts\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (in_subdiv_name && (!tourfname && !masterfname)) {
        fprintf (stderr, "Must specify tour or master for border crossings\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (tourfname && masterfname) {
        fprintf (stderr, "Should not specify both tour and master\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (savfname && !masterfname) {
        fprintf (stderr, "Must specify master if using an LP sav file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (savfname && (!bound_part && !bound_count)) {
        fprintf (stderr, "Must specify bound subproblems to use LP sav file\n");
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
    } else if (tspfname) {
        name = get_problabel (tspfname);
    } else {
        name = get_problabel (masterfname);
    }

    printf ("Name: %s\n", name); fflush (stdout);

    if (masterfname) {
        rval = CCutil_getmaster (masterfname, &ncount, &permdat, &tour);
        CCcheck_rval (rval, "CCutil_getmaster failed")
        CCutil_dat_getnorm (&permdat, &global_norm);

        CC_MALLOC (invtour, ncount, int);
        for (i = 0; i < ncount; i++) invtour[tour[i]] = i;

        /* Subdivision uses unpermuted data, so unpermute the dat struct */

        rval = CCutil_copy_datagroup (ncount, &permdat, &dat);
        CCcheck_rval (rval, "CCutil_copy_datagroup failed");

        rval = CCutil_datagroup_perm (ncount, &dat, invtour); 
        CCcheck_rval (rval, "CCutil_datagroup_perm failed");
    } else {
        if (tsplib_in) {
            rval = CCutil_gettsplib (tspfname, &ncount, &dat);
            CCcheck_rval (rval, "CCutil_gettsplib failed");
            CCutil_dat_getnorm (&dat, &global_norm);
        } else {
            rval = CCutil_getdata (tspfname, binary_in, global_norm, &ncount,
                                   &dat, 0, 0, &rstate);
            CCcheck_rval (rval, "CCutil_getdata failed");
        }
    }

    if ((global_norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE &&
        (global_norm & CC_NORM_SIZE_BITS) != CC_D3_NORM_SIZE) {
        fprintf (stderr, "Only set up for 2D and 3D norms\n");
        rval = 1;  goto CLEANUP;
    }

    if (in_subdiv_name && ((global_norm != CC_EUCLIDEAN) && 
                           (global_norm != CC_MANNORM) &&
                           (global_norm != CC_GEOM))) {
        fprintf (stderr, "Border optimization only Euc/Man/Geom norms\n");
        rval = 1;  goto CLEANUP;
    }

    if (!tour && (in_subdiv_name || lkh_part > 0 || bound_part > 0 ||
                                                    bound_count > 0)) {
        CC_MALLOC (tour, ncount, int);
        if (tourfname) {
            if (simpletour) {
                rval = CCutil_getcycle (ncount, tourfname, tour, 0);
                CCcheck_rval (rval, "CCutil_getcycle failed");
            } else {
                rval = CCutil_getcycle_tsplib (ncount, tourfname, tour);
                CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");
            }
        } else {
            rval = find_tour (ncount, &dat, tour, &rstate);
            CCcheck_rval (rval, "find_tour failed");
        }

        CCutil_cycle_len (ncount, &dat, tour, &val);
        printf ("Tour Length: %.0f\n", val); fflush (stdout);

        /* Write a global master file for use in tcuts */

        rval = CCutil_copy_datagroup (ncount, &dat, &permdat);
        CCcheck_rval (rval, "CCutil_copy_datagroup failed");

        rval = CCutil_datagroup_perm (ncount, &permdat, tour);
        CCcheck_rval (rval, "CCutil_datagroup_perm failed");

        sprintf (buf, "%s.mas", name);
        printf ("Writing global master file: %s\n", buf); fflush (stdout);

        rval = CCutil_putmaster (buf, ncount, &permdat, tour);
        CCcheck_rval (rval, "CCutil_putmaster failed");
    }

    if (tour && !invtour) {
        CC_MALLOC (invtour, ncount, int);
        for (i = 0; i < ncount; i++) invtour[tour[i]] = i;
    }

    if (in_subdiv_name) {
        int incount, ipcount;
        char *iname = (char *) NULL;

        printf ("Build cross-border subproblems\n"); fflush (stdout);

        val = CCutil_read_subdivision_index (in_subdiv_name, &iname, &incount,
                                             &ipcount, &in_trac);
        CCcheck_rval (rval, "CCutil_read_subdivision_index failed");

        printf ("Name: %s\n", iname);
        printf ("ncount = %d, partitions = %d\n", incount, ipcount);
        fflush (stdout);

        rval = border_opt (incount, &dat, iname, ipcount, in_trac, invtour,
                           &rstate);
        CCcheck_rval (rval, "border_opt failed");

        CC_IFFREE (iname, char);
        goto CLEANUP;
    }

    if (bound_count > 0) {
        bound_part = (int) (((double) ncount / (double) bound_count) + 0.5);
    }

    if (lkh_part > 0) {
        if (!tour) {
            fprintf (stderr, "Need to specify a tour for LKH subproblems\n");
            rval = 1;  goto CLEANUP; 
        }
        rval = create_lkh (name, ncount, &dat, tour, lkh_part, first_bucket,
                           val);
        CCcheck_rval (rval, "create_lkh failed");
    } else if (bound_part > 0) {
        if ((global_norm & CC_NORM_BITS) != CC_KD_NORM_TYPE &&
                             global_norm != CC_GEOM &&
                             global_norm != CC_EUCLIDEAN_3D) {
            fprintf (stderr, "Only set for GEOM, KD-tree, and EUC_3D norms\n");
            rval = 1;  goto CLEANUP;
        }
        if (!tour) {
            fprintf (stderr, "Need to specify a tour for bound subproblems\n");
            rval = 1;  goto CLEANUP; 
        }
        rval = CCutil_karp_partition (ncount, &dat, bound_part, &pcount,
                                      &trac, &plist, &rstate, 1);
        CCcheck_rval (rval, "CCutil_karp_partition failed");

        no_overlap (pcount, trac, &lap);
        if (lap == 1) {
            fprintf (stderr, "Error: regions overlap\n");
            rval = 1;  goto CLEANUP;
        }

        if (savfname) {
            CC_MALLOC (R, 1, rcinfo);
            init_rcinfo (R);
            rval = build_rcinfo (R, ncount, &permdat, tour, invtour, savfname,
                                 &rstate, pcount, trac, plist);
            CCcheck_rval (rval, "build_rcinfo failed");
            /* Note: maybe add a depot to capture border conditions?  */
        }

        rval = create_subproblems (ncount, &dat, invtour, name, depotcount,
                  pcount, trac, plist, &rstate, global_norm, set_for_cuts, R);
        CCcheck_rval (rval, "create_subproblems failed");

        rval = CCutil_write_subdivision_index (name, ncount, pcount, trac);
        CCcheck_rval (rval, "CCutil_write_subdivision_index failed");
    } else if (tour_part > 0)  {
        rval = CCutil_karp_partition (ncount, &dat, tour_part, &pcount,
                                      &trac, &plist, &rstate, 1);
        CCcheck_rval (rval, "CCutil_karp_partition failed");

        rval = create_tour_problems (&dat, name, pcount, trac, plist);
        CCcheck_rval (rval, "create_tour_problems");

        rval = CCutil_write_subdivision_index (name, ncount, pcount, trac);
        CCcheck_rval (rval, "CCutil_write_subdivision_index failed");
    } else {
        fprintf (stderr, "No partition type specified\n");
        goto CLEANUP;
    }

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
    if (R) {
        free_rcinfo (R);
        CC_IFFREE (R, rcinfo);
    }
    return rval;
}

static int create_tour_problems (CCdatagroup *dat, char *pname,
       int tcount, CCsubdiv *trac, int **slist)
{
    int i, rval = 0;

    for (i = 0; i < tcount; i++) {
        rval = build_tour_subproblem (dat, pname, trac[i].id, trac[i].cnt,
                                      slist[i]); 
        CCcheck_rval (rval, "build_tour_subproblem failed");
    }

CLEANUP:

    return rval;
}

static int build_tour_subproblem (CCdatagroup *dat, char *pname, int id,
        int scount, int *slist)
{
    int i, inorm, rval = 0;
    char buf[1024];
    FILE *out = (FILE *) NULL;
    CCdatagroup sdat;

    CCutil_init_datagroup (&sdat);
    CCutil_dat_getnorm (dat, &inorm);
    CCutil_dat_setnorm (&sdat, inorm);

    CC_MALLOC (sdat.x, scount, double);
    CC_MALLOC (sdat.y, scount, double);

    for (i = 0; i < scount; i++) {
        sdat.x[i] = dat->x[slist[i]];
        sdat.y[i] = dat->y[slist[i]];
    }

/*
    sprintf (buf, "%s_%d.dat", pname, id);
    printf ("Create %s\n", buf); fflush (stdout);
    rval = CCutil_writedata (buf, 0, scount, &sdat);
    CCcheck_rval (rval, "CCutil_writedata failed");
*/
    sprintf (buf, "%s_%d.tsp", pname, id);
    printf ("Create %s\n", buf); fflush (stdout);
    rval = CCutil_writetsplib (buf, scount, &sdat, "subproblem");
    CCcheck_rval (rval, "CCutil_writetsplib failed");

    sprintf (buf, "%s_%d.nam", pname, id);
    printf ("Create %s\n", buf); fflush (stdout);

    out = fopen (buf, "w");
    if (out == (FILE *) NULL) {
        fprintf (stderr, "Could not open %s for output\n", buf);
        rval = 1;  goto CLEANUP;
    }
    fprintf (out, "%d\n", scount);
    for (i = 0; i < scount; i++) { fprintf (out, "%d\n", slist[i]); }

CLEANUP:
    if (out) fclose (out);
    CCutil_freedatagroup (&sdat);
    return rval;
}

static void no_overlap (int tcount, CCsubdiv *trac, int *yesno)
{
    int i, j;
    CCsubdiv *p, *q;
    double x00, x11, y00, y11;

    *yesno = 0;

    for (i = 0; i < tcount; i++) {
        for (j = i+1; j < tcount; j++) {
            p = &trac[i];
            q = &trac[j];

            x00 = p->xrange[0]; x11 = p->xrange[1];
            y00 = p->yrange[0]; y11 = p->yrange[1];
            in_box (x00, y00, q, yesno, i, j, trac); if (*yesno) return;
            in_box (x00, y11, q, yesno, i, j, trac); if (*yesno) return;
            in_box (x11, y00, q, yesno, i, j, trac); if (*yesno) return;
            in_box (x11, y11, q, yesno, i, j, trac); if (*yesno) return;

            x00 = q->xrange[0]; x11 = q->xrange[1];
            y00 = q->yrange[0]; y11 = q->yrange[1];
            in_box (x00, y00, p, yesno, i, j, trac); if (*yesno) return;
            in_box (x00, y11, p, yesno, i, j, trac); if (*yesno) return;
            in_box (x11, y00, p, yesno, i, j, trac); if (*yesno) return;
            in_box (x11, y11, p, yesno, i, j, trac); if (*yesno) return;
        }
    }
    printf ("Boxes do not overlap -- good\n"); fflush (stdout);
}

static void in_box (double x, double y, CCsubdiv *b, int *yesno, int i, int j,
       CCsubdiv *trac)
{
    if (x > b->xrange[0] && x < b->xrange[1] && 
        y > b->yrange[0] && y < b->yrange[1]) {
        fprintf (stderr, "Box %d (%f, %f, %f, %f) intersects\n",
            i, trac[i].xrange[0], trac[i].xrange[1], trac[i].yrange[0],
            trac[i].yrange[1]);
        fprintf (stderr, "Box %d (%f, %f, %f, %f)\n",
            j, trac[j].xrange[0], trac[j].xrange[1], trac[j].yrange[0],
            trac[j].yrange[1]);
        *yesno = 1;
    } else {
        *yesno = 0;
    }
}

static int create_subproblems (int ncount, CCdatagroup *dat, int *invtour,
        char *pname, int ndepot, int tcount, CCsubdiv *trac, int **slist,
        CCrandstate *rstate, int inorm, int just_cuts, rcinfo *R)
{
    int i, rval = 0;
    char *hit = (char *) NULL;
    CCkdtree kt, *p_kt = (CCkdtree *) NULL;

    if ((inorm != CC_GEOM) && (just_cuts == 0)) {
        rval = CCkdtree_build (&kt, ncount, dat, (double *) NULL, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        p_kt = &kt;
    }

    if (just_cuts) {
        if (ndepot != 0) {
            printf ("Warning:  No depots will be assigned\n"); fflush (stdout);
        }
        ndepot = -1;
    }

    CC_MALLOC (hit, ncount, char);
    for (i = 0; i < ncount; i++) hit[i] = 0;

    for (i = 0; i < tcount; i++) {
        rval = build_subproblem (dat, p_kt, invtour, pname, ndepot,
                           trac[i].id, trac[i].cnt, slist[i], hit, inorm,
                           &trac[i], R); 
        CCcheck_rval (rval, "build_subproblem failed");
    }

    for (i = 0; i < ncount; i++) {
        if (hit[i] == 0) {
            fprintf (stderr, "missed a node in partition\n");
            rval = 1; goto CLEANUP;
        }
    }

CLEANUP:
    if (p_kt) CCkdtree_free (p_kt);
    CC_IFFREE (hit, char);
    return rval;
}

static int build_subproblem (CCdatagroup *dat, CCkdtree *kt, int *invtour,
        char *pname, int ndepot, int id, int scount, int *slist, char *hit,
        int inorm, CCsubdiv *sbox, rcinfo *R)
{
    int rval = 0, i, j, neighbor, len, t;
    int *tpos = (int *) NULL, *perm = (int *) NULL, *p_slist = (int *) NULL;
    char buf[CCutil_FILE_NAME_LEN];
    double xnn, ynn;
    CC_SFILE *out = (CC_SFILE *) NULL;
    FILE *fout = (FILE *) NULL;

    /* Note:  If ndepot == -1 then create cut problems */

    if (ndepot == 0) {
        ndepot = (int) sqrt ((double) scount);
        ndepot = 4*(ndepot + 1);
    }

    for (i = 0; i < scount; i++) {
        xnn = dat->x[slist[i]];
        ynn = dat->y[slist[i]];

        if (xnn < sbox->xrange[0] || xnn > sbox->xrange[1] ||
            ynn < sbox->yrange[0] || ynn > sbox->yrange[1]) {
            fprintf (stderr, "Point (%f,%f) not in Box (%f, %f, %f, %f)\n",
                  xnn, ynn, sbox->xrange[0], sbox->xrange[1], sbox->yrange[0],
                  sbox->yrange[1]);
            rval = 1;  goto CLEANUP;
        }
    }

    sprintf (buf, "%s_%d.mas", pname, id);
    printf ("Create %s, with %d depots, scount = %d\n", buf, ndepot, scount);
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

    if (ndepot != -1) {
        rval = CCutil_swrite_int (out, CC_SUBDIVISION);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        if (R) {
            rval = CCutil_swrite_int (out, CC_SPARSE);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        } else {
            rval = CCutil_swrite_int (out, inorm);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
        rval = CCutil_swrite_int (out, ndepot);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    } else {
        if (R) {
            rval = CCutil_swrite_int (out, CC_SPARSE);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        } else {
            rval = CCutil_swrite_int (out, inorm);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
    }

    if (ndepot != -1) {
        /* The depot costs */

        if (inorm != CC_GEOM) {
            for (i = 0; i < scount; i++) {
                CCkdtree_delete (kt, slist[i]);
            }
        }

        for (i = 0; i < scount; i++) {
            if (inorm == CC_GEOM) {
                xnn = dat->x[p_slist[i]];
                ynn = dat->y[p_slist[i]];
    
                len = CCutil_MAXINT;
                rval = geom_len (xnn, ynn, sbox->xrange[0], ynn, &t);
                CCcheck_rval (rval, "geom_len failed");
                if (t < len) len = t;

                rval = geom_len (xnn, ynn, sbox->xrange[1], ynn, &t);
                CCcheck_rval (rval, "geom_len failed");
                if (t < len) len = t;

                rval = geom_len (xnn, ynn, xnn, sbox->yrange[0], &t);
                CCcheck_rval (rval, "geom_len failed");
                if (t < len) len = t;

                rval = geom_len (xnn, ynn, xnn, sbox->yrange[1], &t);
                CCcheck_rval (rval, "geom_len failed");
                if (t < len) len = t;
            } else {
                CCkdtree_undelete (kt, p_slist[i]);
                neighbor = CCkdtree_node_nearest (kt, p_slist[i], dat,
                                                 (double *) NULL);
                len = CCutil_dat_edgelen (p_slist[i], neighbor, dat);
                len = len / 2;    /* Since edge might get hit in two regions */
                CCkdtree_delete (kt, p_slist[i]);
            }
            rval = CCutil_swrite_int (out, len);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }

        if (inorm != CC_GEOM) {
            for (i = 0; i < scount; i++) {
                CCkdtree_undelete (kt, slist[i]);
            }
        }
    }

    /* The permuted datagroup */

    if (R) {
        /* Gather the rc edges within this part of the partition */
        /* and output a sparse datgroup with this set            */

        int icount = 0, *imarks = (int *) NULL;
        int *ilist = (int *) NULL, *ilen = (int *) NULL;
        int *invp_slist = (int *) NULL;
        int ncount = R->lp->graph.ncount;
        CCdatagroup idat;

        CCutil_init_datagroup (&idat);

        CC_MALLOC (invp_slist, ncount, int);
        for (i = 0; i < ncount; i++) invp_slist[i] = -1;
        for (i = 0; i < scount; i++) invp_slist[p_slist[i]] = i;
       
        CC_MALLOC (imarks, ncount, int);
        CC_MALLOC (ilist, 2*R->ecount, int);
        CC_MALLOC (ilen, R->ecount, int);

        for (i = 0; i < ncount; i++) imarks[i] = 0;
        for (i = 0; i < scount; i++) imarks[slist[i]] = 1;
        for (i = 0; i < R->ecount; i++) {
            if (imarks[R->elist[2*i]] == 1 && imarks[R->elist[2*i+1]] == 1) {
                if (invp_slist[R->elist[2*i]] == -1 ||
                    invp_slist[R->elist[2*i+1]] == -1) {
                    fprintf (stderr, "ERROR: bad index\n");
                    rval = 1; goto CLEANUP;
                }
                ilist[2*icount]   = invp_slist[R->elist[2*i]];
                ilist[2*icount+1] = invp_slist[R->elist[2*i+1]];
                ilen[icount] = R->rlen[i];
                icount++;
            }
        }
        printf ("Local RC Count: %d\n", icount); fflush (stdout);

        rval = CCutil_graph2dat_sparse (scount, icount, ilist, ilen, 10000,
                                        &idat);
        CCcheck_rval (rval, "CCutil_graph2dat_sparse failed");

        /*  Code from UTIL/getdata/writedat_sparse */

        if (CCutil_swrite_int (out, icount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
        }

        for (i = 0; i < scount; i++) {
            if (CCutil_swrite_int (out, idat.degree[i])) {
                fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
            }
            for (j = 0; j < idat.degree[i]; j++) {
                if (CCutil_swrite_int (out, idat.adj[i][j])) {
                    fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
                }
                if (CCutil_swrite_int (out, idat.len[i][j])) {
                    fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
                }
            }
        }

        CC_IFFREE (invp_slist, int);
        CC_IFFREE (imarks, int);
        CC_IFFREE (ilist, int);
        CC_IFFREE (ilen, int);
        CCutil_freedatagroup (&idat);
    } else {
        for (i = 0; i < scount; i++) {
            rval = CCutil_swrite_double (out, dat->x[p_slist[i]]);
            CCcheck_rval (rval, "CCutil_swrite_double failed");
            rval = CCutil_swrite_double (out, dat->y[p_slist[i]]);
            CCcheck_rval (rval, "CCutil_swrite_double failed");
            if (inorm == CC_EUCLIDEAN_3D) {
                rval = CCutil_swrite_double (out, dat->z[p_slist[i]]);
                CCcheck_rval (rval, "CCutil_swrite_double failed");
            }
        }
    }

    /* The permutation */

    for (i = 0; i < scount; i++) {
        rval = CCutil_swrite_int (out, perm[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    /* Original node names (not permuted) */

    if (ndepot != -1) {
        for (i = 0; i < scount; i++) {
            rval = CCutil_swrite_int (out, slist[i]);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
    } else {
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
    }

CLEANUP:
    CC_IFFREE (tpos, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (p_slist, int);
    if (out) CCutil_sclose (out);
    if (fout) fclose (fout);
    return rval;
}

static int geom_len (double x00, double y00, double x11, double y11, int *len)
{
    int rval = 0;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);
    CCutil_dat_setnorm (&dat, CC_GEOM);

    CC_MALLOC (dat.x, 2, double);
    CC_MALLOC (dat.y, 2, double);

    dat.x[0] = x00; dat.y[0] = y00;
    dat.x[1] = x11; dat.y[1] = y11;

    *len = CCutil_dat_edgelen (0, 1, &dat);

CLEANUP:
    CCutil_freedatagroup (&dat);
    return rval;
}

static int create_lkh (char *name, int ncount, CCdatagroup *dat, int *tour,
        int bucketsize, int first, double tourlen)
{
    int rval = 0;
    FILE *out = (FILE *) NULL;
    int *sub = (int *) NULL;
    int i, k, nsub, remain, extra, start;
    CCsubdiv_lkh *plist = (CCsubdiv_lkh *) NULL;

    if (first <= 0 || first >= ncount) first = bucketsize;

    nsub = 1 + ((ncount - first) / bucketsize);
    remain = ncount - first  - ((nsub-1) * bucketsize);
    extra = (remain > bucketsize/16 ? 1 : 0);

    if (first < bucketsize + remain) {
        CC_MALLOC (sub, bucketsize + remain, int);
    } else {
        CC_MALLOC (sub, first, int);
    }

    CC_MALLOC (plist, nsub + extra, CCsubdiv_lkh);
    printf ("Create %d LKH subproblems\n", nsub + extra);
    fflush (stdout);

    start = 0;
    for (i = 0; i < first; i++) sub[i] = tour[start+i];

    rval = lkh_subproblem (0, dat, first, sub, start, plist, name);
    CCcheck_rval (rval, "lkh_subproblem failed");
    start += first;

    for (k = 1; k < nsub-1; k++) {
        for (i = 0; i < bucketsize; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, bucketsize, sub, start, plist, name);
        CCcheck_rval (rval, "lkh_subproblem failed");
        start += bucketsize;
    }
    if (extra == 0) {
        for (i = 0; i < bucketsize+remain; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, bucketsize+remain, sub, start, plist,
                               name);
        CCcheck_rval (rval, "lkh_subproblem failed");
        start += (bucketsize+remain);

    } else {
        for (i = 0; i < bucketsize; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, bucketsize, sub, start, plist, name);
        CCcheck_rval (rval, "lkh_subproblem failed");
        start += bucketsize;
        k++;

        for (i = 0; i < remain; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, remain, sub, start, plist, name);
        CCcheck_rval (rval, "lkh_subproblem failed");
    }

    printf ("Subproblems\n");
    for (i = 0; i < nsub + extra; i++) {
        printf ("%d %d %d %.0f %.0f\n", plist[i].id, plist[i].cnt,
                        plist[i].start, plist[i].origlen, plist[i].newlen);
        fflush (stdout);
    }

    rval = CCutil_write_subdivision_lkh_index (name, ncount, nsub+extra,
                                               plist, tourlen);
    CCcheck_rval (rval, "CCutil_write_subdivision_lkh_index failed");

CLEANUP:
    CC_IFFREE (sub, int);
    CC_IFFREE (plist, CCsubdiv_lkh);
    if (out) fclose (out);
    return rval;
}

static int lkh_subproblem (int id, CCdatagroup *dat, int scount, int *sub,
        int start, CCsubdiv_lkh *plist, char *pname)
{
    int i;
    int rval = 0;
    double total = 0;


    for (i = 1; i < scount; i++) {
        total += ((double) CCutil_dat_edgelen (sub[i-1], sub[i], dat));
    }

    plist[id].id = id; 
    plist[id].cnt = scount;
    plist[id].start = start;
    plist[id].origlen = total;
    plist[id].newlen = -1.0;

    rval =  build_tour_subproblem (dat, pname, id, scount, sub);
    CCcheck_rval (rval, "build_tour_subproblem failed");

CLEANUP:
    return rval;
}

static int find_tour (int ncount, CCdatagroup *dat, int *perm,
        CCrandstate *rstate)
{
    int rval = 0, kicks, istour, silent = 0;
    int ecount, *elist = (int *) NULL;
    int tcount, *tlist = (int *) NULL;
    int *bestcyc = (int *) NULL, *cyc     = (int *) NULL;
    double val, szeit;
    CCedgegengroup plan;

    printf ("Finding a initial global tour...\n"); fflush (stdout);

    szeit = CCutil_zeit ();
    kicks = (ncount > 1000 ? 500 : ncount/2);

    CC_MALLOC (cyc, ncount, int);
    CC_MALLOC (bestcyc, ncount, int);

    CCedgegen_init_edgegengroup (&plan);
    plan.quadnearest = 2;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &ecount,
                            &elist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");
    plan.quadnearest = 0;

    plan.tour.greedy = 1;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &tcount,
                            &tlist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    if (tcount != ncount) {
        fprintf (stderr, "wrong edgeset from CCedgegen_edges\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_edge_to_cycle (ncount, tlist, &istour, cyc);
    CCcheck_rval (rval, "CCutil_edge_to_cycle failed");
    if (istour == 0) {
        fprintf (stderr, "Starting tour has an error\n");
        rval = 1; goto CLEANUP;
    }
    CC_FREE (tlist, int);

    rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
                    cyc, perm, &val, silent, 0.0, 0.0, (char *) NULL,
                    CC_LK_GEOMETRIC_KICK, rstate);
    CCcheck_rval (rval, "CClinkern_tour failed");

    printf ("Time to find global tour: %.2f (seconds)\n", CCutil_zeit()-szeit);
    fflush (stdout);

CLEANUP:
    CC_IFFREE (cyc, int);
    CC_IFFREE (bestcyc, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (tlist, int);
    return rval;
}

static void init_rcinfo (rcinfo *R)
{
    if (R) {
        R->lp = (CCtsp_lp *) NULL;
        R->ecount = 0;
        R->elist = (int *) NULL;
        R->rlen = (int *) NULL;
    }
} 

static void free_rcinfo (rcinfo *R)
{
    if (R) {
        CCtsp_free_tsp_lp_struct (&R->lp);
        R->ecount = 0;
        CC_IFFREE (R->elist, int);
        CC_IFFREE (R->rlen, int);
    }
} 

static int build_rcinfo (rcinfo *R, int ncount, CCdatagroup *permdat,
       int *tour, int *invtour, char *fname, CCrandstate *rstate, int pcount,
       CCsubdiv *trac, int **plist)
{
    int i, j, k, h, ptag, acount, cutcount, rval = 0, infeasible = 0;
    int *cutmarks = (int *) NULL, *pmarks = (int *) NULL;
    int *arr = (int *) NULL;
    CCtsp_lpcut_in C;
    CCtsp_lpclique *q;

    rval = CCtsp_init_lp (&R->lp, (char *) NULL, -1, fname, ncount,
             permdat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
             (int *) NULL, 0, tour, CCtsp_LP_MAXDOUBLE,
             (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 1, rstate,
             &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible");
        rval = 1; goto CLEANUP;
    }

    printf ("Read in the LP file %s\n", fname);
    printf ("Number of Cuts: %d\n", R->lp->cuts.cutcount);
    cutcount = R->lp->cuts.cutcount;


    /* Set pmarks to indicate the nodes in each part of the partition */
    /* The pmarks values correspond to the permuted data, to match LP */

    CC_MALLOC (pmarks, ncount, int);
    for (i = 0; i < ncount; i++) pmarks[i] = -1;
    for (k = 0; k < pcount; k++) {
        for (i = 0; i < trac[k].cnt; i++) {
            pmarks[invtour[plist[k][i]]] = k;
        }
    }
    for (i = 0, k = 0; i < ncount; i++) {
        if (pmarks[i] == -1) k++;
    }
    if (k > 0) {
        printf ("WARNING: %d nodes not in the partition\n", k);
        fflush (stdout);
    }

    /* Run through cuts and mark those with uniform pmarks */

    CC_MALLOC (cutmarks, cutcount, int);
    for (i = 0; i < cutcount; i++) cutmarks[i] = 0;

    for (i = 0; i < cutcount; i++) {
        rval = CCtsp_lpcut_to_lpcut_in (&R->lp->cuts, &R->lp->cuts.cuts[i], &C);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
        ptag = -1;
        for (j = 0; j < C.cliquecount; j++) {
            q = &(C.cliques[j]);
            rval = CCtsp_clique_to_array (q, &arr, &acount);
            CCcheck_rval (rval, "CCtsp_clique_to_array failed");
            if (ptag == -1) ptag = pmarks[arr[0]];
            for (k = 0; k < acount; k++) {
                if (pmarks[arr[k]] != ptag) {
                    goto TAGOUT;
                }
            }
            CC_IFFREE (arr, int);
        }

        for (j = 0; j < C.dominocount; j++) {
            for (h = 0; h < 2; h++) {
                q = &C.dominos[j].sets[h];
                rval = CCtsp_clique_to_array (q, &arr, &acount);
                CCcheck_rval (rval, "CCtsp_clique_to_array failed");
                if (ptag == -1) ptag = pmarks[arr[0]];
                for (k = 0; k < acount; k++) {
                    if (pmarks[arr[k]] != ptag) {
                        goto TAGOUT;
                    }
                }
                CC_IFFREE (arr, int);
            }
        }
        cutmarks[i] = 1;

TAGOUT:
        CCtsp_free_lpcut_in (&C);
    }

    /* cutmarks is used to skip over marked cuts in rc computation */

    rval = build_rc_edges (ncount, permdat, R->lp, rstate, &R->ecount,
                           &R->elist, &R->rlen, cutmarks);
    CCcheck_rval (rval, "build_rc_edges failed");

    /* edge ends are for permuted data, map back to original data */

    for (i = 0; i < R->ecount; i++) {
        R->elist[2*i]   = tour[R->elist[2*i]];
        R->elist[2*i+1] = tour[R->elist[2*i+1]];
    }

CLEANUP:
    CC_IFFREE (cutmarks, int);
    CC_IFFREE (pmarks, int);
    CC_IFFREE (arr, int);
    return rval;
}

static int build_rc_edges (int ncount, CCdatagroup *dat, CCtsp_lp *lp,
       CCrandstate *rstate, int *outcount, int **outlist, int **outlen,
       int *cutmarks)
{
    int i, len, rval = 0;
    int qcount, *qlist = (int *) NULL; 
    int ecount, *elist = (int *) NULL, *elen = (int *) NULL;
    double *rc = (double *) NULL;
    CCedgegengroup plan;
    CCutil_edgehash h;

    /* If lp->full_edges_valid == 1 then use it as the edge set */

    *outcount = 0;
    *outlist = (int *) NULL;
    *outlen = (int *) NULL;

    rval = CCutil_edgehash_init (&h, 2*ncount);
    CCcheck_rval (rval, "CCutil_edgehash_init failed");
    
    CCedgegen_init_edgegengroup (&plan);
    plan.quadnearest = 25;
    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &qcount,
                            &qlist, 0, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    for (i = 0; i < qcount; i++) {
        len = CCutil_dat_edgelen (qlist[2*i], qlist[2*i+1], dat);
        rval = CCutil_edgehash_set (&h, qlist[2*i], qlist[2*i+1], len);
        CCcheck_rval (rval, "CCutil_edgehash_set failed");
    }

    for (i = 0; i < lp->graph.ecount; i++) {
        len = CCutil_dat_edgelen (lp->graph.edges[i].ends[0],
                                  lp->graph.edges[i].ends[1], dat);
        rval = CCutil_edgehash_set (&h, lp->graph.edges[i].ends[0],
                                        lp->graph.edges[i].ends[1], len);
        CCcheck_rval (rval, "CCutil_edgehash_set failed");
    }

    rval = CCutil_edgehash_getall (&h, &ecount, &elist, &elen);
    CCcheck_rval (rval, "CCutil_edgehash_getall failed");

    rval = CCtsp_edgelist_to_genadj (ncount, ecount, elist, elen,
               &lp->fulladj, &lp->fulladjspace);
    CCcheck_rval (rval, "CCtsp_edgelist_to_genadj failed");
    lp->fullcount = ecount;

    rval = CCtsp_reduced_cost_all (lp, outcount, outlist, &rc, cutmarks);
    CCcheck_rval (rval, "CCtsp_reduced_cost_all failed");

    CC_MALLOC (*outlen, *outcount, int);
    for (i = 0; i < *outcount; i++) {
        (*outlen)[i] = (int) (rc[i] * RCMULT);
        if ((*outlen)[i] < 0) {
            /* Might not want negative cost in subproblem */
            (*outlen)[i] = 0;
        }
    }

CLEANUP:
    CC_IFFREE (qlist, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (rc, double);
    CCutil_edgehash_free (&h);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bBCD:f:I:k:j:l:m:M:n:o:s:S:tT:P:N:?", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'B': binary_in = 1; break;
        case 'b': binary_in = 2; break;
        case 'C': set_for_cuts = 1; break;
        case 'D': tspfname = boptarg; break;
        case 'f': first_bucket = atoi(boptarg); break;
        case 'I': in_subdiv_name = boptarg; break;
        case 'l': lkh_part = atoi(boptarg); break;
        case 'j': tour_part = atoi(boptarg); break;
        case 'k': bound_part = atoi(boptarg); break;
        case 'm': depotcount = atoi(boptarg); break;
        case 'M': masterfname  = boptarg; break;
        case 'n': bound_count = atoi(boptarg); break;
        case 'P': probname = boptarg; break;
        case 's': seed = atoi (boptarg); break;
        case 'S': savfname = boptarg; set_for_cuts = 1; break;
        case 't': simpletour = 1; break;
        case 'T': tourfname  = boptarg; break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: global_norm = CC_MAXNORM; break;
            case 1: global_norm = CC_MANNORM; break;
            case 2: global_norm = CC_EUCLIDEAN; break;
            case 3: global_norm = CC_EUCLIDEAN_3D; break;
            case 17: global_norm = CC_GEOM; break;
            case 18: global_norm = CC_EUCLIDEAN_CEIL; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
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
    fprintf (stderr, "   -b    datfile in double binary format\n");
    fprintf (stderr, "   -B    datfile in double integer format\n");
    fprintf (stderr, "   -D f  specify a TSPLIB/DAT file (if no master)\n");
    fprintf (stderr, "   -C    set bound subproblems for cuts (no depots)\n");
    fprintf (stderr, "   -f #  size of first bucket for LKH subproblems\n");
    fprintf (stderr, "   -I f  index file (to create border subprobs)\n");
    fprintf (stderr, "   -j #  create tour subproblems with bucketsize #\n");
    fprintf (stderr, "   -k #  create bound subproblems with bucketsize #\n");
    fprintf (stderr, "   -n #  create bound subproblems with bucketsize ncities/#\n");
    fprintf (stderr, "   -l #  create LKH subproblems with bucketsize #\n");
    fprintf (stderr, "   -m #  number of depots (default 4*sqrt(n)\n");
    fprintf (stderr, "   -M f  specify a master file\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -S f  specify sav file to use reduced-cost distances in subprobs\n");
    fprintf (stderr, "   -T f  specify a tour (needed for LKH, bound, and border, if no master)\n");
    fprintf (stderr, "   -t    tour file in concorde format (default TSPLIB)\n");
    fprintf (stderr, "   -P s  specify a name for the output master files\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 17=GEOM, 18=JOHNSON\n");
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


/**************************************************************************/
/**********************Border crossings from ttour.c  *********************/
/**************************************************************************/

static void subdiv_adj (CCsubdiv *s, CCsubdiv *t, int *len);
static int improve_border (char *name, int ip, int iq, CCdatagroup *dat,
        CCsubdiv *trac, int *invtour, int pnt, CCsubdiv *sbox, int *gotit,
        CCkdtree *kt, int id);
static void find_common_border (CCsubdiv *p, CCsubdiv *q, int *border,
        int *flip);
static int get_part_data (char *name, int ind, int *p_count, int **p_namelist);
static int grab_border_points (int btype, int scount, int *slist,
        CCdatagroup *dat, int *p_bcount, int **p_blist, double xlow,
        double xhi, int maxpoints);


static int border_opt (int ncount, CCdatagroup *dat, char *name, int pcount,
        CCsubdiv *intrac, int *invtour, CCrandstate *rstate)
{
    int rval = 0;
    int i, j, norm, len, gotit, tcount, mcnt = 0, cnt = 0;
    char buf[1024];
    CCsubdiv *trac = (CCsubdiv *) NULL;
    CCkdtree kt;
    CCkdtree *p_kt = (CCkdtree *) NULL;

    CCutil_dat_getnorm (dat, &norm);
    if ((norm != CC_EUCLIDEAN) && (norm != CC_MANNORM) && (norm != CC_GEOM))  {
        fprintf (stderr, "border_opt called with norm %d\b", norm);
        rval = 1;  goto CLEANUP;
    }

    printf ("Border optimization ...\n"); fflush (stdout);

    for (i = 0; i < pcount; i++) {
        if (intrac[i].cnt > mcnt) mcnt = intrac[i].cnt;
        for (j = i+1; j < pcount; j++) {
            subdiv_adj (&intrac[i], &intrac[j], &len);
            if (len == 0) cnt++;
        }
    }

    printf ("Number of adjacent borders: %d\n", cnt);
    printf ("Max Subdivsion: %d\n", mcnt);
    fflush (stdout);

    if (cnt == 0) {
        fprintf (stderr, "No adjacent borders?\n");
        rval = 1;  goto CLEANUP;
    }

    CC_MALLOC (trac, cnt, CCsubdiv);
    for (i = 0; i < cnt; i++) {
        trac[i].id = i;
        trac[i].xrange[0] = 0.0; trac[i].xrange[1] = 1.0;
        trac[i].yrange[0] = 0.0; trac[i].yrange[1] = 1.0;
        trac[i].cnt = 0;
        trac[i].bound = -1.0;
    }

    sprintf (buf, "%s_b", name);

    if (norm != CC_GEOM) {
        rval = CCkdtree_build (&kt, ncount, dat, (double *) NULL, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        p_kt = &kt;
    }

    tcount = 0;
    for (i = 0; i < pcount; i++) {
        for (j = i+1; j < pcount; j++) {
            subdiv_adj (&intrac[i], &intrac[j], &len);
            if (len == 0) {
                rval = improve_border (name, i, j, dat, intrac, invtour,
                          mcnt/2, &(trac[tcount]), &gotit, p_kt, tcount);
                CCcheck_rval (rval, "improve_border failed");
                if (gotit) tcount++;
            }
        }
    }

    rval = CCutil_write_subdivision_index (buf, ncount, tcount, trac);
    CCcheck_rval (rval, "CCutil_write_subdivision_index failed");

CLEANUP:
    if (p_kt) CCkdtree_free (p_kt);
    return rval;
}

#define MIN_BORDER_POINTS 16

static int improve_border (char *name, int ip, int iq, CCdatagroup *dat,
        CCsubdiv *trac, int *invtour, int pnt, CCsubdiv *sbox, int *gotit,
        CCkdtree *kt, int id)
{
    int rval = 0;
    int i, tempi, pcount, qcount, pbcnt, qbcnt, norm, flip = 0, border = 0;
    int *pnamelist = (int *) NULL, *qnamelist = (int *) NULL;
    int *pbord = (int *) NULL, *qbord = (int *) NULL,  *ulist = (int *) NULL;
    double xnn, ynn, xmax, xmin, ymax, ymin;
    char bname[1024];
    CCsubdiv *p = &(trac[ip]), *q = &(trac[iq]), *tempsub;

    CCutil_dat_getnorm (dat, &norm);

    *gotit = 0;

    printf ("Improve border %d-%d\n", ip, iq); fflush (stdout);

    sprintf (bname, "%s_b", name);

    find_common_border (p, q, &border, &flip);
    if (border == -1) {
        fprintf (stderr, "the adjacent pair does not share a border\n");
        rval = 1;  goto CLEANUP;
    }

    if (flip) {
        CC_SWAP (p, q, tempsub);
        CC_SWAP (ip, iq, tempi);
    }

    rval = get_part_data (name, ip, &pcount, &pnamelist);
    CCcheck_rval (rval, "get_part_data failed");

    rval = get_part_data (name, iq, &qcount, &qnamelist);
    CCcheck_rval (rval, "get_part_data failed");

    if (border == 0) {
        rval = grab_border_points (0, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->yrange[0], q->yrange[1], pnt);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (1, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->yrange[0], p->yrange[1], pnt);
        CCcheck_rval (rval, "grab_border_points failed");
    } else {
        rval = grab_border_points (2, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->xrange[0], q->xrange[1], pnt);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (3, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->xrange[0], p->xrange[1], pnt);
        CCcheck_rval (rval, "grab_border_points failed");
    }

    if (pbcnt < MIN_BORDER_POINTS || qbcnt < MIN_BORDER_POINTS) {
        *gotit = 0;
        goto CLEANUP;
    }

    CC_MALLOC (ulist, pbcnt + qbcnt, int);
    for (i = 0; i < pbcnt; i++) ulist[i] = pbord[i];
    for (i = 0; i < qbcnt; i++) ulist[pbcnt+i] = qbord[i];

    /* The data is in ulist, having pbcnt+qbcnt points  */

    sbox->cnt = pbcnt + qbcnt;

    xmin = ymin = CCtsp_LP_MAXDOUBLE;
    xmax = ymax = -CCtsp_LP_MAXDOUBLE;
   
    for (i = 0; i < pbcnt + qbcnt; i++) {
        xnn = dat->x[ulist[i]];
        ynn = dat->y[ulist[i]];
        if (xnn < xmin) xmin = xnn;
        if (xnn > xmax) xmax = xnn;
        if (ynn < ymin) ymin = ynn;
        if (ynn > ymax) ymax = ynn;
    }

    sbox->xrange[0] = xmin;
    sbox->xrange[1] = xmax;
    sbox->yrange[0] = ymin;
    sbox->yrange[1] = ymax;

    *gotit = 1;
    
    /* No need to create depots in a border problem */
    rval = build_subproblem (dat, kt, invtour, bname, -1, id, sbox->cnt, ulist,
                             (char *) NULL, norm, sbox, (rcinfo *) NULL);
    CCcheck_rval (rval, "build_subproblem failed");

CLEANUP:
    CC_IFFREE (pnamelist, int);
    CC_IFFREE (qnamelist, int);
    CC_IFFREE (ulist, int);
    CC_IFFREE (pbord, int);
    CC_IFFREE (qbord, int);
    return rval;
}

static void find_common_border (CCsubdiv *p, CCsubdiv *q, int *border,
        int *flip)
{
    if (p->xrange[0] == q->xrange[1]) {
        *border = 0;  *flip = 0;
    } else if (p->xrange[1] == q->xrange[0]) {
        *border = 0;  *flip = 1;
    } else if (p->yrange[0] == q->yrange[1]) {
        *border = 1;  *flip = 0;
    } else if (p->yrange[1] == q->yrange[0]) {
        *border = 1;  *flip = 1;
    } else {
        *border = -1; *flip = 0;
    }
}

static int grab_border_points (int btype, int scount, int *slist,
        CCdatagroup *dat, int *p_bcount, int **p_blist, double xlow,
        double xhi, int maxpoints)
{
    int rval = 0, bcount = 0, i, k, xcount = 0;
    int *blist = (int *) NULL, *xnames = (int *) NULL, *xperm = (int *) NULL;
    double *x = (double *) NULL;

    *p_bcount = 0;
    *p_blist = (int *) NULL;;

    CC_MALLOC (x, scount, double);
    CC_MALLOC (xnames, scount, int);

    if (btype == 0 || btype == 1) {
        for (i = 0; i < scount; i++) {
            k = slist[i];
            if (dat->y[k] >= xlow && dat->y[k] <= xhi) {
                xnames[xcount] = k;
                if (btype == 0) x[xcount++] = dat->x[k];
                else            x[xcount++] = -dat->x[k];
            }
        }
    } else {
        for (i = 0; i < scount; i++) {
            k = slist[i];
            if (dat->x[k] >= xlow && dat->x[k] <= xhi) {
                xnames[xcount] = k;
                if (btype == 2) x[xcount++] = dat->y[k];
                else            x[xcount++] = -dat->y[k];
            }
        }
    }

    if (xcount == 0) goto CLEANUP;   /* No points to grab */

    CC_MALLOC (xperm, xcount, int);
    for (i = 0; i < xcount; i++) xperm[i] = i;

    CCutil_double_perm_quicksort (xperm, x, xcount);
    if (xcount < maxpoints) bcount = xcount;
    else                    bcount = maxpoints;

    CC_MALLOC (blist, bcount, int);
    for (i = 0; i < bcount; i++) blist[i] = xnames[xperm[i]];

    *p_bcount = bcount;
    *p_blist = blist;

CLEANUP:
    CC_IFFREE (x, double);
    CC_IFFREE (xnames, int);
    CC_IFFREE (xperm, int);
    if (rval) { CC_IFFREE (blist, int); }
    return rval;
}

static int get_part_data (char *name, int ind, int *p_count, int **p_namelist)
{
    int rval = 0, ncount = 0;
    int *namelist = (int *) NULL, *tour = (int *) NULL;
    char local_mastername[1024], local_namesname[1024];;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);
    sprintf (local_mastername, "%s_%d.mas", name, ind);

    rval = CCutil_getmaster (local_mastername, &ncount, &dat, &tour);
    CCcheck_rval (rval, "CCutil_getmaster failed");

    if (dat.orig_names) {
        fprintf (stderr, "border optimization only for cuts (no depots)\n");
        rval = 1; goto CLEANUP;
    }

    if (p_namelist) {
        CC_MALLOC (namelist, ncount, int);
        sprintf (local_namesname,  "%s_%d.nam", name, ind);
        rval = get_names (local_namesname, namelist, ncount);
        CCcheck_rval (rval, "get_names failed")
        *p_namelist = namelist;
    }
    if (p_count) *p_count = ncount;

CLEANUP:
    if (rval) { CC_IFFREE (namelist, int); }
    CC_IFFREE (tour, int);
    CCutil_freedatagroup (&dat);
    return rval;
}

static int get_names (char *fname, int *names, int ncount)
{
    int rval = 0, i, k;
    FILE *in = (FILE *) NULL;

    in = fopen (fname, "r");
    if (!in) {
        fprintf (stderr, "unable to open %s for reading\n", fname);
        rval = 1; goto CLEANUP;
    }

    if (fscanf (in, "%d", &k) != 1) {
        fprintf (stderr, "file %s has wrong format\n", fname);
        rval = 1; goto CLEANUP;
    }

    if (k != ncount) {
        fprintf (stderr, "file %s does not match ncount\n", fname);
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        if (fscanf (in, "%d", &(names[i])) != 1) {
            fprintf (stderr, "file %s has wrong format\n", fname);
            rval = 1; goto CLEANUP;
        }
    }

CLEANUP:
    if (in) fclose (in);
    return rval;
}

static void subdiv_adj (CCsubdiv *s, CCsubdiv *t, int *len)
{
    *len = 1;

    if ((s->xrange[0] == t->xrange[1]) || (s->xrange[1] == t->xrange[0])) {
        if ((t->yrange[0] <= s->yrange[1] && t->yrange[0] >= s->yrange[0]) ||
            (t->yrange[1] <= s->yrange[1] && t->yrange[1] >= s->yrange[0]) ||
            (s->yrange[1] <= t->yrange[1] && s->yrange[1] >= t->yrange[0])) {
                *len = 0;
        }
    } else if ((s->yrange[0]==t->yrange[1]) || (s->yrange[1]==t->yrange[0])) {
        if ((t->xrange[0] <= s->xrange[1] && t->xrange[0] >= s->xrange[0]) ||
            (t->xrange[1] <= s->xrange[1] && t->xrange[1] >= s->xrange[0]) ||
            (s->xrange[1] <= t->xrange[1] && s->xrange[1] >= t->xrange[0])) {
                *len = 0;
        }
    }
}


/********************  End of border crossing code  *********************/

