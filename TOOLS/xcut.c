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
/*                TOOL TO CALL CUTTING ROUTINES FOR AN X-VECTOR             */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: May 15, 2013                                                      */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"
#include "localcut.h"

#define USE_DOMINO 

typedef struct node {
    int *adj;
    int degree;
} node;

typedef struct graph {
    int ncount;
    node *nodelist;
    int *adjspace;
} graph;


static char *xfname = (char *) NULL;
static char *masterfname = (char *) NULL;
static char *targetfname = (char *) NULL;
static int seed = 0;
static double alphamult = 0.5;

int
    main (int ac, char **av);

static int
    check_x (int ncount, int xcount, int *xlist, double *x, int *ok),
    build_graph (graph *G, int ncount, int ecount, int *elist),
    parseargs (int ac, char **av);

static void
    init_graph (graph *G),
    free_graph (graph *G),
    usage (char *fname);


int main (int ac, char **av)
{
    int rval = 0, i, ncount, xcount, tcount, newxcount, cutcount = 0;
    int *ptour = (int *) NULL, *xlist = (int *) NULL, *newxlist = (int *) NULL;
    int *tlist = (int *) NULL;
    int k, ok, n0, n1, temp;
    double *x = (double *) NULL, *tx = (double *) NULL, szeit, z;
    double *newx = (double *) NULL, maxviol;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp lp;
    CCtsp_lpgraph lg;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL, *c, *cnext;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dompool = (CCtsp_lpcuts *) NULL;
    char *probname = (char *) NULL, buf[256]; 
    
    CCtsp_init_lpgraph_struct (&lg);
    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname) {
        fprintf (stderr, "Must specify a master file\n");
        usage (av[0]); goto CLEANUP;
    }

    CCtsp_init_tsp_lp_struct (&lp);

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    probname = CCtsp_problabel (xfname);
    printf ("probname = %s\n", probname); fflush (stdout);

    rval = CCtsp_read_x (ncount, ptour, xfname, &xcount, &xlist, &x);
    CCcheck_rval (rval, "CCtsp_read_x failed");

    rval = CCtsp_build_lpgraph (&lg, ncount, xcount, xlist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
    rval = CCtsp_build_lpadj (&lg, 0, xcount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");

    if (targetfname) {
        rval = CCtsp_read_x (ncount, ptour, targetfname, &tcount, &tlist, &tx);
        CCcheck_rval (rval, "CCtsp_read_x failed");

        /* Build convex combination */ 

        printf ("alpha = %f\n", alphamult); fflush (stdout);

        /* Alloc newx, newxlist, start by putting x into it */

        newx = CC_SAFE_MALLOC (xcount+tcount, double);
        CCcheck_NULL (newx, "out of memory for newx");
        newxlist = CC_SAFE_MALLOC (2*(xcount+tcount), int);
        CCcheck_NULL (newxlist, "out of memory for newxlist");

        for (i = 0; i < xcount; i++) {
            newx[i] = (1.0 - alphamult)*x[i];
            newxlist[2*i] = xlist[2*i];
            newxlist[2*i+1] = xlist[2*i+1];
        }
        newxcount = xcount;

        for (i = 0; i < tcount; i++) {
            n0 = tlist[2*i];  n1 = tlist[2*i+1];
            if (n0 > n1) CC_SWAP (n0, n1, temp);

            k = CCtsp_find_edge (&lg, n0, n1);
            if (k == -1) {
                newx[newxcount] = alphamult*tx[i];
                newxlist[2*newxcount] = n0;
                newxlist[2*newxcount+1] = n1;
                newxcount++;
            } else {
                newx[k] = ((1.0 - alphamult)*x[k]) + (alphamult*tx[i]);
            }
        } 

        CC_IFFREE (x, double);
        CC_IFFREE (xlist, int);
        x = newx;
        xlist = newxlist;
        xcount = newxcount;
    }

    rval = check_x (ncount, xcount, xlist, x, &ok);
    CCcheck_rval (rval, "check_x failed");
    if (!ok) {
        printf ("x vector failed simple tests\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &dompool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    rval = CCtsp_block_combs (&cuts, &cutcount, ncount, xcount, xlist, x, 1,
                              &z, &maxviol);
    CCcheck_rval (rval, "CCtsp_block_combs failed");
    printf ("Found %d block_combs\n", cutcount);
    fflush (stdout);

    rval = CCtsp_edge_comb_grower (&cuts, &cutcount, ncount, xcount, xlist, x,
                                   &lp.stats.extra_tighten_stats, &z, &maxviol);
    CCcheck_rval (rval, "CCtsp_edge_comb_grower failed");
    printf ("Found %d grow_combs\n", cutcount);
    fflush (stdout);

    {
        int  maxchunksize = 28;
        CCchunk_flag flags;
        CCchunk_localcut_timer lc_timer;

        CCchunk_init_localcut_timer (&lc_timer);

        flags.dummy = flags.permute = flags.weighted = 0;
        flags.spheres = flags.uncivilized = flags.noshrink = flags.nolift = 0;
        flags.maxchunksize = maxchunksize;
        flags.spheresize   = maxchunksize - 2;

        rval = CCchunk_localcuts (&cuts, &cutcount, ncount, xcount, xlist, x,
          0.0, flags, &lc_timer, 1, &rstate, &z, &maxviol);
        CCcheck_rval (rval, "CCchunk_localcuts failed");
        printf ("Found %d local cuts, size %d\n", cutcount, maxchunksize);
        fflush (stdout);

        flags.spheres = 1;
        rval = CCchunk_localcuts (&cuts, &cutcount, ncount, xcount,
                     xlist, x, 0.0, flags, &lc_timer, 1, &rstate, &z, &maxviol);
        CCcheck_rval (rval, "CCchunk_localcuts failed");
        printf ("Found %d local cuts, size %d, spheres\n", cutcount,
                     maxchunksize);
        fflush (stdout);
    }

    for (c = cuts; c; c = cnext) {
        cnext = c->next; 
        rval = CCtsp_add_to_cutpool_lpcut_in (pool, c);
        CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
        CCtsp_free_lpcut_in (c);
    }
    cuts = (CCtsp_lpcut_in *) NULL;

    rval = CCtsp_double_decker_lp (pool, &lp.stats.extra_tighten_stats, 
                  &cuts, &cutcount, ncount, xcount, xlist, x, 2.0, 10000,
                  &maxviol, &rstate, &z);
    CCcheck_rval (rval, "CCtsp_double_decker_lp failed");
    printf ("Found %d double deckers from combs\n", cutcount);
    fflush (stdout);

    rval = CCtsp_star_lp (pool, &lp.stats.extra_tighten_stats, &cuts,
                  &cutcount, ncount, xcount, xlist, x, 2.0, 10000, &maxviol,
                  &rstate, &z);
    CCcheck_rval (rval, "CCtsp_star_lp failed");
    printf ("Found %d star inequalities\n", cutcount);
    fflush (stdout);

    rval = CCtsp_teething_lp (pool, &lp.stats.extra_tighten_stats, &cuts,
                  &cutcount, ncount, xcount, xlist, x, 10.0, 10000, &maxviol,
                  &rstate, &z);
    CCcheck_rval (rval, "CCtsp_teething_lp failed");
    printf ("Found %d teething inequalities\n", cutcount);
    fflush (stdout);

    rval = CCtsp_handling_lp (pool, &lp.stats.extra_tighten_stats, &cuts,
                  &cutcount, ncount, xcount, xlist, x, 10.0, 10000, &maxviol,
                  &rstate, &z);
    CCcheck_rval (rval, "CCtsp_handling_lp failed");
    printf ("Found %d handling inequalities\n", cutcount);
    fflush (stdout);

    for (c = cuts; c; c = cnext) {
        cnext = c->next; 
        rval = CCtsp_add_to_cutpool_lpcut_in (pool, c);
        CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
        CCtsp_free_lpcut_in (c);
    }
    cuts = (CCtsp_lpcut_in *) NULL;

    sprintf (buf, "%s.pul", probname);
    rval = CCtsp_write_cutpool (ncount, buf, pool);
    CCcheck_rval (rval, "CCtsp_write_cutpool failed");


#ifdef USE_DOMINO
/*
    rval = CCtsp_shrink_domino (&cuts, &cutcount, ncount, xcount,
                        xlist, x, 2, 5, 1, &rstate, (char *) NULL);
    CCcheck_rval (rval, "CCtsp_shrink_domino failed");
*/
    rval = CCtsp_DP_cuts (&cuts, &cutcount, ncount, xcount,
                         xlist, x, 1, (char *) NULL, &z, &maxviol);
    CCcheck_rval (rval, "CCtsp_DP_cuts failed");

    for (c = cuts; c; c = cnext) {
        cnext = c->next; 
        rval = CCtsp_add_to_cutpool_lpcut_in (dompool, c);
        CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");
        CCtsp_free_lpcut_in (c);
    }
    cuts = (CCtsp_lpcut_in *) NULL;

    sprintf (buf, "%s.dompul", probname);
    rval = CCtsp_write_cutpool (ncount, buf, dompool);
    CCcheck_rval (rval, "CCtsp_write_cutpool failed");
#endif

CLEANUP:
    CC_IFFREE (ptour, int);
    CCutil_freedatagroup (&dat);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    CC_IFFREE (tlist, int);
    CCtsp_free_lpgraph (&lg);
    if (pool) CCtsp_free_cutpool (&pool);
    if (dompool) CCtsp_free_cutpool (&dompool);

    return rval;
}

static int check_x (int ncount, int xcount, int *xlist, double *x, int *ok)
{
    int rval = 0, i, j;
    double y;
    graph G;

    init_graph (&G);
   
    *ok = 0;

    for (i = 0; i < xcount; i++) {
        if (x[i] < -0.0001 || x[i] > 1.0001) {
            printf ("x[%d] outside bounds: %f\n", i, x[i]); fflush (stdout);
            goto CLEANUP;
        }
    }

    rval = build_graph (&G, ncount, xcount, xlist);
    CCcheck_rval (rval, "build_graph failed");

    for (i = 0; i < ncount; i++) {
        y = 0.0;
        for (j = 0; j < G.nodelist[i].degree; j++) {
            y += x[G.nodelist[i].adj[j]];
        }
        if (y < 1.999 || y > 2.001) {
            printf ("deg equation %d violated: %f\n", i, y); fflush (stdout);
            goto CLEANUP;
        }
    }

    *ok = 1;

CLEANUP:
    free_graph (&G);
    return rval;
}

static void init_graph (graph *G)
{
    if (G) {
        G->ncount = 0;
        G->nodelist = (node *) NULL;
        G->adjspace = (int *) NULL;
    }
}

static void free_graph (graph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->adjspace, int);
        init_graph (G);
    }
}

static int build_graph (graph *G, int ncount, int ecount, int *elist)
{
    int rval = 0, i, *p;
    node *n0, *n1;
   
    init_graph (G);
    G->ncount = ncount;
    G->nodelist = CC_SAFE_MALLOC (ncount, node);
    CCcheck_NULL (G->nodelist, "out of memory for nodelist");
    G->adjspace = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (G->adjspace, "out of memory for adjspace");

    for (i = 0; i < ncount; i++) G->nodelist[i].degree = 0;

    for (i = 0; i < ecount; i++) {
        G->nodelist[elist[2*i]].degree++;
        G->nodelist[elist[2*i+1]].degree++;
    }

    p = G->adjspace;
    for (i = 0; i < ncount; i++) {
        G->nodelist[i].adj = p;
        p += G->nodelist[i].degree;
        G->nodelist[i].degree = 0;
    }

    for (i = 0; i < ecount; i++) {
        n0 = &G->nodelist[elist[2*i]];
        n1 = &G->nodelist[elist[2*i+1]];
        n0->adj[n0->degree++] = i;
	n1->adj[n1->degree++] = i;
    }

CLEANUP:
    return rval;
}


static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "a:M:s:V:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'a':
            alphamult = atof(boptarg);
            break;
        case 'M':
            masterfname  = boptarg;
            break;
        case 's':
            seed = atoi(boptarg);
            break;
        case 'V':
            targetfname  = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        xfname = av[boptind++];
    } else{
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] x_file\n", fname);
    fprintf (stderr, "   -a #  alpha multiplier, range (0,1), default 0.5\n");
    fprintf (stderr, "   -M f  specify a master file (required)\n");
    fprintf (stderr, "   -V f  specify a target_x file (usually a tour)\n");
    fprintf (stderr, "   -s #  seed\n");
}

