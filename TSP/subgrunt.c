/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2017 by David Applegate, Robert Bixby,              */
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
/*                  WORKER FOR PARALLEL CUTTING ROUTINES                    */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Updated:     March 4, 2017                                              */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static char *bossname = (char *) NULL;
static int maxchunksize         = 16;
static int multiple_chunker     = 0;
static int usedominos           = 0;

static int process_subproblem (char *probname, int id, int ncount,
        CCdatagroup *dat, int *ptour, double *lbound, CCrandstate *rstate,
        int tmax, int use_mult);
static int build_edges (CCdatagroup *dat, int ncount, int *ecount, int **elist,
        int **elen, int silent, CCrandstate *rstate);
static int parseargs (int ac, char **av);
static void usage (char *fname);
int main (int ac, char **av);

int main (int ac, char **av)
{
    char *bosshost = (char *) NULL;
    double lbound = -1.0, rtime = 0.0, szeit;
    char probname[CCutil_FILE_NAME_LEN];
    int id = -1, rval = 0, seed = 99, ncount;
    int *perm = (int *) NULL;
    CC_SFILE *s = (CC_SFILE *) NULL;
    CCrandstate rstate;
    CCdatagroup dat;

    rval = parseargs (ac, av);
    if (rval) return 1;

    CCutil_printlabel ();
    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed");

    CCutil_sprand (seed, &rstate);

    bosshost = bossname;

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
        rval = CCutil_swrite_double (s, lbound);
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

        CCutil_init_datagroup (&dat);
        perm = (int *) NULL;

        rval = CCutil_readmaster (s, &ncount, &dat, &perm);
        CCcheck_rval (rval, "CCutil_readmaster failed");
        
        CCutil_sclose (s);
        s = (CC_SFILE *) NULL;

        printf ("PROCESSING %s subproblem %d\n", probname, id);
        fflush (stdout);

        szeit = CCutil_zeit ();

        rval = process_subproblem (probname, id, ncount, &dat, perm, &lbound,
                                   &rstate, maxchunksize, multiple_chunker);
        CCcheck_rval (rval, "process_subproblem failed");

        rtime = CCutil_zeit () - szeit;

        CC_IFFREE (perm, int);
        CCutil_freedatagroup (&dat);
    }

DONE:
    printf ("No work available.  Shutting down.\n"); fflush (stdout);

CLEANUP:
    if (s != (CC_SFILE *) NULL) CCutil_sclose (s);
    return rval;
}

static int process_subproblem (char *probname, int id, int ncount,
        CCdatagroup *dat, int *ptour, double *lbound, CCrandstate *rstate,
        int tmax, int use_mult) 
{
    char name[CCutil_FILE_NAME_LEN];
    int i, ecount, yesno, rval = 0, silent = 1, infeasible = 0;
    int *elist = (int *) NULL, *elen = (int *) NULL, *besttour = (int *) NULL;
    double ptour_len, lb;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dompool = (CCtsp_lpcuts *) NULL;
    CCtsp_cutselect sel;
    CCbigguy exbound;

    printf ("ncount = %d, depots = %d\n", ncount, dat->ndepot);
    fflush (stdout);

    CCutil_cycle_len (ncount, dat, ptour, &ptour_len);
    printf ("initial tour: %.2f\n", ptour_len); fflush (stdout);

    rval = build_edges (dat, ncount, &ecount, &elist, &elen, silent, rstate);
    CCcheck_rval (rval, "build_edges failed");

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &dompool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    CC_MALLOC (besttour, ncount, int);
    for (i = 0; i < ncount; i++) besttour[i] = i;

    sprintf (name, "%s_%d", probname, id);

    rval = CCtsp_init_lp (&lp, name,  -1, (char *) NULL, ncount, dat,
                    ecount, elist, elen, 0, (int *) NULL, (int *) NULL,
                    0, ptour, ptour_len, pool, dompool, 0,
                    rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "Initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }
   
    CCtsp_init_cutselect (&sel);
    rval = CCtsp_get_lp_result (lp, &lb, (double *) NULL, (int *) NULL,
              (int **) NULL, (double **) NULL, (double **) NULL,
              (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

    /* The value 2 selects tolerance based on lp value rather than gap */
    rval = CCtsp_cutselect_set_tols (&sel, lp, 2, 0);
    CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
    CCtsp_cutselect_chunksize (&sel, maxchunksize);
    CCtsp_cutselect_dominos (&sel, usedominos, 1);

    if (use_mult) {
        rval = CCtsp_cutting_multiple_loop (lp, &sel, (int *) NULL, 1, tmax, 0,
                                            silent, rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_cutting_multiple_loop failed");
    } else {
        rval = CCtsp_cutting_loop (lp, &sel, (int *) NULL, 1, silent, rstate,
                                   &infeasible);
        CCcheck_rval (rval, "CCtsp_cutting_loop failed");
    }
    if (infeasible) {
        fprintf (stderr, "Infeasible LP after cutting\n");
        rval = 1; goto CLEANUP;
    }

    {
        double tourval;
        CCutil_start_timer (&lp->stats.linkern);
        rval = CCtsp_call_x_heuristic (lp, &tourval, besttour, silent,
                                       rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_call_x_heuristic failed\n");
            goto CLEANUP;
        }
        CCutil_stop_timer (&lp->stats.linkern, 1);
        if (tourval < lp->upperbound) {
            printf ("New upperbound from x-heuristic: %.2f\n", tourval);
            lp->upperbound = tourval;
        }
    }

    printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
            CClp_nrows (lp->lp), CClp_ncols (lp->lp),
            CClp_nnonzeros (lp->lp));
    printf ("Final lower bound %f, upper bound %f\n", lp->lowerbound,
                                                      lp->upperbound);
    fflush (stdout);

{
    int zz, qq = 0;
    for (zz = 0; zz < lp->cuts.cutcount; zz++) {
        if (lp->cuts.cuts[zz].dominocount > 0) qq++;
    }
    printf ("Final LP has %d domino cuts\n", qq); fflush (stdout);
}

    rval = CCtsp_exact_price (lp, &exbound, 1, 0, 0, 0, silent);
    CCcheck_rval (rval, "CCtsp_exact_price failed");

    lp->exact_lowerbound = exbound;
    printf ("Exact lower bound: %.6f\n", CCbigguy_bigguytod (exbound));
    printf ("DIFF: %f\n", lp->lowerbound - CCbigguy_bigguytod (exbound));
    fflush (stdout);

    if (dat->ndepot > 0) {
        rval = CCtsp_depot_valid (lp, dat->ndepot, &yesno);
        CCcheck_rval (rval, "CCtsp_depot_valid failed");
    } else {
        yesno = 1;
    }

    if (yesno == 0) *lbound = -1.0;
    else            *lbound = CCbigguy_bigguytod (exbound);

CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    if (pool) CCtsp_free_cutpool (&pool);
    if (dompool) CCtsp_free_cutpool (&dompool);
    CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int build_edges (CCdatagroup *dat, int ncount, int *ecount, int **elist,
         int **elen, int silent, CCrandstate *rstate)
{
    int i, rval = 0;
    CCedgegengroup plan;

    CCedgegen_init_edgegengroup (&plan);
    plan.linkern.count = 10;
    plan.linkern.quadnearest = 2;
    plan.linkern.greedy_start = 0;
    plan.linkern.nkicks = (ncount / 100) + 1;

    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, ecount,
                            elist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");
    
    CC_MALLOC (*elen, *ecount, int);
    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i], (*elist)[(2*i)+1], dat);
    }

CLEANUP:
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "C:mZ:?", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'C': maxchunksize = atoi(boptarg); break;
        case 'm': multiple_chunker = 1; break;
        case 'Z': usedominos = atoi(boptarg); break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        bossname = av[boptind++];
    } else {
        fprintf (stderr, "Missing tspfile\n");
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] boss\n", fname);
    fprintf (stderr, "   -C #  maximum chunk size in localcuts (default 16)\n");
    fprintf (stderr, "   -m    use multiple passes of cutting loop\n");
#ifdef CCtsp_USE_DOMINO_CUTS
    fprintf (stderr, "   -Z #  dp-cuts (#=1 normal, #=2 shrunk, #=3 both)\n");
#endif
}
