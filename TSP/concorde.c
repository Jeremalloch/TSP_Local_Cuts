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
/*                  THE MAIN PROGRAM FOR CONCORDE                           */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 25, 1995                                                */
/*        November 28, 2003 (bico)                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/*  NOTE:  When using CC_SPARSE edge sets, it is important to specify a     */
/*   a tour or an upperbound if you have one, since our heuristics are      */
/*   not designed for finding Hamilton circuits in sparse graphs.           */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "linkern.h"
#include "heldkarp.h"
#include "bigguy.h"
#include "macrorus.h"

#define CC_JUST_SUBTOUR (1)
#define CC_JUST_BLOSSOM (2)
#define CC_JUST_SUBTOUR_AND_BLOSSOM (3)
#define CC_JUST_FAST_CUTS (4)
#define CC_JUST_WORLD_CUTS (5)
#define CC_JUST_DP (6)
#define CC_JUST_GENERAL (7)

static int norm = CC_EUCLIDEAN;
static char *datfname        = (char *) NULL;
static char *edgegenfname    = (char *) NULL;
static char *problname       = (char *) NULL;
static char *probfname       = (char *) NULL;
static char *edgefname       = (char *) NULL;
static char *fullfname       = (char *) NULL;
static char *tourfname       = (char *) NULL;
static char *masterfname     = (char *) NULL;
static char *poolfname       = (char *) NULL;
static char *dominopoolfname = (char *) NULL;
static char *restartfname    = (char *) NULL;
static char *xfname          = (char *) NULL;
static char *outfname        = (char *) NULL;
static char *filecutname     = (char *) NULL;
static char *moatfname       = (char *) NULL;
static int seed                 = 0;
static int nnodes_want          = 0;
static int binary_in            = 0;
static int tsplib_in            = 1;
static int gridsize             = 0   /* 10000000 */;
static int just_cuts            = 0;
static int dontcutroot          = 0;
static int usetighten           = 0;
static int usedominos           = 0;
static int domino_safeshrink    = 1;
static int usextighten          = 0;
static int usekeepcutting       = 0;
static int maxchunksize         = 16;
static int multiple_chunker     = 0;
static int valid_edges          = 0;
static int dfs_branching        = 0;
static int bfs_branching        = 1;
static int simple_branching     = 0;
static int usebranchcliques     = 1;  
static int tentative_branch_num = 0;
static int complete_price       = 0;
static int want_rcnearest       = 0;
static int output_tour_as_edges = 0;
static int run_silently         = 1;
static int be_nethost           = 0;
static int unlink_files         = 0;
static int dump_root            = 0;
static int dump_the_cliques     = 0;
static double initial_ub = CCtsp_LP_MAXDOUBLE;
static unsigned short hostport = CCtsp_HOST_PORT;
static char *grunthostname  = (char *) NULL;
static char *cutbossname    = (char *) NULL;
static char *domcutbossname = (char *) NULL;
static char *dombossname    = (char *) NULL;
static int eliminate_edges = -1;   /* Set to 1 to force elim, 0 to not elim */
static int longedge_branching = 1; /* Set to 0 to turn off           */
static int save_proof = 0;         /* Set to 1 to save the proof     */
static int standalone_branch = 0;  /* Set to 1 to do a manual branch */

static int dump_cliques (CCtsp_lp *lp);

static void
    adjust_upbound (double *bound, int ncount, CCdatagroup *dat),
    usage (char *f);

static int
    handle_just_cuts (CCtsp_lp *lp, int the_cuts, int tmax, CCrandstate *rstate,
       int silent),
    run_hk (int ncount, CCdatagroup *dat, int *hk_tour),
    build_edges (int *p_ecount, int **p_elist, int **p_elen,
        int ncount, int *ptour, CCdatagroup *dat, char *in_edgefname,
        char *in_edgegenfname, int in_just_cuts, int silent,
        CCrandstate *rstate),
    build_fulledges (int *p_excount, int **p_exlist, int **p_exlen,
        int ncount, int *ptour, char *in_fullfname),
    parseargs (int ac, char **av),
    find_tour (int ncount, CCdatagroup *dat, int *perm, double *ub,
        int trials, int silent, CCrandstate *rstate),
    getedges (CCdatagroup *dat, CCedgegengroup *plan, int ncount, int *ecount,
        int **elist, int **elen, int silent, CCrandstate *rstate),
    dump_rc (CCtsp_lp *lp, int count, char *pname, int usesparse),
    find_moats (int ncount, CCdatagroup *dat, int *perm, CCrandstate *rstate,
        char *fname);



int main (int ac, char **av)
{
    int rval = 0, i, ncount, silent, allow_dups, use_gridsize;
    char *probname = (char *) NULL, buf[1024];
    CCdatagroup dat;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_cutselect sel, tentativesel;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *dominopool = (CCtsp_lpcuts *) NULL;
    int *ptour = (int *) NULL, *snowtour = (int *) NULL;
    int ecount = 0, *elist = (int *) NULL, *elen = (int *) NULL;
    int excount = 0, *exlist = (int *) NULL, *exlen = (int *) NULL;
    int *besttour = (int *) NULL;
    int is_infeasible = 0, bbcount = 0, infeasible = 0;
    double szeit, branchzeit = 0.0, upbound = 0.0;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed");
    
    szeit = CCutil_zeit ();
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    CCutil_printlabel ();
    CCutil_signal_init ();
    CCutil_sprand (seed, &rstate);
    printf ("Using random seed %d\n", seed); fflush (stdout);

    silent = run_silently;

    if (grunthostname) {
#ifdef CC_NETREADY
        rval = CCtsp_grunt (grunthostname, hostport, poolfname,
              dominopoolfname, cutbossname, domcutbossname, problname,
              silent, &rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_grunt failed\n");
        }
        goto CLEANUP;
#else  /* CC_NETREADY */
        fprintf (stderr, "Networking not enabled\n");
        rval = 1; goto CLEANUP;
#endif /* CC_NETREADY */
    }

    if (!be_nethost) { hostport = 0; }

    CCtsp_init_cutselect (&sel);
    CCtsp_init_tentative_cutselect (&tentativesel);
    CCtsp_cutselect_tighten (&sel, usetighten);
    CCtsp_cutselect_tighten (&tentativesel, usetighten);
    CCtsp_cutselect_xtighten (&sel, usextighten);
    CCtsp_cutselect_chunksize (&sel, maxchunksize);
    CCtsp_cutselect_dominos (&sel, usedominos, domino_safeshrink);
    if (usekeepcutting) CCtsp_cutselect_keep_cutting (&sel);
    if (filecutname) CCtsp_cutselect_filecuts (&sel, filecutname);
#ifdef CC_NETREADY
    if (cutbossname != (char *) NULL) {
        CCtsp_cutselect_remotepool (&sel, cutbossname);
        CCtsp_cutselect_remotepool (&tentativesel, cutbossname);
    }
    if (domcutbossname != (char *) NULL) {
        CCtsp_cutselect_remotedompool (&sel, domcutbossname);
        CCtsp_cutselect_remotedompool (&tentativesel, domcutbossname);
    }
    if (dombossname != (char *) NULL) {
        CCtsp_cutselect_domboss (&sel, dombossname);
        CCtsp_cutselect_domboss (&tentativesel, dombossname);
    }
#endif /* CC_NETREADY */

    if (problname)        probname = CCtsp_problabel (problname);
    else if (datfname)    probname = CCtsp_problabel (datfname);
    else if (masterfname) probname = CCtsp_problabel (masterfname);
    else                  probname = CCtsp_problabel ("unnamed");
    if (probname == (char *) NULL) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        rval = 1; goto CLEANUP;
    }
    if (problname == (char *) NULL) { problname = probname; }

    if (masterfname) {
        rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
        CCcheck_rval (rval, "CCutil_getmaster failed");
        if (ncount < 10) {
            fprintf (stderr, "Master file has less than 10 nodes - Abort.\n");
            rval = 1; goto CLEANUP;
        }

        if (tourfname) {
            int *invperm = (int *) NULL;

            CC_MALLOC (snowtour, ncount, int);
            CC_MALLOC (invperm, ncount, int);
            for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;

            rval = CCutil_getcycle (ncount, tourfname, snowtour, 0);
            CCcheck_rval (rval, "CCutil_getcycle failed");
            for (i = 0; i < ncount; i++) {
                snowtour[i] = invperm[snowtour[i]];
            }
            CC_IFFREE (invperm, int);
        }
    } else {
        CCutil_init_datagroup (&dat);
        if (tsplib_in && datfname != (char *) NULL) {
            rval = CCutil_gettsplib (datfname, &ncount, &dat);
            CCcheck_rval (rval, "CCutil_gettsplib failed");
        } else {
            ncount = nnodes_want;
            if (gridsize < 0) {
                use_gridsize = -gridsize; allow_dups = 0;
            } else if (gridsize > 0) {
                use_gridsize = gridsize; allow_dups = 1;
            } else {
                use_gridsize = 10 * nnodes_want; allow_dups = 0;
            }
            rval = CCutil_getdata (datfname, binary_in, norm, &ncount, &dat,
                                   use_gridsize, allow_dups, &rstate);
            CCcheck_rval (rval, "CCutil_getdata failed");
        }

        /* Handle small instances */

        if (ncount < 3) {
            printf ("Only %d nodes -- must have at least 3 nodes in a TSP\n",
                     ncount);
            fflush (stdout);
            rval = 1;  goto CLEANUP;
        } else if (ncount < 10) {
            CC_MALLOC (besttour, ncount, int);
            if (ncount == 3) {
                for (i = 0; i < ncount; i++) besttour[i] = i;
            } else {
                rval = run_hk (ncount, &dat, besttour);
                CCcheck_rval (rval, "run_hk failed");
            }
            CC_MALLOC (ptour, ncount, int);
            for (i = 0; i < ncount; i++) ptour[i] = i;
            rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                                   outfname, output_tour_as_edges, silent);
            CCcheck_rval (rval, "CCtsp_dumptour failed");

            printf ("Total Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
            fflush (stdout);
            goto CLEANUP;
        }

        /***** Get the permutation tour and permute the data  *****/

        CC_MALLOC (ptour, ncount, int);
        if (tourfname) {
            rval = CCutil_getcycle (ncount, tourfname, ptour, 0);
            CCcheck_rval (rval, "CCutil_getcycle failed");
        } else {
            double bnd;
            if (just_cuts > 0) {
                rval = find_tour (ncount, &dat, ptour, &bnd, -1, silent,
                                  &rstate);
            } else if (initial_ub == CCtsp_LP_MAXDOUBLE) {
                rval = find_tour (ncount, &dat, ptour, &bnd, 1, silent,
                                  &rstate);
            } else {
                if (!silent) {
                    printf ("Initial bnd %f - use short LK\n", initial_ub);
                    fflush (stdout);
                }
                rval = find_tour (ncount, &dat, ptour, &bnd, 0, silent,
                                 &rstate);
            }
            CCcheck_rval (rval, "find_tour failed");
        }

        rval = CCutil_datagroup_perm (ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_datagroup_perm failed");

        sprintf (buf, "%s.mas", probname);
        rval = CCutil_putmaster (buf, ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_putmaster failed");
    }

    adjust_upbound (&initial_ub, ncount, &dat);

    if (moatfname) {
        rval = find_moats (ncount, &dat, ptour, &rstate, moatfname);
        CCcheck_rval (rval, "find_moats failed");
        goto CLEANUP;
    }

    if (!probfname && !restartfname) {
        rval = build_edges (&ecount, &elist, &elen, ncount, ptour,
                            &dat, edgefname, edgegenfname, just_cuts,
                            silent, &rstate);
        CCcheck_rval (rval, "build_edges failed");
    }

    rval = build_fulledges (&excount, &exlist, &exlen, ncount, ptour,
                            fullfname);
    CCcheck_rval (rval, "build_fulledges failed");
    
    rval = CCtsp_init_cutpool (&ncount, poolfname, &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
#ifdef CCtsp_USE_DOMINO_CUTS
    rval = CCtsp_init_cutpool (&ncount, dominopoolfname, &dominopool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed for dominos");
#endif

    /***** Initialize besttour to be the permutation tour  ****/

    CC_MALLOC (besttour, ncount, int);
    for (i = 0; i < ncount; i++) { besttour[i] = i; }

    if (restartfname) {
        upbound  = initial_ub;
        bbcount = 0;

        rval = CCtsp_bfs_restart (problname, restartfname, &sel,
                &tentativesel, &upbound, &bbcount, usebranchcliques, &dat,
                ptour, pool, dominopool, ncount, besttour, hostport,
                &branchzeit, save_proof, tentative_branch_num,
                longedge_branching, (double *) NULL, (int *) NULL, silent,
                &rstate);
        CCcheck_rval (rval, "CCtsp_bfs_restart failed");
        goto DONE;
    }

    rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                           (char *) NULL, 0, silent);
    CCcheck_rval (rval, "CCtsp_dumptour failed");

    rval = CCtsp_init_lp (&lp, problname, -1, probfname, ncount, &dat,
                    ecount, elist, elen, excount, exlist, exlen, valid_edges,
                    ptour, initial_ub, pool, dominopool, silent, &rstate,
                    &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        printf ("CCtsp_init_lp reports an infeasible LP\n");
        rval = CCtsp_verify_infeasible_lp (lp, &is_infeasible, silent);
        CCcheck_rval (rval, "CCtsp_verify_infeasible_lp failed");
        if (!is_infeasible) {
            printf ("Couldn't verify infeasible LP\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }
        upbound = CCtsp_LP_MAXDOUBLE;
        bbcount = 1;
        goto DONE;
    } else if (rval) {
        fprintf (stderr, "CCtsp_init_lp failed\n"); goto CLEANUP;
    }

#if 0  /* HACK TO RESET THE LP UPPERBOUND TO LARGER VALUE */
    if (initial_ub != CCtsp_LP_MAXDOUBLE) {
        printf ("Setting LP bound to initial_ub %f\n", initial_ub);
        fflush (stdout);
        lp->upperbound = initial_ub;
    }
#endif

    CCutil_start_timer (&lp->stats.total);
    
    ecount = 0;
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    excount = 0;
    CC_IFFREE (exlist, int);
    CC_IFFREE (exlen, int);

#if 0  /* HACK TO DUMP LP EDGES */
{
    int qq;
    printf ("%d %d\n", lp->graph.ncount, lp->graph.ecount);

    for (qq = 0; qq < lp->graph.ecount; qq++) {
        printf ("%d %d\n", lp->perm[lp->graph.edges[qq].ends[0]],
                           lp->perm[lp->graph.edges[qq].ends[1]]);
    }
    exit (1);
}
#endif

#ifndef CCtsp_WORLD_TSP
    if (lp->full_edges_valid) {   /* For WORLD runs don't test this */
        if (CCtsp_inspect_full_edges (lp)) {
            fprintf (stderr, "full edge set does not contain all LP edges\n");
            rval = 1; goto CLEANUP;
        }
    }
#endif

    if (standalone_branch) {
        rval = CCtsp_do_interactive_branch (lp, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_do_interactive_branch failed");
        printf ("Total Running Time: %.2f (seconds)\n", CCutil_zeit () - szeit);
        goto CLEANUP;
    }

#if 0  /* HACK TO DUMP TSP AS LP (WITH CPLEX) BEFORE CUTTING */
    printf ("Dump the LP problem\n"); fflush (stdout);
    cplex_write_prob (lp->lp);  /* To write the LP */
    printf ("DONE\n"); exit (1);
#endif

    if (just_cuts > 0) {
        rval = handle_just_cuts (lp, just_cuts, maxchunksize, &rstate, silent);
        CCcheck_rval (rval, "handle_just_cuts failed");
        if (want_rcnearest) {
            rval = dump_rc (lp, want_rcnearest, probname, 0);
            CCcheck_rval (rval, "dump_rc failed");
        }
        if (xfname) {
            rval = CCtsp_dump_x (lp, xfname);
            CCcheck_rval (rval, "CCtsp_dump_x failed");
        }
        goto DONE;
    }

    rval = CCtsp_cutselect_set_tols (&sel, lp, 1, silent);
    CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");

    if (dontcutroot == 0) {
        if (multiple_chunker) {
            rval = CCtsp_cutting_multiple_loop (lp, &sel, snowtour, 1,
                           maxchunksize, 1, silent, &rstate, &infeasible);
            CCcheck_rval (rval, "CCtsp_cutting_multiple_loop failed");
        } else {
            rval = CCtsp_cutting_loop (lp, &sel, snowtour, 1, silent, &rstate,
                                       &infeasible);
            CCcheck_rval (rval, "CCtsp_cutting_loop failed");
        }
        if (infeasible) {
            printf ("CCtsp_cutting_loop reports an infeasible LP\n");
            rval = CCtsp_verify_infeasible_lp (lp, &is_infeasible, silent);
            CCcheck_rval (rval, "CCtsp_verify_infeasible_lp failed");
            if (!is_infeasible) {
                printf ("Couldn't verify infeasibile LP\n");
                fflush (stdout);
                rval = 1; goto CLEANUP;
            }
            upbound = CCtsp_LP_MAXDOUBLE;
            bbcount = 1;
            CCutil_stop_timer (&lp->stats.total, 0);
            printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                    CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                    CClp_nnonzeros (lp->lp));
            goto DONE;
        }
    }

    if (ncount < 1000000) {
        double tourval;
        CCutil_start_timer (&lp->stats.linkern);
        rval = CCtsp_call_x_heuristic (lp, &tourval, besttour, silent, &rstate);
        CCcheck_rval (rval, "CCtsp_call_x_heuristic failed");
        CCutil_stop_timer (&lp->stats.linkern, silent);

        if (tourval < lp->upperbound) {
            printf ("New upperbound from x-heuristic: %.2f\n", tourval);
            lp->upperbound = tourval;
            rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                                   (char *) NULL, 0, silent);
            CCcheck_rval (rval, "CCtsp_dumptour failed");
        }
        printf ("Final lower bound %f, upper bound %f\n", lp->lowerbound,
                                                          lp->upperbound);
        fflush (stdout);
    }

    if (xfname) {
        rval = CCtsp_dump_x (lp, xfname);
        CCcheck_rval (rval, "CCtsp_dump_x failed");
    }
    if (want_rcnearest) {
        rval = dump_rc (lp, want_rcnearest, probname, 0);
        CCcheck_rval (rval, "dump_rc failed");
    }
    if (dump_the_cliques) {  /* Used in BONN to get list of cuts July 2013 */
        rval = dump_cliques (lp);
        CCcheck_rval (rval, "dump_cliques failed");
    }

#if 0  /* HACK TO DUMP TSP AS LP (WITH CPLEX) */
    printf ("Dump the LP problem\n"); fflush (stdout);
    cplex_write_prob (lp->lp);  /* To write the LP */
    printf ("DONE\n"); exit (1);
#endif

    if (lp->graph.ncount < 150000 || complete_price) {
        CCbigguy bound;
        CCbigguy bupper;

        if (dat.ndepot > 0) eliminate_edges = 0;
        rval = CCtsp_exact_price (lp, &bound, complete_price, 0,
                                  eliminate_edges, 0, silent);
        CCcheck_rval (rval, "CCtsp_exact_price failed");
        printf ("Exact lower bound: %.6f\n", CCbigguy_bigguytod (bound));
        printf ("DIFF: %f\n", lp->lowerbound - CCbigguy_bigguytod (bound));
        fflush (stdout);

        bupper = CCbigguy_dtobigguy (lp->upperbound);
        CCbigguy_sub (&bupper, CCbigguy_ONE);

        if (CCbigguy_cmp (lp->exact_lowerbound, bupper) > 0) {
            upbound = lp->upperbound;
            bbcount = 1;
            if (!dfs_branching && !bfs_branching) {
                if (!just_cuts) {
                    printf ("Optimal Solution: %.2f\n", upbound);
                    printf ("Number of bbnodes: %d\n", bbcount);
                    fflush (stdout);
                }
            }
            CCutil_stop_timer (&lp->stats.total, silent);
            printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                    CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                    CClp_nnonzeros (lp->lp));

            if (dat.ndepot > 0) {
                rval = CCtsp_depot_valid (lp, dat.ndepot, (int *) NULL);
                CCcheck_rval (rval, "CCtsp_depot_valid failed");
            }
            goto DONE;
        }

        if (dat.ndepot == 0 && eliminate_edges) {
            /* Eliminating with exact bound can remove more edges */
            rval = CCtsp_exact_price (lp, &bound, 0, 0, 1, 1, silent);
            CCcheck_rval (rval, "CCtsp_exact_price failed");
            printf ("Re-Exact lower bound: %.6f\n", CCbigguy_bigguytod (bound));
        }
    } else {
        printf ("During testing, do not exact price large problems\n");
        fflush (stdout);
        CCutil_stop_timer (&lp->stats.total, 0);
        printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
                CClp_nrows (lp->lp), CClp_ncols (lp->lp),
                CClp_nnonzeros (lp->lp));

        goto DONE;
    }

    CCutil_stop_timer (&lp->stats.total, 0);
    printf ("Final LP has %d rows, %d columns, %d nonzeros\n",
            CClp_nrows (lp->lp), CClp_ncols (lp->lp),
            CClp_nnonzeros (lp->lp));
    fflush (stdout);

    if (dat.ndepot > 0) {
        rval = CCtsp_depot_valid (lp, dat.ndepot, (int *) NULL);
        CCcheck_rval (rval, "CCtsp_depot_valid failed");
        goto DONE;
    }
    
    if (dfs_branching) {
        upbound = lp->upperbound;
        bbcount = 0;

        if (simple_branching) CCtsp_init_simple_cutselect (&sel);
        rval = CCtsp_easy_dfs_brancher (lp, &sel, 0, &upbound, &bbcount,
                     usebranchcliques, besttour, longedge_branching,
                     simple_branching, (double *) NULL, (int *) NULL, silent,
                     &rstate);
        CCcheck_rval (rval, "CCtsp_easy_dfs_brancher failed");
    } else if (bfs_branching) {
        double lowbound = lp->lowerbound;
        int id          = lp->id;

        upbound  = lp->upperbound;
        bbcount = 0;

        rval = CCtsp_write_probroot_id (problname, lp);
        CCcheck_rval (rval, "CCtsp_write_probroot_id failed");
        CCtsp_free_tsp_lp_struct (&lp);

        rval = CCtsp_bfs_brancher (problname, id, lowbound, &sel, 
                &tentativesel, &upbound, &bbcount, usebranchcliques, &dat,
                ptour, pool, dominopool, ncount, besttour, hostport,
                &branchzeit, save_proof, tentative_branch_num,
                longedge_branching, (double *) NULL, (int *) NULL, silent,
                &rstate);
        CCcheck_rval (rval, "CCtsp_bfs_brancher failed");
    }

DONE:
    if (dfs_branching || bfs_branching || restartfname) {
        if (!just_cuts) {
            printf ("Optimal Solution: %.2f\n", upbound);
            printf ("Number of bbnodes: %d\n", bbcount);
            fflush (stdout);
            rval = CCtsp_dumptour (ncount, &dat, ptour, probname, besttour,
                                   outfname, output_tour_as_edges, silent);
            CCcheck_rval (rval, "CCtsp_dumptour failed");
        }
    } else {
        if (dump_root) {
/*
            CCtsp_dump_edge_all (lp, "all.edg");
            exit (1);
*/
            rval = CCtsp_write_probroot_id (problname, lp);
            CCcheck_rval (rval, "CCtsp_write_probroot_id failed");
        } else {
            rval = CCtsp_write_probfile_sav (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
        }
    }

    printf ("Total Running Time: %.2f (seconds)", CCutil_zeit () - szeit);
    if (branchzeit != 0.0) {
        printf ("  Branching Time: %.2f (seconds)", branchzeit);
    }
    printf ("\n"); fflush (stdout);

    /*  CCtsp_output_statistics (&lp->stats);  */

    if (pool && pool->cutcount) {
        if (!silent) {
            printf ("Final Pool: %d cuts\n", pool->cutcount); fflush (stdout);
        }
        sprintf (buf, "%s.pul", probname);
        rval = CCtsp_write_cutpool (ncount, buf, pool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }

#ifdef CCtsp_USE_DOMINO_CUTS
    if (dominopool && dominopool->cutcount) {
        if (!silent) {
            printf ("Final Domino Pool: %d cuts\n", dominopool->cutcount);
            fflush (stdout);
        }
        sprintf (buf, "%s.dompul", probname);
        rval = CCtsp_write_cutpool (ncount, buf, dominopool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }
#endif

    if (sel.remotepool && pool && pool->cutcount > pool->savecount) {
        rval = CCtsp_send_newcuts (ncount, pool, sel.remotehost,
                sel.remoteport);
        if (rval) {
            fprintf (stderr, "CCtsp_send_newcuts failed\n");
            rval = 0;
        }
    }
        
    rval = 0;

CLEANUP:
    if (unlink_files) {
        if (!run_silently) {
            printf ("Delete the temporary files: pul sav mas\n");
            fflush (stdout);
        }

        sprintf (buf, "%s.pul", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "O%s.pul", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "%s.sav", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "O%s.sav", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "%s.mas", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }

        sprintf (buf, "O%s.mas", probname);
        rval = unlink (buf);
        if (rval && !run_silently) {
            printf ("CCutil_sdelete_file failed for %s\n", buf);
        }
    }

    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    if (pool) { CCtsp_free_cutpool (&pool); }
    if (dominopool) { CCtsp_free_cutpool (&dominopool); }

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CC_IFFREE (exlist, int);
    CC_IFFREE (exlen, int);
    CC_IFFREE (ptour, int);
    CC_IFFREE (snowtour, int);
    CC_IFFREE (besttour, int);
    CC_IFFREE (probname, char);
    CCutil_freedatagroup (&dat);

    return rval;
}

int CCtsp_test_general_cuts (CCtsp_lp *lp, CCrandstate *rstate);

static int handle_just_cuts (CCtsp_lp *lp, int the_cuts, int tmax,
       CCrandstate *rstate, int silent)
{
    int rval = 0, infeasible = 0;
    CCtsp_cutselect sel;

    if (the_cuts == CC_JUST_FAST_CUTS) {
        printf ("Using only fast-cut routine\n"); fflush (stdout);
        CCtsp_init_fast_cutselect (&sel);
        rval = CCtsp_cutselect_set_tols (&sel, lp, -1, silent);
        CCcheck_rval (rval, "CCtsp_cutselect_set_tols failed");
        rval = CCtsp_cutting_loop (lp, &sel, NULL, 1, silent, rstate,
                                   &infeasible);
        CCcheck_rval (rval, "CCtsp_cutting_loop failed");
    } else if (the_cuts == CC_JUST_SUBTOUR) {
        printf ("Solve over subtour polytope\n"); fflush (stdout);
        rval = CCtsp_subtour_loop (lp, silent, 0.0001, rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_subtour_loop failed");
    } else if (just_cuts == CC_JUST_BLOSSOM) {
        printf ("Solve over two-factor polytope\n"); fflush (stdout);
        rval = CCtsp_blossom_loop (lp, silent, rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_blossom_loop failed");
    } else if (just_cuts == CC_JUST_SUBTOUR_AND_BLOSSOM) {
        printf ("Using only subtours and trivial blossoms\n"); fflush (stdout);
        rval = CCtsp_subtour_and_blossom_loop (lp, silent, rstate, filecutname,
                                               &infeasible);
        CCcheck_rval (rval, "CCtsp_subtour_and_blossom_loop failed");
    } else if (just_cuts == CC_JUST_WORLD_CUTS) {
        printf ("Using only world cuts\n"); fflush (stdout);
        rval = CCtsp_world_loop (lp, 0, tmax, rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_world_loop failed");
#ifdef CCtsp_USE_DOMINO_CUTS
    } else if (just_cuts == CC_JUST_DP) {
        printf ("Using only DP+TP cuts\n"); fflush (stdout);
        rval = CCtsp_domino_loop (lp, silent, rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_domino_loop failed");
#endif
    } else if (just_cuts == CC_JUST_GENERAL) {
        printf ("Testing general-form of cuts\n"); fflush (stdout);
        rval = CCtsp_test_general_cuts (lp, rstate);
        CCcheck_rval (rval, "CCtsp_test_general_cuts failed");
    }

    if (infeasible) {
        printf ("Infeasible LP after cutting loop\n"); fflush (stdout);
        rval = 1; goto CLEANUP;
    }

    printf ("Bound: %f\n", lp->lowerbound); fflush (stdout);
    CCutil_stop_timer (&lp->stats.total, 0);
    printf ("Final Root LP has %d rows, %d columns, %d nonzeros\n",
            CClp_nrows (lp->lp), CClp_ncols (lp->lp), CClp_nnonzeros (lp->lp));

CLEANUP:

    return rval;
}

static int run_hk (int ncount, CCdatagroup *dat, int *hk_tour)
{
    double hk_val;
    int hk_found, hk_yesno;
    int *hk_tlist = (int *) NULL;
    int rval = 0;

    hk_tlist = CC_SAFE_MALLOC (2*ncount, int);
    CCcheck_NULL (hk_tlist, "out of memory for hk_tlist");

    rval = CCheldkarp_small (ncount, dat, (double *) NULL, &hk_val,
                             &hk_found, 0, hk_tlist, 1000000, 2);
    CCcheck_rval (rval, "CCheldkarp_small failed");
    printf ("Optimal Solution: %.2f\n", hk_val); fflush (stdout);

    rval = CCutil_edge_to_cycle (ncount, hk_tlist, &hk_yesno, hk_tour);
    CCcheck_rval (rval, "CCutil_edge_to_cycle failed");

    if (hk_yesno == 0) {
        fprintf (stderr, "Held-Karp returned list that is not a tour\n");
        rval = 1;  goto CLEANUP;
    }

CLEANUP:

     CC_IFFREE (hk_tlist, int);
     return rval;
}

static void adjust_upbound (double *bound, int ncount, CCdatagroup *dat)
{
    double bnd;
    int i;

    bnd = CCutil_dat_edgelen (ncount - 1, 0, dat);
    for (i = 1; i < ncount; i++) {
        bnd += CCutil_dat_edgelen (i-1, i, dat);
    }
    if (bnd < *bound) {
        printf ("Set initial upperbound to %.0f (from tour)\n", bnd);
        fflush (stdout);
        *bound = bnd;
    }
}

static int build_edges (int *p_ecount, int **p_elist, int **p_elen,
        int ncount, int *ptour, CCdatagroup *dat, char *in_edgefname,
        char *in_edgegenfname, int in_just_cuts, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    int *elist = (int *) NULL;
    int ecount;
    int i;
    
    if (in_edgefname) {
        int *invperm = (int *) NULL;

        printf ("Read initial edge set\n"); fflush (stdout);
        
        rval = CCutil_getedgelist (ncount, in_edgefname, p_ecount, p_elist,
                                   p_elen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");
        ecount = *p_ecount;
        elist = *p_elist;
        printf ("Initial edgeset: %d edges (%d nodes)\n", ecount, ncount);
        printf ("Rearrange the edges to match the tour order\n");
        fflush (stdout);

        invperm = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (invperm, "out of memory for invperm");
        for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;
        for (i = 0; i < 2*ecount; i++) elist[i] = invperm[elist[i]];
        CC_FREE (invperm, int);
    } else if (dat) {
        CCedgegengroup plan;
        
        if (in_edgegenfname) {
            rval = CCedgegen_read (in_edgegenfname, &plan);
            CCcheck_rval (rval, "CCedgegen_read failed");
        } else {
            CCedgegen_init_edgegengroup (&plan);
            if (in_just_cuts == CC_JUST_SUBTOUR ||
                in_just_cuts == CC_JUST_BLOSSOM ||
                in_just_cuts == CC_JUST_SUBTOUR_AND_BLOSSOM ||
                in_just_cuts == CC_JUST_WORLD_CUTS) {
                plan.tour.greedy = 1;
                plan.f2match_nearest.number = 4;
            } else {
                plan.linkern.count = 10;
                plan.linkern.quadnearest = 2;
                plan.linkern.greedy_start = 0;
                plan.linkern.nkicks = (ncount / 100) + 1;
            }
        }

        rval = getedges (dat, &plan, ncount, p_ecount, p_elist, p_elen,
                         silent, rstate);
        CCcheck_rval (rval, "getedges failed");
    }

CLEANUP:
    return rval;
}

static int build_fulledges (int *p_excount, int **p_exlist, int **p_exlen,
        int ncount, int *ptour, char *in_fullfname)
{
    int i;
    int rval = 0;
    int *exlist;
    int excount;
    
    if (in_fullfname) {
        int *invperm = (int *) NULL;

        rval = CCutil_getedgelist (ncount, in_fullfname, p_excount, p_exlist,
                                   p_exlen, 0);
        CCcheck_rval (rval, "CCutil_getedgelist failed");

        invperm = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (invperm, "out of memory for invperm");
        for (i = 0; i < ncount; i++) invperm[ptour[i]] = i;
        excount = *p_excount;
        exlist = *p_exlist;
        for (i = 0; i < 2*excount; i++) exlist[i] = invperm[exlist[i]];
        CC_FREE (invperm, int);
    } else {
        *p_excount = 0;
    }

CLEANUP:

    return rval;
}

static int find_tour (int ncount, CCdatagroup *dat, int *perm, double *ub,
        int trials, int silent, CCrandstate *rstate)
{
    int rval = 0;
    CCedgegengroup plan;
    int ecount;
    int *elist = (int *) NULL;
    int tcount;
    int *tlist = (int *) NULL;
    int *bestcyc = (int *) NULL;
    int *cyc     = (int *) NULL;
    int *tmp;
    double val, bestval, szeit;
    int kicks, i, istour;

    szeit = CCutil_zeit ();
    bestval = CCtsp_LP_MAXDOUBLE;

    if (trials == -1) {
        kicks = (ncount > 400 ? 100 : ncount/4);
    } else {
        kicks = (ncount > 1000 ? 500 : ncount/2);
    }

    if (!silent) {
        printf ("Finding a good tour for compression: %d\n", trials);
        fflush (stdout);
    }

    cyc = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (cyc, "out of memory for cyc");
    bestcyc = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (bestcyc, "out of memory for bestcyc");

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
                    cyc, bestcyc, &bestval, silent, 0.0, 0.0,
                    (char *) NULL,
                    CC_LK_GEOMETRIC_KICK, rstate);
    CCcheck_rval (rval, "CClinkern_tour failed");

    for (i = 0; i < trials; i++) {
        rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
                        (int *) NULL, cyc, &val, silent, 0.0, 0.0,
                        (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate);
        CCcheck_rval (rval, "CClinkern_tour failed");
        if (val < bestval) {
            CC_SWAP (cyc, bestcyc, tmp);
            bestval = val;
        }
    }

    if (trials > 0) {
        rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, 2 * kicks,
                        bestcyc, perm, ub, silent, 0.0, 0.0,
                        (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate);
        CCcheck_rval (rval, "CClinkern_tour failed");
    } else {
        for (i = 0; i < ncount; i++) {
            perm[i] = bestcyc[i];
        }
    }

    if (!silent) {
        printf ("Time to find compression tour: %.2f (seconds)\n",
                CCutil_zeit() - szeit);
        fflush (stdout);
    }

CLEANUP:
    CC_IFFREE (cyc, int);
    CC_IFFREE (bestcyc, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (tlist, int);
    return rval;
}

static int getedges (CCdatagroup *dat, CCedgegengroup *plan, int ncount,
        int *ecount, int **elist, int **elen, int silent,
        CCrandstate *rstate)
{
    int i;
    int rval = 0;

    *elist = (int *) NULL;
    *elen = (int *) NULL;

    if (dat == (CCdatagroup *) NULL || plan == (CCedgegengroup *) NULL) {
        fprintf (stderr, "getedges needs CCdatagroup and CCedgegengroup\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCedgegen_edges (plan, ncount, dat, (double *) NULL, ecount,
                           elist, silent, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    *elen = CC_SAFE_MALLOC(*ecount, int);
    CCcheck_NULL (*elen, "out of memory for elen");

    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i], (*elist)[(2*i)+1], dat);
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    return rval;
}

static int dump_rc (CCtsp_lp *lp, int count, char *pname, int usesparse)
{
    int rval = 0;
    char rcnname[1024];

/*
    CCtsp_dump_edge_all (lp, "all.edg");
    exit (1);
*/

    if (count == -1) {
        sprintf (rcnname, "%s.rcall", pname);
        rval = CCtsp_dump_rc_all (lp, rcnname);
        CCcheck_rval (rval, "CCtsp_dump_rc failed");
    } else {
        sprintf (rcnname, "%s.rcn", pname);
        rval = CCtsp_dump_rc_nearest (lp, count, rcnname, usesparse);
        CCcheck_rval (rval, "CCtsp_dump_rc failed");
    }

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    /* remaining: AbcGlnQiVwW */

    while ((c = CCutil_bix_getopt (ac, av, "aBC:dD:e:E:fF:g:HhI:jJ:k:K:L:mM:N:o:O:pP:qr:R:s:S:t:T:u:UvVX:xyY:z:Z:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'a': usextighten = 1;  break;
        case 'B': bfs_branching = 0;  break;
        case 'C': maxchunksize = atoi(boptarg);  break;
        case 'd': dfs_branching = 1; bfs_branching = 0;  break;
        case 'D': edgegenfname = boptarg;  break;
        case 'e': edgefname = boptarg;  break;
        case 'E': fullfname = boptarg; valid_edges = 1;  break;
        case 'f': output_tour_as_edges = 1;  break;
        case 'F': filecutname = boptarg;  break;
        case 'g': grunthostname = boptarg;  break;
        case 'h': be_nethost = 1;  break;
        case 'I': just_cuts = atoi (boptarg);  break;
        case 'j': usekeepcutting = 1;  break;
        case 'J': tentative_branch_num = atoi (boptarg);  break;
        case 'k': nnodes_want = atoi (boptarg);  break;
        case 'K': cutbossname = boptarg;  break;
        case 'L': domcutbossname = boptarg;  break;
        case 'm': multiple_chunker = 1;  break;
        case 'M': masterfname = boptarg;  break;
        case 'o': outfname = boptarg;  break;
        case 'O': moatfname = boptarg;  break;
        case 'p': save_proof = 1;  break;
        case 'P': poolfname = boptarg;  break;
        case 'q': dontcutroot = 1;  break;
        case 'r': gridsize = atoi(boptarg);  break;
        case 'R': restartfname = boptarg;  break;
        case 's': seed = atoi (boptarg);  break;
        case 'S': probfname = boptarg;  break;
        case 't': tourfname =  boptarg;  break;
#ifdef CCtsp_USE_DOMINO_CUTS
        case 'Z': usedominos = atoi(boptarg);  break;
        case 'Y': dominopoolfname = boptarg;  break;
        case 'T': dombossname = boptarg; if (!usedominos) usedominos = 1; break;
        case 'H': domino_safeshrink = 0;  break;
#endif
        case 'u': initial_ub = atof (boptarg);  break;
        case 'U': dump_root = 1;  break;
        case 'v': run_silently = 0;  break;
        case 'X': xfname = boptarg;  break;
        case 'x': unlink_files = 1;  break;
        case 'y': simple_branching = 1;  break;
        case 'z': want_rcnearest = atoi (boptarg);  break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM;  break;
            case 1: norm = CC_MANNORM;  break;
            case 2: norm = CC_EUCLIDEAN;  break;
            case 3: norm = CC_EUCLIDEAN_3D;  break;
            case 4: norm = CC_USER;  break;
            case 5: norm = CC_ATT;  break;
            case 6: norm = CC_GEOGRAPHIC;  break;
            case 7: norm = CC_MATRIXNORM;  break;
            case 8: norm = CC_DSJRANDNORM;  break;
            case 9: norm = CC_CRYSTAL;  break;
            case 10: norm = CC_SPARSE;  break;
            case 11: norm = CC_RHMAP1;  break;
            case 12: norm = CC_RHMAP2;  break;
            case 13: norm = CC_RHMAP3;  break;
            case 14: norm = CC_RHMAP4;  break;
            case 15: norm = CC_RHMAP5;  break;
            case 16: norm = CC_EUCTOROIDAL;  break;
            case 17: norm = CC_GEOM;  break;
            case 18: norm = CC_EUCLIDEAN_CEIL;  break;
            case 20: norm = CC_ROAD;  break;
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
    if (boptind < ac) {
        datfname = av[boptind++];
    }

    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }

    if (datfname == (char *) NULL && nnodes_want == 0 && 
        probfname == (char *) NULL && edgefname == (char *) NULL &&
        masterfname == (char *) NULL && grunthostname == (char *) NULL) {
        usage (av[0]);
        return 1;
    }

    if (masterfname == (char *) NULL && datfname == (char *) NULL &&
          edgefname != (char *) NULL) {
        fprintf (stderr, "cannot give edgefile without a dat or master file\n");
        return 1;
    }

    if (eliminate_edges < 0) {
        eliminate_edges = (bfs_branching || dfs_branching || dump_root);
    }

    /* STAR Branching: turn off elimination 
    eliminate_edges = 0;
    */

    if (dump_root) {
        bfs_branching = 0;
        complete_price = 1;
    }
    
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] [dat_file]\n", f);
    fprintf (stderr, "   -a    use xtighten (mondo tighten)\n");
    fprintf (stderr, "   -B    do not branch\n");
    fprintf (stderr, "   -C #  maximum chunk size in localcuts (default 16)\n");
    fprintf (stderr, "   -d    use dfs branching instead of bfs\n");
    fprintf (stderr, "   -D f  edgegen file for initial edge set\n");
    fprintf (stderr, "   -e f  initial edge file\n");
    fprintf (stderr, "   -E f  full edge file (must contain initial edge set)\n");
    fprintf (stderr, "   -f    write optimal tour as edge file (default is tour file)\n");
    fprintf (stderr, "   -F f  read extra cuts from text file\n");
    fprintf (stderr, "   -g h  be a grunt for boss h\n");
    fprintf (stderr, "   -h    be a boss for the branching\n");
    fprintf (stderr, "   -I #  use only specified class of inequalities\n");
    fprintf (stderr, "         %d=Subtour Polytope, %d=Two-Factor Polytope, %d=Subtours+Blossoms\n", CC_JUST_SUBTOUR, CC_JUST_BLOSSOM, CC_JUST_SUBTOUR_AND_BLOSSOM);
    fprintf (stderr, "         %d=Fast_Cuts, %d=World, %d=DP+TP, %d=General_Test\n", CC_JUST_FAST_CUTS, CC_JUST_WORLD_CUTS, CC_JUST_DP, CC_JUST_GENERAL);
    fprintf (stderr, "   -j    keep cutting the root LP\n");
    fprintf (stderr, "   -J #  number of tentative branches\n");
    fprintf (stderr, "   -k #  number of nodes for random problem\n");
    fprintf (stderr, "   -K h  use cut server h\n");
    fprintf (stderr, "   -L h  use domino cutpool server h\n");
    fprintf (stderr, "   -M f  master file\n");
    fprintf (stderr, "   -m    use multiple passes of cutting loop\n");
    fprintf (stderr, "   -o f  output file name (for optimal tour)\n");
    fprintf (stderr, "   -p    save the branching tree\n");
    fprintf (stderr, "   -P f  cutpool file\n");
    fprintf (stderr, "   -q    do not cut the root lp\n");
    fprintf (stderr, "   -r #  use #x# grid for random points, no dups if #<0\n");
    fprintf (stderr, "   -R f  restart file\n");
    fprintf (stderr, "   -s #  random seed\n");
    fprintf (stderr, "   -S f  problem file\n");
    fprintf (stderr, "   -t f  tour file (in node node node format)\n");
    fprintf (stderr, "   -u v  initial upperbound\n");
    fprintf (stderr, "   -U    just dump a root LP (after elimination)\n");
    fprintf (stderr, "   -v    verbose (turn on lots of messages)\n");
    fprintf (stderr, "   -x    delete files on completion (sav pul mas)\n");
    fprintf (stderr, "   -X f  write the last root fractional solution to f\n");
    fprintf (stderr, "   -y    use simple cutting and branching in DFS\n");
    fprintf (stderr, "   -z #  dump the #-lowest reduced cost edges to file xxx.rcn\n");
#ifdef CCtsp_USE_DOMINO_CUTS
    fprintf (stderr, "   -Z #  dp-cuts (#=1 normal, #=2 shrunk, #=3 both)\n");
    fprintf (stderr, "   -Y f  domino-cutpool file\n");
    fprintf (stderr, "   -T h  use domino server h\n");
    fprintf (stderr, "   -H    do not use safe-shrinking before dominos\n");
#endif
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 3=3D, 4=USER, 5=ATT, 6=GEO, 7=MATRIX,\n");
    fprintf (stderr, "         8=DSJRAND, 9=CRYSTAL, 10=SPARSE, 11-15=RH-norm 1-5, 16=TOROIDAL\n");
    fprintf (stderr, "         17=GEOM, 18=JOHNSON, 20=ROAD\n");
}

/* dump_cliques used in BONN, July 2013, to build a clique graph  */

static int dump_cliques (CCtsp_lp *lp)
{
    int rval = 0, i, j, acount, *a = (int *) NULL;
    int *clique_nums = (int *) NULL;
    double *clique_vals = (double *) NULL;
    int nclique, ncount = lp->graph.ncount, xcount, *xlist = (int *) NULL;
    double *x = (double *) NULL;
    CCtsp_lpclique *c;

    printf ("Dump the cliques\n");

    rval = CCtsp_get_lp_result (lp, NULL, NULL, &xcount, &xlist, &x, NULL,
                                 NULL, NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

    rval = CCtsp_get_clique_prices (lp->pool, &clique_nums, &clique_vals,
                            1.0, 9.99, &nclique, ncount, xcount, xlist, x);
    CCcheck_rval (rval, "CCtsp_get_clique_prices failed");

    printf ("Number of Cliques: %d\n", nclique); fflush (stdout);

    for (i = 0; i < nclique; i++) {
        fflush (stdout);
        if (i == 1000000) {
            printf ("HIT 1000000\n"); return 0;
        }
        rval = CCtsp_get_clique (lp->pool, clique_nums[i], &c);
        CCcheck_rval (rval, "CCtsp_get_clique failed");
        rval = CCtsp_clique_to_array (c, &a, &acount);
        CCcheck_rval (rval, "CCtsp_clique_to_array failed");
        printf ("%f", clique_vals[i]);
        printf (" %d", acount);
        for (j = 0; j < acount; j++) {
            printf (" %d", lp->perm[a[j]]);
        }
        printf ("\n");
        CC_IFFREE (a, int);
    }

CLEANUP:
    return rval;
}

/* quick call to moat building for Hilldale lecture, March 2014 */

typedef struct CCtsp_moat {
    double pi;
    int count;
    int *nodes;
} CCtsp_moat;

typedef struct CCtsp_moatlist {
    int count;
    double *zones;
    CCtsp_moat *moats;
    double bnd;
    int *tour;
    double tourlen;
    int ready;
} CCtsp_moatlist;


void CCtsp_moatlist_init (CCtsp_moatlist *m);
void CCtsp_moatlist_free (CCtsp_moatlist *m);
int CCtsp_geom_dual (int ncount, CCdatagroup *indat, CCtsp_moatlist *moats,
    int boundtype, CCrandstate *rstate);

static int find_moats (int ncount, CCdatagroup *dat, int *perm,
        CCrandstate *rstate, char *fname)
{
    int rval = 0, boundtype = 2, i, j;
    int *invperm = (int *) NULL;
    FILE *out = (FILE *) NULL;
    CCtsp_moatlist moats;

    CCtsp_moatlist_init (&moats);
    rval = CCtsp_geom_dual (ncount, dat, &moats, boundtype, rstate);
    CCcheck_rval (rval, "CCtsp_geom_dual failed"); 

    CC_MALLOC (invperm, ncount, int);
    for (i = 0; i < ncount; i++) invperm[perm[i]] = i;

    out = fopen (fname, "w");
    if (out == (FILE *) NULL) {
        fprintf (stderr, "could not open %s for reading\n", fname);
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "%d\n", ncount);
    for (i = 0; i < ncount; i++) {
        if (moats.zones[invperm[i]] < 0.00001) {
            fprintf (out, "0.000\n");
        } else {
            fprintf (out, "%.3f\n", moats.zones[invperm[i]]);
        }
    }
    fprintf (out, "%d\n", moats.count);
    for (i = 0; i < moats.count; i++) {
        fprintf (out, "%f", moats.moats[i].pi);
        fprintf (out, " %d", moats.moats[i].count);
        for (j = 0; j < moats.moats[i].count; j++) {
            fprintf (out, " %d", perm[moats.moats[i].nodes[j]]);
        }
        fprintf (out, "\n");
    }
   
CLEANUP:
    CCtsp_moatlist_free (&moats);
    CC_IFFREE (invperm, int);
    if (out) fclose (out);
    return rval;
}

