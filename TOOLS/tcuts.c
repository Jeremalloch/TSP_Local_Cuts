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
/*       A PROGRAM TO GATHER CUTS FROM LPS FOR A SUBDIVIDED PROBLEM         */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: September 18, 2006                                                */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

static char *tourfname   = (char *) NULL;
static char *tspfname    = (char *) NULL;
static char *masterfname = (char *) NULL;
static char *indexfname  = (char *) NULL;
static char *borderfname  = (char *) NULL;
static int seed = 0;
static int tsplib_tour = 1;
static int tsplib_in = 1;
static int create_lp = 0;
static int add_subtours = 0;
static int innorm = CC_EUCLIDEAN;

int
    main (int ac, char **av);

static int
    build_edges (int ncount, CCdatagroup *dat, int *ecount, int **elist,
        int **elen, CCrandstate *rstate),
    check_file (char *name, int id, int chk_ncount, int global_ncount,
       int *global_invperm, CCtsp_lpcuts *pool, CCtsp_lpcuts *dompool,
       int *cnt, CCrandstate *rstate, CCtsp_lp *global_lp,
       CCdatagroup *global_dat, CCutil_edgehash *h),
    check_cut_for_depots (CCtsp_lpcut_in *c, int ncount, int *yesno),
    get_names (char *fname, int *names, int ncount),
    lp_value (CCtsp_lp *lp, double *val),
    lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x),
    lp_solve (CCtsp_lp *lp),
    add_subtours_to_global_lp (CCtsp_lp *lp),
    price_edges_in_global_lp (CCtsp_lp *lp, CCrandstate *rstate),
    add_cuts (int count, CCtsp_lpcut_in **cuts, CCtsp_lp *lp, CCtsp_lprow *cr),
    get_blobcuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount, int ecount,
        int *elist, double *x),
    blobcuts (int ncount, int ecount, int *elist, double *x, int *hit,
        int (*doit_fn) (double, int, int *, void *), void *pass_param),
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, ncount, ecount = 0, pcount, cnt, subcount = 0, rval = 0;
    int bncount, bpcount, infeasible = 0, havehash = 0;
    int *ptour = (int *) NULL, *inv_ptour = (int *) NULL;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    char *name = (char *) NULL, *bname = (char *) NULL;
    char pbuf[1024], probname[1024];
    double szeit, rzeit, initial_ub = CCtsp_LP_MAXDOUBLE;
    CCdatagroup dat;
    CCsubdiv *trac = (CCsubdiv *) NULL, *btrac = (CCsubdiv *) NULL;
    CCrandstate rstate;
    CCtsp_lpcuts *gpool = (CCtsp_lpcuts *) NULL;
    CCtsp_lpcuts *gdompool = (CCtsp_lpcuts *) NULL;
    CCtsp_lp *glp = (CCtsp_lp *) NULL;
    CCutil_edgehash h;

    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

#if 0
{
    int numcpu = 0;

    numcpu = sysconf(_SC_NPROCESSORS_ONLN);

    printf ("Number of Available CPUS: %d\n", numcpu);
    fflush (stdout);
}
#endif

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname && (!tourfname || !tspfname)) {
        fprintf (stderr, "Must specify either master or tour+tsp\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (!indexfname) {
        fprintf (stderr, "Must specify an index file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    CCutil_printlabel ();

    szeit = CCutil_zeit ();
    rzeit = CCutil_real_zeit ();
    CCutil_sprand (seed, &rstate);

    if (masterfname) {
        rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
        CCcheck_rval (rval, "CCutil_getmaster failed")
    } else {
        if (tsplib_in) {
            rval = CCutil_gettsplib (tspfname, &ncount, &dat);
            CCcheck_rval (rval, "CCutil_gettsplib failed");
            CCutil_dat_getnorm (&dat, &innorm);
        } else {
            rval = CCutil_getdata (tspfname, 0, innorm, &ncount, &dat,
                                   0, 0, &rstate);
            CCcheck_rval (rval, "CCutil_getdata failed");
        }

        CC_MALLOC (ptour, ncount, int);
        if (tsplib_tour) {
            rval = CCutil_getcycle_tsplib (ncount, tourfname, ptour);
            CCcheck_rval (rval, "CCutil_getcycle failed");
        } else {
            rval = CCutil_getcycle (ncount, tourfname, ptour, 0);
            CCcheck_rval (rval, "CCutil_getcycle failed");
        }

        rval = CCutil_datagroup_perm (ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_datagroup_perm failed");
    }

    CC_MALLOC (inv_ptour, ncount, int);
    for (i = 0; i < ncount; i++) inv_ptour[ptour[i]] = i;

    rval = CCutil_read_subdivision_index (indexfname, &name, &ncount, &pcount,
                                          &trac);
    CCcheck_rval (rval, "CCutil_read_subdivision_index failed");

    printf ("Name: %s\n", name);
    printf ("ncount = %d, partitions = %d\n", ncount, pcount);
    fflush (stdout);

    if (!masterfname) {
        printf ("Write a master file for global problem\n"); fflush (stdout);
        sprintf (pbuf, "%s.mas", name);
        rval = CCutil_putmaster (pbuf, ncount, &dat, ptour);
        CCcheck_rval (rval, "CCutil_putmaster failed")
    }

    if ((innorm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "Only set up for 2D norms\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &gpool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");
    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &gdompool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    if (create_lp) {
        printf ("build initial global LP\n"); fflush (stdout);
        rval = build_edges (ncount, &dat, &ecount, &elist, &elen, &rstate);
        CCcheck_rval (rval, "build_edges failed");

        sprintf (probname, "%s_global", name); 
        rval = CCtsp_init_lp (&glp, probname, -1, (char *) NULL, ncount, &dat,
                  ecount, elist, elen,                 /* elist  */
                  0, (int *) NULL, (int *) NULL, 0,    /* exlist */
                  ptour, initial_ub,
                  (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 1, &rstate,
                  &infeasible);
        CCcheck_rval (rval, "CCtsp_init_lp failed");
        if (infeasible) {
            fprintf (stderr, "initial LP is infeasible\n");
            rval = 1; goto CLEANUP;
        }

        rval = CCutil_edgehash_init (&h, 2*ecount);
        CCcheck_rval (rval, "CCutil_edgehash_init failed");
        havehash = 1;
        for (i = 0; i < ecount; i++) {
            rval = CCutil_edgehash_set (&h, elist[2*i], elist[2*i+1], elen[i]);
            CCcheck_rval (rval, "CCutil_edgehash_set failed");
        }
    }

    for (i = 0; i < pcount; i++) {
        printf ("%d %d %f\n", trac[i].id, trac[i].cnt, trac[i].bound);
        fflush (stdout); 
        rval = check_file (name, trac[i].id, trac[i].cnt, ncount, inv_ptour,
                           gpool, gdompool, &cnt, &rstate, glp, &dat, &h);
        CCcheck_rval (rval, "check_file failed");
        subcount += cnt;
    }

    printf ("Number of Cuts in Pool: %d\n", gpool->cutcount);
    if (gdompool->cutcount > 0) {
        printf ("Number of Cuts in Domino Pool: %d\n", gdompool->cutcount);
    }
    printf ("Number of Subtours: %d\n", subcount);
    fflush (stdout);

    if (borderfname) {
        printf ("Gather cuts from cross-border LPs\n");
        fflush (stdout);

        rval = CCutil_read_subdivision_index (borderfname, &bname, &bncount,
                                              &bpcount, &btrac);
        CCcheck_rval (rval, "CCutil_read_subdivision_index failed");

        if (bncount != ncount) {
            fprintf (stderr, "border index file does not match\n");
            rval = 1; goto CLEANUP;
        }

        for (i = 0; i < bpcount; i++) {
            printf ("%d %d %f\n", btrac[i].id, btrac[i].cnt, btrac[i].bound);
            fflush (stdout); 
            rval = check_file (bname, btrac[i].id, btrac[i].cnt, bncount,
                           inv_ptour, gpool, gdompool, &cnt, &rstate, glp,
                           &dat, &h);
            CCcheck_rval (rval, "check_file failed");
            subcount += cnt;
        }

        printf ("Updated Number of Cuts in Pool: %d\n", gpool->cutcount);
        if (gdompool->cutcount > 0) {
            printf ("Updated Number of Cuts in Pool: %d\n", gdompool->cutcount);
        }
        printf ("Updated Number of Subtours: %d\n", subcount);
        fflush (stdout);
    }

    sprintf (pbuf, "%s_global.pul", name);
    rval = CCtsp_write_cutpool (ncount, pbuf, gpool);
    CCcheck_rval (rval, "CCtsp_write_cutpool failed");

    if (gdompool->cutcount > 0) {
        sprintf (pbuf, "%s_global.dompul", name);
        rval = CCtsp_write_cutpool (ncount, pbuf, gdompool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }

    if (create_lp) {
        double val;
        double lpzeit, lprzeit;

        printf ("Global LP Has %d rows, %d columns, %d nonzeros\n",
                    CClp_nrows (glp->lp), CClp_ncols (glp->lp),
                    CClp_nnonzeros (glp->lp));
        fflush (stdout);

/*
        CClp_dump_lp (glp->lp, "joke");
        exit (1);
*/

        /* printf ("Call dual simplex ....\n"); fflush (stdout); */
        lpzeit = CCutil_zeit ();
        lprzeit = CCutil_real_zeit ();
        rval = CClp_opt (glp->lp, CClp_METHOD_BARRIER, &infeasible); 
        CCcheck_rval (rval, "CClp_opt failed");
        /* rval = CClp_opt (glp->lp, CClp_METHOD_DUAL, &infeasible); */
        if (infeasible) {
            rval = CCtsp_infeas_recover (glp, 1, &rstate, &infeasible);
            CCcheck_rval (rval, "CCtsp_infeas_recover failed");
            if (infeasible) {
                fprintf (stderr, "Could not recover feasible LP\n");
                rval = 1; goto CLEANUP;
            }
        } 
        printf ("LP Time: %.2f seconds\n", CCutil_zeit () - lpzeit);
        printf ("LP Wall Cloc Time: %.2f\n", CCutil_real_zeit () - lprzeit);
        fflush (stdout);
        rval = CCtsp_update_result (glp);
        CCcheck_rval (rval, "CCtsp_update_result failed");
        rval = lp_value (glp, &val);
        CCcheck_rval (rval, "lp_value");
        printf ("Sparse Global LP Value: %.2f\n", val);
        rval = CCtsp_write_probfile_sav (glp);
        CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");

        if (add_subtours) {
            rval = add_subtours_to_global_lp (glp);
            CCcheck_rval (rval, "add_subtours_to_global_lp");
        }

        rval = price_edges_in_global_lp (glp, &rstate);
        CCcheck_rval (rval, "price_edges_in_global_lp"); 
    }

    printf ("Total Running Time: %.2f\n", CCutil_zeit () - szeit);
    printf ("Total Wall Cloc Time: %.2f\n", CCutil_real_zeit () - rzeit);

CLEANUP:
    CC_IFFREE (ptour, int);
    CC_IFFREE (inv_ptour, int);
    CC_IFFREE (name, char);
    CC_IFFREE (bname, char);
    CC_IFFREE (trac, CCsubdiv);
    CC_IFFREE (btrac, CCsubdiv);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    CCutil_freedatagroup (&dat);
    if (gpool) { CCtsp_free_cutpool (&gpool); }
    if (gdompool) { CCtsp_free_cutpool (&gdompool); }
    if (glp) { CCtsp_free_tsp_lp_struct (&glp); }
    if (havehash) CCutil_edgehash_free (&h);
    return rval;
}

static int build_edges (int ncount, CCdatagroup *dat, int *ecount, int **elist,
        int **elen, CCrandstate *rstate)
{
    int i, rval = 0;
    CCedgegengroup plan;

    CCedgegen_init_edgegengroup (&plan);
    plan.linkern.count = 10;
    plan.linkern.quadnearest = 2;
    plan.linkern.greedy_start = 0;
    plan.linkern.nkicks = (ncount / 100) + 1;
    plan.f2match_nearest.number = 10;
/*
    plan.quadnearest = 5;
*/

    rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, ecount,
                            elist, 0 , rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    *elen = CC_SAFE_MALLOC(*ecount, int);
    CCcheck_NULL (*elen, "out of memory for elen");

    for (i = 0; i < *ecount; i++) {
        (*elen)[i] = CCutil_dat_edgelen ((*elist)[2*i], (*elist)[(2*i)+1], dat);    }

CLEANUP:
    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    return rval;
}

static int check_file (char *name, int id, int chk_ncount, int global_ncount,
       int *global_invperm, CCtsp_lpcuts *pool, CCtsp_lpcuts *dompool,
       int *cnt, CCrandstate *rstate, CCtsp_lp *global_lp,
       CCdatagroup *global_dat, CCutil_edgehash *h)
{
    int rval = 0, infeasible = 0;
    int i, scnt, local_ncount, ndepot, ncount, yesno;
    int *local_perm = (int *) NULL, *local_names = (int *) NULL;
    double initial_ub = CCtsp_LP_MAXDOUBLE;
    char local_probname[1024];
    char local_mastername[1024];
    char local_namesname[1024];
    CCdatagroup local_dat;
    CCtsp_lpcut_in new, old;
    CCtsp_lprow cr;
    CCtsp_lpgraph *G;
    CCtsp_predge *pr_list = (CCtsp_predge *) NULL;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;

    CCutil_init_datagroup (&local_dat);
    CCtsp_init_lprow (&cr);

    sprintf (local_mastername, "%s_%d.mas", name, id);
    sprintf (local_probname, "%s_%d.sav", name, id);
    sprintf (local_namesname, "%s_%d.nam", name, id);

    rval = CCutil_getmaster (local_mastername, &local_ncount, &local_dat,
                             &local_perm);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    ndepot = local_dat.ndepot;
    if (ndepot > 0) {
        ncount = local_dat.orig_ncount;
        if (ncount != (local_ncount - ndepot)) {
            fprintf (stderr, "node counts with depots do not agree\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        ncount = local_ncount;
        if (local_dat.orig_names) {
            fprintf (stderr, "orig names but no depots\n");
            rval = 1; goto CLEANUP;
        }
    }

    if (chk_ncount != local_dat.orig_ncount) {
        fprintf (stderr, "ncount if %s does not match index file\n",
                 local_mastername);
        rval = 1;  goto CLEANUP;
    } 


    local_names = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (local_names, "out of memory for local_names");

    if (local_dat.orig_names) {
        for (i = 0; i < ncount; i++) {
            local_names[i] = local_dat.orig_names[i];
        }
    } else {
        rval = get_names (local_namesname, local_names, ncount);
        CCcheck_rval (rval, "get_names failed")
    }

    rval = CCtsp_init_lp (&lp, (char *) NULL, -1, local_probname, local_ncount,
                  &local_dat,
                  0, (int *) NULL, (int *) NULL,       /* elist  */
                  0, (int *) NULL, (int *) NULL, 0,    /* exlist */
                  local_perm, initial_ub,
                  (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 1, rstate,
                  &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    if (global_lp) {
        int end0, end1, len, test, prcnt = 0;
        G = &(lp->graph);
        printf ("Number of local edges: %d\n", G->ecount);
        pr_list = CC_SAFE_MALLOC (G->ecount, CCtsp_predge);
        CCcheck_NULL (pr_list, "out of memory for pr_list");
        for (i = 0; i < G->ecount; i++) {
            end0 = G->edges[i].ends[0];
            end1 = G->edges[i].ends[1];
            if (end0 < ncount && end1 < ncount) {
                end0 = global_invperm[local_names[local_perm[end0]]];
                end1 = global_invperm[local_names[local_perm[end1]]];
                test = CCutil_edgehash_find (h, end0, end1, &len);
                if (test == -1) {
                    len = CCutil_dat_edgelen (end0, end1, global_dat);
                    rval = CCutil_edgehash_add (h, end0, end1, len); 
                    CCcheck_rval (rval, "CCutil_edgehash_add failed");
                    pr_list[prcnt].ends[0] = end0;
                    pr_list[prcnt].ends[1] = end1;
                    pr_list[prcnt].len     = len;
                    pr_list[prcnt].rc      = 0.0;
                    prcnt++;
                }
            }
        }
        if (prcnt > 0) {
            printf ("Add %d edges to LP\n", prcnt); fflush (stdout);
            rval = CCtsp_add_vars_to_lp (global_lp, pr_list, prcnt);
            CCcheck_rval (rval, "CCtsp_add_vars_to_lp failed");
        }
    }

    scnt = 0;
    for (i = 0; i < lp->cuts.cutcount; i++) {
        if (lp->cuts.cuts[i].cliquecount == 1 &&
            lp->cuts.cuts[i].dominocount == 0) {
            scnt++;
            if (global_lp) {
                rval = CCtsp_lpcut_to_lpcut_in (&(lp->cuts),
                                        &(lp->cuts.cuts[i]), &old);
                CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
                if (ndepot > 0) {
                    rval = check_cut_for_depots (&old, ncount, &yesno);
                    CCcheck_rval (rval, "check_cut_for_depots failed");
                } else {
                    yesno = 1;
                }
                if (yesno) {
                    rval = CCtsp_copy_lpcut_in_global (&old, &new, local_names,
                                   local_perm, global_invperm, global_ncount);
                    CCcheck_rval (rval, "CCtsp_copy_lpcut_in_global failed");
                    fflush (stdout);
                    rval = CCtsp_add_cut (global_lp, &new, &cr, 1);
                    CCcheck_rval (rval, "CCtsp_add_cut failed");
                    CCtsp_free_lpcut_in (&new);
                }
                CCtsp_free_lpcut_in (&old);
            }
        } else {
            rval = CCtsp_lpcut_to_lpcut_in (&(lp->cuts), &(lp->cuts.cuts[i]),
                   &old);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            if (ndepot > 0) {
                rval = check_cut_for_depots (&old, ncount, &yesno);
                CCcheck_rval (rval, "check_cut_for_depots failed");
            } else {
                yesno = 1;
            }

            if (yesno) {
                rval = CCtsp_copy_lpcut_in_global (&old, &new, local_names,
                                   local_perm, global_invperm, global_ncount);
                CCcheck_rval (rval, "CCtsp_copy_lpcut_in_global failed");
                if (CCtsp_hypergraph_type_in (&new)) {
                    rval = CCtsp_add_to_cutpool_lpcut_in (pool, &new);
                    CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in");
                } else if (new.dominocount >= 0) {
                    rval = CCtsp_add_to_cutpool_lpcut_in (dompool, &new);
                    CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in");
                }
                if (global_lp) {
                    rval = CCtsp_add_cut (global_lp, &new, &cr, 1);
                    CCcheck_rval (rval, "CCtsp_add_cut failed");
                }
                CCtsp_free_lpcut_in (&new);
            }
            CCtsp_free_lpcut_in (&old);
        }
    }
    printf ("Cuts: %d (with %d subtours)\n", lp->cuts.cutcount, scnt);
    fflush (stdout);
        
    CCtsp_free_tsp_lp_struct (&lp);
    lp = (CCtsp_lp *) NULL;
    *cnt = scnt;

    if (global_lp) {
        if (cr.rowcnt > 0) {
            printf ("Adding %d cuts to global lp\n", cr.rowcnt);
            fflush (stdout);
            rval = CCtsp_add_multiple_rows (global_lp, &cr);
            CCcheck_rval (rval, "CCtsp_add_multiple_rows failed");
        }
    }

CLEANUP:
    CC_IFFREE (local_perm, int);
    CC_IFFREE (local_names, int);
    CC_IFFREE (pr_list, CCtsp_predge);
    CCutil_freedatagroup (&local_dat);
    CCtsp_free_lprow (&cr);
    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int check_cut_for_depots (CCtsp_lpcut_in *c, int ncount, int *yesno)
{
    int i, j, k, acount, rval = 0;
    int *arr = (int *) NULL;
    CCtsp_lpclique *q;
    
    *yesno = 0;

    for (i = 0; i < c->cliquecount; i++) {
        q = &(c->cliques[i]);
        rval = CCtsp_clique_to_array (q, &arr, &acount);
        CCcheck_rval (rval, "CCtsp_clique_to_array failed");
        for (j = 0; j < acount; j++) {
            if (arr[j] >= ncount) {
                /* printf ("X"); */
                goto CLEANUP;
            }
        }
        CC_IFFREE (arr, int);
    }

    for (i = 0; i < c->dominocount; i++) {
        for (k = 0; k < 2; k++) {
            q = &(c->dominos[i].sets[k]);
            rval = CCtsp_clique_to_array (q, &arr, &acount);
            CCcheck_rval (rval, "CCtsp_clique_to_array failed");
            for (j = 0; j < acount; j++) {
                if (arr[j] >= ncount) {
                    goto CLEANUP;
                }
            }
            CC_IFFREE (arr, int);
        }
    }

    *yesno = 1;

CLEANUP:
    CC_IFFREE (arr, int);
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

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval;

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
    return rval;
}

static int lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x)
{
    int rval = 0;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, xcount,
                     xlist, x, (double **) NULL, (double **) NULL,
                     (double **) NULL);
    if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");
    return rval;
}

static int lp_solve (CCtsp_lp *lp)
{
    int rval = 0, infeasible = 0;
    double lpzeit, lprzeit;
    
    lpzeit = CCutil_zeit ();
    lprzeit = CCutil_real_zeit ();

    /* printf ("Call dual simplex ....\n"); fflush (stdout); */
    rval = CClp_opt (lp->lp, CClp_METHOD_BARRIER, &infeasible); 
    /* rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible); */
    CCcheck_rval (rval, "CClp_opt failed");
    if (infeasible) {
        fprintf (stderr, "New LP is infeasible\n"); 
        rval = 1; goto CLEANUP;
    }
    printf ("LP Time: %.2f seconds\n", CCutil_zeit () - lpzeit);
    printf ("LP Wall Cloc Time: %.2f\n", CCutil_real_zeit () - lprzeit);
    fflush (stdout);
    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed");

CLEANUP:
    return rval;
}

#define SUB_ROUNDS 1

static int add_subtours_to_global_lp (CCtsp_lp *lp)
{
    int rval = 0, sround, ncount = lp->graph.ncount;
    int xcount, con_cutcount, blob_cutcount, seg_cutcount, ex_cutcount;
    int *xlist = (int *) NULL;
    double z, val, maxviol = 0.0, *x = (double *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_lprow cr;

    CCtsp_init_lprow (&cr);

    rval = lp_value (lp, &val);
    CCcheck_rval (rval, "lp_value");

    for (sround = 0; sround < SUB_ROUNDS; sround++) {
        CC_IFFREE (xlist, int);
        CC_IFFREE (x, double);
        rval = lp_x (lp, &xcount, &xlist, &x);
        CCcheck_rval (rval, "lp_x");

        rval = CCtsp_connect_cuts (&cuts, &con_cutcount, ncount, xcount, xlist,
                                   x, &z);
        CCcheck_rval (rval, "CCtsp_connect_cuts failed");
        printf ("Found %d connect cuts in %.2f seconds\n", con_cutcount, z);
        fflush (stdout);
        rval = add_cuts (con_cutcount, &cuts, lp, &cr);
        CCcheck_rval (rval, "add_cuts failed");

        CCutil_start_timer (&lp->stats.cuts_connect);
        rval = get_blobcuts (&cuts, &blob_cutcount, ncount, xcount, xlist, x);
        CCcheck_rval (rval, "get_blobcuts failed");
        z = CCutil_stop_timer (&lp->stats.cuts_connect, 1);
        printf ("Found %d blob cuts in %.2f seconds\n", blob_cutcount, z);
        fflush (stdout);
        rval = add_cuts (blob_cutcount, &cuts, lp, &cr);
        CCcheck_rval (rval, "add_cuts failed");

        rval = CCtsp_segment_cuts (&cuts, &seg_cutcount, ncount, xcount, xlist,
                                   x, &z);
        CCcheck_rval (rval, "CCtsp_segment_cuts failed");
        printf ("Found %d segment cuts in %.2f seconds\n", seg_cutcount, z);
        fflush (stdout);
        rval = add_cuts (seg_cutcount, &cuts, lp, &cr);
        CCcheck_rval (rval, "add_cuts failed");

        if (con_cutcount == 0 && 0) {
            CCutil_start_timer (&lp->stats.cuts_exactsubtour);
            rval = CCtsp_exact_subtours (&cuts, &ex_cutcount, ncount, xcount,
                                      xlist, x, (double *) NULL, &z, &maxviol);
            CCcheck_rval (rval, "CCtsp_exact_subtours failed");
            z = CCutil_stop_timer (&lp->stats.cuts_exactsubtour, 1);
            printf ("Found %d exact subtour cuts in %.2f seconds\n",
                               ex_cutcount, z);
            fflush (stdout);
            rval = add_cuts (ex_cutcount, &cuts, lp, &cr);
            CCcheck_rval (rval, "add_cuts failed");
        }

        if (cr.rowcnt > 0) {
            printf ("Adding %d subtours to global lp\n", cr.rowcnt);
            fflush (stdout);
            rval = CCtsp_add_multiple_rows (lp, &cr);
            CCcheck_rval (rval, "CCtsp_add_multiple_rows failed");
            rval = lp_solve (lp);
            CCcheck_rval (rval, "solve_lp failed");
            rval = lp_value (lp, &val);
            CCcheck_rval (rval, "lp_value");
            printf ("Sparse Global LP Value: %.2f\n", val);
            rval = CCtsp_write_probfile_sav (lp);
            CCcheck_rval (rval, "CCtsp_write_probfile_sav failed");
            CCtsp_free_lprow (&cr);
            CCtsp_init_lprow (&cr);
        }
    }

CLEANUP:
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    return rval;
}

static int price_edges_in_global_lp (CCtsp_lp *lp, CCrandstate *rstate)
{
    double penalty, val, rval = 0;
    int nadded;
    CCtsp_edgegenerator eg;

    rval = CCtsp_init_edgegenerator (&eg, lp->graph.ncount, lp->dat,
                (CCtsp_genadj *) NULL, CCtsp_PRICE_COMPLETE_GRAPH, 0,
                rstate);
    CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");

    rval = CCtsp_addbad_variables (lp, &eg, &penalty, &nadded,
            -CCutil_MAXDOUBLE, CCutil_MAXDOUBLE, 0, (int *) NULL, 0, rstate);
    CCcheck_rval (rval, "CCtsp_addbad_variables failed");

    CCtsp_free_edgegenerator (&eg);
    printf ("%d edges added, penalty %f\n", nadded, penalty);
    fflush (stdout);
    rval = lp_value (lp, &val);
    CCcheck_rval (rval, "lp_value");
    printf ("Global LP Bound: %f\n", val + penalty);
    fflush (stdout);

CLEANUP:
    return rval;
}

static int add_cuts (int count, CCtsp_lpcut_in **cuts, CCtsp_lp *lp,
        CCtsp_lprow *cr)
{
    int rval = 0;
    CCtsp_lpcut_in *c, *cnext;

    if (count) {
        for (c = *cuts; c; c = cnext) {
            cnext = c->next;
            rval = CCtsp_add_cut (lp, c, cr, 1);
            CCcheck_rval (rval, "CCtsp_add_cut failed");
            CCtsp_free_lpcut_in (c);
        }
    }
    *cuts = (CCtsp_lpcut_in *) NULL;

CLEANUP:
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm, boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "B:D:I:LM:tT:N:is:S", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 's': seed = atoi(boptarg); break;
        case 'B': borderfname = boptarg; break;
        case 'D': tspfname = boptarg; break;
        case 'I': indexfname = boptarg; break;
        case 'L': create_lp = 1; break;
        case 'S': add_subtours = 1; break;
        case 'M': masterfname  = boptarg; break;
        case 't': tsplib_tour = 0; break;
        case 'T': tourfname  = boptarg; break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: innorm = CC_MAXNORM; break;
            case 1: innorm = CC_MANNORM; break;
            case 2: innorm = CC_EUCLIDEAN; break;
            case 3: innorm = CC_EUCLIDEAN_3D; break;
            case 17: innorm = CC_GEOM; break;
            case 18: innorm = CC_EUCLIDEAN_CEIL; break;
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
    fprintf (stderr, "   -B f  border index file (add cuts from borders)\n");
    fprintf (stderr, "   -D f  specify a tsp file (if no master file)\n");
    fprintf (stderr, "   -I f  index file\n");
    fprintf (stderr, "   -L    build a root LP\n");
    fprintf (stderr, "   -M f  specify a master file\n");
    fprintf (stderr, "   -s #  seed\n");
    fprintf (stderr, "   -S    add subtours to root LP (need -L)\n");
    fprintf (stderr, "   -T f  specify a tour file (if no master file)\n");
    fprintf (stderr, "   -t    tour is in concorde format\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 17=GEOM, 18=JOHNSON\n");
}


/**************************************************************/
/********************                     *********************/
/********************  Old pancake.[c,h]  *********************/
/********************                     *********************/
/**************************************************************/


typedef struct pannode {
        struct panedge **edgelist;
        struct panedge **goodedge;
        int degree;
        struct vaseknode *vptr;
} pannode;

typedef struct panedge {
        double panweight;
        pannode *ends[2];
        int elim;
        int tag;
        struct vaseknode *a[2];
        struct vaseknode *roof;
        struct vaseknode *top;
        struct panedge *next;
        double rc;
        double realrc;
} panedge;

typedef struct vaseknode {
        struct vaseknode *parent;
        struct vaseknode *child;
        struct vaseknode *sibling;
        struct vaseknode *anc;
        struct vaseknode *ptr;
        struct vaseknode *qtr;
        struct vaseknode *listpointer;
        int b;
        int d;
        int n;
        int tag;
        int fringe;
        double y;
        double w;
        double mult;
        struct panedge *tree;
        struct panedge *junk;
        struct triomino *adj;
        struct triomino *scan;
} vaseknode;

typedef struct triomino {
        panedge *edge;
        vaseknode *end;
        struct triomino *next;
} triomino;

typedef struct exactsub_param {
    int             nodecount;
    int             cutcount;
    CCtsp_lpcut_in *cuts;
} exactsub_param;

#define VNODEALLOC(vrequest) {                                          \
        if (vnodestack == NULL) {                                       \
                printf ("Ran out of vnode supply\n");                   \
                exit (1);                                               \
        }                                                               \
        vrequest = vnodestack;                                          \
        vnodestack = vnodestack->ptr;                                   \
    }

#define VNODEFREE(vreturn) {                                            \
        vreturn->ptr = vnodestack;                                      \
        vnodestack = vreturn;                                           \
    }

#define TRIALLOC(trequest) {                                            \
        trequest = tristack;                                            \
        tristack = tristack->next;                                      \
    }

#define TRIFREE(treturn) {                                              \
        treturn->next = tristack;                                       \
        tristack = treturn;                                             \
    }

#define MAGICNODE 0
#define PANEPSILON 0.0001
#define PANFALSE 0
#define PANPOSITIVE 1
#define PANFEW 1
#define PANSHORT 1000
#define PANBIGNEG -10000000000.0
#define PANSWAP(x,y,temp) {temp = x; x = y; y = temp;}

static int panalloc (int ecount, int *elist);
static int buildpanadjlist (int ecount, int *elist);
static int pancakex (int ecount, int *elist, double *x);
static void panfree (void);

static int initpancake (int ecount);
static int buildfirsttree (int ecount);
static int decompositiontree (int ecount);
static int initdecompositiontree (int ecount);
static void drop (panedge *e, vaseknode *x);
static void throw (vaseknode *x, vaseknode *y, panedge *e);
static vaseknode *anc (vaseknode *v);
static void trickledown (int i);
static vaseknode *newcomp (vaseknode *v, double w);
static void hookup (vaseknode *parent, vaseknode *child);
static void distribute (void);
static void initdistribute (void);
static void split (vaseknode *a);
static void bruteforce (panedge *e);
static void update (panedge *e);
static void dealwith (panedge *e, vaseknode **pa);
static void attach (panedge *e);
static void magicrc (void);
static double min2 (panedge **elist);
static double findbound (void);
static double blnode (vaseknode *v, int *hit, int *cutarray, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param);
static void labeler (vaseknode *p, int *cutcount, int *cutarray);
static void freepancake (void);
static int pancakemain (int ecount);
static int blobsviolated (int *hit, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param);
static int add_exact (double val, int count, int *cutarray, void *pass_param);

static int vstep;
static vaseknode *vroot;
static vaseknode *vpannodes, *vnodehit;
static vaseknode *head, *tail;
static vaseknode *vnodestack;
static triomino *trisupply, *tristack;
static panedge *work;
static panedge **vheap;
static int vheapend;
static int vcomponentcount;
static panedge **panedgespace, *panedgelist;
static pannode *pannodelist;
static int nnodes;

static int get_blobcuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
            int ecount, int *elist, double *x)
{
    int rval;
    int hit;
    exactsub_param p;

    *cutcount = 0;

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    rval = blobcuts (ncount, ecount, elist, x, &hit, add_exact, (void *) &p);
    if (rval) {
        fprintf (stderr, "blobcuts failed\n"); goto CLEANUP;
    }

    *cutcount = p.cutcount;
    *cuts = p.cuts;

CLEANUP:

    return rval;
}

static int add_exact (double val, int count, int *cutarray, void *pass_param)
{
    int rval = 0;
    exactsub_param *p = (exactsub_param *) pass_param;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    if (count >= p->nodecount) goto CLEANUP;

    if (val > 2.0) {
        printf ("Warning: Cut of value %f in add_exact\n", val);
        fflush (stdout);
        goto CLEANUP;
    }

    rval = CCtsp_array_to_subtour (&c, cutarray, count, p->nodecount);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_subtour failed\n");
        rval = 1; goto CLEANUP;
    }
    c->next = p->cuts;
    p->cuts = c;
    p->cutcount++;

CLEANUP:
    return rval;
}

static int blobcuts (int ncount, int ecount, int *elist, double *x,
        int *hit, int (*doit_fn) (double, int, int *, void *), void *pass_param)
{
    int rval = 0;

    /* printf ("blobcuts ...\n"); fflush (stdout); */

    nnodes = ncount;

    pancakex (ecount, elist, x);

    rval = blobsviolated (hit, x, doit_fn, pass_param);
    if (rval) {
        fprintf (stderr, "blobsviolated failed\n");  goto CLEANUP;
    }
    freepancake ();
    panfree ();

CLEANUP:
    return rval;
}

static int pancakex (int ecount, int *elist, double *x)
{
    int i, rval = 0;

    /* printf ("pancakex ...\n"); fflush (stdout); */

    rval = panalloc (ecount, elist);
    if (rval) {
        fprintf (stderr, "panalloc failed\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < ecount; i++) {
        panedgelist[i].panweight = -x[i];
    }

    pancakemain (ecount);

CLEANUP:
    return rval;
}

static int panalloc (int ecount, int *elist)
{
    int rval = 0, i;
    panedge *pa;

    /* printf ("panalloc (%d)  ...\n", ecount); fflush (stdout); */

    pannodelist = (pannode *) malloc (nnodes * sizeof (pannode));
    panedgelist = (panedge *) malloc (ecount * sizeof (panedge));
    if (!pannodelist || !panedgelist) {
        fprintf (stderr, "out of memory in panalloc\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, pa = panedgelist; i < ecount; i++, pa++) {
        pa->ends[0] = pannodelist + elist[2*i];
        pa->ends[1] = pannodelist + elist[2*i+1];
    }

    rval = buildpanadjlist (ecount, elist);
    if (rval) {
        fprintf (stderr, "buildpanadjlist failed\n"); goto CLEANUP;
    }

CLEANUP:
    return rval;
}

static void panfree (void)
{
    if (pannodelist) {
        free (pannodelist);
        pannodelist = NULL;
    }
    if (panedgelist) {
        free (panedgelist);
        panedgelist = NULL;
    }
    if (panedgespace) {
        free (panedgespace);
        panedgespace = NULL;
    }
}

static int buildpanadjlist (int ecount, int *elist)
{
    int rval = 0, i, *degrees = (int *) NULL, *pint;
    pannode *pv;
    panedge *pa, **edgespace;

    /* printf ("buildpanadjlist ...\n"); fflush (stdout); */

    edgespace = (panedge **) malloc ((nnodes+(2*ecount)) * sizeof (panedge *));
    degrees = (int *) malloc (nnodes * sizeof (int));
    if (!edgespace || !degrees) {
        fprintf (stderr, "out of memory in buildpanadjlist\n");
        rval = 1;  goto CLEANUP;
    }

    panedgespace = edgespace;
    for (i = 0, pint = degrees; i < nnodes; i++) {
        *pint++ = 0;
    }

    for (i = 0; i < ecount; i++) {
        degrees[elist[2*i]]++;
        degrees[elist[2*i+1]]++;
    }

    for (i = 0, pint = degrees, pv = pannodelist; i < nnodes;
         i++, pint++, pv++) {
        if (*pint) {
            pv->edgelist = pv->goodedge = edgespace;
            edgespace += *pint + 1;
        } else
            pv->edgelist = pv->goodedge = (panedge **) NULL;
    }

    for (i = 0, pa = panedgelist; i < ecount; i++, pa++) {
        *(pa->ends[0]->goodedge++) = pa;
        *(pa->ends[1]->goodedge++) = pa;
    }

    for (i = 0, pv = pannodelist; i < nnodes; i++, pv++) {
        *(pv->goodedge) = (panedge *) NULL;
    }

CLEANUP:
    if (degrees) { free (degrees); }
    return rval;
}

/**************************************************************************/
/*                         Old dtree.c                                    */
/**************************************************************************/

static int pancakemain (int ecount)
{
    int rval = 0;

    rval = initpancake (ecount);
    if (rval) {
        fprintf (stderr, "initpancake failed\n"); goto CLEANUP;
    }
    rval = buildfirsttree (ecount);
    if (rval) {
        fprintf (stderr, "buildfirsttree failed\n"); goto CLEANUP;
    }

    findbound ();

CLEANUP:
    return rval;
}

static int initpancake (int ecount)
{
    int rval = 0, i;
    pannode *pn;
    vaseknode *pv;
    triomino *pt;

    vpannodes = (vaseknode *) malloc ((2*nnodes-1) * sizeof (vaseknode));
    if (!vpannodes) {
        fprintf (stderr, "out of memory in initpancake\n");
        rval = 1; goto CLEANUP;
    } 

    for (i = nnodes, pn = pannodelist, pv = vpannodes; i; i--, pn++, pv++) {
        pv->parent = (vaseknode *) NULL;
        pv->child = (vaseknode *) NULL;
        pv->sibling = (vaseknode *) NULL;
        pv->adj = (triomino *) NULL;
        pv->n = 0;
        pv->b = 1;
        pv->anc = pv;
        pv->junk = (panedge *) NULL;
        pv->w = PANBIGNEG;
        pv->y = 0.0;
        pn->vptr = pv;
    }

    vnodestack = vpannodes + nnodes;
    for (i = nnodes - 2, pv = vnodestack; i; i--, pv++)
        pv->ptr = pv + 1;
    pv->ptr = (vaseknode *) NULL;

    trisupply = (triomino *) malloc ((2*ecount) * sizeof (triomino));
    if (!trisupply) {
        fprintf (stderr, "out of memory for trisupply\n");
        rval = 1;  goto CLEANUP;
    }
    tristack = trisupply;
    for (i = 2 * ecount - 1, pt = tristack; i; i--, pt++)
        pt->next = pt + 1;

    VNODEALLOC (head);

CLEANUP:
    return rval;
}

static void freepancake (void)
{
    if (vpannodes) {
        free (vpannodes);
        vpannodes = NULL;
    }
    if (trisupply) {
        free (trisupply);
        trisupply = NULL;
    }
    if (vheap) {
        free (vheap);
        vheap = NULL;
    }
}

static int buildfirsttree (int ecount)
{
    int rval = 0;
    vaseknode *p, *q;
    double parw;

    rval = decompositiontree (ecount);
    if (rval) {
        fprintf (stderr, "decompositiontree failed\n"); goto CLEANUP;
    }
    distribute ();


    head->ptr = tail = vroot;
    p = head;
    do {
        p = p->ptr;
        parw = p->w;
        for (q = p->child; q != (vaseknode *) NULL; q = q->sibling) {
            q->mult = parw - q->w;
            if (q->child != (vaseknode *) NULL) {
                tail->ptr = q;
                tail = q;
            }
        }
    } while (p != tail);
    vroot->mult = -vroot->w;

    magicrc ();

CLEANUP:
    return rval;
}

static int decompositiontree (int ecount)
{
    int rval = 0, i;
    panedge *e;
    double w, ub;
    vaseknode *x, *y;

    /* printf ("decompositiontree ...\n"); fflush (stdout); */

    rval = initdecompositiontree (ecount);
    if (rval) {
        fprintf (stderr, "initdecompositiontree failed\n"); goto CLEANUP;
    }

    for (vcomponentcount = nnodes - 2; vcomponentcount;) {
        for (;;) {
            e = *vheap;
            *vheap = vheap[vheapend--];
            trickledown (0);
            x = anc (e->ends[0]->vptr);
            y = anc (e->ends[1]->vptr);
            if (x != y)
                break;
            drop (e, x);
        }

        w = e->panweight;
        throw (x, y, e);

        ub = w + 0.01;

        while (vheapend >= 0 && (e = *vheap)->panweight < ub) {
            *vheap = vheap[vheapend--];
            trickledown (0);
            x = anc (e->ends[0]->vptr);
            y = anc (e->ends[1]->vptr);
            if (x != y)
                throw (x, y, e);
            else
                drop (e, x);
        }

        for (; vnodehit != (vaseknode *) NULL; vnodehit = vnodehit->ptr) {
            if (vnodehit->n)
                vroot = newcomp (vnodehit, w);
        }
    }
    i = vheapend + 1;
    while (i) {
	drop (vheap[--i], vroot);
    }

CLEANUP:
    return rval;
}

static int initdecompositiontree (int ecount)
{
    int rval = 0, i;
    panedge *pe, **ph;
    pannode *pm;

    vnodehit = (vaseknode *) NULL;
    work = (panedge *) NULL;

    for (i = ecount, pe = panedgelist; i; i--, pe++)
        pe->tag = PANFALSE;

    vheap = (panedge **) malloc (ecount * sizeof (panedge *));
    if (!vheap) {
        fprintf (stderr, "out of memory for vheap\n");
        rval = 1;  goto CLEANUP;
    }

    pm = pannodelist + MAGICNODE;
    for (i = ecount, pe = panedgelist, ph = vheap; i; i--, pe++)
        if (pe->ends[0] != pm && pe->ends[1] != pm)
            *(ph++) = pe;
        else
            pe->rc = pe->panweight;

    vheapend = (ph - vheap) - 1;
    for (i = vheapend / 2; i >= 0; i--) {
        trickledown (i);
    }

CLEANUP:
    return rval;
}

static void drop (panedge *e, vaseknode *x)
{
    e->a[0] = e->ends[0]->vptr;
    e->a[1] = e->ends[1]->vptr;
    e->roof = x;
    e->rc = PANPOSITIVE;
    e->next = work;
    work = e;
}

static void throw (vaseknode *x, vaseknode *y, panedge *e)
{
    e->a[0] = x;
    e->a[1] = y;

    e->rc = 0.0;
    attach (e);

    if (!(x->n)) {
        x->n = 1;
        x->ptr = vnodehit;
        vnodehit = x;
    }
    if (!(y->n)) {
        y->n = 1;
        y->ptr = vnodehit;
        vnodehit = y;
    }
    e->next = work;
    work = e;
}

static vaseknode *anc (vaseknode *v)
{
    vaseknode *hand, *va;

    hand = v;
    while (hand != hand->anc)
	hand = hand->anc;

    va = hand;
    for (hand = v; hand != va; hand = hand->anc)
        hand->anc = va;

    return va;
}

static void trickledown (int i)
{
    panedge *memo;
    int k, minchild;

    memo = vheap[i];

    while ((k = (2 * i) + 2) <= vheapend) {
        minchild = (vheap[k - 1]->panweight <= vheap[k]->panweight ? k - 1 : k);

        if (memo->panweight > vheap[minchild]->panweight) {
            vheap[i] = vheap[minchild];
            i = minchild;
        } else {
            vheap[i] = memo;
            return;
        }
    }
    if (k - 1 == vheapend && memo->panweight > vheap[vheapend]->panweight) {
        vheap[i] = vheap[vheapend];
        i = vheapend;
    }
    vheap[i] = memo;
}

static vaseknode *newcomp (vaseknode *v, double w)
{
    vaseknode *new, *stack;
    triomino *t;

    VNODEALLOC (new);
    new->parent = (vaseknode *) NULL;
    new->child = (vaseknode *) NULL;
    new->sibling = (vaseknode *) NULL;
    new->anc = new;
    new->n = 0;
    new->b = 0;
    new->w = w;
    new->adj = (triomino *) NULL;
    new->junk = (panedge *) NULL;
    new->tag = PANFALSE;

    hookup (new, v);
    v->qtr = (vaseknode *) NULL;

    do {
        stack = v->qtr;
        for (t = v->adj; t != (triomino *) NULL; t = t->next) {
            t->edge->roof = new;
            v = t->end;
            if (v->n) {
                v->qtr = stack;
                stack = v;
                hookup (new, v);
                vcomponentcount--;
            }
        }
    } while ((v = stack) != (vaseknode *) NULL);

    return new;
}

static void hookup (vaseknode *parent, vaseknode *child)
{
    child->n = 0;
    child->parent = parent;
    child->anc = parent;
    child->sibling = parent->child;
    parent->child = child;
}

static void distribute (void)
{
    vaseknode *active;
    panedge *e, *f;

    /* printf ("distribute ...\n"); fflush (stdout); */

    initdistribute ();
    active = (vaseknode *) NULL;
    vroot->n = 0;

    for (e = work, work = (panedge *) NULL; e != (panedge *) NULL; e = f) {
        f = e->next;
        e->top = vroot;
        dealwith (e, &active);
    }

    while (work) {
        for (; active != (vaseknode *) NULL; active = active->ptr)
            if (active->n < PANFEW)
                active->n = -1;
            else {
                active->n = 0;
                split (active);
            }
        for (e = work, work = (panedge *) NULL; e != (panedge *) NULL;
                                                                    e = f) {
            f = e->next;
            if (e->top->n >= 0) {
                update (e);
                dealwith (e, &active);
            } else
                bruteforce (e);
        }
        vstep = vstep / 2;
    }
}

static void initdistribute (void)
{
    vaseknode *stack, *finger, *x;
    int maxd, twice, d;

    maxd = 0;

    vroot->d = 0;
    vroot->ptr = (vaseknode *) NULL;

    for (finger = vroot; finger != (vaseknode *) NULL; finger = stack) {
        stack = finger->ptr;
        finger->anc = vroot;
        if ((x = finger->child) != (vaseknode *) NULL) {
            d = finger->d + 1;
            do {
                x->d = d;
                x->ptr = stack;
                stack = x;
                x = x->sibling;
            } while (x != (vaseknode *) NULL);
            if (d > maxd)
                maxd = d;
        }
    }
    vstep = 1;
    twice = 2;
    while (twice < maxd) {
	vstep = twice;
	twice = vstep + vstep;
    }
}

static void split (vaseknode *a)
{
    int mid, bot;
    vaseknode *stack, *hand, *foot, *memo, *x;

    mid = vstep + a->d;
    bot = vstep + mid;

    a->qtr = (vaseknode *) NULL;
    for (hand = a; hand != (vaseknode *) NULL; hand = stack) {
        stack = hand->qtr;
        if (hand->d == mid) {
            memo = hand->qtr;
            hand->qtr = (vaseknode *) NULL;
            for (foot = hand; foot != (vaseknode *) NULL; foot = stack) {
                stack = foot->qtr;
                foot->anc = hand;
                if (foot->d != bot) {
                    for (x = foot->child; x != (vaseknode *) NULL;
                         x = x->sibling) {
                        x->qtr = stack;
                        stack = x;
                    }
                }
            }
            hand->qtr = memo;
        } else
            for (x = hand->child; x != (vaseknode *) NULL; x = x->sibling) {
                x->qtr = stack;
                stack = x;
            }
    }
}

static void bruteforce (panedge *e)
{
    vaseknode *x, *y, *nx, *ny;
    int dx, dy;

    x = e->a[0];
    y = e->a[1];

    if (x == y) {
        printf ("Tough luck Pal 1.\n");
        exit (1);
    }
    dx = x->d;
    dy = y->d;

    while (dx > dy) {
        x = x->parent;
        dx--;
    }
    if (x == y) {
        printf ("Tough luck Pal 2.\n");
        exit (1);
    }
    while (dy > dx) {
        y = y->parent;
        dy--;
    }
    if (x == y) {
        printf ("Tough luck Pal 3.\n");
        exit (1);
    }
    nx = x->parent;
    ny = y->parent;
    while (nx != ny) {
	x = nx;
	y = ny;
	nx = x->parent;
	ny = y->parent;
    }

    e->a[0] = x;
    e->a[1] = y;

    e->roof = nx;
    /* if (e->rc > 0.0) { e->next = nx->junk; nx->junk = e; } else { e->tag =
     * PANFALSE; attach (e); } */
    e->next = nx->junk;
    nx->junk = e;
}

static void update (panedge *e)
{
    vaseknode *x, *y, *v;

    x = e->a[0]->anc;
    y = e->a[1]->anc;
    v = e->top;

    if (x == v) {
        if (y != v)
            e->a[1] = y;
    } else if (y == v)
        e->a[0] = x;
    else if (x != y) {
        e->a[0] = x;
        e->a[1] = y;
    } else {
        e->top = x;
        if (x->d > e->roof->d)
            e->roof = x;
    }
}

static void dealwith (panedge *e, vaseknode **pa)
{
    if ((e->roof->d) - (e->a[0]->d) < PANSHORT &&
        (e->roof->d) - (e->a[1]->d) < PANSHORT) {
        bruteforce (e);
    } else {
        e->next = work;
        work = e;
        if (!e->top->n) {
            e->top->ptr = *pa;
            *pa = e->top;
        }
        (e->top->n)++;
    }
}

static void attach (panedge *e)
{
    triomino *cell;

    TRIALLOC (cell);
    cell->edge = e;
    cell->end = e->a[1];
    cell->next = e->a[0]->adj;
    e->a[0]->adj = cell;

    TRIALLOC (cell);
    cell->edge = e;
    cell->end = e->a[0];
    cell->next = e->a[1]->adj;
    e->a[1]->adj = cell;
}

static void magicrc (void)
{
    double a;
    panedge **pee, *pe;

    a = min2 ((pannodelist + MAGICNODE)->edgelist);

    for (pee = (pannodelist + MAGICNODE)->edgelist;
                (pe = *pee) != (panedge *) NULL; pee++)
        pe->rc -= a;

    (pannodelist + MAGICNODE)->vptr->y += a;
}

static double min2 (panedge **elist)
{
    double minweight, minweight2, td;
    panedge *e;

    if (elist == (panedge **) NULL ||
     elist[0] == (panedge *) NULL || elist[1] == (panedge *) NULL) {
        fprintf (stderr, "Vertex has degree < two\n");
        exit (1);
    }
    minweight = elist[0]->rc;
    minweight2 = elist[1]->rc;
    if (minweight > minweight2) {
        PANSWAP (minweight, minweight2, td);
    }
    for (elist += 2; (e = *elist) != (panedge *) NULL; elist++) {
        if (e->rc < minweight2) {
            minweight2 = e->rc;
            if (minweight > minweight2) {
                PANSWAP (minweight, minweight2, td);
            }
        }
    }
    return minweight2;
}

static double findbound (void)
{
    vaseknode *p, *q, *stack;
    triomino *tri;
    panedge **pee, *pe;
    double tree_bound = 0.0;
    double star_bound = 0.0;
    double edge_bound = 0.0;

    vroot->ptr = (vaseknode *) NULL;
    vroot->n = 1;
    vroot->b = 0;


    for (p = stack = vroot; p; p = stack) {
        if (p->n) {
            p->n = 0;
            q = p->child;
            if (q)
                for (; q; q = q->sibling) {
                    q->ptr = stack;
                    stack = q;
                    q->n = 1;
                    q->b = 0;
                }
            else {
                stack = p->ptr;
                (p->parent->b)++;
                star_bound += p->y;
            }
            for (tri = p->adj; tri; tri = tri->next)
		edge_bound += tri->edge->rc;
        } else {
            stack = p->ptr;
            if (stack)
                (p->parent->b) += p->b;
            tree_bound -= (p->mult) * ((p->b) - 1);
        }
    }
    star_bound *= 2.0;
    edge_bound /= 2.0;

    for (pee = ((pannodelist + MAGICNODE)->edgelist);
                 (pe = *pee) != (panedge *) NULL; pee++)
        if (pe->rc < 0.0)
            edge_bound += pe->rc;
    star_bound += (pannodelist + MAGICNODE)->vptr->y * 2;

    return tree_bound + star_bound + edge_bound;
}

static int blobsviolated (int *hit, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param)
{
    int *cutarray = (int *) NULL;
    int rval = 0;

    /* printf ("blobsviolated ....\n"); fflush (stdout); */

    if (hit) *hit = 0;

    cutarray = (int *) malloc (nnodes * sizeof (int));
    if (!cutarray) {
        fprintf (stderr, "out of memory for cutarray\n");
        rval = 1; goto CLEANUP;
    }

    blnode (vroot, hit, cutarray, x, doit_fn, pass_param);

CLEANUP:
  
    if (cutarray) free (cutarray);
    return rval;
}

#define CUTTOLERANCE 0.01

static double blnode (vaseknode *v, int *hit, int *cutarray, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param)
{
    double w = 0.0;
    double t;
    panedge *e;
    vaseknode *c;
    int rval = 0;
    int cutcount;

    if (!v->child)
        return 0.0;
    else {
        for (e = v->junk; e; e = e->next)
            w += (x[e - panedgelist]);
        for (c = v->child; c; c = c->sibling)
            w += blnode (c, hit, cutarray, x, doit_fn, pass_param);
        t = v->b;
        if (w > t - 1.0 + CUTTOLERANCE) {
            cutcount = 0;
            labeler (v, &cutcount, cutarray);
            if (doit_fn) {
                rval = doit_fn (1.9, cutcount, cutarray, pass_param);
                if (rval) {
                    fprintf (stderr, "doit_fn failed\n");
                    exit (1);
                }
            }
            (*hit)++;
        }
        return w;
    }
}

static void labeler (vaseknode *p, int *cutcount, int *cutarray)
{
    vaseknode *c;

    if (!p->child) {
        cutarray[*cutcount] = p - vpannodes;
        (*cutcount)++;
    } else {
        for (c = p->child; c; c = c->sibling) {
            labeler (c, cutcount, cutarray);
        }
    }
}

