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
/*                TOOL TO GATHER STATISTICS ON CUTTING PLANES               */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: July 16, 2012                                                     */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

static char *lpfname    = (char *) NULL;
static char *masterfname = (char *) NULL;
static int show_dd_toothweights = 0;
static int show_dp_stats = 0;
static int seed = 0;

int
    main (int ac, char **av);

static int
    dp_stats (CCtsp_lp *lp, CCtsp_lpgraph *lg, double *x),
    dp_look (CCtsp_lpcut_in *c, CCtsp_lpgraph *lg, double *x, int *marks),
    ddecker_stats (CCtsp_lp *lp, CCtsp_lpgraph *lg, double *x, int ddcount,
        int *ddlist),
    ddecker_look (CCtsp_lpcut_in *c, CCtsp_lpgraph *lg, double *x),
    cross_weight (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *T,
        CCtsp_lpclique *H, double *r),
    lp_value (CCtsp_lp *lp, double *val),
    lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x),
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int i, ncount, rval = 0, yesno;
    int subtour_count = 0, comb_count = 0, dd_count = 0, dp_count = 0;
    int twop_count = 0, xcount;
    int *ptour = (int *) NULL, *ddlist = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_lpcut_in c;
    double val, szeit;
    double *x = (double *) NULL;
    int *xlist = (int *) NULL;
    CCtsp_lpgraph lg;

    CCtsp_init_lpgraph_struct (&lg);
    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed")

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!masterfname) {
        fprintf (stderr, "Must specify a master file\n");
        usage (av[0]);
        goto CLEANUP;
    }

    szeit = CCutil_zeit ();
    CCutil_sprand (seed, &rstate);

    rval = CCutil_getmaster (masterfname, &ncount, &dat, &ptour);
    CCcheck_rval (rval, "CCutil_getmaster failed")

    rval = CCtsp_init_lp (&lp, (char *) NULL, -1, lpfname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 0, &rstate);
    CCcheck_rval (rval, "CCtsp_init_lp failed");

    rval = lp_value (lp, &val);
    CCcheck_rval (rval, "lp_value failed");
    printf ("LP Value: %f\n", val); fflush (stdout);

    rval = lp_x (lp, &xcount, &xlist, &x);
    CCcheck_rval (rval, "lp_x failed");

    rval = CCtsp_build_lpgraph (&lg, ncount, xcount, xlist, (int *) NULL);
    CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
    rval = CCtsp_build_lpadj (&lg, 0, xcount);
    CCcheck_rval (rval, "CCtsp_build_lpadj failed");


    ddlist = CC_SAFE_MALLOC (lp->cuts.cutcount, int);
    CCcheck_NULL (ddlist, "out of memory for ddlist");

    for (i = 0; i < lp->cuts.cutcount; i++) {
        if (lp->cuts.cuts[i].branch) continue;
        if (lp->cuts.cuts[i].twodom_cliquecount > 0) {
            twop_count++;
        } else if (lp->cuts.cuts[i].dominocount > 0) {
            dp_count++;
        } else if (lp->cuts.cuts[i].cliquecount == 1) {
            subtour_count++; 
        } else if (lp->cuts.cuts[i].cliquecount > 3) {
            rval = CCtsp_lpcut_to_lpcut_in (&lp->cuts, &lp->cuts.cuts[i], &c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

            rval = CCtsp_test_pure_comb (ncount, &c, &yesno, (int *) NULL);
            CCcheck_rval (rval, "CCtsp_test_pure_comb failed");
            if (yesno) {
                comb_count++;
            } else {
                rval = CCtsp_test_pure_double_decker (&c, &yesno,
                              (int *) NULL, (int *) NULL, (int *) NULL,
                              (int **) NULL, 0);
                CCcheck_rval (rval, "CCtsp_test_pure_double_decker failed");
                if (yesno) {
                    ddlist[dd_count++] = i;
                }
            }
            CCtsp_free_lpcut_in (&c);
        }
    }

    printf ("Number of Cuts: %3d\n", lp->cuts.cutcount);
    printf ("Subtours:       %3d\n", subtour_count);
    printf ("Combs:          %3d\n", comb_count);
    printf ("DDeckers:       %3d\n", dd_count);
    if (dp_count > 0) printf ("DP-inequalties: %3d\n", dp_count);
    if (twop_count > 0) printf ("2P-inequalties: %3d\n", twop_count);
    printf ("\n");
    fflush (stdout);

    if (show_dd_toothweights) {
        rval = ddecker_stats (lp, &lg, x, dd_count, ddlist);
        CCcheck_rval (rval, "ddecker_stats failed");
    }

    if (show_dp_stats) {
        rval = dp_stats (lp, &lg, x);
        CCcheck_rval (rval, "dp_stats failed");
    }

CLEANUP:
    CC_IFFREE (ptour, int);
    CC_IFFREE (ddlist, int);
    CCutil_freedatagroup (&dat);
    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    CCtsp_free_lpgraph (&lg);

    return rval;
}

static int dp_stats (CCtsp_lp *lp, CCtsp_lpgraph *lg, double *x)
{
    int rval = 0, i;
    CCtsp_lpcut *cuts = lp->cuts.cuts;
    int cutcount = lp->cuts.cutcount;
    CCtsp_lpcut_in cin;
    int *marks = (int *) NULL;

    marks = CC_SAFE_MALLOC (lg->ncount, int);
    CCcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < lg->ncount; i++) marks[i] = 0;

    for (i = 0; i < cutcount; i++) {
        if (cuts[i].dominocount > 0) {
            rval = CCtsp_lpcut_to_lpcut_in (&lp->cuts, &cuts[i], &cin);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
            rval = dp_look (&cin, lg, x, marks);
            CCcheck_rval (rval, "dp_look failed");
            CCtsp_free_lpcut_in (&cin);
        }
    }

CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

static int dp_look (CCtsp_lpcut_in *c, CCtsp_lpgraph *lg, double *x, int *marks)
{
    int rval = 0, i, j, cnt, hcnt;
    CCtsp_lpdomino *dominos = c->dominos;

    CCtsp_mark_clique (&c->cliques[0], marks, 1);
    
    for (i = 0; i < c->dominocount; i++) {
        printf ("[");
        for (j = 0; j < 2; j++) {
            if (j == 1) printf (", ");
            CCtsp_clique_count (&dominos[i].sets[j], &cnt);
            CCtsp_clique_marked_count (&dominos[i].sets[j], marks, 1, &hcnt);
            printf ("%d(%d)", cnt, hcnt);
        }
        printf ("] ");
    }
    printf ("\n");

    CCtsp_mark_clique (&c->cliques[0], marks, 0);

CLEANUP:
    return rval;
}

static int ddecker_stats (CCtsp_lp *lp, CCtsp_lpgraph *lg, double *x,
        int ddcount, int *ddlist)
{
    int j, d, rval = 0;
    CCtsp_lpcut *cuts = lp->cuts.cuts;
    CCtsp_lpcut_in cin;

    for (j = 0; j < ddcount; j++) {
        d = ddlist[j];
        rval = CCtsp_lpcut_to_lpcut_in (&lp->cuts, &cuts[d], &cin);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");
        rval = ddecker_look (&cin, lg, x);
        CCcheck_rval (rval, "ddecker_look failed");
        CCtsp_free_lpcut_in (&cin);
    }

CLEANUP:
    return rval;
}

static int ddecker_look (CCtsp_lpcut_in *c, CCtsp_lpgraph *lg, double *x)
{
    int rval = 0, yesno, i, j, k;
    int hand1, hand2, dt_count;
    int *dt_list = (int *) NULL, *dt_ind = (int *) NULL;
    double *w = (double *) NULL, q;

    /* teeth with beta=2 appear twice in the tooth list */
    /* dt_list will give the teeth having beta = 2 */
    /* dt_ind is 0 if beta=1; 1 if beta=2; 2 if second copy of tooth */

    rval = CCtsp_test_pure_double_decker (c, &yesno, &hand1, &hand2, 
                                          &dt_count, &dt_list, 0);
    CCcheck_rval (rval, "CCtsp_test_pure_double_decker failed");
    if (!yesno) { rval = 1; goto CLEANUP; }

    dt_ind = CC_SAFE_MALLOC (c->cliquecount, int);
    CCcheck_rval (rval, "out of memory for dt_ind");
    for (i = 0; i < c->cliquecount; i++) dt_ind[i] = 0;
    for (i = 0; i < dt_count; i++) {
        k = dt_list[i];
        dt_ind[k] = 1;
        for (j = 0; j < c->cliquecount; j++) {
            if (j != k && j != hand1 && j != hand2) {
                CCtsp_clique_eq (&c->cliques[k], &c->cliques[j], &yesno);
                if (yesno) {
                    dt_ind[j] = 2; break;
                }
            }
        }
    }

    /* w(T) = x(\delta(T))+x(E(T\setminus H_1:H_1))+x(E(T\setminus H_2:H_2)) */

    w = CC_SAFE_MALLOC (c->cliquecount, double);
    CCcheck_NULL (w, "out of memory for w");
    for (i = 0; i < c->cliquecount; i++) w[i] = 0.0;

    /* Add x(\delta(T) */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != hand1 && i != hand2 && dt_ind[i] != 2) {
            rval = CCtsp_clique_delta (lg, x, &c->cliques[i], &q);
            CCcheck_rval (rval, "CCtsp_clique_delta failed");
            w[i] += q;
        }
    }

    /* Add x(E(T\setminus H_1:H_1)) */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != hand1 && i != hand2 && dt_ind[i] != 2) {
            rval = cross_weight (lg, x, &c->cliques[i], &c->cliques[hand1], &q);
            CCcheck_rval (rval, "cross_weight failed");
            w[i] += q;
        }
    }
         
    /* Add x(E(T\setminus H_2:H_2)) */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != hand1 && i != hand2 && dt_ind[i] != 2) {
            rval = cross_weight (lg, x, &c->cliques[i], &c->cliques[hand2], &q);
            CCcheck_rval (rval, "cross_weight failed");
            w[i] += q;
        }
    }

    printf ("DD: ");
    for (i = 0; i < c->cliquecount; i++) {
        if (i != hand1 && i != hand2 && dt_ind[i] != 2) {
            if (dt_ind[i] == 1) printf ("*");
            else                printf (" ");
            printf ("%.2f  ", w[i]);
        }
    }
    printf ("\n");

CLEANUP:
    CC_IFFREE (dt_list, int);
    CC_IFFREE (dt_ind, int);
    CC_IFFREE (w, double);
    return rval;
}

static int cross_weight (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *T,
        CCtsp_lpclique *H, double *r)
{
    /* Sets q to x(E(T\setminus H:H)) */
    int rval = 0;
    int j, k, tmp;
    int *marks = (int *) NULL;
    CCtsp_lpnode *n;

    *r = 0;

    marks = CC_SAFE_MALLOC (g->ncount, int);
    CCcheck_NULL (marks, "out of memory for marks");

    CCtsp_mark_clique_and_neighbors (g, H, marks, 0);
    CCtsp_mark_clique (T, marks, 1);
    CCtsp_mark_clique (H, marks, 0);

    CC_FOREACH_NODE_IN_CLIQUE (j, *H, tmp) {
        n = &g->nodes[j];
        for (k = 0; k < n->deg; k++) {
            if (marks[n->adj[k].to]) {
                (*r) += x[n->adj[k].edge];
            }
        }
    }

CLEANUP:
    CC_IFFREE (marks, int);
    return rval;
}

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval;

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}

static int lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x)
{
    int rval;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, xcount,
                     xlist, x, (double **) NULL, (double **) NULL,
                     (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}


static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "M:ps:t", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'M':
            masterfname  = boptarg;
            break;
        case 'p':
            show_dp_stats = 1;
            break;
        case 's':
            seed = atoi(boptarg);
            break;
        case 't':
            show_dd_toothweights = 1;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        lpfname = av[boptind++];
    } else{
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] savfile_root\n", fname);
    fprintf (stderr, "   -M f  specify a master file (required)\n");
    fprintf (stderr, "   -s #  seed\n");
    fprintf (stderr, "   -p    display stats for DP inequalities\n");
    fprintf (stderr, "   -t    display weights of double-decker teeth\n");
}


