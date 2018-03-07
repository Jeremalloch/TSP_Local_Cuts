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
/*                      MIP MODEL FOR COMB SEPARATION                       */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 9, 2011                                                     */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"
#include "macrorus.h"

typedef struct varset {
    int    count;
    int    *varindex;
    int    *sol;
} varset;

typedef struct graph {
    int ncount;
    int ecount;
    int *degree;
    int **adj;
    int *adjspace;
} graph;

static char *rootfname    = (char *) NULL;
static char *masterfname = (char *) NULL;
static int seed = 0;

int
    main (int ac, char **av);

static int
    mip_combs (int ncount, int ecount, int *elist, double *x),
    varset_addcols (CClp *lp, varset *v, double *obj, int varupper),
    varset_add_gamma (CClp *lp, varset *s, varset *gamma, int ecount, 
        int *elist),
    varset_add_double_gamma (CClp *lp, varset *s1, varset *s2, 
        varset *gamma, int ncount, int ecount, int *elist),
    varset_add_delta (CClp *lp, varset *s, varset *delta, int ecount,
        int *elist),
    varset_disjoint_sets (CClp *lp, int cnt, varset **vlist, varset *extra),
    varset_disjoint_edgeteeth (CClp *lp, varset *edgeteeth, int nbig,
        varset **tin, varset **tout, graph *G),
    varset_teeth_parity (CClp *lp, int ecount, varset *edgeteeth,
        int nbig, varset *k),
    varset_nonempty (CClp *lp, varset *v),
    lp_value (CCtsp_lp *lp, double *val),
    lp_x (CCtsp_lp *lp, int *xcount, int **xlist, double **x),
    shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand),
    varset_alloc (varset **v, int n),
    graph_build (graph *G, int ncount, int ecount, int *elist),
    parseargs (int ac, char **av);

static void
    graph_init (graph *G),
    graph_free (graph *G),
    varset_free (varset **v),
    usage (char *fname);

void CClp_changeobjectivesense (CClp *lp, int sense);
int CClp_integervariables (CClp *lp);
int CClp_mipopt (CClp *lp);

int main (int ac, char **av)
{
    int xcount, sncount, secount, ncount, rval = 0, infeasible = 0;
    int *ptour = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    double val, szeit;
    double *x = (double *) NULL, *sx = (double *) NULL;
    int *xlist = (int *) NULL, *selist = (int *) NULL;
    CC_SRKexpinfo expand;

    CCcut_SRK_init_expinfo (&expand);
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

    rval = CCtsp_init_lp (&rootlp, (char *) NULL, -1, rootfname, 0,
               &dat, 0, (int *) NULL, (int *) NULL, 0, (int *) NULL,
               (int *) NULL, 0, ptour, CCtsp_LP_MAXDOUBLE,
               (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, 0, &rstate,
               &infeasible);
    CCcheck_rval (rval, "CCtsp_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "initial LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Run Exact Subtour Loop ...\n");  fflush (stdout);
    rval = CCtsp_subtour_loop (rootlp, 0, 0.0001, &rstate);
    CCcheck_rval (rval, "CCtsp_subtour_loop failed");

    rval = lp_value (rootlp, &val);
    CCcheck_rval (rval, "lp_value failed");
    printf ("LP Value: %f\n", val); fflush (stdout);

    rval = lp_x (rootlp, &xcount, &xlist, &x);
    CCcheck_rval (rval, "lp_x failed");
    printf ("X Count: %d\n", xcount); fflush (stdout);

    rval = shrink_ones (ncount, xcount, xlist, x, &sncount, &secount, &selist,
                        &sx, &expand);
    CCcheck_rval (rval, "shrink_ones failed");
    printf ("Shrunk Graph: %d nodes, %d edges\n", sncount, secount);
    fflush (stdout);

    rval = mip_combs (sncount, secount, selist, sx);
    CCcheck_rval (rval, "mip_combs failed");

#if 0
{
    int dncount = 6, decount = 9;
    int delist[18];
    double dx[9];
    int i;

    i = 0;
    delist[2*i] = 0; delist[2*i+1] = 1; dx[i] = .5;
    i++;
    delist[2*i] = 0; delist[2*i+1] = 2; dx[i] = .5;
    i++;
    delist[2*i] = 1; delist[2*i+1] = 2; dx[i] = .5;
    i++;
    delist[2*i] = 2; delist[2*i+1] = 3; dx[i] = 1.0;
    i++;
    delist[2*i] = 0; delist[2*i+1] = 4; dx[i] = 1.0;
    i++;
    delist[2*i] = 1; delist[2*i+1] = 5; dx[i] = 1.0;
    i++;
    delist[2*i] = 3; delist[2*i+1] = 4; dx[i] = 0.5;
    i++;
    delist[2*i] = 3; delist[2*i+1] = 5; dx[i] = 0.5;
    i++;
    delist[2*i] = 4; delist[2*i+1] = 5; dx[i] = 0.5;
   
    rval = mip_combs (dncount, decount, delist, dx);
    CCcheck_rval (rval, "mip_combs failed");
}
#endif


CLEANUP:

    CC_IFFREE (ptour, int);
    CC_IFFREE (xlist, int);
    CC_IFFREE (x, double);
    CC_IFFREE (selist, int);
    CC_IFFREE (sx, double);
    CCutil_freedatagroup (&dat);
    CCcut_SRK_free_expinfo (&expand);
    if (rootlp) CCtsp_free_tsp_lp_struct (&rootlp);

    return rval;
}

#define NBIG 3

static int mip_combs (int ncount, int ecount, int *elist, double *x)
{
    int i, j, rval = 0;
    CClp *lp = (CClp *) NULL;
    varset *h = (varset *) NULL, *hgamma = (varset *) NULL;
    varset *tin[NBIG], *tout[NBIG], *tgamma[NBIG];
    varset *edgeteeth = (varset *) NULL;
    varset *k = (varset *) NULL;
    varset *vlist[NBIG+1];
    double *negobj = (double *) NULL;
    double *eobj = (double *) NULL;
    double o[1];
    double val;
    graph G;

    printf ("mip_combs(%d, %d) ...\n", ncount, ecount); fflush (stdout);

    for (i = 0; i < NBIG; i++) {
        tin[i] = (varset *) NULL;
        tout[i] = (varset *) NULL;
    }
    graph_init (&G);

    rval = CClp_init (&lp);
    CCcheck_rval (rval, "CClp_init failed");

    rval = CClp_create (lp, "combmip");
    CCcheck_rval (rval, "CClp_create");

    negobj = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (negobj, "out of memory for negobj");
    for (i = 0; i < ncount; i++) negobj[i] = -1;

    rval = varset_alloc (&h, ncount);
    CCcheck_rval (rval, "varset_alloc failed");
    rval = varset_addcols (lp, h, negobj, 1);
    CCcheck_rval (rval, "varset_addcols failed");

    rval = varset_alloc (&hgamma, ecount);
    CCcheck_rval (rval, "varset_alloc failed");
    rval = varset_addcols (lp, hgamma, x, 1);
    CCcheck_rval (rval, "varset_addcols failed");

    for (i = 0; i < NBIG; i++) {
        rval = varset_alloc (&(tin[i]), ncount);
        CCcheck_rval (rval, "varset_alloc failed");
        rval = varset_addcols (lp, tin[i], negobj, 1);
        CCcheck_rval (rval, "varset_addcols failed");

        rval = varset_alloc (&(tout[i]), ncount);
        CCcheck_rval (rval, "varset_alloc failed");
        rval = varset_addcols (lp, tout[i], negobj, 1);
        CCcheck_rval (rval, "varset_addcols failed");

        rval = varset_alloc (&(tgamma[i]), ecount);
        CCcheck_rval (rval, "varset_alloc failed");
        rval = varset_addcols (lp, tgamma[i], x, 1);
        CCcheck_rval (rval, "varset_addcols failed");
    }

    eobj = CC_SAFE_MALLOC (ecount, double);
    CCcheck_NULL (eobj, "out of memory for eobj");
    for (i = 0; i < ecount; i++) eobj[i] = x[i] - 1.0;

    rval = varset_alloc (&edgeteeth, ecount);
    CCcheck_rval (rval, "varset_alloc failed");
    rval = varset_addcols (lp, edgeteeth, eobj, 1);
    CCcheck_rval (rval, "varset_addcols failed");

    rval = varset_alloc (&k, 1);
    CCcheck_rval (rval, "varset_alloc failed");
    o[0] = 1.0;
    rval = varset_addcols (lp, k, o, ncount/2);
    CCcheck_rval (rval, "varset_addcols failed");
    printf ("Number of teeth is stored in x%d\n", k->varindex[0] + 1);
    fflush (stdout);

    rval = varset_add_gamma (lp, h, hgamma, ecount, elist);
    CCcheck_rval (rval, "varset_add_gamma failed");

    for (i = 0; i < NBIG; i++) {
        rval = varset_add_double_gamma (lp, tin[i], tout[i], tgamma[i], ncount,
                                        ecount, elist);
        CCcheck_rval (rval, "varset_add_double_gamma failed");
    }

    for (i = 0; i < NBIG; i++) {
        vlist[i] = tout[i];
    }
    vlist[NBIG] = h;
    
    rval = varset_disjoint_sets (lp, NBIG+1, vlist, (varset *) NULL);
    CCcheck_rval (rval, "varset_disjoint_sets failed");

    rval = varset_disjoint_sets (lp, NBIG, tin, h);
    CCcheck_rval (rval, "varset_disjoint_sets failed");

    for (i = 0; i < NBIG; i++) {
        rval = varset_nonempty (lp, tin[i]);
        CCcheck_rval (rval, "varset_nonempty failed");
        rval = varset_nonempty (lp, tout[i]);
        CCcheck_rval (rval, "varset_nonempty failed");
    }

    rval = varset_add_delta (lp, h, edgeteeth, ecount, elist);
    CCcheck_rval (rval, "varset_add_delta failed");

    rval = graph_build (&G, ncount, ecount, elist);
    CCcheck_rval (rval, "graph_build failed");

    rval = varset_disjoint_edgeteeth (lp, edgeteeth, NBIG, tin, tout, &G);
    CCcheck_rval (rval, "varset_disjoint_edgeteeth failed");

    rval = varset_teeth_parity (lp, ecount, edgeteeth, NBIG, k);
    CCcheck_rval (rval, "varset_teeth_parity failed");

    CClp_changeobjectivesense (lp, 1);
    rval = CClp_integervariables (lp);
    CCcheck_rval (rval, "CClp_integervariables failed");

    
    rval = CClp_dump_lp (lp, "looky.lp");
    CCcheck_rval (rval, "CClp_dump_lp failed");

/*
    rval = CClp_opt (lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed");
*/

    rval = CClp_mipopt (lp);
    CCcheck_rval (rval, "CClp_mipopt failed");

    rval = CClp_objval (lp, &val);
    CCcheck_rval (rval, "CClp_objval failed");

    printf ("Comb LP: %lf\n", val); fflush (stdout);

    {
        double *y = (double *) NULL;
        int ncols;

        ncols = CClp_ncols (lp);
        y = CC_SAFE_MALLOC (ncols, double);
        CCcheck_NULL (y, "out of memory for y");

        rval = CClp_x (lp, y);
        CCcheck_rval (rval, "CClp_x failed");

        printf ("Handle:");
        for (i = 0; i < ncount; i++) { 
            if (y[h->varindex[i]] > 0.5) printf (" %d", i);
        }
        printf ("\n"); fflush (stdout);
        printf ("Number of Teeth: %.0lf\n", y[k->varindex[0]]);
        for (i = 0; i < NBIG; i++) {
            printf ("Tooth %d:", i);
            for (j = 0; j < ncount; j++) {
                if (y[tin[i]->varindex[j]] > 0.5) printf (" %d (in)", j);
            } 
            for (j = 0; j < ncount; j++) {
                if (y[tout[i]->varindex[j]] > 0.5) printf (" %d (out)", j);
            } 
            printf ("\n"); fflush (stdout);
        }
        printf ("Edgeteeth:");
        for (i = 0; i < ecount; i++) {
            if (y[edgeteeth->varindex[i]] > 0.5) printf (" %d", i);
        }
        printf ("\n");

        val = 0.0;
        printf ("gamma(h):");
        for (i = 0; i < ecount; i++) {
            if (y[hgamma->varindex[i]] > 0.5) {
                printf (" %d", i);
                val += x[i];
            }
        }
        printf ("  Val = %.2f\n", val); fflush (stdout);

        for (j = 0; j < NBIG; j++) {
            val = 0.0;
            printf ("gamma(tooth(%d)):", j);
            for (i = 0; i < ecount; i++) {
                if (y[tgamma[j]->varindex[i]] > 0.5) {
                    printf (" %d", i);
                    val += x[i];
                }
            }
            printf ("  Val = %.2f\n", val); fflush (stdout);
        }
/*
        for (j = 0; j < NBIG; j++) {
            printf ("Variables for Gamma Tooth %d:", j);
            for (i = 0; i < ecount; i++) {
                printf (" x%d", tgamma[j]->varindex[i] + 1);
            }
            printf ("\n"); fflush (stdout);
        }
*/

        CC_IFFREE (y, double);
    }
    

CLEANUP:

    varset_free (&h);
    varset_free (&hgamma);
    for (i = 0; i < NBIG; i++) {
        varset_free (&(tin[i]));
        varset_free (&(tout[i]));
        varset_free (&(tgamma[i]));
    }
    varset_free (&edgeteeth);
    varset_free (&k);
    CClp_free (&lp);
    CC_IFFREE (negobj, double);
    CC_IFFREE (eobj, double);
    graph_free (&G);
    return rval;
}

static int varset_addcols (CClp *lp, varset *v, double *obj, int varupper)
{
    int start, finish, i, rval = 0;
    double *lower = (double *) NULL, *upper = (double *) NULL;

    lower = CC_SAFE_MALLOC (v->count, double);
    CCcheck_NULL (lower, "out of memory for lower");
    upper = CC_SAFE_MALLOC (v->count, double);
    CCcheck_NULL (lower, "out of memory for lower");

    for (i = 0; i < v->count; i++) {
        lower[i] = 0.0;
        upper[i] = (double) varupper;
    }

    start = CClp_ncols (lp);
    printf ("CClp_addcols (%d) ...\n", v->count);  fflush (stdout);
    rval = CClp_addcols (lp, v->count, 0, obj, (int *) NULL, (int *) NULL,
                         (double *) NULL, lower, upper);
    CCcheck_rval (rval, "CClp_addcols failed");

    finish = CClp_ncols (lp);
    if (finish - start != v->count) {
        printf ("Some columns did not enter LP\n"); fflush (stdout);
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < v->count; i++) {
        v->varindex[i] = i+start;
    }
    printf ("Number of Columns = %d\n", finish); fflush (stdout);

CLEANUP:

    CC_IFFREE (lower, double);
    CC_IFFREE (upper, double);
    return rval;
}

static int varset_add_gamma (CClp *lp, varset *s, varset *gamma, int ecount,
        int *elist)
{
    int i, start, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;

    /* h_e - h_u <= 0  and  h_e - h_v <= 0 for all edges e = uv */

    start = CClp_nrows (lp);
    printf ("c%d-c%d:  h_e - h_u <= 0, h_e - h_v <= 0 for all e\n",
               start+1, start+2*ecount);
    fflush (stdout);

    rhs = CC_SAFE_MALLOC (2*ecount, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    for (i = 0; i < 2*ecount; i++) rhs[i] = 0.0;

    sense = CC_SAFE_MALLOC (2*ecount, char);
    CCcheck_NULL (sense, "out of memory for sense");
    for (i = 0; i < 2*ecount; i++) sense[i] = 'L';

    rmatbeg = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (4*ecount, int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (4*ecount, double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    for (i = 0; i < ecount; i++) {
        rmatbeg[2*i] = 4*i;
        rmatind[4*i]   = gamma->varindex[i];
        rmatval[4*i]   = 1.0;
        rmatind[4*i+1] = s->varindex[elist[2*i]];
        rmatval[4*i+1] = -1.0;
        rmatbeg[2*i+1] = 4*i+2;
        rmatind[4*i+2] = gamma->varindex[i];
        rmatval[4*i+2] = 1.0;
        rmatind[4*i+3] = s->varindex[elist[2*i+1]];
        rmatval[4*i+3] = -1.0;
    }

    rval = CClp_addrows (lp, 2*ecount, 4*ecount, rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

    /* h_e - h_u - h_v >= -1 for all edges e = uv */

    start = CClp_nrows (lp);
    printf ("c%d-c%d:  h_e - h_u - h_v >= -1 for all e\n",
               start+1, start+ecount);
    fflush (stdout);

    for (i = 0; i < ecount; i++) rhs[i] = -1.0;
    for (i = 0; i < ecount; i++) sense[i] = 'G';

    for (i = 0; i < ecount; i++) {
        rmatbeg[i] = 3*i;
        rmatind[3*i]   = gamma->varindex[i];
        rmatval[3*i]   = 1.0;
        rmatind[3*i+1] = s->varindex[elist[2*i]];
        rmatval[3*i+1]   = -1.0;
        rmatind[3*i+2] = s->varindex[elist[2*i+1]];
        rmatval[3*i+2]   = -1.0;
    }

    rval = CClp_addrows (lp, ecount, 3*ecount, rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
    return rval;
}

static int varset_add_double_gamma (CClp *lp, varset *s1, varset *s2, 
        varset *gamma, int ncount, int ecount, int *elist)
{
    int i, start, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;

    /* t_e - tin_u - tout_u <= 0 and t_e - tin_v - tin_u <= 0 */

    start = CClp_nrows (lp);
    printf ("c%d-c%d: t_e - tin_u - tout_u <= 0, t_e - tin_v - tin_u <= 0\n", 
               start+1, start+2*ecount);
    fflush (stdout);

    rhs = CC_SAFE_MALLOC (2*ecount, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    for (i = 0; i < 2*ecount; i++) rhs[i] = 0.0;

    sense = CC_SAFE_MALLOC (2*ecount, char);
    CCcheck_NULL (sense, "out of memory for sense");
    for (i = 0; i < 2*ecount; i++) sense[i] = 'L';

    rmatbeg = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (6*ecount, int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (6*ecount, double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    for (i = 0; i < ecount; i++) {
        rmatbeg[2*i] = 6*i;
        rmatind[6*i]   = gamma->varindex[i];
        rmatval[6*i]   = 1.0;
        rmatind[6*i+1] = s1->varindex[elist[2*i]];
        rmatval[6*i+1] = -1.0;
        rmatind[6*i+2] = s2->varindex[elist[2*i]];
        rmatval[6*i+2] = -1.0;
        rmatbeg[2*i+1] = 6*i+3;
        rmatind[6*i+3] = gamma->varindex[i];
        rmatval[6*i+3] = 1.0;
        rmatind[6*i+4] = s1->varindex[elist[2*i+1]];
        rmatval[6*i+4] = -1.0;
        rmatind[6*i+5] = s2->varindex[elist[2*i+1]];
        rmatval[6*i+5] = -1.0;
    }

    rval = CClp_addrows (lp, 2*ecount, 6*ecount, rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

    /* t_e - tin_u - tout_u - tin_v - tout_v >= -1 for all edges e = uv */

    start = CClp_nrows (lp);
    printf ("c%d-c%d: t_e - tin_u - tout_u - tin_v - tout_v >= -1 for all e\n",
               start+1, start+ecount);
    fflush (stdout);

    for (i = 0; i < ecount; i++) rhs[i] = -1.0;
    for (i = 0; i < ecount; i++) sense[i] = 'G';

    for (i = 0; i < ecount; i++) {
        rmatbeg[i] = 5*i;
        rmatind[5*i]   = gamma->varindex[i];
        rmatval[5*i]   = 1.0;
        rmatind[5*i+1] = s1->varindex[elist[2*i]];
        rmatval[5*i+1]   = -1.0;
        rmatind[5*i+2] = s2->varindex[elist[2*i]];
        rmatval[5*i+2]   = -1.0;
        rmatind[5*i+3] = s1->varindex[elist[2*i+1]];
        rmatval[5*i+3]   = -1.0;
        rmatind[5*i+4] = s2->varindex[elist[2*i+1]];
        rmatval[5*i+4]   = -1.0;
    }

    rval = CClp_addrows (lp, ecount, 5*ecount, rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

    /* tin_v + tout_v <= 1 for all v */

    start = CClp_nrows (lp);
    printf ("c%d-c%d: tin_v + tout_v <= 1 for all v\n",
               start+1, start+ecount);
    fflush (stdout);

    for (i = 0; i < ncount; i++) rhs[i] = 1.0;
    for (i = 0; i < ncount; i++) sense[i] = 'L';

    for (i = 0; i < ncount; i++) {
        rmatbeg[i] = 2*i;
        rmatind[2*i]   = s1->varindex[i];
        rmatval[2*i]   = 1.0;
        rmatind[2*i+1] = s2->varindex[i];
        rmatval[2*i+1] = 1.0;
    }

    rval = CClp_addrows (lp, ncount, 2*ncount, rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");


CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
    return rval;
}

static int varset_add_delta (CClp *lp, varset *s, varset *delta, int ecount,
        int *elist)
{
    int i, start, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;

    /* y_e + h_u + h_v <= 2, y_e - h_u - h_v <= 0 for all edges e = uv */
    /* The y variables are in the delta struct, the h are in s struct  */

    start = CClp_nrows (lp);
    printf ("c%d-c%d:  y_e + h_u + h_v <= 2, y_e - h_u - h_v <= 0 for all e\n",
               start+1, start+2*ecount);
    fflush (stdout);

    rhs = CC_SAFE_MALLOC (2*ecount, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    for (i = 0; i < ecount; i++) {
        rhs[2*i] = 2.0;
        rhs[2*i+1] = 0.0;
    }

    sense = CC_SAFE_MALLOC (2*ecount, char);
    CCcheck_NULL (sense, "out of memory for sense");
    for (i = 0; i < ecount; i++) {
        sense[2*i] = 'L';
        sense[2*i+1] = 'L';
    }

    rmatbeg = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (6*ecount, int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (6*ecount, double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    for (i = 0; i < ecount; i++) {
        rmatbeg[2*i] = 6*i;
        rmatind[6*i]   = delta->varindex[i];
        rmatval[6*i]   = 1.0;
        rmatind[6*i+1] = s->varindex[elist[2*i]];
        rmatval[6*i+1] = 1.0;
        rmatind[6*i+2] = s->varindex[elist[2*i+1]];
        rmatval[6*i+2] = 1.0;
        rmatbeg[2*i+1] = 6*i+3;
        rmatind[6*i+3] = delta->varindex[i];
        rmatval[6*i+3] = 1.0;
        rmatind[6*i+4] = s->varindex[elist[2*i]];
        rmatval[6*i+4] = -1.0;
        rmatind[6*i+5] = s->varindex[elist[2*i+1]];
        rmatval[6*i+5] = -1.0;
    }

    rval = CClp_addrows (lp, 2*ecount, 6*ecount, rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
    return rval;
}

static int varset_disjoint_sets (CClp *lp, int cnt, varset **vlist,
        varset *extra)
{
    int i, j, t, start, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;
    int ncount;

    t = (extra ? 1 : 0);

    ncount = vlist[0]->count;
    for (i = 1; i < cnt; i++) {
	if (vlist[i]->count != ncount) {
	    printf ("Mismatch in vlist.\n"); fflush (stdout);
            rval = 1;  goto CLEANUP;
	}
    }
    if (extra && extra->count != ncount) {
        printf ("Mismatch in extra.\n"); fflush (stdout);
        rval = 1;  goto CLEANUP;

    }

    /* extra == null:  sum var_v <= 1 for all vertices v       */
    /* extra != null:  sum var_v <= extra_v for all vertices v */

    start = CClp_nrows (lp);
    if (!extra) {
        printf ("c%d-c%d: sum tout_v + h_v <= 1\n", start+1, start+ncount);
    } else {
        printf ("c%d-c%d: sum tout_v <= h_v\n", start+1, start+ncount);
    }
    fflush (stdout);

    rhs = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    for (i = 0; i < ncount; i++) {
        rhs[i] = (extra ? 0.0 : 1.0);
    }

    sense = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (sense, "out of memory for sense");
    for (i = 0; i < ncount; i++) sense[i] = 'L';

    rmatbeg = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (ncount*(cnt+t), int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (ncount*(cnt+t), double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    for (i = 0; i < ncount; i++) {
        rmatbeg[i] = i*(cnt+t);
        for (j = 0; j < cnt; j++) {
            rmatind[i*(cnt+t) + j] = vlist[j]->varindex[i];
            rmatval[i*(cnt+t) + j] = 1.0;
        }
        if (extra) {
            rmatind[i*(cnt+t) + j] = extra->varindex[i];
            rmatval[i*(cnt+t) + j] = -1.0;
        }
    }

    rval = CClp_addrows (lp, ncount, ncount*(cnt+t), rhs, sense, rmatbeg,
                         rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
    return rval;
}

static int varset_disjoint_edgeteeth (CClp *lp, varset *edgeteeth, int nbig,
        varset **tin, varset **tout, graph *G)
{
    int i, j, k, start, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;
    int ncount = G->ncount;
    int ecount = G->ecount;

    /* sum(y_e: e in delta(v) + sum(tin(v)) + sum(toutv) <= 1 for all v */
    /* y_e are stored in edgeteeth struct, nbig is number of big teeth  */

    start = CClp_nrows (lp);
    printf ("c%d-c%d: sum (y_e : e meets v) + sum (tin_v + tout_v) <= 1\n",
                      start+1, start+ncount);
    fflush (stdout);

    rhs = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    for (i = 0; i < ncount; i++) rhs[i] = 1.0;

    sense = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (sense, "out of memory for sense");
    for (i = 0; i < ncount; i++) sense[i] = 'L';

    rmatbeg = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (2*ncount*nbig + 2*ecount, int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (2*ncount*nbig + 2*ecount, double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    k = 0;
    for (i = 0; i < ncount; i++) {
        rmatbeg[i] = k;
        for (j = 0; j < G->degree[i]; j++) {
            rmatind[k+j] = edgeteeth->varindex[G->adj[i][j]];
            rmatval[k+j] = 1.0;
        }
        k += G->degree[i];
        for (j = 0; j < nbig; j++) {
            rmatind[k+(2*j)]   = tin[j]->varindex[i];
            rmatval[k+(2*j)]   = 1.0;
            rmatind[k+(2*j+1)] = tout[j]->varindex[i];
            rmatval[k+(2*j+1)] = 1.0;
        }
        k += (2*nbig);
    }

    rval = CClp_addrows (lp, ncount, 2*ncount*nbig + 2*ecount, rhs, sense,
                         rmatbeg, rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
    return rval;
}

static int varset_teeth_parity (CClp *lp, int ecount, varset *edgeteeth,
        int nbig, varset *k)
{
    int i, start, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;

    /* sum(y_e: e in E) + nbig = 2k+1 */
    /* y_e are stored in edgeteeth struct, nbig is number of big teeth  */

    start = CClp_nrows (lp);
    printf ("c%d: parity constraint for number of teeth\n", start+1);
    fflush (stdout);

    rhs = CC_SAFE_MALLOC (1, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    rhs[0] = 1.0 - (double) nbig;

    sense = CC_SAFE_MALLOC (1, char);
    CCcheck_NULL (sense, "out of memory for sense");
    sense[0] = 'E';

    rmatbeg = CC_SAFE_MALLOC (1, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (ecount+1, int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (ecount+1, double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    rmatbeg[0] = 0;
    for (i = 0; i < ecount; i++) {
        rmatind[i] = edgeteeth->varindex[i];
        rmatval[i] = 1.0;
    }
    rmatind[ecount] = k->varindex[0];
    rmatval[ecount] = -2.0;

    rval = CClp_addrows (lp, 1, ecount+1, rhs, sense, rmatbeg, rmatind,
                         rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");


CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
    return rval;
}

static int varset_nonempty (CClp *lp, varset *v)
{
    int i, rval = 0;
    double *rhs = (double *) NULL;
    char *sense = (char *) NULL;
    int *rmatbeg = (int *) NULL, *rmatind = (int *) NULL;
    double *rmatval = (double *) NULL;
    int ncount = v->count;

    rhs = CC_SAFE_MALLOC (1, double);
    CCcheck_NULL (rhs, "out of memory for rhs");
    rhs[0] = 1.0;

    sense = CC_SAFE_MALLOC (1, char);
    CCcheck_NULL (sense, "out of memory for sense");
    sense[0] = 'G';

    rmatbeg = CC_SAFE_MALLOC (1, int);
    CCcheck_NULL (rmatbeg, "out of memory for rmatbeg");
    rmatind = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (rmatind, "out of memory for rmatind");
    rmatval = CC_SAFE_MALLOC (ncount, double);
    CCcheck_NULL (rmatval, "out of memory for rmatval");

    rmatbeg[0] = 0;
    for (i = 0; i < ncount; i++) {
        rmatind[i] = v->varindex[i];
        rmatval[i] = 1.0;
    }

    rval = CClp_addrows (lp, 1, ncount, rhs, sense, rmatbeg, rmatind, rmatval);
    CCcheck_rval (rval, "CClp_addrows failed");

CLEANUP:

    CC_IFFREE (rhs, double);
    CC_IFFREE (sense, char);
    CC_IFFREE (rmatbeg, int);
    CC_IFFREE (rmatind, int);
    CC_IFFREE (rmatval, double);
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
    int rval;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL, xcount,
                     xlist, x, (double **) NULL, (double **) NULL,
                     (double **) NULL);
    if (rval) fprintf (stderr, "CCtsp_get_lp_result failed\n");

    return rval;
}

static int shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand)
{
    CC_SRKgraph G;
    int k, rval = 0;

    rval = CCcut_SRK_buildgraph (&G, ncount, ecount, elist, dlen);
    CCcheck_rval (rval, "CCcut_SRK_buildgraph failed");

    rval = CCcut_SRK_defluff (&G);
    CCcheck_rval (rval, "CCcut_SRK_defluff failed");

    CCcut_SRK_identify_paths_to_edges (&G, &k, 0);
    CCcut_SRK_identify_one_triangles (&G, &k, (CC_SRKnode *) NULL, 0.001, 2.0,
                                      0);

    rval = CCcut_SRK_grab_edges (&G, oncount, oecount, olist, olen, expand);
    CCcheck_rval (rval, "CCcut_SRK_grab_edges failed");

CLEANUP:

    CCcut_SRK_free_graph (&G);
    return 0;
}

static void graph_init (graph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->degree = (int *) NULL;
        G->adj = (int **) NULL;
        G->adjspace = (int *) NULL;
    }
}

static void graph_free (graph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        CC_IFFREE (G->degree, int);
        CC_IFFREE (G->adj, int *);
        CC_IFFREE (G->adjspace, int);
    }
}

static int graph_build (graph *G, int ncount, int ecount, int *elist)
{
    int rval = 0;
    int i, j;
    int *p;

    G->ncount = ncount;
    G->ecount = ecount;
    G->degree = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (G->degree, "out of memory for G->degree");
    G->adj = CC_SAFE_MALLOC (ncount, int *);
    CCcheck_NULL (G->adj, "out of memory for G->adj");
    G->adjspace = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (G->adjspace, "out of memory for G->adjspace");

    for (i = 0; i < ncount; i++) G->degree[i] = 0;
    for (i = 0; i < ecount; i++) {
        G->degree[elist[2*i]]++;
        G->degree[elist[2*i+1]]++;
    }

    p = G->adjspace;
    for (i = 0; i < ncount; i++) {
        G->adj[i] = p;
        p += G->degree[i];
        G->degree[i] = 0;
    }

    for (i = 0; i < ecount; i++) {
        j = elist[2*i];
        G->adj[j][G->degree[j]] = i;
        G->degree[j]++;
        j = elist[2*i+1];
        G->adj[j][G->degree[j]] = i;
        G->degree[j]++;
    }

CLEANUP:

    if (rval) graph_free (G);
    return rval;
}

static int varset_alloc (varset **v, int n)
{
    int rval = 0;

    *v = CC_SAFE_MALLOC (1, varset);
    CCcheck_NULL (*v, "out of memory for v");

    (*v)->varindex = CC_SAFE_MALLOC (n, int);
    (*v)->sol = CC_SAFE_MALLOC (n, int);
    CCcheck_NULL ((*v)->varindex, "out of memory for varindex");
    CCcheck_NULL ((*v)->sol, "out of memory for sol");
    (*v)->count = n;

CLEANUP:
   
    if (rval) {
        if (*v) varset_free (v);
    }
    return rval;
}

static void varset_free (varset **v)
{
    if (*v) {
        CC_IFFREE ((*v)->varindex, int);
        CC_IFFREE ((*v)->sol, int);
        CC_IFFREE (*v, varset);
        *v = (varset *) NULL;
    }
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "M:s:", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'M':
            masterfname  = boptarg;
            break;
        case 's':
            seed = atoi(boptarg);
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        rootfname = av[boptind++];
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
    fprintf (stderr, "   -u #  upperbound on optimal tour length\n");
}


