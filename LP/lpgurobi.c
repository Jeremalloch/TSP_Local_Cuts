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
/*                                                                          */ /****************************************************************************/ /****************************************************************************/
/*                                                                          */
/*              Interface Routines to the Gurobi Solver                     */
/*                                                                          */
/*  NOTE: Use this code in place of lp_none.c to access the Gurobi 0.91     */
/*   libraries. You will also need to link the gurobi library via the       */
/*   makefile.                                                              */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 31, 2009                                                    */
/*                                                                          */
/*  NOTES:                                                                  */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "lp.h"
#include "tsp.h"
#include <gurobi_c.h>

#undef  CC_GUROBI_DISPLAY

struct CClp {
    GRBenv    *gurobi_env;
    GRBmodel  *gurobi_lp;
};

struct CClp_warmstart {
    int      rcount;
    int      ccount;
    int     *rstat;
    int     *cstat;
    double  *dnorm;
};

struct CClp_info {
    int  rcount;
    int  ccount;
    int *rstat;
    int *cstat;
};

#define SOLVER_WARMSTART_NAME "GRB1"


static int
    primalopt (CClp *lp),
    dualopt (CClp *lp),
    baropt (CClp *lp);
static void
    printerrorcode (int c);


int CClp_init (CClp **lp)
{
    int rval = 0;

    CClp_free (lp);

    (*lp) = CC_SAFE_MALLOC (1, CClp);
    if ((*lp) == (CClp *) NULL) {
        fprintf (stderr, "Out of memory in CClp_init\n");
        rval = 1; goto CLEANUP;
    }

    (*lp)->gurobi_env = (GRBenv *) NULL;
    (*lp)->gurobi_lp = (GRBmodel *) NULL;

    rval = GRBloadenv (&((*lp)->gurobi_env), NULL);
    CCcheck_rval (rval, "GRBloadenv failed");


#ifdef CC_GUROBI_DISPLAY
    rval = GRBsetintparam ((*lp)->gurobi_env, GRB_INT_PAR_OUTPUTFLAG, 1);
#else
    rval = GRBsetintparam ((*lp)->gurobi_env, GRB_INT_PAR_OUTPUTFLAG, 0);
#endif
    CCcheck_rval (rval, "GRBsetintparam OUTPUTFLAG failed");

    /* Set to dual steepest edge */
    rval = GRBsetintparam ((*lp)->gurobi_env, GRB_INT_PAR_METHOD, 1);
    CCcheck_rval (rval, "GRBsetintparam failed for METHOD");
    rval = GRBsetintparam ((*lp)->gurobi_env, GRB_INT_PAR_SIMPLEXPRICING, 1);
    CCcheck_rval (rval, "GRBsetintparam failed for SIMPLEXPRICING");

    /* Following REB, 14 October 1997 for CPLEX settings */

    rval = GRBsetdblparam ((*lp)->gurobi_env,GRB_DBL_PAR_PERTURBVALUE,1.0E-6);
    CCcheck_rval (rval, "GRBsetdblparam GRB_DLP_PAR_PERTURBVALUE failed");

    rval = GRBsetdblparam ((*lp)->gurobi_env,GRB_DBL_PAR_OPTIMALITYTOL,1.0E-9);
    CCcheck_rval (rval, "GRBsetdblparam GRB_DLP_PAR_OPTIMALITYTOL failed");

    rval = GRBsetdblparam ((*lp)->gurobi_env,GRB_DBL_PAR_FEASIBILITYTOL,1.0E-9);
    CCcheck_rval (rval, "GRBsetdblparam GRB_DLP_PAR_FEASIBILITYTOL failed");

CLEANUP:

    if (rval) {
        CClp_free (lp);
    }

    return rval;
}

int CClp_force_perturb (CClp *lp)
{
    int rval;

    /* Not implemented in Gurobi */

    if (!lp) {
        fprintf (stderr, "CClp_force_perturb called without an lp\n");
        rval = 1;
    }

    return 0;
}

int CClp_tune_small (CClp *lp)
{
    int rval;

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp),
                           GRB_INT_PAR_SIMPLEXPRICING, -1);
    CCcheck_rval (rval, "GRBsetintparam failed for SIMPLEXPRICING");

CLEANUP:
    return rval;
}

int CClp_disable_presolve (CClp *lp)
{
    int rval = 0;

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_PRESOLVE, 0);
    CCcheck_rval (rval, "GRBsetintparam PRESOLVE failed");

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_AGGREGATE, 0);
    CCcheck_rval (rval, "GRBsetintparam AFFREGATE failed");

CLEANUP:
    return rval;
}

void CClp_free (CClp **lp)
{
    if (*lp) {
        if ((*lp)->gurobi_env) {
            if ((*lp)->gurobi_lp) {
                GRBfreemodel ((*lp)->gurobi_lp);
            }
            GRBfreeenv ((*lp)->gurobi_env);
        }
        CC_FREE (*lp, CClp);
    }
}

void CClp_freelp (CClp **lp)
{
    if (*lp) {
        if ((*lp)->gurobi_lp) {
            GRBfreemodel ((*lp)->gurobi_lp);
        }
    }
}

int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,
        int objsense, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub)
{
    int rval = 0;

    rval = GRBloadmodel (lp->gurobi_env, &(lp->gurobi_lp), name, ncols, nrows,
                     objsense, 0.0, obj, sense, rhs, matbeg, matcnt,
                     matind, matval, lb, ub, (char *) NULL, NULL, NULL);
    if (rval) {
        fprintf (stderr, "GRBloadmodel failed, code: %d\n", rval);
        goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int CClp_create (CClp *lp, const char *name)
{
    int rval = 0;

    rval = GRBnewmodel (lp->gurobi_env, &lp->gurobi_lp, name, 0, 
                  (double *) NULL, (double *) NULL, (double *) NULL,
                  (char *) NULL, NULL);
    if (rval) {
       fprintf (stderr, "GRBnewmodel failed, return code %d\n", rval);
       goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int CClp_new_row (CClp *lp, char sense, double rhs)
{
    int rval = 0;

    rval = GRBaddconstr (lp->gurobi_lp, 0, (int *) NULL, (double *) NULL,
                         sense, rhs, NULL);
    CCcheck_rval (rval, "GRBaddconstr failed");
 
    rval = GRBupdatemodel (lp->gurobi_lp);
    CCcheck_rval (rval, "GRBupdatemodel failed");
 

CLEANUP:
    return rval;
}

int CClp_change_sense (CClp *lp, int row, char sense)
{
    int rval = 0;
    char isense;

    if (sense == 'L')      isense = '<';
    else if (sense == 'E') isense = '=';
    else                   isense = '>';

    CCcheck_rval (rval, "CPXchgsense failed");
    rval = GRBsetcharattrelement (lp->gurobi_lp, GRB_CHAR_ATTR_SENSE,
                                  row, isense);
    CCcheck_rval (rval, "GRBsetcharattrelement failed");

CLEANUP:
    return rval;
}

int CClp_opt (CClp *lp, int method)
{
    int  rval = 0;

    switch (method) {
        case CClp_METHOD_PRIMAL:
            rval = primalopt (lp);
            break;
        case CClp_METHOD_DUAL:
            rval = dualopt (lp);
            break;
        case CClp_METHOD_BARRIER:
            rval = baropt (lp);
            break;
        default:
            rval = 1;
            fprintf (stderr, "Nonexistent method in CClp_opt\n");
            goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static int primalopt (CClp *lp)
{
    int rval;
    int solstat;

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_METHOD, 0);
    CCcheck_rval (rval, "GRBsetintparam failed for METHOD_0");

    rval = GRBoptimize (lp->gurobi_lp);
    if (rval) {
        fprintf (stderr, "GRBoptimize failed, return code %d\n", rval);
        goto CLEANUP;
    }
    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_STATUS, &solstat);
    CCcheck_rval (rval, "GRBgetintattr STATUS");

    if (solstat == GRB_INFEASIBLE) {
        return 2;
    } else if (solstat != GRB_OPTIMAL) {
        fprintf (stderr, "Gurobi optimization status %d\n", solstat);
        return 1;
    }

CLEANUP:

    if (rval) return 1;
    else      return 0;
}

static int dualopt (CClp *lp)
{
    int rval;
    int solstat;

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_METHOD, 1);
    CCcheck_rval (rval, "GRBsetintparam failed for METHOD_1");

    rval = GRBoptimize (lp->gurobi_lp);
    if (rval) {
        fprintf (stderr, "GRBoptimize failed, return code %d\n", rval);
        goto CLEANUP;
    }
    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_STATUS, &solstat);
    CCcheck_rval (rval, "GRBgetintattr STATUS");

    if (solstat == GRB_INFEASIBLE) {
        return 2;
    } else if (solstat != GRB_OPTIMAL) {
        fprintf (stderr, "Gurobi optimization status %d\n", solstat);
        return 1;
    }

CLEANUP:

    if (rval) return 1;
    else      return 0;
}

static int baropt (CClp *lp)
{
    int rval;
    int solstat;

    printf ("baropt ...\n"); fflush (stdout);

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_OUTPUTFLAG, 1);
    CCcheck_rval (rval, "GRBsetintparam failed for OUTPUTFLAG");

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_METHOD, 2);
    CCcheck_rval (rval, "GRBsetintparam failed for METHOD_2");

    /* For the big-instance tests, turn off crossover */

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_CROSSOVER, 0);
    CCcheck_rval (rval, "GRBsetintparam failed for CROSSOVER");

    /* Ed Rothberg suggests BARORDER=0 to try to avoid big-memory */

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_BARORDER, 0);
    CCcheck_rval (rval, "GRBsetintparam failed for CROSSOVER");

    rval = GRBoptimize (lp->gurobi_lp);
    if (rval) {
        fprintf (stderr, "GRBoptimize failed, return code %d\n", rval);
        goto CLEANUP;
    }
    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_STATUS, &solstat);
    CCcheck_rval (rval, "GRBgetintattr STATUS");

    if (solstat == GRB_INFEASIBLE) {
        return 2;
    } else if (solstat != GRB_OPTIMAL) {
        fprintf (stderr, "Gurobi optimization status %d\n", solstat);
        return 1;
    }

CLEANUP:

    rval = GRBsetintparam (GRBgetenv(lp->gurobi_lp), GRB_INT_PAR_OUTPUTFLAG, 0);
    CCcheck_rval (rval, "GRBsetintparam failed for OUTPUTFLAG");

    if (rval) return 1;
    else      return 0;
}


int CClp_limited_dualopt (CClp *lp, int iterationlim, int *status,
                          double *objupperlim)
{
    int rval = 0;
    int sval = 0;
    int solstat;

    /* Set for dual */

    double old_iterationlim;
    int    got_iterationlim = 0;
    double old_objupperlim;
    int    got_objupperlim = 0;

    rval = GRBgetdblparam (GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_ITERATIONLIMIT,
                           &old_iterationlim);
    CCcheck_rval (rval, "GRBgetintparam INTERATION_LIMIT failed\n");
    got_iterationlim = 1;

    rval = GRBgetdblparam (GRBgetenv(lp->gurobi_lp),GRB_DBL_PAR_CUTOFF,&old_objupperlim);
    CCcheck_rval (rval, "GRBgetdblparam CUTOFF failed\n");
    got_objupperlim = 1;

    rval = GRBsetdblparam (GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_ITERATIONLIMIT, 
                          (double) iterationlim);
    CCcheck_rval (rval, "GRBsetintparam ITERATIONLIMIT failed\n");

    if (objupperlim) {
        rval = GRBsetdblparam (GRBgetenv(lp->gurobi_lp),GRB_DBL_PAR_CUTOFF,*objupperlim);
        CCcheck_rval (rval, "GRBsetdblparam CUTOFF failed\n");
    }

    rval = GRBoptimize (lp->gurobi_lp);
    if (rval) {
        fprintf (stderr, "GRBoptimize failed, return code %d\n", rval);
        return 1;
    }
    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_STATUS, &solstat);
    if (rval) {
        fprintf (stderr, "GRBgetintattr failed, return code %d\n", rval);
        return 1;
    }

    if (solstat == GRB_INFEASIBLE) {
        printf ("Infeasible in dual\n"); fflush (stdout);
        if (status) *status = CClp_INFEASIBLE;
    } else if (solstat == GRB_UNBOUNDED) {
        if (status) *status = CClp_UNBOUNDED;
    } else if (solstat != GRB_OPTIMAL          &&
               solstat != GRB_CUTOFF    &&
               solstat != GRB_ITERATION_LIMIT) {
        fprintf (stderr, "Gurobi optimization status %d\n", solstat);
        if (status) *status = CClp_FAILURE;
    } else {
        if (status) *status = CClp_SUCCESS;
    }

CLEANUP:

    if (got_iterationlim == 1) {
        sval = GRBsetdblparam (GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_ITERATIONLIMIT,
                               old_iterationlim);
        if (sval) {
            fprintf (stderr, "GRBsetdblparam ITERATIONLIMIT failed\n");
            rval = 1;
        }
    }

    if (got_objupperlim == 1) {
        sval = GRBsetdblparam (GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_CUTOFF,
                               old_objupperlim);
        if (sval) {
            fprintf (stderr, "GRBsetdblparam CUTOFF failed\n");
            rval = 1;
        }
    }

    if (rval && status) {
        *status = CClp_FAILURE;
    }

    return rval;
}

int CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs, char *sense,
                  int *rmatbeg, int *rmatind, double *rmatval)
{
    int rval = 0;

    rval = GRBaddconstrs (lp->gurobi_lp, newrows, newnz, rmatbeg, rmatind,
                        rmatval, sense, rhs, NULL);
    if (rval) {
        printf ("GRBaddconstrs returned %d\n", rval); fflush (stdout);
        printerrorcode (rval);
    }
    CCcheck_rval (rval, "GRBaddconstrs failed");

    rval = GRBupdatemodel (lp->gurobi_lp);
    CCcheck_rval (rval, "GRBupdatemodel failed");

CLEANUP:

    return rval;
}

int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,
                  int *cmatbeg, int *cmatind, double *cmatval,
                  double *lb, double *ub)
{
    int rval = 0;

    rval = GRBaddvars (lp->gurobi_lp, newcols, newnz, cmatbeg, cmatind,
                       cmatval, obj, lb, ub, (char *) NULL, NULL);
    if (rval) {
        printf ("GRBaddvars returned %d\n", rval); fflush (stdout);
        printerrorcode (rval);
    }
    CCcheck_rval (rval, "GRBaddvars failed");

    rval = GRBupdatemodel (lp->gurobi_lp);
    CCcheck_rval (rval, "GRBupdatemodel failed");

CLEANUP:
    return rval;
}

int CClp_delete_row (CClp *lp, int i)
{
    int rval = 0;
    int locali[1];

    locali[0] = i;
    rval = GRBdelconstrs (lp->gurobi_lp, 1, locali);
    CCcheck_rval (rval, "GRBdelconstrs failed");

CLEANUP:

    return rval;
}

int CClp_delete_set_of_rows (CClp *lp, int *delstat)
{
    int i, j, rval = 0;
    int *dellist = (int *) NULL;
    int delcnt = 0;
    int rcnt = CClp_nrows (lp);

    for (i=0; i<rcnt; i++) {
        if (delstat[i]) delcnt++;
    }
    if (delcnt == 0) {
        fprintf (stderr, "delete_set_of_rows with no deleted rows\n");
        return 0;
    }
    dellist = CC_SAFE_MALLOC (delcnt, int);
    CCcheck_NULL (dellist, "out of memory for dellist");

    for (i=0, j=0; i<rcnt; i++) {
        if (delstat[i]) {
            dellist[j++] = i;
        }
    }
    if (j != delcnt) {
        fprintf (stderr, "Lost some deleted rows\n");
        CC_FREE (dellist, int);
        return 1;
    }

    rval = GRBdelconstrs (lp->gurobi_lp, delcnt, dellist);
    CCcheck_rval (rval, "GRBdelconstrs failed");

CLEANUP:

    CC_IFFREE (dellist, int);
    return rval;
}

int CClp_delete_column (CClp *lp, int i)
{
    int rval = 0;
    int locali[1];

    locali[0] = i;
    rval = GRBdelvars (lp->gurobi_lp, 1, locali);
    CCcheck_rval (rval, "GRBdelvars failed");

CLEANUP:

    return rval;
}

int CClp_delete_set_of_columns (CClp *lp, int *delstat)
{
    int i, j, rval = 0;
    int *dellist = (int *) NULL;
    int delcnt = 0;
    int ccnt = CClp_ncols (lp);

    for (i=0; i<ccnt; i++) {
        if (delstat[i]) delcnt++;
    }
    if (delcnt == 0) {
        fprintf (stderr, "delete_set_of_columns with no deleted columns\n");
        rval = 0; goto CLEANUP;
    }
    dellist = CC_SAFE_MALLOC (delcnt, int);
    CCcheck_NULL (dellist, "out of memory for dellist");
    for (i=0, j=0; i<ccnt; i++) {
        if (delstat[i]) {
            dellist[j++] = i;
        }
    }
    if (j != delcnt) {
        fprintf (stderr, "Lost some deleted columns\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBdelvars (lp->gurobi_lp, delcnt, dellist);
    CCcheck_rval (rval, "GRBdelvars failed");

CLEANUP:

    CC_IFFREE (dellist, int);
    return rval;
}

int CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    if (lower_or_upper == 'L') {
        rval = GRBsetdblattrelement (lp->gurobi_lp, GRB_DBL_ATTR_LB, col, bnd);
    } else {
        rval = GRBsetdblattrelement (lp->gurobi_lp, GRB_DBL_ATTR_UB, col, bnd);
    }
    CCcheck_rval (rval, "GRBsetdblattr LB or UB failed");

CLEANUP:

    return rval;
}

int CClp_get_warmstart (CClp *lp, CClp_warmstart **w)
{
    int nrows, ncols, rval = 0;

    CClp_free_warmstart (w);

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->ccount = 0;
    (*w)->rcount = 0;
    (*w)->cstat  = (int *) NULL;
    (*w)->rstat  = (int *) NULL;
    (*w)->dnorm  = (double *) NULL;

    ncols = CClp_ncols (lp);
    nrows = CClp_nrows (lp);
    (*w)->ccount = ncols;
    (*w)->rcount = nrows;

    (*w)->cstat = CC_SAFE_MALLOC ((*w)->ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC ((*w)->rcount, int);
    (*w)->dnorm = CC_SAFE_MALLOC ((*w)->rcount, double);
    if (!(*w)->cstat || !(*w)->rstat || !(*w)->dnorm) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    CC_IFFREE ((*w)->dnorm, double);   /* No dnorm function thus far */

    rval = GRBgetintattrarray(lp->gurobi_lp, GRB_INT_ATTR_VBASIS, 0, ncols,
                                                                (*w)->cstat);
    CCcheck_rval (rval, "GRBgetintattrarray GRB_INT_ATTR_VBASIS failed");

    rval = GRBgetintattrarray(lp->gurobi_lp, GRB_INT_ATTR_CBASIS, 0, nrows,
                                                                (*w)->rstat);
    CCcheck_rval (rval, "GRBgetintattrarray GRB_INT_ATTR_CBASIS failed");


    return 0;

CLEANUP:

    CClp_free_warmstart (w);
    return rval;
}

int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)
{
    int rval = 0;

    if (w->cstat && w->rstat && w->dnorm) {
        fprintf (stderr, "No Gurobi call for dnorms\n");
        rval = 1;  goto CLEANUP;
    } else if (w->cstat && w->rstat) {
        rval = GRBsetintattrarray(lp->gurobi_lp, GRB_INT_ATTR_VBASIS, 0,
                                   w->ccount, w->cstat);
        CCcheck_rval (rval, "GRBsetintattrarray VBASIS failed");

        rval = GRBsetintattrarray(lp->gurobi_lp, GRB_INT_ATTR_CBASIS, 0,
                                   w->rcount, w->rstat);
        CCcheck_rval (rval, "GRBsetintattrarray CBASIS failed");
    } else {
        printf ("WARNING: No basis in call to load_warmstart\n");
        fflush (stdout);
    }

CLEANUP:

    return rval;
}

int CClp_build_warmstart (CClp_warmstart **w, CClp_info *i)
{
    int rval = 0;
    int j;

    CClp_free_warmstart (w);

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->ccount = 0;
    (*w)->rcount = 0;
    (*w)->cstat = (int *) NULL;
    (*w)->rstat = (int *) NULL;
    (*w)->dnorm = (double *) NULL;

    (*w)->ccount = i->ccount;
    if ((*w)->ccount == 0) {
        fprintf (stderr, "No columns in CClp_info\n");
        rval = 1; goto CLEANUP;
    }
    (*w)->rcount = i->rcount;
    if ((*w)->rcount == 0) {
        fprintf (stderr, "No rows in CClp_info\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->cstat = CC_SAFE_MALLOC ((*w)->ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC ((*w)->rcount, int);
    if (!(*w)->cstat || !(*w)->rstat) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    for (j=0; j<(*w)->ccount; j++) {
        (*w)->cstat[j] = i->cstat[j];
    }
    for (j=0; j<(*w)->rcount; j++) {
        (*w)->rstat[j] = i->rstat[j];
    }

    return 0;

CLEANUP:

    CClp_free_warmstart (w);
    return rval;
}

void CClp_free_warmstart (CClp_warmstart **w)
{
    if (*w != (CClp_warmstart *) NULL) {
        CC_IFFREE ((*w)->cstat, int);
        CC_IFFREE ((*w)->rstat, int);
        CC_IFFREE ((*w)->dnorm, double);
        CC_FREE (*w, CClp_warmstart);
    }
}

int CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w)
{
    char name[5];
    int i;
    int ccount;
    int rcount;
    int has_dnorms;

    CClp_free_warmstart (w);

    for (i=0; i<4; i++) {
        if (CCutil_sread_char (f, &name[i])) goto CLEANUP;
    }
    name[4] = '\0';

    if (strncmp (name, SOLVER_WARMSTART_NAME, 4)) {
        fprintf (stderr, "warmstart for another solver (%s) ignored\n", name);
        return 0;
    }

    if (CCutil_sread_int (f, &ccount)) goto CLEANUP;
    if (CCutil_sread_int (f, &rcount)) goto CLEANUP;

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_sread_warmstart\n");
        goto CLEANUP;
    }

    (*w)->ccount = 0;
    (*w)->rcount = 0;
    (*w)->cstat  = (int *) NULL;
    (*w)->rstat  = (int *) NULL;
    (*w)->dnorm  = (double *) NULL;

    (*w)->cstat = CC_SAFE_MALLOC (ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC (rcount, int);
    if ((*w)->cstat == (int *) NULL ||
        (*w)->rstat == (int *) NULL) {
        fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
        goto CLEANUP;
    }
    for (i = 0; i < ccount; i++) {
        if (CCutil_sread_bits (f, &(((*w)->cstat)[i]), 2))
            goto CLEANUP;
        (*w)->cstat[i] *= -1;
    }
    for (i = 0; i < rcount; i++) {
        if (CCutil_sread_bits (f, &(((*w)->rstat)[i]), 1))
            goto CLEANUP;
        (*w)->rstat[i] *= -1;
    }

    if (CCutil_sread_int (f, &has_dnorms)) goto CLEANUP;

    if (has_dnorms) {
        (*w)->dnorm = CC_SAFE_MALLOC (rcount, double);
        if ((*w)->dnorm == (double *) NULL) {
            fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
            goto CLEANUP;
        }
        for (i = 0; i < rcount; i++) {
            if (CCutil_sread_double (f, &(((*w)->dnorm)[i]))) goto CLEANUP;
        }
    }

    (*w)->ccount = ccount;
    (*w)->rcount = rcount;
    
    return 0;

CLEANUP:

    CClp_free_warmstart (w);
    return 1;
}

int CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w)
{
    int i, k;
    const char *name = SOLVER_WARMSTART_NAME;

    for (i=0; i<4; i++) {
        if (CCutil_swrite_char (f, name[i])) return 1;
    }

    if (CCutil_swrite_int (f, w->ccount)) return 1;
    if (CCutil_swrite_int (f, w->rcount)) return 1;

    for (i = 0; i < w->ccount; i++) {
        if (w->cstat[i] < 0) k = -w->cstat[i];
        else                 k =  w->cstat[i];
        if (CCutil_swrite_bits (f, k, 2)) return 1;
    }

    for (i = 0; i < w->rcount; i++) {
        if (w->rstat[i] < 0) k = -w->rstat[i];
        else                 k =  w->rstat[i];
        if (CCutil_swrite_bits (f, k, 1)) return 1;
    }

    if (w->dnorm == (double *) NULL) {
        if (CCutil_swrite_int (f, 0)) return 1;
    } else {
        if (CCutil_swrite_int (f, 1)) return 1;
        for (i = 0; i < w->rcount; i++) {
            if (CCutil_swrite_double (f, w->dnorm[i])) return 1;
        }
    }

    return 0;
}

int CClp_get_info (CClp *lp, CClp_info **i)
{
    int ncols, nrows, rval = 0;

    CClp_free_info (i);

    (*i) = CC_SAFE_MALLOC (1, CClp_info);
    if ((*i) == (CClp_info *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->ccount = 0;
    (*i)->rcount = 0;
    (*i)->cstat = (int *) NULL;
    (*i)->rstat = (int *) NULL;

    ncols = CClp_ncols (lp);
    nrows = CClp_nrows (lp);
    (*i)->ccount = ncols;
    (*i)->rcount = nrows;

    (*i)->cstat = CC_SAFE_MALLOC (ncols, int);
    (*i)->rstat = CC_SAFE_MALLOC (nrows, int);
    if (!(*i)->cstat || !(*i)->rstat) {
        fprintf (stderr, "out of memory in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    rval = GRBgetintattrarray(lp->gurobi_lp, GRB_INT_ATTR_VBASIS, 0, ncols,
                                                                (*i)->cstat);
    CCcheck_rval (rval, "GRBgetintattrarray GRB_INT_ATTR_VBASIS failed");

    rval = GRBgetintattrarray(lp->gurobi_lp, GRB_INT_ATTR_CBASIS, 0, nrows,
                                                                (*i)->rstat);
    CCcheck_rval (rval, "GRBgetintattrarray GRB_INT_ATTR_CBASIS failed");

    return 0;

CLEANUP:

    CClp_free_info (i);
    return rval;
}

int CClp_create_info (CClp_info **i, int rcount, int ccount)
{
    int rval = 0;
    int j;

    CClp_free_info (i);

    (*i) = CC_SAFE_MALLOC (1, CClp_info);
    if ((*i) == (CClp_info *) NULL) {
        fprintf (stderr, "Out of memory in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->ccount = 0;
    (*i)->rcount = 0;
    (*i)->cstat = (int *) NULL;
    (*i)->rstat = (int *) NULL;

    (*i)->ccount = ccount;
    if (ccount == 0) {
        fprintf (stderr, "No columns in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }
    (*i)->rcount = rcount;
    if (rcount == 0) {
        fprintf (stderr, "No rows in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->cstat = CC_SAFE_MALLOC ((*i)->ccount, int);
    (*i)->rstat = CC_SAFE_MALLOC ((*i)->rcount, int);
    if (!(*i)->cstat || !(*i)->rstat) {
        fprintf (stderr, "out of memory in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    for (j=0; j<ccount; j++) {
        (*i)->cstat[j] = 0;
    }
    for (j=0; j<rcount; j++) {
        (*i)->rstat[j] = 0;
    }

    return 0;

CLEANUP:

    CClp_free_info (i);
    return rval;
}

int CClp_is_col_active (CClp_info *i, int c)
{
    if (c < 0 || c >= i->ccount) return 0;
    return i->cstat[c] == 0 || i->cstat[c] == -2;
}

int CClp_is_row_active (CClp_info *i, int r)
{
    if (r < 0 || r >= i->rcount) return 0;
    return i->rstat[r] == -1;
}

void CClp_set_col_active (CClp_info *i, int c)
{
    if (c >= 0 && c < i->ccount) i->cstat[c] = 0;
}

void CClp_set_col_inactive (CClp_info *i, int c)
{
    if (c >= 0 && c < i->ccount) i->cstat[c] = -1;
}

void CClp_set_col_upper (CClp_info *i, int c)
{
    if (c >= 0 && c < i->ccount) i->cstat[c] = -2;
}

void CClp_set_row_active (CClp_info *i, int r)
{
    if (r >= 0 && r < i->rcount) i->rstat[r] = -1;
}

void CClp_set_row_inactive (CClp_info *i, int r)
{
    if (r >= 0 && r < i->rcount) i->rstat[r] = 0;
}

void CClp_free_info (CClp_info **i)
{
    if ((*i) != (CClp_info *) NULL) {
        CC_IFFREE ((*i)->cstat, int);
        CC_IFFREE ((*i)->rstat, int);
        CC_FREE (*i, CClp_info);
    }
}

int CClp_x (CClp *lp, double *x)
{
    int rval = 0;
    int ncols;

    ncols = CClp_ncols (lp);
    rval = GRBgetdblattrarray(lp->gurobi_lp, GRB_DBL_ATTR_X, 0, ncols, x);
    CCcheck_rval (rval, "GRBgetdblattrarray GRB_DBL_ATTR_X failed");

CLEANUP:
    return rval; 
}

int CClp_rc (CClp *lp, double *rc)
{
    int rval = 0;
    int ncols;

    ncols = CClp_ncols (lp);
    rval = GRBgetdblattrarray(lp->gurobi_lp, GRB_DBL_ATTR_RC, 0, ncols, rc);
    CCcheck_rval (rval, "GRBgetdblattrarray GRB_DBL_ATTR_RC failed");

CLEANUP:
    return rval;
}

int CClp_pi (CClp *lp, double *pi)
{
    int rval = 0;
    int nrows;
    int solstat;

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_STATUS, &solstat);
    CCcheck_rval (rval, "GRBgetintattr GRB_INT_ATTR_STATUS failed");
    if (solstat == GRB_INFEASIBLE) {
        fprintf (stderr, "In Gurobi not sure how to get infeas array\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_NUMCONSTRS, &nrows);
    if (nrows == 0) {
        fprintf (stderr, "No rows in LP\n"); 
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetdblattrarray(lp->gurobi_lp, GRB_DBL_ATTR_PI, 0, nrows, pi);
    CCcheck_rval (rval, "GRBgetdblattrarray GRB_DBL_ATTR_PI failed");

CLEANUP:

    return rval;
}

int CClp_objval (CClp *lp, double *obj)
{
    int rval = 0;

    rval = GRBgetdblattr (lp->gurobi_lp, GRB_DBL_ATTR_OBJVAL, obj);
    CCcheck_rval (rval, "CPXgetobjval failed");

CLEANUP:
    return rval;
}

int CClp_nrows (CClp *lp)
{
    int rval, nrows;

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_NUMCONSTRS, &nrows);
    if (rval) {
        fprintf (stderr, "GRBgetintattr GRB_INT_ATTR_NUMCONSTRS failed\n");
        return 0;
    }

    return nrows;
}

int CClp_ncols (CClp *lp)
{
    int rval, ncols;

    rval = GRBupdatemodel (lp->gurobi_lp);
    if (rval) {
        fprintf (stderr, "GRBupdatemodel failed\n");
        return 0;
    }

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_NUMVARS, &ncols);
    if (rval) {
        fprintf (stderr, "GRBgetintattr GRB_INT_ATTR_NUMVARS failed\n");
        return 0;
    }
    return ncols;
}

int CClp_nnonzeros (CClp *lp)
{
    int rval, nzs;

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_NUMNZS, &nzs);
    if (rval) {
        fprintf (stderr, "GRBgetintattr GRB_INT_ATTR_NUMNZS failed\n");
        return 0;
    }
    return nzs;
}

int CClp_status (CClp *lp, int *status)
{
    int solstat;
    int rval;

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_STATUS, &solstat);
    CCcheck_rval (rval, "GRBgetintattr GRB_INT_ATTR_STATUS failed");

    if (solstat == GRB_OPTIMAL) {
        *status = 0;
        return 0;
    } else if (solstat == GRB_INFEASIBLE) {
        *status = 1;
        return 0;
    } else {
        fprintf (stderr, "lp in an unknown state: %d\n",  solstat);
        *status = -1;
        return 1;
    }

CLEANUP:
    return rval;
}

int CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,
                    double *rmatval, double *weight)
{
    printf ("CClp_getweight ...\n"); fflush (stdout);

    if (lp || nrows || rmatbeg || rmatind || rmatval || weight) {
        fprintf (stderr, "getweight not set up yet\n"); return 1;
    } else {
        return 0;
    }
}

int CClp_dump_lp (CClp *lp, const char *fname)
{
    int rval = 0;
    char nambuf[128];

    sprintf (nambuf, "%s.lp", fname);
    rval = GRBwrite (lp->gurobi_lp, nambuf);
    CCcheck_rval (rval, "GRBwrite failed");

CLEANUP:
    return rval;
}

int CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,
                      double *downpen, double *uppen)
{
    /* Replace by something looking for non-degenerate pivots */

    int  rval = 0;
    int  ncols, i, j;
    double *x = (double *) NULL;

    printf ("QQQ CClp_getgoodlist ...\n"); fflush (stdout);

    *goodlen_p = 0;

    rval = GRBgetintattr(lp->gurobi_lp, GRB_INT_ATTR_NUMVARS, &ncols);
    CCcheck_rval (rval, "GRBgetintattr GRB_INT_ATTR_NUMVARS");
    if (ncols == 0) {
        fprintf (stderr, "No columns in LP\n");
        rval = 1; goto CLEANUP;
    }

    x = CC_SAFE_MALLOC (ncols, double);
    CCcheck_NULL (x, "out of memory for x");

    rval = GRBgetdblattrarray(lp->gurobi_lp, GRB_DBL_ATTR_X, 0, ncols, x);
    CCcheck_rval (rval, "GRBgetdblattrarray GRB_DBL_ATTR_X failed");

    for (i = 0, j = 0; i < ncols; i++) {
        if (x[i] >= CCtsp_INTTOL && x[i] <= 1.0 - CCtsp_INTTOL) {
            goodlist[j] = i;
            downpen[j]  = x[i];
            uppen[j]    = 1.0 - x[i];
            j++;
        }
    }

    *goodlen_p = j;

CLEANUP:

    CC_IFFREE (x, double);
    return rval;
}

int CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,
                       double *downpen, double *uppen, int iterations,
                       double upperbound)
{
    double oldcutoff, olditlim;
    int i, status, rval = 0;

    printf ("strongbranch (%d) ...\n", ncand); fflush (stdout);

#define USE_GRBstrongbranch
#ifdef  USE_GRBstrongbranch
    if (!lp || !candidatelist || ncand == 0 || !downpen || !uppen ||
         iterations == 0 || upperbound < 0) {
        fprintf (stderr, "CClp_strongbranch called with bad input\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetdblparam(GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_CUTOFF, 
                          &oldcutoff);
    CCcheck_rval (rval, "GRBgetdblparam failed");
   
    rval = GRBgetdblparam(GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_ITERATIONLIMIT,
                          &olditlim);
    CCcheck_rval (rval, "GRBgetdblparam failed");

    rval = GRBsetdblparam(GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_CUTOFF,
                          upperbound);
    CCcheck_rval (rval, "GRBsetdblparam failed");

    rval = GRBsetdblparam(GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_ITERATIONLIMIT,
                          (double) iterations);
    CCcheck_rval (rval, "GRBsetdblparam failed");
   

    /* status =  0, OK; 
              = -1, lp isn't a min problem; 
              =  1  basis changed, downpen = upperpen = current obj. value */

    rval = GRBstrongbranch(lp->gurobi_lp, ncand, candidatelist, downpen, 
                           uppen, &status);
    CCcheck_rval (rval, "GRBstrongbranch failed");

    rval = GRBsetdblparam(GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_CUTOFF,
                          oldcutoff);
    CCcheck_rval (rval, "GRBsetdblparam failed");

    rval = GRBsetdblparam(GRBgetenv(lp->gurobi_lp), GRB_DBL_PAR_ITERATIONLIMIT,
                          olditlim);
    CCcheck_rval (rval, "GRBsetdblparam failed");

    for (i=0; i<ncand; i++) {
        if (downpen[i] > upperbound) downpen[i] = upperbound;
        if (uppen[i] > upperbound) uppen[i] = upperbound;
    }
#else
    {
        double objval;

        rval = CClp_objval (lp, &objval);
        CCcheck_rval (rval, "CClp_objval failed");

        for (i=0; i<ncand; i++) {
            downpen[i] = objval;
            uppen[i]   = objval;
        }
    }
#endif

CLEANUP:
    return rval;
}

static void printerrorcode (int c)
{
    switch (c) {
    case GRB_ERROR_OUT_OF_MEMORY:
        printf ("Available memory was exhausted\n");
        break;
    case GRB_ERROR_NULL_ARGUMENT:
        printf ("NULL input value provided for a required argument\n");
        break;
    case GRB_ERROR_INVALID_ARGUMENT:
        printf ("An invalid value was provided for a routine argument\n");
        break;
    case GRB_ERROR_UNKNOWN_ATTRIBUTE:
        printf ("Tried to query or set an unknown attribute\n");
        break;
    case GRB_ERROR_DATA_NOT_AVAILABLE:
        printf ("Attempted to query or set an attribute that could\n");
        printf ("not be accessed at that time\n");
        break;
    case GRB_ERROR_INDEX_OUT_OF_RANGE:
        printf ("Tried to query or set an attribute, but one or more\n");
        printf ("of the provided indices (e.g., constraint index, variable \n");
        printf ("index) was outside the range of valid values\n");
        break;
    case GRB_ERROR_UNKNOWN_PARAMETER:
        printf ("Tried to query or set an unknown parameter\n");
        break;
    case GRB_ERROR_VALUE_OUT_OF_RANGE:
        printf ("Tried to set a parameter to a value that is outside\n");
        printf ("the parameter's valid range\n");
        break;
    case GRB_ERROR_NO_LICENSE:
        printf ("Failed to obtain a valid license\n");
        break;
    case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
        printf ("Attempted to solve a model that is larger than the\n");
        printf ("limit for a demo license\n");
        break;
    case GRB_ERROR_CALLBACK:
        printf ("Problem in callback\n");
        break;
    case GRB_ERROR_FILE_READ:
        printf ("Failed to read the requested file\n");
        break;
    case GRB_ERROR_FILE_WRITE:
        printf ("Failed to write the requested file\n");
        break;
    case GRB_ERROR_NUMERIC:
        printf ("Numerical error during requested operation\n");
        break;
    case GRB_ERROR_IIS_NOT_INFEASIBLE:
        printf ("Attempted to perform infeasibility analysis on a\n");
        printf ("feasible model\n"); 
        break;
    case GRB_ERROR_NOT_FOR_MIP:
        printf ("Requested operation not valid for a MIP model\n");
        break;
    case GRB_ERROR_OPTIMIZATION_IN_PROGRESS:
        printf ("Tried to query or modify a model while optimization\n");
        printf ("was in progress\n");
        break;
     default:
        printf ("Unknown error code\n");
    }
    fflush (stdout);
}
