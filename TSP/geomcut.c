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
/*                     Cuts from cloud of x-vectors                         */
/*                                                                          */
/*                               TSP CODE                                   */
/*                                                                          */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: July 7, 2016                                                      */
/*                                                                          */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int CCtsp_geom_cuts (CCtsp_lpcut_in **cuts, int *cutcount,              */
/*      int ncount, int ecount, int *elist, double *x, double *fullzeit)    */
/*    FINDS violated subtour inequalities via connectivity.                 */
/*     -cuts will return any new cuts found (they will be added to the      */
/*      head of the linked list)                                            */
/*     -cutcount will return the number of new cuts added                   */
/*     -ncount is the number of nodes                                       */
/*     -ecount is the number of edges                                       */
/*     -elist contains the LP edges in node node format                     */
/*     -x is an LP solution                                                 */
/*     -fullzeit if not NULL returns the total run time in seconds          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "tsp.h"
#include "cut.h"
#include "combs.h"
#include "localcut.h"

static int lp_value (CCtsp_lp *lp, double *val);

int CCtsp_geom_cuts (CCtsp_lpcut_in **cuts, int *cutcount, CCtsp_lp *lp, 
        int xcount, int *xlist, double *x, int *snowtour, double *fullzeit,
        CCrandstate *rstate)
{
    int rval = 0, infeasible, i;
    CClp *newlp = (CClp *) NULL, *currentlp;
    CClp_info *binfo = (CClp_info *) NULL;
    double szeit, val, xval;

    printf ("CCtsp_geom_cuts ...\n"); fflush (stdout);

    if (!snowtour) {
        printf ("really should have a snowtour!\n");
        rval = 1; goto CLEANUP;
    }

    /* build a new LP with perturbed costs */

    rval = CClp_init (&newlp);
    CCcheck_rval (rval, "CClp_init failed");

    /* we will change the LP and steepest-edge norms will not be valid */
    /* so build just a basis via a call to CClp_get_info               */

    rval = CClp_get_info (lp->lp, &binfo);
    CCcheck_rval (rval, "CClp_get_info failed");

    rval = CClp_build_warmstart (&lp->warmstart, binfo);
    CCcheck_rval (rval, "CClp_build_warmstart failed");
    CClp_free_info (&binfo);

    currentlp = lp->lp;
    lp->lp = newlp;
    rval = CCtsp_load_lp (lp, rstate, 1, 0);
    CCcheck_rval (rval, "CCtsp_load_lp failed");

    szeit = CCutil_zeit();
    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed");
    if (infeasible) {
        printf ("new LP is infeasible!\n");
        rval = 1; goto CLEANUP;
    }
    printf ("Optimization Time: %.2f\n", CCutil_zeit() - szeit);
    fflush (stdout);

    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed");

    rval = lp_value (lp, &val);
    CCcheck_rval (rval, "lp_value failed");

    printf ("LP Value: %f\n", val);
    fflush (stdout);

    /* grab x from perturbed LP and use it to find cuts */

    /* if x was objective value greater than snowtour, reject the vector */

    for (i = 0; i < lp->graph.ecount; i++) {

    }

    /* restore the actual lp */

    lp->lp = currentlp;

    /* Note: should save old warmstart and reload it into lp */

CLEANUP:
    if (newlp) CClp_freelp (&newlp);
    CClp_free_info (&binfo);
    return rval;
}

static int lp_value (CCtsp_lp *lp, double *val)
{
    int rval = 0;

    rval = CCtsp_get_lp_result (lp, val, (double *) NULL, (int *) NULL,
                 (int **) NULL, (double **) NULL, (double **) NULL,
                 (double **) NULL, (double **) NULL);
    CCcheck_rval (rval, "CCtsp_get_lp_result failed");

CLEANUP:
    return rval;
}


