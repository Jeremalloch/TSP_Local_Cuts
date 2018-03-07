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
/*               A PROGRAM TO CHECK A SEARCH TREE                          */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 19, 2007                                                 */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/*  This file currently has a copy of the restart reading code from         */
/*  TSP/bcontrol.c.  In the long run, this code should be shared instead    */
/*  of copied.                                                              */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "lp.h"

typedef struct tsp_bbnode {
    int id;
    int status;
    int workstatus;
    int numtentative;
    struct tsp_bbnode *prev;
    struct tsp_bbnode *next;
    struct tsp_bbnode *parent;
    struct tsp_bbnode *child0;
    struct tsp_bbnode *child1;
    struct tsp_tnode *tentative_nodes;
    struct tsp_tnode *tparent;
    double lowerbound;
    double t;
    double cputime;
    int number;
    CCtsp_branchobj *branch;
} tsp_bbnode;

typedef struct tsp_tnode {
    tsp_bbnode *parent;
    tsp_bbnode *child0;
    tsp_bbnode *child1;
} tsp_tnode;

typedef struct proofnode {
    int number; 
    int parent;
    int side;
    int child0;
    int child1;
    CCtsp_branchobj *branch;
    CCtsp_bigdual *exact_dual;
} proofnode;

#define BB_NEEDS_CUTTING             (1)
#define BB_NEEDS_TENTATIVE_CUTTING   (2)
#define BB_NEEDS_BRANCHING           (3)
#define BB_NEEDS_TENTATIVE_BRANCHING (4)
#define BB_DONE                      (5)
#define BB_IDLE              (1)
#define BB_WORKING           (2)
#define BB_PRUNED            (3)

#define TASK_WAIT_SECONDS  (60)


CC_PTRWORLD_ROUTINES (tsp_bbnode, tsp_bbnode_alloc, tsp_bbnode_bulkalloc,
        tsp_bbnode_free)
CC_PTRWORLD_LEAKS_ROUTINE (tsp_bbnode, tsp_bbnode_check_leaks, id, int)

static char *resfname   =  (char *) NULL;
static char *prooffname = (char *) NULL;
static int run_silent = 1;
static int check_files = 0;
static int printit = 0;
static int buildroot = 0;
static char *masterfname = (char *) NULL;
static char *copyfname = (char *) NULL;

static int
    check_tree (char *probname, int ncount, tsp_bbnode *bbnode,
        tsp_bbnode *parent, int  depth, int side, CCdatagroup *dat, int *ptour,
        double upbound, CCtsp_lpcuts *pool, char *copyit, CCtsp_lp **rootlp,
        CCtsp_lp *prooflp, proofnode *lproof, CCrandstate *rstate, int silent),
    dump_tree (tsp_bbnode *b),
    add_branch (CCtsp_lp *lp, CCtsp_branchobj *b, int side, int depth,
        CCrandstate *rstate),
    remove_branch (CCtsp_lp *lp, CCtsp_branchobj *b, int side, int depth),
    read_restart (char *restart_name, char **p_probname,
        tsp_bbnode **p_rootbbnode, double *p_upbound, int *p_ncount,
        int *p_bbcount, double *p_branchzeit, CCptrworld *bbnode_world),
    read_bbtree (FILE *f, tsp_bbnode **p_b, CCptrworld *bbnode_world),
    read_tentative_nodes (FILE *f, int count, tsp_tnode **list,
        tsp_bbnode *parent, CCptrworld *bbnode_world),
    save_proof (char *probname, int ncount, int bbcount, proofnode *lproof),
    copy_branchobj (CCtsp_branchobj *a, CCtsp_branchobj *b),
    prob_name (char *buf, int buflen, const char *f, int n),
    lp_value (CCtsp_lp *lp, double *val),
    addbad (CCtsp_lp *lp, CCrandstate *rstate),
    parseargs (int ac, char **av);

static void
    number_tree (tsp_bbnode *bbnode, int *cnt),
    print_tree (tsp_bbnode *bbnode, int depth),
    init_bbnode (tsp_bbnode *bbnode),
    free_tree (tsp_bbnode **bbnode, CCptrworld *bbnode_world),
    free_tree_branch (tsp_bbnode *bbnode),
    init_proofnode (proofnode *p),
    free_proofnode (proofnode *p),
    print_branchobj (CCtsp_branchobj *b),
    print_branch (CCtsp_lpclique *c, char sense, int rhs, int depth),
    usage (char *fname);


int main (int ac, char **av)
{
    int rval = 0, infeasible = 0;
    tsp_bbnode *rootbbnode  = (tsp_bbnode *) NULL;
    char *probname = (char *) NULL;
    int *ptour = (int *) NULL;
    double restart_upbound = 0.0;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    int ncount = 0;
    int master_ncount = 0;
    int i, cnt, bbcount = 0;
    double branchzeit = 0.0;
    CCptrworld bbnode_world;
    CCrandstate rstate;
    CCdatagroup dat;
    char buf[1024];
    int seed = 99;
    proofnode *lproof = (proofnode *) NULL;
    CCtsp_lp *rootlp = (CCtsp_lp *) NULL;
    CCtsp_lp *prooflp = (CCtsp_lp *) NULL;
    CCtsp_lp **plp;

    CCptrworld_init (&bbnode_world);
    CCutil_init_datagroup (&dat);
    
    rval = parseargs (ac, av);
    if (rval) return 1;

    if (check_files && masterfname == (char *) NULL) {
        fprintf (stderr, "Need masterfile for file check\n");
        rval = 1; goto CLEANUP;
    }

    if (prooffname && masterfname == (char *) NULL) {
        fprintf (stderr, "Need masterfile for to build proof from root LP\n");
        rval = 1; goto CLEANUP;
    }

    if (buildroot) plp = &rootlp;
    else           plp = (CCtsp_lp **) NULL;

    CCutil_sprand (seed, &rstate);

    rval = read_restart (resfname, &probname, &rootbbnode, &restart_upbound,
                         &ncount, &bbcount, &branchzeit, &bbnode_world);
    if (rval) {
        fprintf (stderr, "read_restart failed\n");
        goto CLEANUP;
    }

    printf ("Problem Name: %s\n", probname); fflush (stdout);
    printf ("Number of Nodes: %d\n", ncount); fflush (stdout);
    printf ("Upperbound: %.0f\n", restart_upbound); fflush (stdout);

    if (masterfname) {
        rval = CCutil_getmaster (masterfname, &master_ncount, &dat, &ptour);
        CCcheck_rval (rval, "CCutil_getmaster failed");

        if (master_ncount != ncount) {
            fprintf (stderr, "master does not match problem\n");
            rval = 1;  goto CLEANUP;
        }
    }

    if (prooffname) {
        printf ("Call CCtsp_init_lp ...\n"); fflush (stdout);
        rval = CCtsp_init_lp (&prooflp, probname, -1, prooffname,
                  ncount, &dat,
                  0, (int *) NULL, (int *) NULL,       /* elist  */
                  0, (int *) NULL, (int *) NULL, 0,    /* exlist */
                  ptour, restart_upbound,
                  (CCtsp_lpcuts *) NULL, (CCtsp_lpcuts *) NULL, run_silent,
                  &rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_init_lp failed");
        if (infeasible) {
            fprintf (stderr, "Initial LP is infeasible\n");
            rval = 1; goto CLEANUP;
        }
        printf ("DONE CCtsp_init_lp\n"); fflush (stdout);

        lproof = CC_SAFE_MALLOC (bbcount, proofnode);
        for (i = 0; i < bbcount; i++) init_proofnode (&lproof[i]);
    }

    rootbbnode->number = 0;
    cnt = 1;
    number_tree (rootbbnode, &cnt);

    if (printit) {
        print_tree (rootbbnode, 0);
    }

    if (check_files) {
        rval = CCtsp_init_cutpool (&master_ncount, (char *) NULL, &pool);
        CCcheck_rval (rval, "CCtsp_init_cutpool failed");

        rval =  check_tree (probname, ncount, rootbbnode, (tsp_bbnode *) NULL,
                          0, -1, &dat, ptour, restart_upbound, pool, copyfname,
                          plp, prooflp, lproof, &rstate, run_silent);
        CCcheck_rval (rval, "check_tree failed");

        printf ("Number of cuts: %d\n", pool->cutcount);
        fflush (stdout);
        sprintf (buf, "%s_total.pul", probname);
        rval = CCtsp_write_cutpool (ncount, buf, pool);
        CCcheck_rval (rval, "CCtsp_write_cutpool failed");
    }

    if (buildroot) {
        if (rootlp) {
            double objval, penalty;
            int nadded;
            CCtsp_edgegenerator eg;

            rval = CClp_opt (rootlp->lp, CClp_METHOD_DUAL, &infeasible);
            CCcheck_rval (rval, "CClp_opt failed\n");
            if (infeasible) {
                fprintf (stderr, "Infeasible root LP\n");
                rval = 1; goto CLEANUP;
            }
            rval = CCtsp_update_result (rootlp);
            CCcheck_rval (rval, "CCtsp_update_result failed\n");

            printf ("have the root\n"); fflush (stdout);
            rval = lp_value (rootlp, &objval);
            CCcheck_rval (rval, "lp_value failed");
            printf ("lp value: %f\n", objval); fflush (stdout);

            eg.ncount = 0;
            rval = CCtsp_init_edgegenerator (&eg, rootlp->graph.ncount,
                         rootlp->dat, rootlp->fulladj, 0, run_silent, &rstate);
            CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");

            rval = CCtsp_addbad_variables (rootlp, &eg, &penalty, &nadded,
                     CCtsp_PRICE_RCTHRESH, CCtsp_PRICE_MAXPENALTY, 0,
                     (int *) NULL, run_silent, &rstate);
            CCcheck_rval (rval, "CCtsp_addbad_variables failed");
            rval = lp_value (rootlp, &objval);
            CCcheck_rval (rval, "lp_value failed");

            printf ("added %d edges to root, penalty %f\n", nadded, penalty);
            printf ("Root LP Value: %f\n", objval); fflush (stdout);
            fflush (stdout);
            rootlp->lowerbound = objval + penalty;

            rval = CCtsp_update_result (rootlp);
            CCcheck_rval (rval, "CCtsp_update_result failed\n");

            sprintf (buf, "new_%s", probname);
            rval  = CCtsp_write_probroot_id (buf, rootlp);
            CCcheck_rval (rval, "CCtsp_write_probroot_id failed");

            if (eg.ncount) CCtsp_free_edgegenerator (&eg);
        } else {
            printf ("Warning: Did not pick up the root LP\n");
            fflush (stdout);
        }
    }

    rval = dump_tree (rootbbnode);
    CCcheck_rval (rval, "dump_tree failed");

    if (lproof) {
        rval = save_proof (probname, ncount, bbcount, lproof);
        CCcheck_rval (rval, "save_proof failed");
        sprintf (buf, "plp_%s", probname);
        rval  = CCtsp_write_probroot_id (buf, prooflp);
        CCcheck_rval (rval, "CCtsp_write_probroot_id failed");
    }

CLEANUP:

    CC_IFFREE (ptour, int);
    CC_IFFREE (probname, char);
    if (pool) {CCtsp_free_cutpool (&pool); }
    CCutil_freedatagroup (&dat);
    if (rootlp) CCtsp_free_tsp_lp_struct (&rootlp);
    if (prooflp) CCtsp_free_tsp_lp_struct (&prooflp);
    free_tree_branch (rootbbnode);
    free_tree (&rootbbnode, &bbnode_world);
    if (lproof) {
        for (i = 0; i < bbcount; i++) {
            free_proofnode (&lproof[i]);
        }
        CC_FREE (lproof, proofnode);
    }
    CCptrworld_delete (&bbnode_world);
    return rval;
}

static void number_tree (tsp_bbnode *bbnode, int *cnt)
{
    if (bbnode->child0 != (tsp_bbnode *) NULL) {
        bbnode->child0->number = *cnt;
        (*cnt)++;
    }
    if (bbnode->child1 != (tsp_bbnode *) NULL) {
        bbnode->child1->number = *cnt;
        (*cnt)++;
    }

    if (bbnode->child0 != (tsp_bbnode *) NULL) {
        number_tree (bbnode->child0, cnt);
    }

    if (bbnode->child1 != (tsp_bbnode *) NULL) {
        number_tree (bbnode->child1, cnt);
    }
}
    
static void print_tree (tsp_bbnode *bbnode, int depth)
{
    printf ("Node %d: Depth %d", bbnode->number, depth);
    if (bbnode->child0 == (tsp_bbnode *) NULL &&
        bbnode->child1 == (tsp_bbnode *) NULL) {
        printf (" Leaf\n");
    } else {
        if (bbnode->child0 != (tsp_bbnode *) NULL) {
            printf (" Child0 = %d", bbnode->child0->number);
        } else {
            printf (" Child0 = X");
        }
        
        if (bbnode->child1 != (tsp_bbnode *) NULL) {
            printf (" Child1 = %d", bbnode->child1->number);
        } else {
            printf (" Child1 = X");
        }
        printf ("\n");
    }


    if (bbnode->child0 != (tsp_bbnode *) NULL) {
        print_tree (bbnode->child0, depth+1);
    }

    if (bbnode->child1 != (tsp_bbnode *) NULL) {
        print_tree (bbnode->child1, depth+1);
    }
}

static int check_tree (char *probname, int ncount, tsp_bbnode *bbnode,
        tsp_bbnode *parent, int  depth, int side, CCdatagroup *dat, int *ptour,
        double upbound, CCtsp_lpcuts *pool, char *copyit, CCtsp_lp **rootlp,
        CCtsp_lp *prooflp, proofnode *lproof, CCrandstate *rstate, int silent)
{
    int rval = 0, infeasible = 0, i, k, kcnt, ncuts;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_lpcut_in c;
    char buf[1024];
    char copybuf[2048];
    double *pi = (double *) NULL;
    double objval;
    CCtsp_lprow cr;    

    CCtsp_init_lprow (&cr);  

    if ((bbnode->child0 == (tsp_bbnode *) NULL &&
         bbnode->child1 != (tsp_bbnode *) NULL)  ||
        (bbnode->child1 == (tsp_bbnode *) NULL &&
         bbnode->child0 != (tsp_bbnode *) NULL)) {
        printf ("Error: BBnode with only one child\n");
        rval = 1; goto CLEANUP;
    }

    if (depth == 0 && bbnode->number != 0) {
        printf ("Eroor: Root of tree does not have number 0\n");
        rval = 1; goto CLEANUP;
    }
      
    if (lproof) {
        k = bbnode->number;
        lproof[k].number = k;
        if (parent) lproof[k].parent = parent->number;
        else        lproof[k].parent = -1;
        if (bbnode->child0) lproof[k].child0 = bbnode->child0->number;
        else                lproof[k].child0 = -1;
        if (bbnode->child1) lproof[k].child1 = bbnode->child1->number;
        else                lproof[k].child1 = -1;
        lproof[k].side = side;
    }
 
    if (depth == 0)                                    printf ("Root ");
    else if (bbnode->child0 == (tsp_bbnode *) NULL)    printf ("Leaf ");
    else                                               printf ("Node ");

    printf ("%2d: ID=%5d ", bbnode->number, bbnode->id);

    rval = prob_name (buf, 1024, probname, bbnode->id);
    CCcheck_rval (rval, "prob_name failed");
    printf (" %s ", buf); fflush (stdout);

    if (copyit) {
        sprintf (copybuf, "cp -p %s %s/%s", buf, copyit, buf);
        i = system (copybuf);
        if (i == -1) {
            fprintf (stderr, "\nprob copy failed\n");
            rval = 1;  goto CLEANUP;
        }
    }

    rval = CCtsp_bb_init_lp (&lp, probname, bbnode->id, ncount, dat, ptour,
                         upbound, (CCtsp_lpcuts *) NULL,
                         (CCtsp_lpcuts *) NULL, silent, rstate, &infeasible);
    CCcheck_rval (rval, "CCtsp_bb_init_lp failed");
    if (infeasible) {
        fprintf (stderr, "init LP is infeasible\n");
        rval = 1; goto CLEANUP;
    }
    if (prooflp && (depth == 0 || (bbnode->child0 == (tsp_bbnode *) NULL))) {
        rval = lp_value (lp, &objval);
        CCcheck_rval (rval, "lp_value failed");
        printf ("\nNode LP Value: %f\n", objval); fflush (stdout);
    }

    if (lp->branchdepth != depth) {
        printf ("ERROR: Depth does not match branchdepth\n"); 
        rval = 1; goto CLEANUP;
    }

    if (lp->branchdepth >= 1) {
        print_branchobj (&lp->branchhistory[lp->branchdepth-1]);
    } else {
        printf ("Root "); fflush (stdout);
    }
    printf ("\n"); 

    if (parent && (parent->branch == (CCtsp_branchobj *) NULL)) {
        parent->branch = CC_SAFE_MALLOC (1, CCtsp_branchobj);
        CCcheck_NULL (parent->branch, "out of memory for parent->branch");
        rval = copy_branchobj (&(lp->branchhistory[lp->branchdepth-1]),
                               parent->branch);
        CCcheck_rval (rval, "copy_branchobj failed");
        if (lproof) {
            if (lproof[parent->number].branch) {
                printf ("Error: Already recorded a branch\n");
                rval = 1; goto CLEANUP;
            }
            lproof[parent->number].branch = CC_SAFE_MALLOC (1, CCtsp_branchobj);
            CCcheck_NULL (lproof[parent->number].branch,
                          "out of memory for parent->branch");
            rval = copy_branchobj (&(lp->branchhistory[lp->branchdepth-1]),
                                   lproof[parent->number].branch);
        }
    }

    if (depth == 0 || bbnode->child0 == (tsp_bbnode *) NULL) {
        ncuts = lp->cuts.cutcount;
        pi = CC_SAFE_MALLOC (ncount+ncuts, double);
        CCcheck_NULL (pi, "out of memory for pi");
        rval = CClp_pi (lp->lp, pi);
        CCcheck_rval (rval, "CClp_pi failed");
        for (i = 0; i < lp->cuts.cutcount; i++) {
            if (lp->cuts.cuts[i].branch == 0 && pi[i+ncount] > 0.0) {
                rval = CCtsp_lpcut_to_lpcut_in (&(lp->cuts),
                                                &(lp->cuts.cuts[i]), &c);
                CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

                kcnt = pool->cutcount;
                rval = CCtsp_add_to_cutpool_lpcut_in (pool, &c);
                CCcheck_rval (rval, "CCtsp_add_to_cutpool_lpcut_in failed");

                if (depth > 0 && rootlp != (CCtsp_lp **) NULL &&
                         pool->cutcount > kcnt) {
                    rval = CCtsp_add_cut (*rootlp, &c, &cr, 1);  
                    CCcheck_rval (rval, "CCtsp_add_cut failed");
                }

                CCtsp_free_lpcut_in (&c);
            }
        }
        CC_IFFREE (pi, double);
    }

    if (cr.rowcnt > 0) {
        printf ("Adding %d cuts to root lp\n", cr.rowcnt);
        fflush (stdout);
        rval = CCtsp_add_multiple_rows (*rootlp, &cr);
        CCcheck_rval (rval, "CCtsp_add_multiple_rows failed");
    }

    if (depth == 0 && rootlp != (CCtsp_lp **) NULL) {
        *rootlp = lp;
    } else {
        CCtsp_free_tsp_lp_struct (&lp);
    }
    lp = (CCtsp_lp *) NULL;


    if (prooflp && depth > 0) {
        if (!parent || !parent->branch || side == -1) {
            printf ("Error: No branchobj to add to LP\n");
            rval = 1;  goto CLEANUP;
        }

        rval = add_branch (prooflp, parent->branch, side, depth, rstate);
        CCcheck_rval (rval, "add_branch failed");
    }

    if (prooflp && (depth == 0 || (bbnode->child0 == (tsp_bbnode *) NULL))) {
        CCbigguy bound;

        rval = addbad (prooflp, rstate);
        CCcheck_rval (rval, "addbad failed");


        rval = lp_value (prooflp, &objval);
        CCcheck_rval (rval, "lp_value failed");
        printf ("Proof LP Value: %f\n", objval); fflush (stdout);

        CCtsp_exact_dual (prooflp);
        lproof[bbnode->number].exact_dual =
                      CC_SAFE_MALLOC (1, CCtsp_bigdual);
        CCcheck_NULL (lproof[bbnode->number].exact_dual, 
                      "out of memory for exact_dual");
        lproof[bbnode->number].exact_dual->cutcount =
                      prooflp->exact_dual->cutcount;
        lproof[bbnode->number].exact_dual->node_pi =
                      prooflp->exact_dual->node_pi;
        lproof[bbnode->number].exact_dual->cut_pi =
                      prooflp->exact_dual->cut_pi;

        rval = CCtsp_exact_price (prooflp, &bound, 0, 0, 0, 0, 1);
        CCcheck_rval (rval, "CCtsp_exact_price failed");
        printf ("Exactbound: %f\n", CCbigguy_bigguytod (bound));
        fflush (stdout);


        prooflp->exact_dual->node_pi = (CCbigguy *) NULL;
        prooflp->exact_dual->cut_pi = (CCbigguy *) NULL;
        CC_FREE (prooflp->exact_dual, CCtsp_bigdual);
    }


    if (bbnode->child0 != (tsp_bbnode *) NULL) {
        rval = check_tree (probname, ncount, bbnode->child0, bbnode, depth+1,
                        0, dat, ptour, upbound, pool, copyit, rootlp, prooflp,
                        lproof, rstate, silent);
        if (rval) goto CLEANUP;
    }

    if (bbnode->child1 != (tsp_bbnode *) NULL) {
        rval = check_tree (probname, ncount, bbnode->child1, bbnode, depth+1,
                        1, dat, ptour, upbound, pool, copyit, rootlp, prooflp,
                        lproof, rstate, silent);
        if (rval) goto CLEANUP;
    }

    if (prooflp && depth > 0) {
        rval = remove_branch (prooflp, parent->branch, side, depth);
        CCcheck_rval (rval, "remove_branch failed");
    }


CLEANUP:

    CC_IFFREE (pi, double);
    if (lp) CCtsp_free_tsp_lp_struct (&lp);
    CCtsp_free_lprow (&cr);

    return rval;
}

static int add_branch (CCtsp_lp *lp, CCtsp_branchobj *b, int side, int depth,
        CCrandstate *rstate)
{
    int rval = 0, infeasible = 0;
    int ar[2];
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;
    CCtsp_lprow cr;
    CCtsp_lpcut_in d;

    CCtsp_init_lpcut_in (&d);

    printf ("ADD (%d): ", depth); fflush (stdout);
    print_branchobj (b);
    printf ("\n"); fflush (stdout);

    c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    CCcheck_NULL (c, "out of memory for c");

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        rval = CCtsp_array_to_lpclique (ar, 2, c);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
        if (side == 1) {
            d.rhs = 2.0;
            d.sense = 'L';
        } else {
            d.rhs = 4.0;
            d.sense = 'G';
        }
    } else {
        if (!b->clique) {
            fprintf (stderr, "CCtsp_branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_copy_lpclique (b->clique, c);
        CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
        if (side == 0) {
            d.rhs = 2.0;
            d.sense = 'L';
        } else {
            d.rhs = 4.0;
            d.sense = 'G';
        }
    }

    d.cliquecount = 1;
    d.branch = 1;
    d.cliques = c;

    print_branch (c, d.sense, d.rhs, depth);

    rval = CCtsp_construct_skeleton (&d, lp->graph.ncount);
    CCcheck_rval (rval, "CCtsp_construct_skeleton failed");

    CCtsp_init_lprow (&cr);
    rval = CCtsp_add_cut (lp, &d, &cr, 0);
    CCcheck_rval (rval, "CCtsp_add_cut failed");

    rval = CCtsp_add_multiple_rows (lp, &cr);
    CCcheck_rval (rval, "CCtsp_add_multiple_rows failed");

    CCtsp_free_lprow (&cr);
    CCtsp_free_lpcut_in (&d);

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed");
    if (infeasible) {
        printf ("LP is infeasible, try to recover\n"); fflush (stdout);
        rval = CCtsp_infeas_recover (lp, 1, rstate, &infeasible);
        CCcheck_rval (rval, "CCtsp_infeas_recover failed");
        if (infeasible) {
            printf ("LP is really infeasible, try to recover\n");
            fflush (stdout);
            rval = 1; goto CLEANUP;
        }
    }

    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed\n");

CLEANUP:
    return rval;
}

static int remove_branch (CCtsp_lp *lp, CCtsp_branchobj *b, int side,
        int depth)
{
    int rval = 0, infeasible = 0;
    int ar[2];
    int k;
    char sense;
    int rhs, num;
    CCtsp_lpcut *cu;
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;

    printf ("REMOVE (%d): ", depth); fflush (stdout);
    print_branchobj (b);
    printf ("\n"); fflush (stdout);

    if (b->ends[0] != -1) {
        ar[0] = b->ends[0];
        ar[1] = b->ends[1];

        c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        CCcheck_NULL (c, "out of memory for c");

        rval = CCtsp_array_to_lpclique (ar, 2, c);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

        if (side == 1) {
            sense = 'L';
            rhs = 2;
        } else {
            sense = 'G';
            rhs = 4;
        }
    } else {
        if (!b->clique) {
            printf ("Error: Un-branchobj has no edge or clique\n");
            rval = 1; goto CLEANUP;
        }
        c = b->clique;
        if (side == 0) {
            sense = 'L';
            rhs = 2;
        } else {
            sense = 'G';
            rhs = 4;
        }
    }


    print_branch (c, sense, rhs, depth);

    num = lp->cuts.cutcount - 1;

    /* Check that last cut in the LP is really the branching cut */

    cu = &(lp->cuts.cuts[num]);
    if (cu->cliquecount != 1 || cu->sense != sense || cu->rhs != rhs) {
        printf ("Error: Last LP row does not match branch\n");
        printf ("Count = %d, Sense = %c, RHS = %d\n", cu->cliquecount,
                 cu->sense, cu->rhs);
        rval = 1;  goto CLEANUP;
    }

    CCtsp_lpclique_compare (&(lp->cuts.cliques[cu->cliques[0]]),c, &k);
    if (k != 0) {
        printf ("Error: Last LP row clique does not match branch\n");
        rval = 1;  goto CLEANUP;
    }
        
    /* Remove the last cut from the LP */

    rval = CCtsp_delete_cut (lp, num);
    CCcheck_rval (rval, "CCtsp_delete_cut failed");
    CCtsp_delete_cut_from_cutlist (&lp->cuts, num);

    rval = CClp_opt (lp->lp, CClp_METHOD_DUAL, &infeasible);
    CCcheck_rval (rval, "CClp_opt failed\n");
    if (infeasible) {
        fprintf (stderr, "Infeasible LP\n");
        rval = 0;  goto CLEANUP;
    }

    rval = CCtsp_update_result (lp);
    CCcheck_rval (rval, "CCtsp_update_result failed\n");

    if (b->ends[0] != -1) {
        CCtsp_free_lpclique (c);
        CC_FREE (c, CCtsp_lpclique);
    }

CLEANUP:
    return rval;
}

static int read_restart (char *restart_name, char **p_probname,
        tsp_bbnode **p_rootbbnode, double *p_upbound, int *p_ncount,
        int *p_bbcount, double *p_branchzeit, CCptrworld *bbnode_world)
{
    FILE *f = (FILE *) NULL;
    char *probname = (char *) NULL;
    int rval = 0;

    f = fopen (restart_name, "r");
    if (f == (FILE*) NULL) {
        perror (restart_name);
        fprintf (stderr, "Unable to open %s for input in read_restart\n",
                 restart_name);
        rval = 1; goto CLEANUP;
    }

    probname = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    if (probname == (char *) NULL) {
        fprintf (stderr, "Out of memory in read_restart\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_readstr (f, probname, CCtsp_PROB_FILE_NAME_LEN);
    
    rval = fscanf (f, "%d%lf%d%lf\n", p_ncount, p_upbound, p_bbcount,
            p_branchzeit);
    if (rval <= 0) {
        perror (restart_name);
        fprintf (stderr, "fscanf from %s failed\n", restart_name);
        rval = 1; goto CLEANUP;
    }
    rval = read_bbtree (f, p_rootbbnode, bbnode_world);
    if (rval) {
        fprintf (stderr, "read_bbtree failed\n"); goto CLEANUP;
    }
    
    rval = fclose (f);
    if (rval) {
        perror (restart_name);
        fprintf (stderr, "fclose %s failed\n", restart_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    *p_probname = probname;
    
    rval = 0;

  CLEANUP:
    if (f != (FILE *) NULL) {
        fclose (f);
    }
    if (rval) {
        CC_IFFREE (probname, char);
        free_tree (p_rootbbnode, bbnode_world);
    }
    return rval;
}

static int read_bbtree (FILE *f, tsp_bbnode **p_b, CCptrworld *bbnode_world)
{
    int rval = 0;
    int child0, child1;
    tsp_bbnode *b = (tsp_bbnode *) NULL;

    b = tsp_bbnode_alloc (bbnode_world);
    if (b == (tsp_bbnode *) NULL) {
        fprintf (stderr, "tsp_bbnode_alloc failed\n");
        rval = 1; goto CLEANUP;
    }
    init_bbnode (b);
    
    rval = fscanf (f, "%d%d%d%d%d%lf%lf\n", &b->status, &b->id,
                    &child0, &child1, &b->numtentative, &b->lowerbound,
                    &b->cputime);
    if (rval <= 0) {
        perror ("restart_file");
        fprintf (stderr, "fscanf failed reading restart file\n");
        rval = 1; goto CLEANUP;
    }
    b->workstatus = BB_IDLE;

    if (b->numtentative > 0) {
        rval = read_tentative_nodes (f, b->numtentative, &b->tentative_nodes,
                                     b, bbnode_world);
        if (rval) {
            fprintf (stderr, "read_tentative_nodes failed\n");
            goto CLEANUP;
        }
    }

    if (child0 != -1) {
        rval = read_bbtree (f, &(b->child0), bbnode_world);
        if (rval) goto CLEANUP;
        if (b->child0->id != child0) {
            fprintf (stderr, "syntax error in restart file\n");
            rval = 1; goto CLEANUP;
        }
        b->child0->parent = b;
    }
    if (child1 != -1) {
        rval = read_bbtree (f, &(b->child1), bbnode_world);
        if (rval) goto CLEANUP;
        if (b->child1->id != child1) {
            fprintf (stderr, "syntax error in restart file\n");
            rval = 1; goto CLEANUP;
        }
        b->child1->parent = b;
    }

    *p_b = b;
    rval = 0;

  CLEANUP:
    if (rval) {
        tsp_bbnode_free (bbnode_world, b);
    }
    return rval;
}

static int read_tentative_nodes (FILE *f, int count, tsp_tnode **list,
        tsp_bbnode *parent, CCptrworld *bbnode_world)
{
    int rval = 0;
    int i;
    int obtained = 0;
    tsp_tnode *s;
    tsp_bbnode *child0 = (tsp_bbnode *) NULL;
    tsp_bbnode *child1 = (tsp_bbnode *) NULL;

    *list = CC_SAFE_MALLOC (count, tsp_tnode);
    if (*list == (tsp_tnode *) NULL) {
        fprintf (stderr, "out of memory in read_tentative_nodes\n");
        rval = 1; goto CLEANUP;
    }

    for (obtained = 0; obtained < count; obtained++) {
        s = &((*list)[obtained]);
        child0 = tsp_bbnode_alloc (bbnode_world);
        if (child0 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tsp_bbnode_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        init_bbnode (child0);
        child1 = tsp_bbnode_alloc (bbnode_world);
        if (child1 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tsp_bbnode_alloc failed\n");
            tsp_bbnode_free (bbnode_world, child0);
            rval = 1; goto CLEANUP;
        }
        init_bbnode (child1);

        rval = fscanf (f, "%d%d%lf%lf %d%d%lf%lf\n",
                    &child0->status, &child0->id,
                    &child0->lowerbound, &child0->cputime,
                    &child1->status, &child1->id,
                    &child1->lowerbound, &child1->cputime);
        if (rval <= 0) {
            perror ("restart_file");
            fprintf (stderr, "fscanf failed reading tentative line\n");
            rval = 1; goto CLEANUP;
        }
        child0->tparent = s;
        child1->tparent = s;
        s->child0 = child0;
        s->child1 = child1;
        s->parent = parent;
    }
    rval = 0;

CLEANUP:

    if (rval) {
        for (i = 0; i < obtained; i++) {
            tsp_bbnode_free (bbnode_world, (*list)[i].child0);
            tsp_bbnode_free (bbnode_world, (*list)[i].child1);
        }
        CC_IFFREE (*list, tsp_tnode);
    }
    return rval;
}

static int dump_tree (tsp_bbnode *b)
{
    int rval = 0;

    printf ("Node %d: ", b->number);
    if ((b->child0 != (tsp_bbnode *) NULL &&
         b->child1 == (tsp_bbnode *) NULL) ||
        (b->child1 != (tsp_bbnode *) NULL &&
         b->child0 == (tsp_bbnode *) NULL)) {
        printf ("ERROR: Node with only one child\n");
        rval = 1; goto CLEANUP;
    }

    if (b->child0 != (tsp_bbnode *) NULL) {
        printf ("%d ",  b->child0->number);
    } else {
        printf ("-1 ");
    }
    if (b->child1 != (tsp_bbnode *) NULL) {
        printf ("%d ",  b->child1->number);
    } else {
        printf ("-1 ");
    }

    if (b->child1 != (tsp_bbnode *) NULL) {
        print_branchobj (b->branch);
    }
    printf ("\n"); fflush (stdout);


    if (b->child0 != (tsp_bbnode *) NULL) {
        rval = dump_tree (b->child0);
        if (rval) goto CLEANUP;
    }

    if (b->child1 != (tsp_bbnode *) NULL) {
        rval = dump_tree (b->child1);
        if (rval) goto CLEANUP;
    }

    
CLEANUP:

    return rval;
}

static int prob_name (char *buf, int buflen, const char *f, int n)
{
    int l = (int) strlen(f);
    int lastslash;
    int i;
    int d;

    if (l + 5 > buflen || n < 0) {
        fprintf (stderr, "Cannot generate filename for %s node %d\n",
                 f, n);
        return -1;
    }

    for (i = 0, lastslash = -1; i < l; i++) {
        if (f[i] == '/') lastslash = i;
        buf[i] = f[i];
    }
    if (l > lastslash + 9) l = lastslash + 9;
    for (i = lastslash+1; i < l; i++) {
        if (buf[i] == '.') buf[i] = '_';
    }
    if (n < 1000) {
        buf[l++] = '.';
        d = n/100;
        buf[l++] = '0' + ((unsigned int) d);
        n -= d*100;
        d = n/10;
        buf[l++] = '0' + d;
        n -= d*10;
        d = n;
        buf[l++] = '0' + ((unsigned int) d);
    } else if (n < 1000 + (26*36*36 - 5)) {
        buf[l++] = '.';
#define NAMESTRNUM(xa,xb,xc) (((xa)-'a') * 1296 + ((xb)-'a'+10) * 36 + \
                              ((xc)-'a'+10))
        n -= 1000;
        if (n >= NAMESTRNUM('m','a','s')) n++;
        if (n >= NAMESTRNUM('p','u','l')) n++;
        if (n >= NAMESTRNUM('r','e','s')) n++;
        if (n >= NAMESTRNUM('s','a','v')) n++;
        if (n >= NAMESTRNUM('s','o','l')) n++;
        d = n/1296;
        buf[l++] = 'a' + ((unsigned int) d);
        n -= d*1296;
        d = n/36;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*36;
        d = n;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
    } else if (n < 1000 + (26*36*36 - 5) + 36*36*36*36) {
        n -= 1000;
        n -= 26*36*36 - 5;
        d = n/(36*36*36);
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*36*36*36;
        d = n/(36*36);
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*(36*36);
        d = n/36;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
        n -= d*36;
        d = n;
        buf[l++] = (d < 10) ? '0' + ((unsigned int) d)
                            : 'a' + ((unsigned int) (d-10));
    } else {
        fprintf (stderr, "Node number %d too large\n", n);
        return -1;
    }
        
    buf[l] = '\0';
    return 0;
}

static void init_bbnode (tsp_bbnode *bbnode)
{
    bbnode->id         = 0;
    bbnode->lowerbound = 0.0;
    bbnode->status     = BB_NEEDS_CUTTING;
    bbnode->workstatus = BB_IDLE;
    bbnode->prev       = (tsp_bbnode *) NULL;
    bbnode->next       = (tsp_bbnode *) NULL;
    bbnode->parent     = (tsp_bbnode *) NULL;
    bbnode->child0     = (tsp_bbnode *) NULL;
    bbnode->child1     = (tsp_bbnode *) NULL;
    bbnode->cputime    = 0.0;
    bbnode->branch     = (CCtsp_branchobj *) NULL;
}

static void free_tree (tsp_bbnode **bbnode, CCptrworld *bbnode_world)
{
    if (!(*bbnode))  return;
    free_tree (&((*bbnode)->child0), bbnode_world);
    free_tree (&((*bbnode)->child1), bbnode_world);
    tsp_bbnode_free (bbnode_world, *bbnode);
    *bbnode = (tsp_bbnode *) NULL;
}

static void free_tree_branch (tsp_bbnode *bbnode)
{
    if (!bbnode)  return;
    free_tree_branch (bbnode->child0);
    free_tree_branch (bbnode->child1);
    if (bbnode->branch) {
        CCtsp_free_branchobj (bbnode->branch);
        CC_FREE (bbnode->branch, CCtsp_branchobj);
    }
}

static void init_proofnode (proofnode *p) 
{
    if (p) {
        p->number = -1;
        p->parent = -1;
        p->child0 = -1;
        p->child1 = -1;
        p->branch = (CCtsp_branchobj *) NULL;
        p->exact_dual = (CCtsp_bigdual *) NULL;
    }
}

static void free_proofnode (proofnode *p) 
{
    if (p) {
        if (p->branch) {
            CCtsp_free_branchobj (p->branch);
            CC_FREE (p->branch, CCtsp_branchobj);
        }
        if (p->exact_dual) {
            CCtsp_free_bigdual (&p->exact_dual);
        }
        init_proofnode (p);
    }
}

static void print_branchobj (CCtsp_branchobj *b)
{
    int i;

    printf ("Depth=%d ", b->depth);
    if (b->ends[0] != -1) {
        printf ("Edge (%d,%d)  = %d", b->ends[0], b->ends[1], b->rhs);
    } else {
        printf ("Clique ");
        for (i = 0; i < b->clique->segcount; i++) {
            printf ("%d->%d ", b->clique->nodes[i].lo, b->clique->nodes[i].hi);
        }
        if (b->sense == 'L') {
            printf ("<= %d", b->rhs);
        } else {
            printf (">= %d", b->rhs);
        }
    }
    fflush (stdout);
}

static void print_branch (CCtsp_lpclique *c, char sense, int rhs, int depth)
{
    int i;

    printf ("Depth=%d ", depth);
    printf ("Clique ");
    for (i = 0; i < c->segcount; i++) {
        printf ("%d->%d ", c->nodes[i].lo, c->nodes[i].hi);
    }
    if (sense == 'L') {
        printf ("<= %d\n", rhs);
    } else {
        printf (">= %d\n",rhs);
    }
    fflush (stdout);
}

static int copy_branchobj (CCtsp_branchobj *a, CCtsp_branchobj *b)
{
    int rval = 0;

    CCtsp_init_branchobj (b);
    b->depth = a->depth;
    b->rhs   = a->rhs;
    b->ends[0] = a->ends[0];
    b->ends[1] = a->ends[1];
    b->sense   = a->sense;
    if (a->clique) {
        b->clique = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        CCcheck_NULL (b->clique, "out of memory for b->clique");
        rval = CCtsp_copy_lpclique (a->clique, b->clique);
        CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
    }

CLEANUP:

    return rval;

}

static int addbad (CCtsp_lp *lp, CCrandstate *rstate)
{
    int rval = 0;
    int nadded;
    double penalty;
    CCtsp_edgegenerator eg;

    eg.ncount = 0;
    rval = CCtsp_init_edgegenerator (&eg, lp->graph.ncount,
                     lp->dat, lp->fulladj, 0, 1, rstate);
    CCcheck_rval (rval, "CCtsp_init_edgegenerator failed");

    rval = CCtsp_addbad_variables (lp, &eg, &penalty, &nadded,
             CCtsp_PRICE_RCTHRESH, CCtsp_PRICE_MAXPENALTY, 0,
             (int *) NULL, 1, rstate);
    CCcheck_rval (rval, "CCtsp_addbad_variables failed");

    printf ("add %d edges\n", nadded); fflush (stdout);

CLEANUP:

    if (eg.ncount) CCtsp_free_edgegenerator (&eg);
    return rval;
}


static int save_proof (char *probname, int ncount, int bbcount,
        proofnode *lproof)
{
    int rval = 0;
    int i, j;
    CC_SFILE *f = (CC_SFILE *) NULL;
    CCtsp_lpclique *c;
    char buf[1024];
    CCtsp_bigdual *d;


    sprintf (buf, "%s.proof", probname);
    f = CCutil_sopen (buf, "w");
    if (!f) {
        fprintf (stderr, "could not safe open %s for writing\n", buf);
        rval =1; goto CLEANUP;
    }  

    rval = CCutil_swrite_int (f, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (f, bbcount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < bbcount; i++) {
        rval = CCutil_swrite_int (f, lproof[i].number);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].parent);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].side);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].child0);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (f, lproof[i].child1);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        if (lproof[i].branch) {
            rval = CCutil_swrite_int (f, 1);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            rval = CCutil_swrite_int (f, lproof[i].branch->ends[0]);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            rval = CCutil_swrite_int (f, lproof[i].branch->ends[1]);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            if (lproof[i].branch->clique) {
                c = lproof[i].branch->clique;
                rval = CCutil_swrite_int (f, c->segcount);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
                for (j = 0; j < c->segcount; j++) {
                    rval = CCutil_swrite_int (f, c->nodes[j].lo);
                    CCcheck_rval (rval, "CCutil_swrite_int failed");
                    rval = CCutil_swrite_int (f, c->nodes[j].hi);
                    CCcheck_rval (rval, "CCutil_swrite_int failed");
                }
            } else {
                rval = CCutil_swrite_int (f, 0);
                CCcheck_rval (rval, "CCutil_swrite_int failed");
            }
        } else {
            rval = CCutil_swrite_int (f, 0);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
        if (lproof[i].exact_dual) {
            d = lproof[i].exact_dual;
            rval = CCutil_swrite_int (f, d->cutcount);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
            for (j = 0; j < ncount; j++) {
                rval = CCbigguy_swrite (f, d->node_pi[j]);
                CCcheck_rval (rval, "CCbigguy_swrite failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = CCbigguy_swrite (f, d->cut_pi[j]);
                CCcheck_rval (rval, "CCbigguy_swrite failed");
            }
        } else {
            rval = CCutil_swrite_int (f, 0);
            CCcheck_rval (rval, "CCutil_swrite_int failed");
        }
    }

CLEANUP:

    if (f) CCutil_sclose (f);

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

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "C:fM:prR:v?", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'C':
            copyfname = boptarg;
            check_files = 1;
            break;
        case 'f':
            check_files = 1;
            break;
        case 'M':
            masterfname = boptarg;
            break;
        case 'p':
            printit = 1;
            break;
        case 'r':
            buildroot = 1;
            check_files = 1;
            break;
        case 'R':
            prooffname = boptarg;
            break;
        case 'v':
            run_silent = 0;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }

    resfname = av[boptind++];

    if (boptind < ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] restart_fname\n", fname);
    fprintf (stderr, "   -C f  copy files to subdirectory f [needs -M]\n");
    fprintf (stderr, "   -f    check the files [needs -M]\n");
    fprintf (stderr, "   -M f  master file\n");
    fprintf (stderr, "   -p    print the tree nodes\n");
    fprintf (stderr, "   -r    build a root LP [needs -M]\n");
    fprintf (stderr, "   -R f  root LP to build proof [needs -M]\n");
    fprintf (stderr, "   -v    verbose\n");
}

