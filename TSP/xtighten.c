/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2005 by David Applegate, Robert Bixby,              */
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
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_xtighten_lpcut_in (CCtsp_lpgraph *g, CCtsp_lpcut_in *c,       */
/*      double *x, CCtsp_lpcut_in *d, double *pimprove)                     */
/*    Use flow methods to carry extreme-tighten procedure.                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "verify.h"

typedef struct clique {
    int mark;
    int atomcount;
    int *atomlist;
} clique;

typedef struct atom {
    int mark;
    int cliquecount;
    int *cliquelist;
    int nodecount;
    int *nodes;
} atom;

typedef struct cutinfo {
    int ncount;
    int atomcount;
    atom *atoms;
    int cliquecount;
    clique *cliques;
    int **smatrix;
} cutinfo;


static void atom_init (atom *a);
static void atom_free (atom *a);
static void cutinfo_init (cutinfo *A);
static void cutinfo_free (cutinfo *A);
static int cutinfo_build (cutinfo *A, CCtsp_lpcut_in *c, CCtsp_lpgraph *lg);
static int find_atoms (cutinfo *A, CCtsp_lpcut_in *c);
static int build_smatrix (cutinfo *A);
static int xtighten_cut (cutinfo *A, CCtsp_lpgraph *g, double *x, double *win);
static int xtighten_border (cutinfo *A, CCtsp_lpgraph *g, double *x, int atom1,
    int atom2, int *labels, double *improve);
static int nonempty_cut (int ncount, int ecount, int *elist, double *ecap,
        int s, int t, double *cval, int **cut, int *cutcount);
static int gather_cut (cutinfo *A, CCtsp_lpcut_in *old, CCtsp_lpcut_in *ctry);


#define XTIGHT_EPS 0.0001

static void clique_init (clique *c)
{
    if (c) {
        c->mark = 0;
        c->atomcount = 0;
        c->atomlist = (int *) NULL;
    }
}

static void clique_free (clique *c)
{
    if (c) {
        CC_IFFREE (c->atomlist, int);
        clique_init (c);
    }
}

static void atom_init (atom *a)
{
    if (a) {
        a->mark = 0;
        a->cliquecount = 0;
        a->cliquelist = (int *) NULL;
        a->nodecount = 0;
        a->nodes = (int *) NULL;
    }
}

static void atom_free (atom *a)
{
    if (a) {
        CC_IFFREE (a->cliquelist, int);
        CC_IFFREE (a->nodes, int); 
        atom_init (a);
    }
}

static void cutinfo_init (cutinfo *A)
{
    if (A) {
        A->atomcount = 0;
        A->atoms = (atom *) NULL;
        A->cliquecount = 0;
        A->cliques = (clique *) NULL;
        A->smatrix = (int **) NULL;
    }
}

static void cutinfo_free (cutinfo *A)
{
    int i;

    if (A) {
        if (A->atoms) {
            for (i = 0; i < A->atomcount; i++) {
                atom_free (&(A->atoms[i]));
            }
            CC_IFFREE (A->atoms, atom);
        }
        if (A->cliques) {
            for (i = 0; i < A->cliquecount; i++) {
                clique_free (&(A->cliques[i]));
            }
            CC_IFFREE (A->cliques, clique);
        }
        if (A->smatrix) {
            for (i = 0; i < A->atomcount; i++) {
                CC_IFFREE (A->smatrix[i], int);
            }
            CC_IFFREE (A->smatrix, int *);
        }
        cutinfo_init (A);
    }
}

static int cutinfo_build (cutinfo *A, CCtsp_lpcut_in *cut, CCtsp_lpgraph *lg)
{
    int rval = 0;

    if (cut->skel.atomcount == 0) {
        fprintf (stderr, "no atoms\n");  rval = 1;  goto CLEANUP;
    }

    A->ncount = lg->ncount;
    rval = find_atoms (A, cut);
    CCcheck_rval (rval, "find_atoms failed");

    rval = build_smatrix (A);
    CCcheck_rval (rval, "build_smatrix failed");

CLEANUP:

    return rval;
}

static int find_atoms (cutinfo *A, CCtsp_lpcut_in *c)
{
    /* Note: c->skeleton may be missing some atoms due to tighten.  */
    /* To get a clean implementation of xtighten we will use a full */
    /* set of atoms (adopting code from CCtsp_construct_skeleton).  */

    int rval = 0;
    int cliquecount = c->cliquecount;
    CCtsp_lpclique *cliques = c->cliques;
    int *label = (int *) NULL;
    int atomcount;
    int atomcount_save;
    int *atomsize = (int *) NULL;
    int *atomnew = (int *) NULL;
    int *atomwork = (int *) NULL;
    int *cnodes = (int *) NULL;
    int *marks = (int *) NULL;
    int i, j, tmp, ccount, countsav, cnt;
    int ncount = A->ncount;
    atom *atoms;

    /* List the nodes in each atom */

    label = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (label, "out of memory for label");
    for (i = 0; i < ncount; i++) label[i] = 0;

    atomsize = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (atomsize, "out of memory for atomsize");
    atomnew = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (atomnew, "out of memory for atomnew");
    atomwork = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (atomwork, "out of memory for atomwork");

    ccount = 0;
    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            if (label[j] == 0) {
                label[j] = 1;
                ccount++;
            }
        }
    }
    cnodes = CC_SAFE_MALLOC (ccount, int);
    CCcheck_NULL (cnodes, "out of memory for cnodes");
    ccount = 0;
    for (i=0; i<cliquecount; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            if (label[j] == 1) {
                label[j] = 0;
                cnodes[ccount++] = j;
            }
        }
    }


    /* refine atoms */
    atomsize[0] = ccount;
    atomcount = 1;
    for (i=0; i<cliquecount; i++) {
        for (j=0; j<atomcount; j++) {
            atomwork[j] = 0;
        }
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            atomwork[label[j]]++;
        }
        atomcount_save = atomcount;
        for (j=0; j<atomcount_save; j++) {
            if (atomwork[j] == 0) {
                atomnew[j] = -1;
            } else if (atomwork[j] == atomsize[j]) {
                atomnew[j] = j;
            } else {
                atomsize[atomcount] = atomwork[j];
                atomsize[j] -= atomwork[j];
                atomnew[j] = atomcount;
                atomcount++;
            }
        }
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            label[j] = atomnew[label[j]];
        }
    }


    countsav = atomcount;
    if (ccount < ncount) {   /* Have an outside node */
        atomcount += 1;
    }

    A->atoms = CC_SAFE_MALLOC (atomcount, atom);
    CCcheck_NULL (A->atoms, "out of memory for atoms");
    for (i = 0; i < atomcount; i++) {
        atom_init (&(A->atoms[i]));
    }
    atoms = A->atoms;

    for (i = 0; i < ccount; i++) {
        atoms[label[cnodes[i]]].nodecount++;
    }

    for (j = 0; j < countsav; j++) {
        if (atoms[j].nodecount <= 0) {
            fprintf (stderr, "Empty atom: %d\n", j);
            rval = 1;  goto CLEANUP;
        }
        atoms[j].nodes = CC_SAFE_MALLOC (atoms[j].nodecount, int);
        CCcheck_NULL (atoms[j].nodes, "out of memory for atom nodes");
        atoms[j].nodecount = 0;
    }
    for (i = 0; i < ccount; i++) {
        j = label[cnodes[i]];
        atoms[j].nodes[atoms[j].nodecount++] = cnodes[i];
    }

    if (atomcount > countsav) {
/*
        printf ("Outside Atom: %d nodes\n", ncount - ccount);
        fflush (stdout);
*/
        marks = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (marks, "out of memory for marks");

        for (i = 0; i < ncount; i++) marks[i] = 0;
        for (i = 0; i < ccount; i++) {
            marks[cnodes[i]] = 1;
        }
        atoms[atomcount-1].nodes = CC_SAFE_MALLOC (ncount-ccount, int);
        j = 0;
        for (i = 0; i < ncount; i++) {
            if (marks[i] == 0) {
                atoms[atomcount-1].nodes[j++] = i;
            }
        }
        atoms[atomcount-1].nodecount = ncount - ccount;
    }
    A->atomcount = atomcount;

/*
    if (atomcount > countsav) {
        printf ("Outside Atom: %d nodes\n", atoms[atomcount-1].nodecount);
        fflush (stdout);
    }
*/


    /* Now list the atoms in each clique */

    A->cliques = CC_SAFE_MALLOC (cliquecount, clique);
    CCcheck_NULL (A->cliques, "out of memory for cliques");
    for (i = 0; i < cliquecount; i++) {
        clique_init (&(A->cliques[i]));
    }


    for (i = 0; i < cliquecount; i++) {
        for (j = 0; j < atomcount; j++) {
            atomwork[j] = 0;
        }
        cnt = 0;
        CC_FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp) {
            if (atomwork[label[j]] == 0) {
                cnt++;
                atomwork[label[j]] = 1;
            }
        }
        if (cnt == 0) {
            fprintf (stderr, "clique has no atoms");
            rval = 1; goto CLEANUP;
        }
        A->cliques[i].atomlist = CC_SAFE_MALLOC (cnt, int);
        CCcheck_NULL (A->cliques[i].atomlist, "out of memory for clique atoms");
        A->cliques[i].atomcount = cnt;
        cnt = 0;
        for (j = 0; j < atomcount; j++) {
            if (atomwork[j] == 1) {
                A->cliques[i].atomlist[cnt++] = j;
            }
        }
        if (cnt != A->cliques[i].atomcount) {
            fprintf (stderr, "dropped an atom from a clique");
            rval = 1; goto CLEANUP;
        }
    }
    A->cliquecount = cliquecount;

/*
    for (i = 0; i < A->cliquecount; i++) {
        printf ("Clique %d:", i);
        for (j = 0; j < A->cliques[i].atomcount; j++) {
            printf (" %d", A->cliques[i].atomlist[j]);
        }
        printf ("\n");
    }
*/

CLEANUP:

    CC_IFFREE (atomwork, int);
    CC_IFFREE (atomnew, int);
    CC_IFFREE (atomsize, int);
    CC_IFFREE (cnodes, int);
    CC_IFFREE (label, int)
    CC_IFFREE (marks, int)

    return rval;
}

static int build_smatrix (cutinfo *A)
{
    int rval = 0;
    int atomcount = A->atomcount;
    int i, j, k, ato;
    int **smatrix;
    int *marks = (int *) NULL;

    A->smatrix = CC_SAFE_MALLOC (atomcount, int *);
    CCcheck_NULL (A->smatrix, "out of memory for smatrix");
    smatrix = A->smatrix;

    for (i = 0; i < atomcount; i++) {
        smatrix[i] = (int *) NULL;
    }
    for (i = 0; i < atomcount; i++) {
        smatrix[i] = CC_SAFE_MALLOC (atomcount, int);
        CCcheck_NULL (A->smatrix[i], "out of memory for smatrix index");
    }
    for (i = 0; i < atomcount; i++) {
        for (j = 0; j < atomcount; j++) {
            smatrix[i][j] = 0;
        }
    }

    marks = CC_SAFE_MALLOC (atomcount, int);
    CCcheck_NULL (marks, "out of memory for marks");
    for (j = 0; j < atomcount; j++) {
        marks[j] = 0;
    }

    for (i = 0; i < A->cliquecount; i++) {
        for (j = 0; j < A->cliques[i].atomcount; j++) {
            marks[A->cliques[i].atomlist[j]] = 1;
        }
        for (j = 0; j < A->cliques[i].atomcount; j++) {
            ato = A->cliques[i].atomlist[j];
            for (k = 0; k < atomcount; k++) {
                if (marks[k] == 0) {
                    smatrix[ato][k]++;
                    smatrix[k][ato]++;
                }
            }
        }
        for (j = 0; j < A->cliques[i].atomcount; j++) {
            marks[A->cliques[i].atomlist[j]] = 0;
        }
    }

/*
    printf ("Smatrix\n");
    for (i = 0; i < atomcount; i++) {
        for (j = 0; j < atomcount; j++) {
            printf ("%d ", smatrix[i][j]);
        }
        printf ("\n");
    }
*/
    

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

static int xtighten_cut (cutinfo *A, CCtsp_lpgraph *g, double *x, double *win)
{
    int rval = 0;
    int atomcount = A->atomcount;
    int ncount = g->ncount;
    int i, j;
    int *labels = (int *) NULL;
    double improve;

    *win = 0.0;

    if (atomcount < 3) {
        fprintf (stderr, "trying to xtighten a cut with only two atoms\n");
        rval = 1;  goto CLEANUP;
    }

    labels = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (labels, "out of memory for labels");
    for (i = 0; i < ncount; i++) labels[i] = -1;

    for (i = 0; i < atomcount; i++) {
        for (j = 0; j < A->atoms[i].nodecount; j++) {
            if (labels[A->atoms[i].nodes[j]] != -1) {
                fprintf (stderr, "atoms not disjoint\n");
                rval = 1;  goto CLEANUP;
            }
            labels[A->atoms[i].nodes[j]] = i;
        }
    }
    for (i = 0; i < ncount; i++) {
        if (labels[i] == -1) {
            fprintf (stderr, "atoms do not cover all nodes\n");
            rval = 1;  goto CLEANUP;
        }
    }

    for (i = 0; i < atomcount; i++) {
        for (j = i+1; j < atomcount; j++) {
            if (A->smatrix[i][j] == 1) {
                rval = xtighten_border (A, g, x, i, j, labels, &improve);
                CCcheck_rval (rval, "xtighten_border failed");
                (*win) += improve;
                if (improve > 0.0) goto CLEANUP;
            }
        }
    }

/*
    if (*win > 0) {
        printf ("Improved cut by %f\n", *win); fflush (stdout);
    }
*/

CLEANUP:

    CC_IFFREE (labels, int);
    return rval;

}

static int xtighten_border (cutinfo *A, CCtsp_lpgraph *g, double *x, int atom0,
        int atom1, int *labels, double *improve)
{
    int rval = 0;
    int ncount = g->ncount;
    int i, j, k, n, m, s, t, scnt, tcnt, fcount, fncount;
    int *marks = (int *) NULL;
    double cap0, cap1, cval, val = 0.0;
    int **smatrix = A->smatrix;
    atom *ato[2];
    int *numbers = (int *) NULL;
    int *invnumbers = (int *) NULL;
    int *flist = (int *) NULL;
    double *fcap = (double *) NULL;
    int cutcount;
    int *cut = (int *) NULL;
    int *tcut = (int *) NULL;
    int *scut = (int *) NULL;
    int *side0 = (int *) NULL, *side1 = (int *) NULL;
    int side0count = 0, side1count = 0;

    *improve = 0.0;

    ato[0] = &(A->atoms[atom0]);
    ato[1] = &(A->atoms[atom1]);

    marks = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    numbers = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (marks, "out of memory for numbers");
    invnumbers = CC_SAFE_MALLOC (ncount+2, int);
    CCcheck_NULL (marks, "out of memory for invnumbers");

    fncount = 2;
    invnumbers[0] = -2;
    invnumbers[1] = -1;
    for (k = 0; k < 2; k++) {
        for (i = 0; i < ato[k]->nodecount; i++) {
            invnumbers[fncount] = ato[k]->nodes[i];
            numbers[ato[k]->nodes[i]] = fncount++;
        }
    }

    if (fncount == 4) goto CLEANUP;

    for (i = 0; i < ato[1]->nodecount; i++) {
        marks[ato[1]->nodes[i]] = 1;
    }

    for (i = 0; i < ato[0]->nodecount; i++) {
        n = ato[0]->nodes[i];
        for (j = 0; j < g->nodes[n].deg; j++) {
            if (marks[g->nodes[n].adj[j].to] == 1) {
                val += x[g->nodes[n].adj[j].edge];
            }
        }
    }

    for (i = 0; i < ato[0]->nodecount; i++) {
        marks[ato[0]->nodes[i]] = 1;
    }
    for (k = 0; k < 2; k++) {
        if (k == 0) t = atom0;
        else        t = atom1;
        for (i = 0; i < ato[k]->nodecount; i++) {
            n = ato[k]->nodes[i];
            for (j = 0; j < g->nodes[n].deg; j++) {
                m = g->nodes[n].adj[j].to;
                if (marks[m] == 0) {
                    val += (smatrix[t][labels[m]] * x[g->nodes[n].adj[j].edge]);
                }
            }
        }
    }

    fcount = 0; 
    for (k = 0; k < 2; k++) {
        for (i = 0; i < ato[k]->nodecount; i++) {
            n = ato[k]->nodes[i];
            for (j = 0; j < g->nodes[n].deg; j++) {
                m = g->nodes[n].adj[j].to;
                if ((n < m) && (marks[m] == 1)) {
                    fcount++;
                }
            }
        }
    }
    fcount += (2 * (ato[0]->nodecount + ato[1]->nodecount));

    flist = CC_SAFE_MALLOC (2*fcount, int);
    CCcheck_NULL (flist, "out of memory for flist");
    fcap = CC_SAFE_MALLOC (fcount, double);
    CCcheck_NULL (fcap, "out of memory for fcap");

    fcount = 0;
    for (k = 0; k < 2; k++) {
        for (i = 0; i < ato[k]->nodecount; i++) {
            n = ato[k]->nodes[i];
            for (j = 0; j < g->nodes[n].deg; j++) {
                m = g->nodes[n].adj[j].to;
                if ((n < m) && (marks[m] == 1)) {
                    flist[2*fcount] = numbers[n];
                    flist[2*fcount+1] = numbers[m];
                    fcap[fcount] = x[g->nodes[n].adj[j].edge];
                    fcount++;
                }
            }
        }
    }
    for (k = 0; k < 2; k++) {
        for (i = 0; i < ato[k]->nodecount; i++) {
            n = ato[k]->nodes[i];
            cap0 = cap1 = 0.0;
            for (j = 0; j < g->nodes[n].deg; j++) {
                m = g->nodes[n].adj[j].to;
                if (marks[m] == 0) {
                    cap0 += (smatrix[atom0][labels[m]] *
                            x[g->nodes[n].adj[j].edge]);
                    cap1 += (smatrix[atom1][labels[m]] *
                            x[g->nodes[n].adj[j].edge]);
                }
            }
            flist[2*fcount] = numbers[n];
            flist[2*fcount+1] = 0;
            fcap[fcount] = cap1;  /* on 0 side of cut pay cost to node 1 */
            fcount++;
            flist[2*fcount] = numbers[n];
            flist[2*fcount+1] = 1;
            fcap[fcount] = cap0;  /* on 1 side of cut pay cost to node 0 */
            fcount++;
        }
    }
   
    if (ato[0]->nodecount > ato[1]->nodecount) {
        s = 0;  t = 1;
    } else {
        s = 1;  t = 0;
    }

    rval = nonempty_cut (fncount, fcount, flist, fcap, s, t, &cval, &cut,
                         &cutcount);
    CCcheck_rval (rval, "nonempty_cut failed");

    if (val - cval > XTIGHT_EPS && cutcount > 1 && cutcount < fncount-1) {
        scut = CC_SAFE_MALLOC (ato[0]->nodecount + ato[1]->nodecount, int);
        CCcheck_NULL (scut, "out of memory for scut");
        scnt = 0;
        tcut = CC_SAFE_MALLOC (ato[0]->nodecount + ato[1]->nodecount, int);
        CCcheck_NULL (tcut, "out of memory for tcut");
        tcnt = 0;
        for (i = 0; i < ncount; i++) marks[i] = 0;
        for (i = 0; i < cutcount; i++) {
            if (cut[i] > 1) {
                scut[scnt++] = invnumbers[cut[i]];
                marks[invnumbers[cut[i]]] = 1;
            }
        }
        for (k = 0; k < 2; k++) {
            for (i = 0; i < ato[k]->nodecount; i++) {
                if (marks[ato[k]->nodes[i]] == 0) {
                    tcut[tcnt++] = ato[k]->nodes[i];
                }
            }
        }
        for (i = 0; i < cutcount; i++) {
            if (cut[i] < 2) {
                if (cut[i] == 0) {
                    side0 = scut;  side1 = tcut;
                    side0count = scnt;  side1count = tcnt;
                } else {
                    side1 = scut;  side0 = tcut;
                    side1count = scnt;  side0count = tcnt;
                }          
                break;
            }
        }
        if (i == cutcount) {
            fprintf (stderr, "no side node in cut\n");
            fflush (stdout);
        }

/*
        printf ("Old Atom 0:");
        for (i = 0; i < ato[0]->nodecount; i++) {
            printf (" %d", ato[0]->nodes[i]);
        }
        printf ("\n");
        printf ("New Atom 0:");
        for (i = 0; i < side0count; i++) {
            printf (" %d", side0[i]);
        }
        printf ("\n"); fflush (stdout);
        printf ("Old Atom 1:");
        for (i = 0; i < ato[1]->nodecount; i++) {
            printf (" %d", ato[1]->nodes[i]);
        }
        printf ("\n");
        printf ("New Atom 1:");
        for (i = 0; i < side1count; i++) {
            printf (" %d", side1[i]);
        }
        printf ("\n"); fflush (stdout);
*/

        CC_IFFREE (ato[0]->nodes, int);
        ato[0]->nodes = side0;
        ato[0]->nodecount = side0count;
        CC_IFFREE (ato[1]->nodes, int);
        ato[1]->nodes = side1;
        ato[1]->nodecount = side1count;
        scut = tcut = (int *) NULL;

        *improve = val - cval;
    }
   
CLEANUP:

    CC_IFFREE (marks, int);
    CC_IFFREE (numbers, int);
    CC_IFFREE (invnumbers, int);
    CC_IFFREE (flist, int);
    CC_IFFREE (fcap, double);
    CC_IFFREE (cut, int);
    CC_IFFREE (tcut, int);
    return rval;
}

static int nonempty_cut (int ncount, int ecount, int *elist, double *ecap,
        int s, int t, double *cval, int **cut, int *cutcount)
{
    int rval = 0;
    int *tcut = (int *) NULL;
    int *bestcut = (int *) NULL;
    int *exlist = (int *) NULL;
    int *trials = (int *) NULL;
    int trialcount = 0;
    double *excap = (double *) NULL;
    int tcutcount, bestcutcount;
    double tval, bestval, smax, tmax;
    int i, si, ti, cnt, u, v;

/*
    for (i = 0; i < ecount; i++) {
        printf ("%5.3f ", ecap[i]);
        if (i % 10 == 9) printf ("\n");
    }
    printf ("\n");
*/

    /* Find min cut and check if both sides are nonempty */

    if (cval) *cval = CCtsp_LP_MAXDOUBLE;
    if (cut) *cut = (int *) NULL;
    if (cutcount) *cutcount = 0;

    bestval = CCtsp_LP_MAXDOUBLE;
    bestcutcount = 0;
    bestcut = (int *) NULL;

    rval = CCcut_mincut_st (ncount, ecount, elist, ecap, s, t, &tval,
                            &tcut, &tcutcount);
    CCcheck_rval (rval, "CCcut_mincut_st failed");

    if (tcutcount > 1 && tcutcount < ncount - 1) {
        bestval = tval;
        bestcut = tcut;
        bestcutcount = tcutcount;
        tcut = (int *) NULL;
        goto CLEANUP;
    } else {
        CC_IFFREE (tcut, int);
    }

    /* Run through each edge (u,v), with best cut having u in S, v in T, */
    /* and also best cut having v in S, u in T.                          */

    exlist = CC_SAFE_MALLOC (2*ecount + 4, int);
    CCcheck_NULL (exlist, "out of memory for exlist");
    excap = CC_SAFE_MALLOC (ecount + 2, double);
    CCcheck_NULL (exlist, "out of memory for excap");
    for (i = 0; i < ecount; i++) {
        exlist[2*i] = elist[2*i];
        exlist[2*i+1] = elist[2*i+1];
        excap[i] = ecap[i];
    }
    exlist[2*ecount] = s;
    exlist[2*ecount+2] = t;
    excap[ecount] = 1000.0;
    excap[ecount+1] = 1000.0;

    cnt = 0;
    for (i = 0; i < ecount; i++) {
        u = elist[2*i];
        v = elist[2*i+1];
        if (u != s &&  u != t &&  v != s && v != t) {
            cnt++;
        }
    }

    if (cnt > 0 && cnt <= 20) {
        trials = CC_SAFE_MALLOC (2*cnt, int);
        CCcheck_NULL (trials, "out of memory for trials");
        trialcount = 0;
        for (i = 0; i < ecount; i++) {
            u = elist[2*i];
            v = elist[2*i+1];
            if (u != s &&  u != t &&  v != s && v != t) {
                 trials[2*trialcount] = u;
                 trials[2*trialcount+1] = v;
                 trialcount++;
            }
        }
    } else {
        smax = tmax = -1.0;
        si = ti = -1;
        for (i = 0; i < ecount; i++) {
            u = elist[2*i];
            v = elist[2*i+1];
            if (u == s || v == s) {
                if (ecap[i] > smax) {
                    smax = ecap[i];
                    si = (u == s ? v : u);
                }
            } else if (u == t || v == t) {
                if (ecap[i] > tmax) {
                    tmax = ecap[i];
                    ti = (u == t ? v : u);
                }
            }
        }
        if (si == -1 || ti == -1) {
            fprintf (stderr, "didn't find si or ti\n");
            rval = 1;  goto CLEANUP;
        }
        cnt = 1;
        trials = CC_SAFE_MALLOC (2*cnt, int);
        CCcheck_NULL (trials, "out of memory for trials");
        trialcount = 1;
        trials[0] = ti;
        trials[1] = si;
    }

    for (i = 0; i < trialcount; i++) {
        u = trials[2*i];
        v = trials[2*i+1];
        if (u != s &&  u != t &&  v != s && v != t) {
            exlist[2*ecount+1] = u;
            exlist[2*ecount+3] = v;
            rval = CCcut_mincut_st (ncount, ecount+2, exlist, excap, s, t,
                                    &tval, &tcut, &tcutcount);
            CCcheck_rval (rval, "CCcut_mincut_st failed");
            if (tcutcount > 1 && tcutcount < ncount - 1 && tval < bestval &&
                tval < 1000.0) {
                CC_IFFREE (bestcut, int);
                bestval = tval;
                bestcut = tcut;
                bestcutcount = tcutcount;
                tcut = (int *) NULL;
            } else {
                CC_IFFREE (tcut, int);
            }

            exlist[2*ecount+1] = v;
            exlist[2*ecount+3] = u;
            rval = CCcut_mincut_st (ncount, ecount+2, exlist, excap, s, t,
                                    &tval, &tcut, &tcutcount);
            CCcheck_rval (rval, "CCcut_mincut_st failed");
            if (tcutcount > 1 && tcutcount < ncount - 1 && tval < bestval &&
                tval < 1000.0) {
                CC_IFFREE (bestcut, int);
                bestval = tval;
                bestcut = tcut;
                bestcutcount = tcutcount;
                tcut = (int *) NULL;
            } else {
                CC_IFFREE (tcut, int);
            }
        }
    }

CLEANUP:

    if (!rval) {
        if (cut) *cut = bestcut;
        if (cutcount) *cutcount = bestcutcount;
        if (cval) *cval = bestval;
    }
    CC_IFFREE (tcut, int);
    CC_IFFREE (trials, int);
    CC_IFFREE (exlist, int);
    CC_IFFREE (excap, double);

    return rval;
}

static int gather_cut (cutinfo *A, CCtsp_lpcut_in *cold, CCtsp_lpcut_in *cnew)
{
    int rval = 0;
    int *ar = (int *) NULL;
    int i, j, k, a, count;

    CCtsp_init_lpcut_in (cnew);
    ar = CC_SAFE_MALLOC (A->ncount, int);
    CCcheck_NULL (ar, "out of memory for ar");

    cnew->rhs = cold->rhs;
    cnew->sense = 'G';
    cnew->cliquecount = A->cliquecount;
    cnew->cliques = CC_SAFE_MALLOC (cnew->cliquecount, CCtsp_lpclique);
    CCcheck_NULL (cnew->cliques, "out of memory for cliques");

    for (i = 0; i < A->cliquecount; i++) {
        count = 0;
        for (j = 0; j < A->cliques[i].atomcount; j++) {
            a = A->cliques[i].atomlist[j];
            for (k = 0; k < A->atoms[a].nodecount; k++) {
                ar[count++] = A->atoms[a].nodes[k];
            }
        }
        if (count <= 0 || count >= A->ncount) {
            fprintf (stderr, "formed bad clique\n");
            rval = 1;  goto CLEANUP; 
        }
        rval = CCtsp_array_to_lpclique (ar, count, &(cnew->cliques[i]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
    }

    rval = CCtsp_construct_skeleton (cnew, A->ncount);
    CCcheck_rval (rval, "CCtsp_construct_skeleton failed");

CLEANUP:

    CC_IFFREE (ar, int);
    return rval;
}

int CCtsp_xtighten_lpcut_in (CCtsp_lpgraph *g, CCtsp_lpcut_in *c, double *x,
        CCtsp_lpcut_in *cout, double *pimprove)
{
    cutinfo A;
    int rval = 0;
    double oldslack, newslack, win = 0.0;
    int i;

    if (pimprove) *pimprove = 0.0;

    cutinfo_init (&A);
    CCtsp_init_lpcut_in (cout);

    if (c->branch != 0 || !CCtsp_hypergraph_type_in (c)) {
        fprintf (stderr, "try to xtighten a impropercut\n");
        rval = 1; goto CLEANUP;
    }
    if (!cout) {
        fprintf (stderr, "xtighten called without an output cut\n");
        rval = 1; goto CLEANUP;
    }
        
/*
    printf ("lpcut_in skeleton:");
    for (i=0; i<c->skel.atomcount; i++) {
        printf (" %d", c->skel.atoms[i]);
    }
    printf ("\n");
*/

/*
    CCtsp_print_lpcut_in (c);
    fflush (stdout);
*/

    rval = cutinfo_build (&A, c, g);
    CCcheck_rval (rval, "cutinfo_build failed");
    
    rval = xtighten_cut (&A, g, x, &win);
    CCcheck_rval (rval, "xtighten_cut failed");

    if (win > 0.0) {
        rval = gather_cut (&A, c, cout);
        CCcheck_rval (rval, "gather_cut failed");

        oldslack = CCtsp_cutprice (g, c, x);
        newslack = CCtsp_cutprice (g, cout, x);
        if (newslack > oldslack - XTIGHT_EPS) {
            printf ("OOPS: win = %f while oldslack = %f and newslack = %f\n",
                     win, oldslack, newslack);
            fflush (stdout);

            rval = CCverify_cut (cout, g->ncount, CC_TYPE_ALL, &i, 0,
                             (int *) NULL, (char *) NULL);
            if (rval) printf ("Cut is not even valid!\n");
            CCtsp_print_lpcut_in (cout);
            CCtsp_print_lpcut_in (c);
            exit (1);

            CCtsp_free_lpcut_in (cout);
        } else {
            rval = CCverify_cut (cout, g->ncount, CC_TYPE_ALL, &i, 0,
                             (int *) NULL, (char *) NULL);
            if (rval) {
                fprintf (stderr, "Discard cut, could not verify\n");
                CCtsp_free_lpcut_in (cout);
                rval = 0; goto CLEANUP;
            }
            printf ("Xtighten %f\n", oldslack - newslack);
            fflush (stdout);
            if (pimprove) *pimprove = oldslack - newslack;
        }
    }

CLEANUP:

    cutinfo_free (&A);
    return rval;
}

