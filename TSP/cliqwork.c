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
/*            SOME ROUTINES FOR HANDLING CLIQUES AND CUTS                   */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: July 15, 1997                                                     */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCtsp_mark_clique (CCtsp_lpclique *c, int *marks, int marker)      */
/*    MARKS the nodes in clique c                                           */
/*     -marks an array of length at least ncount                            */
/*     -marker an int that is used to mark the clique entries in marks      */
/*                                                                          */
/*  void CCtsp_mark_domino (CCtsp_lpdomino *c, int *marks, int marker)      */
/*    MARKS the nodes in domino c                                           */
/*                                                                          */
/*  void CCtsp_mark_clique_and_neighbors (CCtsp_lpgraph *g,                 */
/*      CCtsp_lpclique *c, int *marks, int marker)                          */
/*    MARKS the clique and the neighbors of the clique                      */
/*                                                                          */
/*  void CCtsp_mark_domino_and_neighbors (CCtsp_lpgraph *g,                 */
/*      CCtsp_lpdomino *c, int *marks, int marker)                          */
/*    MARKS the domino and the neighbors of the domino                      */
/*                                                                          */
/*  void CCtsp_mark_clique_and_neighbors_double (CCtsp_lpgraph *g,          */
/*      CCtsp_lpclique *c, double *marks, double marker)                    */
/*    MARKS the clique and the neighbors of the clique in a double array.   */
/*                                                                          */
/*  void CCtsp_mark_cut (CCtsp_lpcut_in *c, int *marks, int marker)         */
/*    MARKS the nodes in the cut.                                           */
/*                                                                          */
/*  void CCtsp_mark_cut_and_neighbors (CCtsp_lpgraph *g,                    */
/*      CCtsp_lpcut_in *c, int *marks, int marker)                          */
/*    MARKS the nodes in the cut and their neighbors                        */
/*                                                                          */
/*  void CCtsp_is_clique_marked (CCtsp_lpclique *c, int *marks,             */
/*      int marker, int *yes_no)                                            */
/*    CHECKS if a node in the clique is marked with the value marker.       */
/*     -yesno returns the result (1 is yes and 0 is no)                     */
/*                                                                          */
/*  void CCtsp_clique_count (CCtsp_lpclique *c, int *count)                 */
/*    COUNTS the nodes in the clique.                                       */
/*                                                                          */
/*  void CCtsp_clique_marked_count (CCtsp_lpclique *c, int *marks,          */
/*      int marker, int *count)                                             */
/*    COUNTS the nodes in the clique have mark equal to marker.             */
/*                                                                          */
/*  int CCtsp_clique_to_array (CCtsp_lpclique *c, int **ar, int *count)     */
/*    EXPANDS a clique into an array of integers.                           */
/*                                                                          */
/*  int CCtsp_clique_delta (CCtsp_lpgraph *g, double *x,                    */
/*      CCtsp_lpclique *c, double *delta)                                   */
/*    COMPUTES the sum of the x-edges in the coboundary of the clique,      */
/*     that is, x(delta(c)).                                                */
/*     -delta returns the sum                                               */
/*                                                                          */
/*  int CCtsp_segment_to_subtour (CCtsp_lpcut_in **cut, int a, int b,       */
/*      int ncount)                                                         */
/*    BUILDS a subtour CCtsp_lpcut_in from an the segment.                  */
/*     -cut will return the subtour (it will be allocated).                 */
/*                                                                          */
/*  int CCtsp_array_to_subtour (CCtsp_lpcut_in **cut, int *ar,              */
/*      int acount, int ncount)                                             */
/*    BUILDS a subtour CCtsp_lpcut_in from an array.                        */
/*     -cut will return the subtour (it will be allocated).                 */
/*                                                                          */
/*  void CCtsp_init_lpcut_in (CCtsp_lpcut_in *c)                            */
/*    INITIALIZE the fields of the CCtsp_lpcut_in structure                 */
/*                                                                          */
/*  void CCtsp_init_lpcut (CCtsp_lpcut *c)                                  */
/*    INITIALIZE the fields of the CCtsp_lpcut structure                    */
/*                                                                          */
/*  void CCtsp_init_lpclique (CCtsp_lpclique *c)                            */
/*    INITIALIZE the fields of the CCtsp_lpclique structure                 */
/*                                                                          */
/*  void CCtsp_init_lpdomino (CCtsp_lpdomino *c)                            */
/*    INITIALIZE the fields of the CCtsp_lpdomino structure                 */
/*                                                                          */
/*  int CCtsp_array_to_lpclique (int *ar, int acount,                       */
/*      CCtsp_lpclique *cliq)                                               */
/*    BUILDS an CCtsp_lpclique represented the nodes in an array.           */
/*     -ar is an array of node numbers                                      */
/*     -acount is the length of ar                                          */
/*     -cliq's segcount and nodes will be filled with the compressed        */
/*      version of the nodes in ar                                          */
/*                                                                          */
/*  int CCtsp_seglist_to_lpclique (int nseg, int *list,                     */
/*      CCtsp_lpclique *cliq)                                               */
/*    BUILDS the CCtsp_lpclique represented by a list of segments (it       */
/*     will sort the segments before it puts them into the CCtsp_segment    */
/*     structures)                                                          */
/*     -list is an array of segments in lo-hi-lo-hi format                  */
/*     -clig's segcount and nodes will be filled in (nodes will be          */
/*      allocated)                                                          */
/*                                                                          */
/*  int CCtsp_shrunk_set_to_lpclique (int cnt, int *set, int *wset,         */
/*      CC_SRKexpinfo *expand, CCtsp_lpclique *cliq)                        */
/*    BUILDS an lpclique by exanding a shrunk set of nodes.                 */
/*     -cnt is the number of nodes in set                                   */
/*     -set is an array of the nodes                                        */
/*     -wset is a working array, it should be ncount long and the values    */
/*      may be changed by this function                                     */
/*     -expand contains the info to expand the nodes                        */
/*     -cliq returns the clique                                             */
/*                                                                          */
/*  int CCtsp_add_nodes_to_lpclique (CCtsp_lpclique *cin,                   */
/*      CCtsp_lpclique *cout, int addcount, int *adda)                      */
/*    ADDS nodes to clique cin, and returns the new clique in cout          */
/*     -addcount has number of nodes to be added.                           */
/*     -adda has the indices of the addcount nodes to be added.             */
/*                                                                          */
/*  int CCtsp_delete_nodes_from_lpclique (CCtsp_lpclique *cin,              */
/*      CCtsp_lpclique *cout, int delcount, int *del)                       */
/*    DELETES nodes from clique cin, and returns the new clique in cout     */
/*     -delcount has number of nodes to be deleted.                         */
/*     -del has the indices of the delcount nodes to be deleted.            */
/*                                                                          */
/*  void CCtsp_print_lpcut_in (CCtsp_lpcut_in *c)                           */
/*    PRINTS the CCtsp_lpcut_in                                             */
/*                                                                          */
/*  void CCtsp_print_lpclique (CCtsp_lpclique *c)                           */
/*    PRINTS the segments in the clique to stdout.                          */
/*                                                                          */
/*  void CCtsp_print_lpdomino (CCtsp_lpdomino *d)                           */
/*    PRINTS the cliques of the domino.                                     */
/*                                                                          */
/*  int CCtsp_copy_lpcut_in (CCtsp_lpcut_in *c, CCtsp_lpcut_in *new)        */
/*    COPIES an CCtsp_lpcut_in                                              */
/*     -c is a pointer to an CCtsp_lpcut_in                                 */
/*     -new returns the copied CCtsp_lpcut_in                               */
/*                                                                          */
/*  int CCtsp_lpcut_to_lpcut_in (CCtsp_lpcuts *cuts, CCtsp_lpcut *c,        */
/*      CCtsp_lpcut_in *new)                                                */
/*    COPIES an CCtsp_lpcut to an CCtsp_lpcut_in                            */
/*     -cuts is a pointer to the structure holding the set of cuts          */
/*     -c is the cut to be copied                                           */
/*     -new returns the copied cut                                          */
/*                                                                          */
/*  void CCtsp_lpclique_compare (CCtsp_lpclique *a, CCtsp_lpclique *b,      */
/*      int *diff)                                                          */
/*    COMPARES two CCtsp_lpcliques.                                         */
/*     -diff is set to 1 if they differ and 0 if they are the same          */
/*      NOTES: Assumes segments are ordered.                                */
/*                                                                          */
/*  int CCtsp_copy_lpclique (CCtsp_lpclique *c, CCtsp_lpclique *new)        */
/*    COPIES an CCtsp_lpclique                                              */
/*     -c is a pointer to an CCtsp_lpclique                                 */
/*     -new returns the copied clique                                       */
/*                                                                          */
/*  int CCtsp_copy_lpdomino (CCtsp_lpdomino *c, CCtsp_lpdomino *new)        */
/*    COPIES an CCtsp_lpdomino                                              */
/*                                                                          */
/*  int CCtsp_create_lpcliques (CCtsp_lpcut_in *c, int cliquecount)         */
/*    ALLOCATES and INITIALIZES the cliques space for c.                    */
/*                                                                          */
/*  int CCtsp_max_node (CCtsp_lpcut_in *c)                                  */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_build_dp_cut (CCtsp_lpcut_in **cut, int ndomino, int *Acount, */
/*      int **A, int *Bcount, int **B, int handlecount, int *handle)        */
/*    BUILD the lpcut_in for a domino parity cut.                           */
/*     -ndomino number of dominos                                           */
/*     -Acount[i] specifies the number of nodes in one side of ith domino   */
/*     -Bcount[i] specifies the number of nodes in other side of  domino    */
/*     -A[i] is an array of the nodes in one side of ith domino             */
/*     -B[i] is an array of the nodes in other side of ith domino           */
/*     -handlecount is number of nodes in the handle                        */
/*     -handle is an array of the nodes in the handle                       */
/*                                                                          */
/* int CCtsp_extract_dp_cut (CCtsp_lpcut_in *c, int *pndomino,              */
/*     int **pAcount, int ***pA, int **pBcount, int ***pB,                  */
/*     int *phandlecount, int **phandle)                                    */
/*    EXTRACTS the set form of a DP-cut from the lpcut_in form.             */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "tsp.h"
#include "cut.h"

static int build_dominos (CCtsp_lpdomino **pdom, int ndomino, int *Acount,
        int **A, int *Bcount, int **B);
static int copy_lpclique_global (CCtsp_lpclique *c, CCtsp_lpclique *new,
        int *global_names, int *local_perm, int *global_invperm);

void CCtsp_mark_clique (CCtsp_lpclique *c, int *marks, int marker)
{
    int j, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        marks[j] = marker;
    }
}

void CCtsp_mark_domino (CCtsp_lpdomino *c, int *marks, int marker)
{
    CCtsp_mark_clique (&(c->sets[0]), marks, marker);
    CCtsp_mark_clique (&(c->sets[1]), marks, marker);
}

void CCtsp_mark_clique_and_neighbors (CCtsp_lpgraph *g, CCtsp_lpclique *c,
        int *marks, int marker)
{
    int j, k, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        marks[j] = marker;
        for (k = 0; k < g->nodes[j].deg; k++) {
            marks[g->nodes[j].adj[k].to] = marker;
        }
    }
}

void CCtsp_mark_domino_and_neighbors (CCtsp_lpgraph *g, CCtsp_lpdomino *c,
        int *marks, int marker)
{
    CCtsp_mark_clique_and_neighbors (g, &(c->sets[0]), marks, marker);
    CCtsp_mark_clique_and_neighbors (g, &(c->sets[1]), marks, marker);
}

void CCtsp_mark_cut (CCtsp_lpcut_in *c, int *marks, int marker)
{
    int i;

    for (i = 0; i < c->cliquecount; i++) {
        CCtsp_mark_clique (&(c->cliques[i]), marks, marker);
    }
    for (i = 0; i < c->dominocount; i++) {
        CCtsp_mark_domino (&(c->dominos[i]), marks, marker);
    }
    if (c->TP_handles != (CCtsp_lpclique *) NULL) {
        for (i = 0; i < 2; i++) {
            CCtsp_mark_clique (&(c->TP_handles[i]), marks, marker);
        }
        for (i = 0; i < c->TP_domcount0; i++) {
            CCtsp_mark_domino (&(c->TP_dominos0[i]), marks, marker);
        }
        for (i = 0; i < c->TP_domcount1; i++) {
            CCtsp_mark_domino (&(c->TP_dominos1[i]), marks, marker);
        }
        for (i = 0; i < c->TP_tricount; i++) {
            CCtsp_mark_clique (&(c->TP_tsets[i]), marks, marker);
            CCtsp_mark_domino (&(c->TP_semicuts0[i]), marks, marker);
            CCtsp_mark_domino (&(c->TP_semicuts1[i]), marks, marker);
        }
    }
}

void CCtsp_mark_cut_and_neighbors (CCtsp_lpgraph *g, CCtsp_lpcut_in *c,
        int *marks, int marker)
{
    int i;

    for (i = 0; i < c->cliquecount; i++) {
        CCtsp_mark_clique_and_neighbors (g, &(c->cliques[i]), marks, marker);
    }
    for (i = 0; i < c->dominocount; i++) {
        CCtsp_mark_domino_and_neighbors (g, &(c->dominos[i]), marks, marker);
    }
    if (c->TP_handles != (CCtsp_lpclique *) NULL) {
        for (i = 0; i < 2; i++) {
            CCtsp_mark_clique_and_neighbors (g, &(c->TP_handles[i]),
                                             marks, marker);
        }
        for (i = 0; i < c->TP_domcount0; i++) {
            CCtsp_mark_domino_and_neighbors (g, &(c->TP_dominos0[i]),
                                             marks, marker);
        }
        for (i = 0; i < c->TP_domcount1; i++) {
            CCtsp_mark_domino_and_neighbors (g, &(c->TP_dominos1[i]),
                                             marks, marker);
        }
        for (i = 0; i < c->TP_tricount; i++) {
            CCtsp_mark_clique_and_neighbors (g, &(c->TP_tsets[i]),
                                             marks, marker);
            CCtsp_mark_domino_and_neighbors (g, &(c->TP_semicuts0[i]),
                                             marks, marker);
            CCtsp_mark_domino_and_neighbors (g, &(c->TP_semicuts1[i]),
                                             marks, marker);
        }
    }
}

void CCtsp_mark_clique_and_neighbors_double (CCtsp_lpgraph *g,
        CCtsp_lpclique *c, double *marks, double marker)
{
    int j, k, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        marks[j] = marker;
        for (k = 0; k < g->nodes[j].deg; k++) {
            marks[g->nodes[j].adj[k].to] = marker;
        }
    }
}

void CCtsp_is_clique_marked (CCtsp_lpclique *c, int *marks, int marker,
        int *yes_no)
{
    int j, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        if (marks[j] == marker) {
            *yes_no = 1;
            return;
        }
    }
    *yes_no = 0;
}

void CCtsp_clique_count (CCtsp_lpclique *c, int *count)
{
    int i;
    CCtsp_segment *nodes = c->nodes;

    *count = 0;
    for (i = 0; i < c->segcount; i++) {
        (*count) += (nodes[i].hi - nodes[i].lo + 1);
    }
}

void CCtsp_clique_marked_count (CCtsp_lpclique *c, int *marks, int marker,
        int *count)
{
    int j, tmp;

    *count = 0;
    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        if (marks[j] == marker) {
            (*count)++;
        }
    }
}

int CCtsp_clique_to_array (CCtsp_lpclique *c, int **ar, int *count)
{
    int rval = 0;
    int j, tmp;
    int k = 0;

    *ar = (int *) NULL;

    CCtsp_clique_count (c, count);
    if (*count) {
        *ar = CC_SAFE_MALLOC (*count, int);
        if (!(*ar)) {
            fprintf (stderr, "out of memory in CCtsp_clique_to_array\n");
            rval = 1; goto CLEANUP;
        }
        CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
            (*ar)[k++] = j;
        }
    }

CLEANUP:

    return rval;
}

int CCtsp_clique_delta (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *c,
        double *delta)
{
    int rval = 0;
    int j, k, tmp;
    int *marks = (int *) NULL;
    CCtsp_lpnode *n;

    *delta = 0.0;

    marks = CC_SAFE_MALLOC (g->ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in CCtsp_clique_delta\n");
        rval = 1; goto CLEANUP;
    }

    CCtsp_mark_clique_and_neighbors (g, c, marks, 0);
    CCtsp_mark_clique (c, marks, 1);

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        n = &g->nodes[j];
        for (k = 0; k < n->deg; k++) {
            if (!marks[n->adj[k].to]) {
                (*delta) += x[n->adj[k].edge];
            }
        }
    }

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

void CCtsp_init_lpcut_in (CCtsp_lpcut_in *c)
{
    if (c) {
        c->cliquecount = 0;
        c->cliquemult  = (int *) NULL;
        c->dominocount = 0;
        c->semicount   = 0;
        c->semimult    = (int *) NULL;
        c->coefcount   = 0;
        c->coef        = (int *) NULL;
        c->coefedges   = (int *) NULL;
        c->rhs         = 0;
        c->rhs         = 0;
        c->sense       = 'X';
        c->branch      = 0;
        c->cliques     = (CCtsp_lpclique *) NULL;
        c->dominos     = (CCtsp_lpdomino *) NULL;
        c->semicuts    = (CCtsp_lpdomino *) NULL;

        c->TP_handles   = (CCtsp_lpclique *) NULL;
        c->TP_tricount  = 0;
        c->TP_domcount0 = 0;
        c->TP_domcount1 = 0;
        c->TP_dominos0  = (CCtsp_lpdomino *) NULL;
        c->TP_dominos1  = (CCtsp_lpdomino *) NULL;
        c->TP_tsets     = (CCtsp_lpclique *) NULL;
        c->TP_semicuts0 = (CCtsp_lpdomino *) NULL;
        c->TP_semicuts1 = (CCtsp_lpdomino *) NULL;

        CCtsp_init_skeleton (&c->skel);
        c->next        = (CCtsp_lpcut_in *) NULL;
        c->prev        = (CCtsp_lpcut_in *) NULL;
    }
}

void CCtsp_init_lpcut (CCtsp_lpcut *c)
{
    if (c) {
        c->cliquecount = 0;
        c->cliquemult  = (int *) NULL;
        c->dominocount = 0;
        c->semicount   = 0;
        c->semimult    = (int *) NULL;
        c->coefcount   = 0;
        c->coef        = (int *) NULL;
        c->coefedges   = (int *) NULL;
        c->modcount    = 0;
        c->age         = 0;
        c->rhs         = 0;
        c->sense       = 'X';
        c->branch      = 0;
        c->cliques     = (int *) NULL;
        c->dominos     = (int *) NULL;
        c->semicuts    = (int *) NULL;

        c->TP_handles   = (int *) NULL;
        c->TP_tricount  = 0;
        c->TP_domcount0 = 0;
        c->TP_domcount1 = 0;
        c->TP_dominos0  = (int *) NULL;
        c->TP_dominos1  = (int *) NULL;
        c->TP_tsets     = (int *) NULL;
        c->TP_semicuts0 = (int *) NULL;
        c->TP_semicuts1 = (int *) NULL;

        c->mods = (CCtsp_sparser *) NULL;
        CCtsp_init_skeleton (&c->skel);
    }
}

int CCtsp_copy_lpcut_in (CCtsp_lpcut_in *c, CCtsp_lpcut_in *new)
{
    int i, rval = 0;

    CCtsp_init_lpcut_in (new);
    new->cliquecount    = c->cliquecount;
    new->dominocount    = c->dominocount;
    new->semicount      = c->semicount;
    new->coefcount      = c->coefcount;
    new->rhs            = c->rhs;
    new->sense          = c->sense;

    if (c->cliquemult) {
        if (!c->cliquecount) {
            fprintf (stderr, "multipliers but no cliques\n");
            rval = 1; goto CLEANUP;
        }
        CC_MALLOC (new->cliquemult, c->cliquecount, int);
        for (i = 0; i < c->cliquecount; i++) {
            new->cliquemult[i] = c->cliquemult[i];
        }
    }

    if (c->semimult) {
        if (!c->semicount) {
            fprintf (stderr, "multipliers but no cliques\n");
            rval = 1; goto CLEANUP;
        }
        CC_MALLOC (new->semimult, c->semicount, int);
        for (i = 0; i < c->semicount; i++) {
            new->semimult[i] = c->semimult[i];
        }
    }

    if (c->coefcount) {
        CC_MALLOC (new->coef, c->coefcount, int);
        CC_MALLOC (new->coefedges, 2*c->coefcount, int);
        for (i = 0; i < c->coefcount; i++) {
            new->coef[i] = c->coef[i];
            new->coefedges[2*i] = c->coefedges[2*i];
            new->coefedges[2*i+1] = c->coefedges[2*i+1];
        }
    }

    if (c->cliquecount) {
        CC_MALLOC (new->cliques, c->cliquecount, CCtsp_lpclique);
        for (i = 0; i < c->cliquecount; i++) {
            rval = CCtsp_copy_lpclique (&c->cliques[i], &new->cliques[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
        }
    }

    if (c->dominocount) {
        CC_MALLOC (new->dominos, c->dominocount, CCtsp_lpdomino);
        for (i = 0; i < c->dominocount; i++) {
            rval = CCtsp_copy_lpdomino (&c->dominos[i], &new->dominos[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpdomino failed");
        }
    }

    if (c->semicount) {
        CC_MALLOC (new->semicuts, c->semicount, CCtsp_lpdomino);
        for (i = 0; i < c->semicount; i++) {
            rval = CCtsp_copy_lpdomino (&c->semicuts[i], &new->semicuts[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpdomino failed");
        }
    }

    if (c->TP_handles != (CCtsp_lpclique *) NULL) {
        CC_MALLOC (new->TP_handles, 2, CCtsp_lpclique);
        for (i = 0; i < 2; i++) {
            rval = CCtsp_copy_lpclique (&c->TP_handles[i], &new->TP_handles[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
        }
    }

    if (c->TP_tricount) {
        CC_MALLOC (new->TP_tsets, c->TP_tricount, CCtsp_lpclique);
        CC_MALLOC (new->TP_semicuts0, c->TP_tricount, CCtsp_lpdomino);
        CC_MALLOC (new->TP_semicuts1, c->TP_tricount, CCtsp_lpdomino);
        for (i = 0; i < c->TP_tricount; i++) {
            rval = CCtsp_copy_lpclique (&c->TP_tsets[i], &new->TP_tsets[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
            rval = CCtsp_copy_lpdomino (&c->TP_semicuts0[i],
                                      &new->TP_semicuts0[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpdomino failed");
            rval = CCtsp_copy_lpdomino (&c->TP_semicuts1[i],
                                      &new->TP_semicuts1[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpdomino failed");
        }
        new->TP_tricount = c->TP_tricount;
    }

    if (c->TP_domcount0) {
        CC_MALLOC (new->TP_dominos0, c->TP_domcount0, CCtsp_lpdomino);
        for (i = 0; i < c->TP_domcount0; i++) {
            rval = CCtsp_copy_lpdomino (&c->TP_dominos0[i],
                                      &new->TP_dominos0[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpdomino failed");
        }
        new->TP_domcount0 = c->TP_domcount0;
    }

    if (c->TP_domcount1) {
        CC_MALLOC (new->TP_dominos1, c->TP_domcount1, CCtsp_lpdomino);
        for (i = 0; i < c->TP_domcount1; i++) {
            rval = CCtsp_copy_lpdomino (&c->TP_dominos1[i],
                                      &new->TP_dominos1[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpdomino failed");
        }
        new->TP_domcount1 = c->TP_domcount1;
    }

    rval = CCtsp_copy_skeleton (&c->skel, &new->skel);
    CCcheck_rval (rval, "CCtsp_copy_skeleton failed");

CLEANUP:
    if (rval) { CCtsp_free_lpcut_in (new); }
    return rval;
}

int CCtsp_copy_lpcut_in_global (CCtsp_lpcut_in *c, CCtsp_lpcut_in *new,
        int *global_names, int *local_perm, int *global_invperm,
        int global_ncount)
{
    int rval = 0;
    int i, k;
    CCtsp_lpdomino *dom, *newdom;

    CCtsp_init_lpcut_in (new);

    new->cliquecount    = c->cliquecount;
    new->dominocount    = c->dominocount;
    new->rhs            = c->rhs;
    new->sense          = c->sense;

    if (c->cliquecount) {
        CC_MALLOC (new->cliques, c->cliquecount, CCtsp_lpclique);
        for (i = 0; i < c->cliquecount; i++) {
            rval = copy_lpclique_global (&c->cliques[i],
                       &new->cliques[i], global_names, local_perm,
                       global_invperm);
            CCcheck_rval (rval, "copy_lpclique_global failed");
        }
    }
    if (c->dominocount > 0) {
        CC_MALLOC (new->dominos ,c->dominocount, CCtsp_lpdomino);
        for (i = 0; i < c->dominocount; i++) {
            dom = &c->dominos[i];
            newdom = &new->dominos[i];
            CCtsp_init_lpdomino (newdom);
            for (k = 0; k < 2; k++) {
                rval = copy_lpclique_global (&dom->sets[k], &newdom->sets[k],
                                  global_names, local_perm, global_invperm);
                CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
            }
        }
    }

    if (c->TP_handles != (CCtsp_lpclique *) NULL) {
        CC_MALLOC (new->TP_handles, 2, CCtsp_lpclique);
        for (i = 0; i < 2; i++) {
            rval = copy_lpclique_global (&c->TP_handles[i], &new->TP_handles[i],
                                 global_names, local_perm, global_invperm);
            CCcheck_rval (rval, "copy_lpclique_global failed");
        }
    }

    if (c->TP_tricount) {
        CC_MALLOC (new->TP_tsets, c->TP_tricount, CCtsp_lpclique);
        CC_MALLOC (new->TP_semicuts0, c->TP_tricount, CCtsp_lpdomino);
        CC_MALLOC (new->TP_semicuts1, c->TP_tricount, CCtsp_lpdomino);
        for (i = 0; i < c->TP_tricount; i++) {
            rval = copy_lpclique_global (&c->TP_tsets[i], &new->TP_tsets[i],
                                 global_names, local_perm, global_invperm);
            CCcheck_rval (rval, "copy_lpclique_global failed");

            dom = &c->TP_semicuts0[i];
            newdom = &new->TP_semicuts0[i];
            CCtsp_init_lpdomino (newdom);
            for (k = 0; k < 2; k++) {
                rval = copy_lpclique_global (&dom->sets[k], &newdom->sets[k],
                                     global_names, local_perm, global_invperm);
                CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
            }

            dom = &c->TP_semicuts1[i];
            newdom = &new->TP_semicuts1[i];
            CCtsp_init_lpdomino (newdom);
            for (k = 0; k < 2; k++) {
                rval = copy_lpclique_global (&dom->sets[k], &newdom->sets[k],
                                     global_names, local_perm, global_invperm);
                CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
            }
        }
        new->TP_tricount = c->TP_tricount;
    }

    if (c->TP_domcount0) {
        CC_MALLOC (new->TP_dominos0, c->TP_domcount0, CCtsp_lpdomino);
        for (i = 0; i < c->TP_domcount0; i++) {
            dom = &c->TP_dominos0[i];
            newdom = &new->TP_dominos0[i];
            CCtsp_init_lpdomino (newdom);
            for (k = 0; k < 2; k++) {
                rval = copy_lpclique_global (&dom->sets[k], &newdom->sets[k],
                                     global_names, local_perm, global_invperm);
                CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
            }
        }
        new->TP_domcount0 = c->TP_domcount0;
    }

    if (c->TP_domcount1) {
        CC_MALLOC (new->TP_dominos1, c->TP_domcount1, CCtsp_lpdomino);
        for (i = 0; i < c->TP_domcount1; i++) {
            dom = &c->TP_dominos1[i];
            newdom = &new->TP_dominos1[i];
            CCtsp_init_lpdomino (newdom);
            for (k = 0; k < 2; k++) {
                rval = copy_lpclique_global (&dom->sets[k], &newdom->sets[k],
                                     global_names, local_perm, global_invperm);
                CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
            }

        }
        new->TP_domcount1 = c->TP_domcount1;
    }

    if (c->dominocount == 0 && c->TP_handles == (CCtsp_lpclique *) NULL) {
        rval = CCtsp_construct_skeleton (new, global_ncount);
        CCcheck_rval (rval, "CCtsp_construct_skeleton failed");
    }

CLEANUP:

    if (rval) { CCtsp_free_lpcut_in (new); }
    return rval;
}

int CCtsp_segment_to_subtour (CCtsp_lpcut_in **cut, int a, int b, int ncount)
{
    int rval = 0;
    int list[2];
    int t;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    *cut = (CCtsp_lpcut_in *) NULL;
    if (a > b) CC_SWAP (a, b, t);

    c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (c, "out of memory in CCtsp_seqment_to_subtour");
    CCtsp_init_lpcut_in (c);

    c->cliquecount = 1;
    c->cliques = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    if (!c->cliques) {
        fprintf (stderr, "out of memory in CCtsp_segment_to_subtour\n");
        rval = 1; goto CLEANUP;
    }

    list[0] = a;
    list[1] = b;
    rval = CCtsp_seglist_to_lpclique (1, list, &(c->cliques[0]));
    if (rval) {
        goto CLEANUP;
    }
    c->rhs = 2;
    c->sense = 'G';
    c->branch = 0;

    rval = CCtsp_construct_skeleton (c, ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }

    *cut = c;
    rval = 0;

CLEANUP:

    if (rval) {
        if (c) {
            CCtsp_free_lpcut_in (c);
            CC_FREE (c, CCtsp_lpcut_in);
        }
    }
    return rval;
}

int CCtsp_array_to_subtour (CCtsp_lpcut_in **cut, int *ar, int acount,
        int ncount)
{
    int rval = 0;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    *cut = (CCtsp_lpcut_in *) NULL;

    c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (c, "out of memory in CCtsp_array_to_subtour");
    CCtsp_init_lpcut_in (c);

    c->cliquecount = 1;
    c->cliques = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    if (!c->cliques) {
        fprintf (stderr, "out of memory in CCtsp_array_to_subtour\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_array_to_lpclique (ar, acount, &(c->cliques[0]));
    if (rval) {
        goto CLEANUP;
    }
    c->rhs = 2;
    c->sense = 'G';
    c->branch = 0;

    rval = CCtsp_construct_skeleton (c, ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }

    *cut = c;
    rval = 0;

CLEANUP:

    if (rval) {
        if (c) {
            CCtsp_free_lpcut_in (c);
            CC_FREE (c, CCtsp_lpcut_in);
        }
    }
    return rval;
}

void CCtsp_init_lpclique (CCtsp_lpclique *c)
{
    if (c) {
        c->segcount = 0;
        c->nodes = (CCtsp_segment *) NULL;
        c->hashnext = 0;
        c->refcount = 0;
    }
}

void CCtsp_init_lpdomino (CCtsp_lpdomino *c)
{
    if (c) {
        CCtsp_init_lpclique (&(c->sets[0]));
        CCtsp_init_lpclique (&(c->sets[1]));
        c->hashnext = 0;
        c->refcount = 0;
    }
}

int CCtsp_array_to_lpclique (int *ar, int acount, CCtsp_lpclique *cliq)
{
    int i, nseg;

    /* Function will alter the order on the array */

    CCutil_int_array_quicksort (ar, acount);
    nseg = 0;
    i = 0;
    while (i < acount) {
        while (i < (acount - 1) && ar[i + 1] == (ar[i] + 1))
            i++;
        i++;
        nseg++;
    }

    cliq->nodes = CC_SAFE_MALLOC (nseg, CCtsp_segment);
    if (!cliq->nodes) {
        fprintf (stderr, "out of memory in CCtsp_array_to_lpclique\n");
        return 1;
    }
    cliq->segcount = nseg;

    nseg = 0;
    i = 0;
    while (i < acount) {
        cliq->nodes[nseg].lo = ar[i];
        while (i < (acount - 1) && ar[i + 1] == (ar[i] + 1))
            i++;
        cliq->nodes[nseg].hi = ar[i++];
        nseg++;
    }
    return 0;
}

int CCtsp_seglist_to_lpclique (int nseg, int *list, CCtsp_lpclique *cliq)
{
    int i;
    int *perm = (int *) NULL;
    int *len  = (int *) NULL;
    int rval = 0;

    perm = CC_SAFE_MALLOC (nseg, int);
    len  = CC_SAFE_MALLOC (nseg, int);
    if (!perm || !len) {
        fprintf (stderr, "out of memory in CCtsp_seglist_to_lpclique\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nseg; i++) {
        perm[i] = i;
        len[i] = list[2*i];
    }
    CCutil_int_perm_quicksort (perm, len, nseg);

    cliq->nodes = CC_SAFE_MALLOC (nseg, CCtsp_segment);
    if (!cliq->nodes) {
        fprintf (stderr, "out of memory in CCtsp_seglist_to_lpclique\n");
        rval = 1; goto CLEANUP;
    }
    cliq->segcount = nseg;

    for (i = 0; i < nseg; i++) {
        cliq->nodes[i].lo = list[2*perm[i]];
        cliq->nodes[i].hi = list[2*perm[i]+1];
    }

    nseg = 0;

CLEANUP:

    CC_IFFREE (perm, int);
    CC_IFFREE (len, int);

    return rval;
}

int CCtsp_shrunk_set_to_lpclique (int cnt, int *set, int *wset,
        CC_SRKexpinfo *expand, CCtsp_lpclique *cliq)
{
    int rval = 0;
    int wcount = 0;
    int i, j;

    CCtsp_init_lpclique (cliq);

    for (i = 0; i < cnt; i++) {
        for (j = expand->memindex[set[i]];
             j < expand->memindex[set[i]+1]; j++) {
            wset[wcount++] = expand->members[j];
        }
    }
    rval = CCtsp_array_to_lpclique (wset, wcount, cliq);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CCtsp_add_nodes_to_lpclique (CCtsp_lpclique *cin, CCtsp_lpclique *cout,
        int addcount, int *adda)
{
    int rval = 0, count = 0, maxn = 0, i, tmp, *ar = (int *) NULL;
    char *marks = (char *) NULL;

    if (addcount == 0) { return CCtsp_copy_lpclique (cin, cout); }

    CCtsp_init_lpclique (cout);
    for (i = 0; i < cin->segcount; i++) {
        if (cin->nodes[i].lo > cin->nodes[i].hi) {
            /* not needed, but test makes Valgrind happy -- 140226 */
            printf ("No way -- lo is bigger than hi.\n");
            rval = 1; goto CLEANUP;
        }
        if (cin->nodes[i].hi > maxn) maxn = cin->nodes[i].hi;
    }
    for (i = 0; i < addcount; i++) {
        if (adda[i] > maxn) maxn = adda[i];
    }

    CC_MALLOC (marks, maxn + 1, char);
    CC_MALLOC (ar, maxn + 1, int);
    for (i = 0; i < addcount; i++) { marks[adda[i]] = 0; }

    CC_FOREACH_NODE_IN_CLIQUE (i, *cin, tmp) {
        ar[count++] = i;
        marks[i] = 1;
    }

    for (i = 0; i < addcount; i++) {
        if (marks[adda[i]] == 0) {
            ar[count++] = adda[i];
            marks[adda[i]] = 1;
        }
    }
    rval = CCtsp_array_to_lpclique (ar, count, cout);
    CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

CLEANUP:

    CC_IFFREE (ar, int);
    CC_IFFREE (marks, char);
    return rval;
}

int CCtsp_delete_nodes_from_lpclique (CCtsp_lpclique *cin,
        CCtsp_lpclique *cout, int delcount, int *del)
{
    int count = 0;
    int rval  = 0;
    int *ar   = (int *) NULL;
    int i, tmp;
    char *marks = (char *) NULL;
    int maxn = 0;

    CCtsp_init_lpclique (cout);

    for (i = 0; i < cin->segcount; i++) {
        if (cin->nodes[i].hi > maxn) maxn = cin->nodes[i].hi;
    }
    for (i = 0; i < delcount; i++) {
        if (del[i] > maxn) maxn = del[i];
    }
    marks = CC_SAFE_MALLOC (maxn + 1, char);
    ar = CC_SAFE_MALLOC (maxn + 1, int);
    if (!marks || !ar) {
        fprintf (stderr, "out of memory in CCtsp_delete_nodes_from_lpclique\n");
        rval = 1; goto CLEANUP;
    }
    CC_FOREACH_NODE_IN_CLIQUE (i, *cin, tmp) {
        marks[i] = 0;
    }
    for (i = 0; i < delcount; i++) {
        marks[del[i]] = 1;
    }

    CC_FOREACH_NODE_IN_CLIQUE (i, *cin, tmp) {
        if (!marks[i]) {
            ar[count++] = i;
        }
    }
    rval = CCtsp_array_to_lpclique (ar, count, cout);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (ar, int);
    CC_IFFREE (marks, char);
    return rval;
}

void CCtsp_print_lpcut_in (CCtsp_lpcut_in *c)
{
    int i, j, k;

    printf ("%d %d\n", c->cliquecount, c->rhs);
    for (i = 0; i < c->cliquecount; i++) {
        for (j = 0; j < c->cliques[i].segcount; j++) {
            for (k = c->cliques[i].nodes[j].lo;
                k <= c->cliques[i].nodes[j].hi; k++) {
                printf ("%d ", k);
            } 
        }
        printf ("-1\n");
    }
    printf ("\n");

    if (c->TP_handles != (CCtsp_lpclique *) NULL) {
        if (c->cliquecount != 0) {
            printf ("Bad Domino, should not have any cliques\n");
        } else {
            printf ("Triomino Inequality\n");
            for (i = 0; i < 2; i++) {
                printf ("      Handle %d\n", i);
                printf ("      ");
                CCtsp_print_lpclique (&(c->TP_handles[i]));
            }
            printf ("      %d dominos for handle 0\n", c->TP_domcount0);
            for (i = 0; i < c->TP_domcount0; i++) {
                printf ("      ");
                CCtsp_print_lpdomino (&(c->TP_dominos0[i]));
            }
            printf ("      %d dominos for handle 1\n", c->TP_domcount1);
            for (i = 0; i < c->TP_domcount1; i++) {
                printf ("      ");
                CCtsp_print_lpdomino (&(c->TP_dominos1[i]));
            }
            printf ("      %d triominos\n", c->TP_tricount);
            for (i = 0; i < c->TP_tricount; i++) {
                printf ("      ");
                CCtsp_print_lpclique (&(c->TP_tsets[i]));
                printf ("      ");
                CCtsp_print_lpdomino (&(c->TP_semicuts0[i]));
                printf ("      ");
                CCtsp_print_lpdomino (&(c->TP_semicuts1[i]));
            }
        }

    } else if (c->dominocount > 0) {
        if (c->cliquecount != 1) {
            printf ("Bad Domino, more than one handle\n");
        } else {
            printf ("Domino Inequality\n");
            printf ("      ");
            CCtsp_print_lpclique (&(c->cliques[0]));
            for (i = 0; i < c->dominocount; i++) {
                printf ("      ");
                CCtsp_print_lpdomino (&(c->dominos[i]));
            }
        }
    } else {
        if (c->cliquecount == 1) {
            printf ("Subtour\n");
            printf ("      ");
            CCtsp_print_lpclique (&(c->cliques[0]));
        } else {
            printf ("Comb, Clique Tree or Wild Thing (rhs %d)\n", c->rhs);
            for (i = 0; i < c->cliquecount; i++) {
                printf ("      ");
                CCtsp_print_lpclique (&(c->cliques[i]));
            }
        }
    }
    printf ("\n"); fflush (stdout);
}

void CCtsp_print_lpclique (CCtsp_lpclique *c)
{
    int i;

    if (c->segcount == 0) {
        printf ("Empty Clique\n"); fflush (stdout);
    } else {
        for (i = 0; i < c->segcount; i++) {
            printf ("%d->%d ", c->nodes[i].lo, c->nodes[i].hi);
        }
        printf ("\n"); fflush (stdout);
    }
}

void CCtsp_print_lpdomino (CCtsp_lpdomino *d)
{
    int i, k;
    CCtsp_lpclique *c;

    for (k = 0; k < 2; k++) { 
        c = &(d->sets[k]);
        if (c->segcount == 0) {
            printf ("Empty Clique "); fflush (stdout);
        } else {
            for (i = 0; i < c->segcount; i++) {
                printf ("%d->%d ", c->nodes[i].lo, c->nodes[i].hi);
            }
        }
        if (k == 0) printf (" |  ");
    }
   printf ("\n"); fflush (stdout);
}

int CCtsp_lpcut_to_lpcut_in (CCtsp_lpcuts *cuts, CCtsp_lpcut *c,
        CCtsp_lpcut_in *new)
{
    int i;
    CCtsp_lpclique *cl;
    CCtsp_lpdomino *dom;
    int rval = 0;

    CCtsp_init_lpcut_in (new);
    
    new->rhs = c->rhs;
    new->sense = c->sense;
    new->branch = c->branch;
    new->next =  (CCtsp_lpcut_in *) NULL;
    new->prev = (CCtsp_lpcut_in *) NULL;

    new->cliques = CC_SAFE_MALLOC (c->cliquecount, CCtsp_lpclique);
    CCcheck_NULL (new->cliques, "out of memory for cliques");
    if (c->dominocount > 0) {
        new->dominos = CC_SAFE_MALLOC (c->dominocount, CCtsp_lpdomino);
        CCcheck_NULL (new->dominos, "out of memory for dominos");
    }

    for (i = 0; i < c->cliquecount; i++) {
        cl = &(cuts->cliques[c->cliques[i]]);
        rval = CCtsp_copy_lpclique (cl, &new->cliques[i]);
        CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
        new->cliquecount++;
    }

    if (c->dominocount > 0) {
        for (i = 0; i < c->dominocount; i++) {
            dom = &(cuts->dominos[c->dominos[i]]);
            rval = CCtsp_copy_lpdomino (dom, &new->dominos[i]);
            CCcheck_rval (rval , "CCtsp_copy_lpdomino failed");
            new->dominocount++;
        }
    }

    if (c->TP_handles != (int *) NULL) {
        CC_MALLOC (new->TP_handles, 2, CCtsp_lpclique);
        for (i = 0; i < 2; i++) {
            cl = &(cuts->cliques[c->TP_handles[i]]);
            rval = CCtsp_copy_lpclique (cl, &new->TP_handles[i]);
            CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
        }
        if (c->TP_domcount0 > 0) {
            CC_MALLOC (new->TP_dominos0, c->TP_domcount0, CCtsp_lpdomino);
            for (i = 0; i < c->TP_domcount0; i++) {
                dom = &(cuts->dominos[c->TP_dominos0[i]]);
                rval = CCtsp_copy_lpdomino (dom, &new->TP_dominos0[i]);
                CCcheck_rval (rval , "CCtsp_copy_lpdomino failed");
                new->TP_domcount0++;
            }
        }
        if (c->TP_domcount1 > 0) {
            CC_MALLOC (new->TP_dominos1, c->TP_domcount1, CCtsp_lpdomino);
            for (i = 0; i < c->TP_domcount1; i++) {
                dom = &(cuts->dominos[c->TP_dominos1[i]]);
                rval = CCtsp_copy_lpdomino (dom, &new->TP_dominos1[i]);
                CCcheck_rval (rval , "CCtsp_copy_lpdomino failed");
                new->TP_domcount1++;
            }
        }
        if (c->TP_tricount > 0) {
            CC_MALLOC (new->TP_tsets, c->TP_tricount, CCtsp_lpclique);
            CC_MALLOC (new->TP_semicuts0, c->TP_tricount, CCtsp_lpdomino);
            CC_MALLOC (new->TP_semicuts1, c->TP_tricount, CCtsp_lpdomino);
            for (i = 0; i < c->TP_tricount; i++) {
                cl = &(cuts->cliques[c->TP_tsets[i]]);
                rval = CCtsp_copy_lpclique (cl, &new->TP_tsets[i]);
                CCcheck_rval (rval, "CCtsp_copy_lpclique failed");
                dom = &(cuts->dominos[c->TP_semicuts0[i]]);
                rval = CCtsp_copy_lpdomino (dom, &new->TP_semicuts0[i]);
                CCcheck_rval (rval , "CCtsp_copy_lpdomino failed");
                dom = &(cuts->dominos[c->TP_semicuts1[i]]);
                rval = CCtsp_copy_lpdomino (dom, &new->TP_semicuts1[i]);
                CCcheck_rval (rval , "CCtsp_copy_lpdomino failed");
                new->TP_tricount++;
            }
        }
    }

    rval = CCtsp_copy_skeleton (&c->skel, &new->skel);
    CCcheck_rval (rval, "CCtsp_copy_skeleton failed");

    if (new->cliquecount != c->cliquecount ||
        new->dominocount != c->dominocount) {
        fprintf (stderr, "error in counts\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    if (rval) CCtsp_free_lpcut_in (new);
    return rval;
}

int CCtsp_copy_lpclique (CCtsp_lpclique *c, CCtsp_lpclique *new)
{
    int rval = 0, k;
    CCtsp_segment *s = (CCtsp_segment *) NULL;

    CCtsp_init_lpclique (new);
    if (c->segcount) {
        CC_MALLOC (s, c->segcount, CCtsp_segment);
        for (k = 0; k < c->segcount; k++) {
            s[k].lo = c->nodes[k].lo;
            s[k].hi = c->nodes[k].hi;
        }
    }
    new->segcount = c->segcount;
    new->nodes = s;
CLEANUP:
    return rval;
}

static int copy_lpclique_global (CCtsp_lpclique *c, CCtsp_lpclique *new,
        int *global_names, int *local_perm, int *global_invperm)
{
    int rval = 0;
    int i, acount;
    int *arr = (int *) NULL;

    CCtsp_init_lpclique (new);

    rval = CCtsp_clique_to_array (c, &arr, &acount);
    CCcheck_rval (rval, "CCtsp_clique_to_array failed");

    for (i = 0; i < acount; i++) {
        arr[i] = global_invperm[global_names[local_perm[arr[i]]]];
    }

    rval = CCtsp_array_to_lpclique (arr, acount, new);
    CCcheck_rval (rval, "CCtsp_array_lpclique failed");

CLEANUP:

    CC_IFFREE (arr, int);
    return rval;
}

int CCtsp_copy_lpdomino (CCtsp_lpdomino *c, CCtsp_lpdomino *new)
{
    int k;
    int rval = 0;
     
    CCtsp_init_lpdomino (new);
    for (k = 0; k < 2; k++) {
        rval = CCtsp_copy_lpclique (&(c->sets[k]), &(new->sets[k]));
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpclique failed\n");
            CCtsp_free_lpdomino (new);
            goto CLEANUP;
        }
    }

CLEANUP:
    return rval;
}

void CCtsp_lpclique_compare (CCtsp_lpclique *a, CCtsp_lpclique *b, int *diff)
{
    int i;

    if (a->segcount != b->segcount) {
        *diff = 1; return;
    } else {
        for (i = 0; i < a->segcount; i++) {
            if (a->nodes[i].lo != b->nodes[i].lo) {
                *diff = 1; return;
            }
            if (a->nodes[i].hi != b->nodes[i].hi) {
                *diff = 1; return;
            }
        }
    }
    *diff = 0; return;
}

int CCtsp_create_lpcliques (CCtsp_lpcut_in *c, int cliquecount)
{
    int i;

    c->cliques = CC_SAFE_MALLOC (cliquecount, CCtsp_lpclique);
    if (c->cliques == (CCtsp_lpclique *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_create_lpcliques\n");
        return 1;
    }
    for (i=0; i<cliquecount; i++) {
        CCtsp_init_lpclique (&c->cliques[i]);
    }
    c->cliquecount = cliquecount;
    return 0;
}

int CCtsp_max_node (CCtsp_lpcut_in *c)
{
    int i, j, k;
    int maxn;
    int cliquecount = c->cliquecount;
    CCtsp_lpclique *cliques = c->cliques;
    int dominocount = c->dominocount;
    CCtsp_lpdomino *dominos = c->dominos;

    maxn = 0;

    for (i = 0; i < cliquecount; i++) {
        for (j = 0; j < cliques[i].segcount; j++) {
            if (cliques[i].nodes[j].hi > maxn) {
                maxn = cliques[i].nodes[j].hi;
            }
        }
    }

    for (i = 0; i < dominocount; i++) {
        for (k = 0; k < 2; k++) {
            for (j = 0; j < dominos[i].sets[k].segcount; j++) {
                if (dominos[i].sets[k].nodes[j].hi > maxn) {
                    maxn = dominos[i].sets[k].nodes[j].hi;
                }
            }
        }
    }

    return maxn;
}

int CCtsp_build_dp_cut (CCtsp_lpcut_in **cut, int ndomino, int *Acount,
        int **A, int *Bcount, int **B, int handlecount, int *handle)
{
    int rval = 0;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    *cut = (CCtsp_lpcut_in *) NULL;

    if (ndomino == 0) {
        printf ("Cannot build DP cut with no dominos\n");
        rval = 1;  goto CLEANUP;
    }
    if (handlecount == 0) {
        printf ("Skip empty handle DP cut\n"); fflush (stdout);
        goto CLEANUP;
    }

    c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (c, "out of memory in CCtsp_build_dp_cut");
    CCtsp_init_lpcut_in (c);

    c->cliquecount = 1;
    c->cliques = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    CCcheck_NULL (c->cliques, "out of memory in CCtsp_build_dp_cut");
    CCtsp_init_lpclique (&(c->cliques[0]));

    rval = CCtsp_array_to_lpclique (handle, handlecount, &(c->cliques[0]));
    CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

    c->dominocount = ndomino;
    rval = build_dominos (&c->dominos, ndomino, Acount, A, Bcount, B);
    CCcheck_rval (rval, "build_dominos failed");

    c->rhs = CCtsp_DOMINORHS(c);
    c->sense = 'G';
    c->branch = 0;

    *cut = c;

CLEANUP:

    if (rval) {
        if (c) {
            CCtsp_free_lpcut_in (c);
            CC_FREE (c, CCtsp_lpcut_in);
        }
    }
    return rval;
}

static int build_dominos (CCtsp_lpdomino **pdom, int ndomino, int *Acount,
        int **A, int *Bcount, int **B)
{
    int i, j, Amin, Bmin, tag, rval = 0;
    CCtsp_lpdomino *dominos = (CCtsp_lpdomino *) NULL;

    *pdom = (CCtsp_lpdomino *) NULL;

    if (ndomino == 0) {
        fprintf (stderr, "ndomino = 0 in build_dominos\n");
        rval = 1;  goto CLEANUP;
    }

    dominos = CC_SAFE_MALLOC (ndomino, CCtsp_lpdomino);
    CCcheck_NULL (dominos, "out of memory in build_dominos");
    for (i = 0; i < ndomino; i++) {
        CCtsp_init_lpdomino (&dominos[i]);
    }

    for (i = 0; i < ndomino; i++) {
        Amin = Bmin = CCutil_MAXINT;
        for (j = 0; j < Acount[i]; j++) {
            if (A[i][j] < Amin) Amin = A[i][j];
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (B[i][j] < Bmin) Bmin = B[i][j];
        }
        if (Amin < Bmin) tag = 0;
        else             tag = 1;

        rval = CCtsp_array_to_lpclique (A[i], Acount[i],
                                        &dominos[i].sets[tag]);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");

        rval = CCtsp_array_to_lpclique (B[i], Bcount[i],
                                        &dominos[i].sets[1-tag]);
        CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
    }
    *pdom = dominos;

CLEANUP:

    if (rval) {
        CC_IFFREE (dominos, CCtsp_lpdomino);  /* Misses the contents */
    }
    return rval;
}

#if 0
static int get_complement_array (int *marks, int fullcount, int *fulllist,
       int acount, int *alist, int *newcount, int **newlist)
{
    int rval = 0;
    int k = 0;
    int i;
    int *nlist = (int *) NULL;

    *newcount = 0;
    *newlist  = (int *) NULL;

    for (i = 0; i < fullcount; i++) {
        marks[fulllist[i]] = 1;
    }
    for (i = 0; i < acount; i++) {
        if (marks[alist[i]] != 1) {
            printf ("duplicate or outside node in alist\n"); fflush (stdout);
            rval = 1;  goto CLEANUP;
        }
        marks[alist[i]] = 0;
    }
    for (i = 0; i < fullcount; i++) {
        if (marks[fulllist[i]] == 1) k++;
    }

    if (k == 0) goto CLEANUP;

    nlist = CC_SAFE_MALLOC (k, int);
    CCcheck_NULL (nlist, "out of memory for nlist");

    k = 0;
    for (i = 0; i < fullcount; i++) {
        if (marks[fulllist[i]] == 1) {
            nlist[k] = fulllist[i];
            k++;
        }
    }

    *newcount = k;
    *newlist = nlist;

CLEANUP:

    for (i = 0; i < fullcount; i++) {
        marks[fulllist[i]] = 0;
    }

    return rval;
}
#endif

int CCtsp_extract_dp_cut (CCtsp_lpcut_in *c, int *pndomino, int **pAcount,
        int ***pA, int **pBcount, int ***pB, int *phandlecount, int **phandle)
{
    int rval = 0;
    int ndomino, handlecount;
    int *Acount = (int *) NULL;
    int *Bcount = (int *) NULL;
    int **A = (int **) NULL;
    int **B = (int **) NULL;
    int *handle = (int *) NULL;
    int i;

    if (c->dominocount <= 0) {
        fprintf (stderr, "no dominos to extract\n");
        rval = 1;  goto CLEANUP;
    }
    if (c->cliquecount != 1) {
        fprintf (stderr, "no domino-handle to extract\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCtsp_clique_to_array (&c->cliques[0], &handle, &handlecount);
    CCcheck_rval (rval, "CCtsp_clique_to_array failed");

    ndomino = c->dominocount;
    Acount = CC_SAFE_MALLOC (ndomino, int);
    CCcheck_NULL (Acount, "out of memory for Acount");
    Bcount = CC_SAFE_MALLOC (ndomino, int);
    CCcheck_NULL (Acount, "out of memory for Bcount");
    A = CC_SAFE_MALLOC (ndomino, int *);
    CCcheck_NULL (Acount, "out of memory for A");
    B = CC_SAFE_MALLOC (ndomino, int *);
    CCcheck_NULL (Acount, "out of memory for B");

    for (i = 0; i < ndomino; i++) {
        rval = CCtsp_clique_to_array (&c->dominos[i].sets[0],
                     &A[i], &Acount[i]);
        CCcheck_rval (rval, "CCtsp_clique_to_array failed");
        rval = CCtsp_clique_to_array (&c->dominos[i].sets[1],
                     &B[i], &Bcount[i]);
        CCcheck_rval (rval, "CCtsp_clique_to_array failed");
    }

    *pndomino     = ndomino;
    *pAcount      = Acount;
    *pA           = A;
    *pBcount      = Bcount;
    *pB           = B;
    *phandlecount = handlecount;
    *phandle      = handle;

CLEANUP:

    if (rval) {
        CC_IFFREE (Acount, int);
        CC_IFFREE (Bcount, int);
        CC_IFFREE (A, int *);  /* Note:  Misses the sets */
        CC_IFFREE (B, int *);  /* Note:  Misses the sets */
        CC_IFFREE (handle, int);
    }

    return rval;
}

#if 0
static void array_compare (int acount, int *a, int bcount, int *b, int *diff)
{
    int i;

    if (acount != bcount)  {
        *diff = 1; return;
    }

    for (i = 0; i < acount; i++) {
        if (a[i] != b[i]) {
            *diff = 1; return;
        }
    }

    *diff = 0;
}

static int array_union (int acount, int *a, int bcount, int *b, int *newcount,
        int **new)
{
    int i;
    int rval = 0;

    *newcount = acount+bcount;
    *new = CC_SAFE_MALLOC (acount+bcount, int);
    CCcheck_NULL ((*new), "out of memory for new");

    for (i = 0; i < acount; i++) {
        (*new)[i] = a[i];
    }
    for (i = acount; i < acount+bcount; i++) {
        (*new)[i] = b[i-acount];
    }

CLEANUP:

    if (rval) {
        *newcount = 0;
        *new = (int *) NULL;
    }

    return rval;
}
#endif
