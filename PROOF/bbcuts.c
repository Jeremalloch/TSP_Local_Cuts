/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2007 by David Applegate, Robert Bixby,              */
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
/*  Code for checking validity of cuts (modified version of verify.c)       */ 
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 19, 2007                                                    */ 
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int BBcuts_verify (BBtsp_lpcut_in *cut, int check_types, int *type,     */
/*      BBdatagoup *getdat)                                                 */
/*    getdata -- just return a dat group for the TSP cut (can be NULL)      */
/*                                                                          */
/*  int BBcuts_backbone (BBtsp_lpcut_in *cut, BBtsp_lpcut_in *new,          */
/*      int *no_outside)                                                    */
/*                                                                          */
/****************************************************************************/

#include "bbproof.h"

#define BB_HELDKARP_LIMIT -1

typedef struct atom {
    int mark;
    int representative;
    int cliquecount;
    int *cliquelist;
} atom;

typedef struct vclique {
    int mark;
    int atomcount;
    int *atomlist;
    int crosscount;
    int *crosslist;
    int color;
    int flipped;
    struct vclique *parent;
    struct vclique *child;
    struct vclique *sibling;
    struct vclique **prevsib;
} vclique;

#define CHILD_ADD(n,c) {                                       \
    (c)->sibling = (n)->child;                                 \
    if ((c)->sibling) (c)->sibling->prevsib = &((c)->sibling); \
    (c)->prevsib = &((n)->child);                              \
    (n)->child = (c);                                          \
    (c)->parent = (n);                                         \
}

#define CHILD_DEL(c) {                                         \
    *((c)->prevsib) = (c)->sibling;                            \
    if ((c)->sibling) (c)->sibling->prevsib = (c)->prevsib;    \
}

typedef struct node {
    int atom;
} node;

typedef struct graph {
    int nodecount;
    int edgecount;
    int *elist;
    int *elen;
} graph;

typedef struct atom_info {
    int atomcount;
    int cliquecount;
    int nodespace;
    char sense;
    int rhs;
    atom *atomlist;
    node *nodelist;
    vclique *cliquelist;
    vclique family[2]; /* dummy roots for the nested families of cliques */
} atom_info;

static int
    build_atom_info (BBtsp_lpcut_in *cut, atom_info *atoms),
    build_atoms (BBtsp_lpcut_in *cut, atom_info *atoms),
    build_clique_atoms (BBtsp_lpclique *clique, atom_info *atoms, int cnum),
    build_atom_cliques (atom_info *atoms),
    try_numbering (atom_info *atoms, int *phand, int *invperm),
    verify_atom_info (BBtsp_lpcut_in *cut, atom_info *atoms),
    verify_nodes (BBtsp_lpcut_in *cut, atom_info *atoms),
    verify_atom_clique (BBtsp_lpclique *clique, atom_info *atoms, int cnum),
    verify_atom_lists (atom_info *atoms),
    find_clique_atom (atom_info *atoms, int cnum, int anum),
    find_atom_clique (atom_info *atoms, int anum, int cnum),
    verify_subtour (atom_info *atoms),
    atom_singleton (atom *a, vclique *cliquelist, int c, int nflipped),
    atom_doubleton (atom *a, vclique *cliquelist, int c1, int c2,
        int nflipped),
    verify_chvatal_comb (atom_info *atoms),
    verify_clean_tooth (atom_info *atoms, int tooth, int handle, int nflipped),
    build_binesting (atom_info *atoms),
    build_cross_cliques (atom_info *atoms, int cliq),
    binesting_dfs (atom_info *atoms, int cliq, int color),
    nest_atom_color (atom_info *atoms, int color),
    nest_clique (atom_info *atoms, int cliqnum, vclique *family),
    clique_complement (atom_info *atoms, vclique *cliq, int *catomcount,
        int **catomlist),
    verify_star (atom_info *atoms, int handnum),
    verify_bipartition (atom_info *atoms, int handnum),
    flip_disjoint (vclique *family),
    flip_disjoint_one (vclique *family),
    flip_concentric (vclique *family),
    child_count (vclique *root),
    family_count (vclique *root, int *nflipped),
    check_mark_concentric (vclique *root, int mark),
    collect_atom_mark_counts (atom_info *atoms, vclique *c, int *alist,
        int *p_acnt, int nflipped),
    verify_other (atom_info *atoms, BBdatagroup *getdat),
    build_complete_graph (int nodecount, graph **p_g),
    verify_rhs (atom_info *atoms, graph *g),
    tooth_has_cavity (atom_info *atoms, vclique *c, int nflipped),
    max_node (BBtsp_lpcut_in *c),
    test_dp_cut (int ncount, BBtsp_lpcut_in *c, int *yes_no),
    extract_dp_cut (BBtsp_lpcut_in *c, int *pndomino, int **pAcount,
        int ***pA, int **pBcount, int ***pB, int *phandlecount, int **phandle),
    check_raw_domino (int ncount, int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B, int *valid, int *comb);

static void
    free_atom_info (atom_info *atoms),
    init_atom_info (atom_info *atoms),
    flip_marked_cliques (vclique *family),
    flip_clique (vclique *family, vclique *cliq),
    free_graph (graph **p_g),
    compute_lhs (atom_info *atoms, graph *g),
    cleanup_atom_info (BBtsp_lpcut_in *cut, atom_info *atoms),
    clear_atom_marks (atom_info *atoms),
    clear_clique_marks (atom_info *atoms);

static vclique
   *concentric_innermost (vclique *root);


int BBcuts_verify (BBtsp_lpcut_in *cut, int ncount, int check_types, int *type,
        BBdatagroup *getdat)
{
    atom_info atoms;
    int rval = 0;

    init_atom_info (&atoms);
 
    if (type != (int *) NULL) *type = -1;
    if (getdat && cut->dominocount > 0) {
        fprintf (stderr, "cannot get a cut-TSP for a dp-inequality\n");
        rval = -1; goto CLEANUP;
    }

    if (cut->dominocount > 0) {
        int dvalid;
        rval = test_dp_cut (ncount, cut, &dvalid);
        if (dvalid) {
            if (type != (int *) NULL) *type = BB_TYPE_DOMINO;
            goto CLEANUP;
        } else {
            fprintf (stderr, "Unable to verify domino parity cut\n");
            rval = -1;
            goto CLEANUP;
        }
    }
    
    if (cut->skel.atomcount == 0 || cut->skel.atoms == (int *) NULL) {
        fprintf (stderr, "Cut to verify has no skeleton\n");
        rval = 1; goto CLEANUP;
    }
    
    rval = build_atom_info (cut, &atoms);
    if (rval) {
        fprintf (stderr, "build_atom_info failed\n");
        goto CLEANUP;
    }

    rval = verify_atom_info (cut, &atoms);
    if (rval) {
        fprintf (stderr, "atom_info failed verification\n");
        goto CLEANUP;
    }

    if (getdat) {
        rval = verify_other (&atoms, getdat);
        BBcheck_rval (rval, "verify_other for getdat failed");
        goto CLEANUP;
    }

    if (check_types & BB_TYPE_SUBTOUR) {
        rval = verify_subtour (&atoms);
        if (rval == 0) {
            if (type != (int *) NULL) *type = BB_TYPE_SUBTOUR;
            goto CLEANUP;
        }
    }

    if (check_types & BB_TYPE_COMB) {
        rval = verify_chvatal_comb (&atoms);
        if (rval == 0) {
            if (type != (int *) NULL) *type = BB_TYPE_COMB;
            goto CLEANUP;
        }
    }

    if (check_types & (BB_TYPE_STAR | BB_TYPE_BIPARTITION)) {
        rval = build_binesting (&atoms);
        if (rval == 0) {
            if (check_types & BB_TYPE_STAR) {
                rval = verify_star (&atoms, 0);
                if (rval == 0) {
                    if (type != (int *) NULL) *type = BB_TYPE_STAR;
                    goto CLEANUP;
                }
                
                rval = verify_star (&atoms, 1);
                if (rval == 0) {
                    if (type != (int *) NULL) *type = BB_TYPE_STAR;
                    goto CLEANUP;
                }
            }
            if (check_types & BB_TYPE_BIPARTITION) {
                rval = verify_bipartition (&atoms, 0);
                if (rval == 0) {
                    if (type != (int *) NULL) *type = BB_TYPE_BIPARTITION;
                    goto CLEANUP;
                }
                
                rval = verify_bipartition (&atoms, 1);
                if (rval == 0) {
                    if (type != (int *) NULL) *type = BB_TYPE_BIPARTITION;
                    goto CLEANUP;
                }
            }
        }
    }

    if (check_types & BB_TYPE_OTHER) {
        rval = verify_other (&atoms, (BBdatagroup *) NULL);
        if (rval == 0) {
            if (type != (int *) NULL) *type = BB_TYPE_OTHER;
            goto CLEANUP;
        }
        fprintf (stderr, "Unable to verify cut\n");
    }
    rval = -1;

CLEANUP:
    free_atom_info (&atoms);
    return rval;
}

int BBcuts_backbone (BBtsp_lpcut_in *cut, BBtsp_lpcut_in *new,
        int *no_outside)
{
    int rval = 0;
    int *hit = (int *) NULL;
    int *perm = (int *) NULL;
    int *invperm = (int *) NULL;
    int *degree = (int *) NULL;
    atom_info atoms;
    int x, i, j, hand = -1;
    int *hlist = (int *) NULL, *hperm = (int *) NULL;

    if (no_outside) *no_outside = 0;
    init_atom_info (&atoms);
    BButil_init_lpcut_in (new);

    if (cut->skel.atomcount == 0 || cut->skel.atoms == (int *) NULL) {
        fprintf (stderr, "Cut to verify has no skeleton\n");
        rval = 1; goto CLEANUP;
    }
    
    rval = build_atom_info (cut, &atoms);
    BBcheck_rval (rval, "build_atom_info failed");

    rval = verify_atom_info (cut, &atoms);
    BBcheck_rval (rval, "verify_atom_info failed");

    new->cliques = BB_SAFE_MALLOC (atoms.cliquecount, BBtsp_lpclique);
    BBcheck_NULL (new->cliques, "out of memory for new->cliques");
    for (i = 0; i < atoms.cliquecount; i++) {
        BButil_init_lpclique (&new->cliques[i]);
    }
    new->cliquecount = atoms.cliquecount;

    invperm = BB_SAFE_MALLOC (atoms.atomcount, int);
    BBcheck_NULL (invperm, "out of memory for invperm");

    rval = try_numbering (&atoms, &hand, invperm);
    BBcheck_rval (rval, "try_numbering failed");

    if (hand == -1) {
        perm = BB_SAFE_MALLOC (atoms.atomcount, int);
        BBcheck_NULL (perm, "out of memory for perm");
        for (i = 0; i < atoms.atomcount; i++) perm[i] = i;

        degree = BB_SAFE_MALLOC (atoms.atomcount, int);
        BBcheck_NULL (degree, "out of memory for degree");
        for (i = 0; i < atoms.atomcount; i++) {
            degree[i] = atoms.atomlist[i].cliquecount;
        }
        BButil_int_perm_quicksort (perm, degree, atoms.atomcount);

        for (i = 0; i < atoms.atomcount; i++) {
            invperm[perm[i]] = i;
        }
    }

    hit = BB_SAFE_MALLOC (atoms.atomcount, int);
    BBcheck_NULL (hit, "out of memory for hit");
    for (i = 0; i < atoms.atomcount; i++) hit[i] = 0;
    for (i = 0; i < atoms.atomcount; i++) {
        if (invperm[i]<0 || invperm[i] >= atoms.atomcount || hit[invperm[i]]) {
            fprintf (stderr, "Error in invperm array\n");
            rval = 1;  goto CLEANUP;
        }
        hit[invperm[i]] = 1;
    }

    hlist = BB_SAFE_MALLOC (atoms.cliquecount, int);
    BBcheck_NULL (hlist, "out of memory for hlist");
    hperm = BB_SAFE_MALLOC (atoms.cliquecount, int);
    BBcheck_NULL (hperm, "out of memory for hperm");

    for (i = 0; i < atoms.cliquecount; i++) {
        for (j = 0; j < atoms.cliquelist[i].atomcount; j++) {
            atoms.cliquelist[i].atomlist[j] =
               invperm[atoms.cliquelist[i].atomlist[j]];
        }
        BButil_int_array_quicksort (atoms.cliquelist[i].atomlist, 
                                    atoms.cliquelist[i].atomcount);
        x = atoms.cliquelist[i].atomcount;
        for (j = 0; j < atoms.cliquelist[i].atomcount; j++) {
            x = x * 4099 + atoms.cliquelist[i].atomlist[j];
        }
        hlist[i] = x;
    }

    /* Sort the cliques to enable isomorphism checking later */

    for (i = 0; i < atoms.cliquecount; i++) hperm[i] = i;
    BButil_int_perm_quicksort (hperm, hlist, atoms.cliquecount);
 
    for (i = 0; i < atoms.cliquecount; i++) {
        rval = BButil_array_to_lpclique (atoms.cliquelist[hperm[i]].atomlist,
                                        atoms.cliquelist[hperm[i]].atomcount,
                                        &(new->cliques[i]));
        BBcheck_rval (rval, "BButil_array_to_lpclique failed");
    }

    new->rhs = cut->rhs;
    new->sense = cut->sense;
    new->branch = 0;

    if (no_outside) {
        *no_outside = 1;
        for (i = 0; i < atoms.atomcount; i++) {
            if (atoms.atomlist[i].cliquecount == 0) {
                *no_outside = 0; break;
            }
        }
    }

CLEANUP:
    BB_IFFREE (hit, int);
    BB_IFFREE (perm, int);
    BB_IFFREE (degree, int);
    BB_IFFREE (invperm, int);
    BB_IFFREE (hlist, int);
    BB_IFFREE (hperm, int);
    free_atom_info (&atoms);
    return rval;
}

static int try_numbering (atom_info *atoms, int *phand, int *invperm)
{
    int acount = atoms->atomcount, ccount = atoms->cliquecount;
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int i, j, k, hand = -1, done = 0, a, c, a1, maxc, maxi, rval = 0;
    int *hit = (int *) NULL;
    int *perm = (int *) NULL;
    int *list = (int *) NULL;
    int *degree = (int *) NULL;

    *phand = -1;
    hit = BB_SAFE_MALLOC (acount, int);
    BBcheck_NULL (hit, "out of memory for hit");
    for (i = 0; i < acount; i++) hit[i] = 0;

    perm = BB_SAFE_MALLOC (acount, int);
    BBcheck_NULL (perm, "out of memory for perm");
    degree = BB_SAFE_MALLOC (acount, int);
    BBcheck_NULL (perm, "out of memory for degree");
    list = BB_SAFE_MALLOC (acount, int);
    BBcheck_NULL (list, "out of memory for list");

    maxc = 0;
    maxi = -1;
    for (i = 0; i < ccount; i++) {
        if (cliquelist[i].atomcount > maxc) {
            maxc = cliquelist[i].atomcount;
            maxi = i;
        }
    }

    if (maxc <= 3) {
        goto CLEANUP;
    } else {
        hand = maxi;
    }

    for (i = 0; i < cliquelist[hand].atomcount; i++) {
        degree[i] = atomlist[cliquelist[hand].atomlist[i]].cliquecount;
        perm[i] = i;
    }
    BButil_int_perm_quicksort (perm, degree, cliquelist[hand].atomcount);

    for (i = 0; i < cliquelist[hand].atomcount; i++) {
        list[i] = cliquelist[hand].atomlist[perm[i]];
        hit[list[i]] = 1;
        done++;
    }

    for (i = 0; i < cliquelist[hand].atomcount; i++) {
        a = list[i];
        for (j = 0; j < atomlist[a].cliquecount; j++) {
            c = atomlist[a].cliquelist[j];
            if (c != hand) {
                for (k = 0; k < cliquelist[c].atomcount; k++) {
                    a1 = cliquelist[c].atomlist[k];
                    if (hit[a1] == 0) {
                        list[done++] = a1;
                        hit[a1] = 1;
                    }
                }
            }
        }
    }

    if (done < acount) {
        for (i = 0; i < acount; i++) {
            if (hit[i] == 0) {
                list[done++] = i;
                hit[i] = 1;
            }
        }
    }

    if (done != acount) {
        printf ("Bad done count in try_numbering\n");
        rval = 1;   goto CLEANUP;
    }

    for (i = 0; i < acount; i++) hit[i] = 0;
    for (i = 0; i < acount; i++) {
        if (list[i] < 0 || list[i] >= acount || hit[list[i]] == 1) {
            printf ("Bad list in try_numbering\n");
            rval = 1;   goto CLEANUP;
        }
        hit[list[i]] = 1;
    }
    for (i = 0; i < acount; i++) invperm[list[i]] = i;
    
     *phand = hand;

CLEANUP:
    BB_IFFREE (hit, int);
    BB_IFFREE (perm, int);
    BB_IFFREE (list, int);
    BB_IFFREE (degree, int);
    return rval;
}

static int build_atom_info (BBtsp_lpcut_in *cut, atom_info *atoms)
{
    int rval,i, maxn;
    BBtsp_skeleton *skel = &cut->skel;
    int cliquecount = cut->cliquecount;

    atoms->atomlist = (atom *) NULL;
    atoms->atomcount = 0;
    atoms->nodelist = (node *) NULL;
    atoms->cliquelist = (vclique *) NULL;

    maxn = max_node (cut);
    for (i=0; i<skel->atomcount; i++) {
        if (skel->atoms[i] > maxn) maxn = skel->atoms[i];
    }
    atoms->nodespace = maxn + 1;

    atoms->nodelist = BB_SAFE_MALLOC (atoms->nodespace, node);
    if (atoms->nodelist == (node *) NULL) {
        fprintf (stderr, "Out of memory in build_atom_info\n");
        goto FAILURE;
    }

    cleanup_atom_info (cut, atoms);

    rval = build_atoms (cut, atoms);
    if (rval) {
        fprintf (stderr, "build_atoms failed\n");
        goto FAILURE;
    }
    
    atoms->cliquecount = cliquecount;
    atoms->sense = cut->sense;
    atoms->rhs = cut->rhs;

    atoms->cliquelist = BB_SAFE_MALLOC (cliquecount, vclique);
    if (atoms->cliquelist == (vclique *) NULL) {
        fprintf (stderr, "Out of memory in build_atom_info\n");
        rval = -1; goto FAILURE;
    }
    for (i=0; i<cliquecount; i++) {
        atoms->cliquelist[i].atomlist = (int *) NULL;
        atoms->cliquelist[i].atomcount = 0;
        atoms->cliquelist[i].crosslist = (int *) NULL;
        atoms->cliquelist[i].crosscount = 0;
    }

    for (i=0; i<cliquecount; i++) {
        rval = build_clique_atoms (&cut->cliques[i], atoms, i);
        if (rval) {
            fprintf (stderr, "build_clique_atoms failed\n"); goto FAILURE;
        }
    }

    rval = build_atom_cliques (atoms);
    if (rval) {
        fprintf (stderr, "build_atom_cliques failed\n"); goto FAILURE;
    }

    return 0;

FAILURE:
    free_atom_info (atoms);
    return -1;
}

static int build_atoms (BBtsp_lpcut_in *cut, atom_info *atoms)
{
    node *nodelist = atoms->nodelist;
    BBtsp_skeleton *skel = &cut->skel;
    atom *atomlist;
    int atomcount, i;

    atomcount = skel->atomcount;
    atomlist = BB_SAFE_MALLOC (atomcount, atom);
    if (atomlist == (atom *) NULL) {
        goto FAILURE;
    }

    atoms->atomcount = atomcount;
    atoms->atomlist = atomlist;

    for (i=0; i<atomcount; i++) {
        atomlist[i].mark = 0;
        atomlist[i].representative = skel->atoms[i];
        atomlist[i].cliquecount = 0;
        atomlist[i].cliquelist = (int *) NULL;
        nodelist[atomlist[i].representative].atom = i;
    }

    return 0;

FAILURE:
    BB_IFFREE (atoms->atomlist, atom);
    atoms->atomcount = 0;
    return -1;
}

static int build_clique_atoms (BBtsp_lpclique *clique, atom_info *atoms,
        int cnum)
{
    node *nodelist = atoms->nodelist;
    atom *atomlist = atoms->atomlist;
    vclique *newclique = &(atoms->cliquelist[cnum]);
    int i, j, cnt;

    clear_atom_marks (atoms);
    cnt = 0;
    BB_FOREACH_NODE_IN_CLIQUE (i, *clique, j) {
        if (nodelist[i].atom >= 0) {
            cnt++;
            atomlist[nodelist[i].atom].mark = 1;
        }
    }
    
    newclique->atomcount = cnt;
    newclique->atomlist = BB_SAFE_MALLOC (cnt, int);
    if (newclique->atomlist == (int *) NULL) {
        fprintf (stderr, "Out of memory in build_clique_atoms\n");
        return -1;
    }
    for (i=0, cnt=0; i<atoms->atomcount; i++) {
        if (atomlist[i].mark) {
            newclique->atomlist[cnt++] = i;
        }
    }
    if (cnt != newclique->atomcount) {
        fprintf (stderr, "SURPRISE - lost some atoms\n");
        BB_IFFREE (newclique->atomlist, int);
        return -1;
    }

    return 0;
}

static int build_atom_cliques (atom_info *atoms)
{
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    vclique *cliquelist = atoms->cliquelist;
    int cliquecount = atoms->cliquecount;
    int i, j, a;

    for (i=0; i<atomcount; i++) {
        atomlist[i].cliquecount = 0;
    }

    for (i=0; i<cliquecount; i++) {
        for (j=0; j<cliquelist[i].atomcount; j++) {
            atomlist[cliquelist[i].atomlist[j]].cliquecount++;
        }
    }

    for (i=0; i<atomcount; i++) {
        if (atomlist[i].cliquecount > 0) {
            atomlist[i].cliquelist = BB_SAFE_MALLOC (atomlist[i].cliquecount,
                                                     int);
            if (atomlist[i].cliquelist == (int *) NULL) {
                fprintf (stderr, "Out of memory in build_atom_cliques\n");
                goto FAILURE;
            }
        } else {
            atomlist[i].cliquelist = (int *) NULL;
        }
        atomlist[i].cliquecount = 0;
    }

    for (i=0; i<cliquecount; i++) {
        for (j=0; j<cliquelist[i].atomcount; j++) {
            a = cliquelist[i].atomlist[j];
            atomlist[a].cliquelist[atomlist[a].cliquecount++] = i;
        }
    }

    return 0;

FAILURE:
    for (i=0; i<atomcount; i++) {
        BB_IFFREE (atomlist[i].cliquelist, int);
    }
    return -1;
}

static void init_atom_info (atom_info *atoms)
{
    if (atoms) {
        atoms->atomlist = (atom *) NULL;
        atoms->cliquelist = (vclique *) NULL;
        atoms->nodelist = (node *) NULL;
        atoms->nodespace = 0;
        atoms->atomcount = 0;
        atoms->cliquecount = 0;
    }
}

static void free_atom_info (atom_info *atoms)
{
    int i;

    if (atoms->atomlist) {
        for (i=0; i<atoms->atomcount; i++) {
            BB_IFFREE (atoms->atomlist[i].cliquelist, int);
        }
        BB_FREE (atoms->atomlist, atom);
    }
    if (atoms->cliquelist) {
        for (i=0; i<atoms->cliquecount; i++) {
            BB_IFFREE (atoms->cliquelist[i].atomlist, int);
            BB_IFFREE (atoms->cliquelist[i].crosslist, int);
        }
        BB_FREE (atoms->cliquelist, vclique);
    }

    BB_IFFREE (atoms->nodelist, node);

    atoms->nodespace = 0;
    atoms->atomcount = 0;
    atoms->cliquecount = 0;
}

/*    verify_atom_info verifies that the atom information in atoms
 *    corresponds to a subset of the atoms of cut.  In particular, the atom
 *    structure has the following properties:
 * (0) Every node of every clique is < nodespace (the size of nodelist)
 * (1) Every atom has a representative node (trivial by listing the
 *     representative for each atom)
 * (2) Every node that is a representative (ie, atom >= 0) is the
 *     representative for that atom.
 * (3) The representative node for every atom listed in a clique's
 *     atomlist is inside the clique
 * (4) The representative node for every atom not listed in a clique's
 *     atomlist is outside the clique
 * (5) Every atom listed in a clique's atomlist lists the clique in its
 *     cliquelist
 * (6) Every clique listed in an atom's cliquelist lists the atom in its
 *     atomlist
 * (7) The rhs, sense, and cliquecount are the same in the cut and the atom */

static int verify_atom_info (BBtsp_lpcut_in *cut, atom_info *atoms)
{
    int i, rval = 0;

    /* (0) and (2) */
    rval = verify_nodes (cut, atoms);
    if (rval) goto CLEANUP;

    /* (3) and (4) */
    for (i=0; i<cut->cliquecount; i++) {
        rval = verify_atom_clique (&(cut->cliques[i]), atoms, i);
        if (rval) goto CLEANUP;
    }

    /* (5) and (6) */
    rval = verify_atom_lists (atoms);
    if (rval) goto CLEANUP;

    /* (7) */
    if (cut->rhs != atoms->rhs) {
        fprintf (stderr, "cut rhs %d != atom rhs %d\n", cut->rhs,
                 atoms->rhs);
        rval = -1; goto CLEANUP;
    }
    if (cut->sense != atoms->sense) {
        fprintf (stderr, "cut sense %c != atom sense %c\n", cut->sense,
                 atoms->sense);
        rval = -1; goto CLEANUP;
    }
    if (cut->cliquecount != atoms->cliquecount) {
        fprintf (stderr, "cut cliquecount %d != atom cliquecount %d\n",
                 cut->cliquecount, atoms->cliquecount);
        rval = -1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

static int verify_nodes (BBtsp_lpcut_in *cut, atom_info *atoms)
{
    node *nodelist = atoms->nodelist;
    int nodespace = atoms->nodespace;
    BBtsp_lpclique *cliques = cut->cliques;
    int cliquecount = cut->cliquecount;
    int atomcount = atoms->atomcount;
    atom *atomlist = atoms->atomlist;
    int i, k, j;

    for (i=0; i<cliquecount; i++) {
        BB_FOREACH_NODE_IN_CLIQUE (j, cliques[i], k) {
            if (j < 0 || j >= nodespace) {
                fprintf (stderr, "clique %d contains bogus node %d\n",
                         i, j);
                return -1;
            }
            if (nodelist[j].atom < -1 || nodelist[j].atom >= atomcount) {
                fprintf (stderr, "node %d is in bogus atom %d\n", j,
                         nodelist[j].atom);
                return -1;
            }
            if (nodelist[j].atom >= 0 &&
                atomlist[nodelist[j].atom].representative != j) {
                fprintf (stderr, "node %d is not the rep of its atom %d\n",
                         j, nodelist[j].atom);
                return -1;
            }
        }
    }
    return 0;
}

static int verify_atom_clique (BBtsp_lpclique *clique, atom_info *atoms,
        int cnum)
{
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    node *nodelist = atoms->nodelist;
    vclique *cliqueinfo = &(atoms->cliquelist[cnum]);
    int i, j, a, rval = 0;

    clear_atom_marks (atoms);

    BB_FOREACH_NODE_IN_CLIQUE (i, *clique, j) {
        a = nodelist[i].atom;
        if (a >= 0) {
            if (atomlist[a].mark) {
                fprintf (stderr, "Duplicate atom %d in clique\n", a);
                rval = -1; goto CLEANUP;
            }
            atomlist[a].mark = 1;
        }
    }

    for (i=0; i<cliqueinfo->atomcount; i++) {
        j = cliqueinfo->atomlist[i];
        if (atomlist[j].mark != 1) {
            fprintf (stderr, "Atom %d not in clique %d\n", j, cnum);
            rval = -1; goto CLEANUP;
        }
        atomlist[j].mark = 0;
    }

    for (i=0; i<atomcount; i++) {
        if (atomlist[i].mark != 0) {
            fprintf (stderr, "Atom %d not in cliquelist %d\n", j, cnum);
            rval = -1; goto CLEANUP;
        }
    }

CLEANUP:
    return rval;
}

static int verify_atom_lists (atom_info *atoms)
{
    int atomcount = atoms->atomcount;
    atom *atomlist = atoms->atomlist;
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    int i, j, rval;

    for (i=0; i<atomcount; i++) {
        for (j=0; j<atomlist[i].cliquecount; j++) {
            rval = find_clique_atom (atoms, atomlist[i].cliquelist[j], i);
            if (rval) {
                fprintf (stderr, "clique %d missing atom %d\n",
                         atomlist[i].cliquelist[j], i);
                return rval;
            }
        }
    }

    for (i=0; i<cliquecount; i++) {
        for (j=0; j<cliquelist[i].atomcount; j++) {
            rval = find_atom_clique (atoms, cliquelist[i].atomlist[j], i);
            if (rval) {
                fprintf (stderr, "atom %d missing clique %d\n",
                         cliquelist[i].atomlist[j], i);
                return rval;
            }
        }
    }
    return 0;
}

static int find_clique_atom (atom_info *atoms, int cnum, int anum)
{
    int i;
    vclique *cliquelist = atoms->cliquelist;

    for (i=0; i<cliquelist[cnum].atomcount; i++) {
        if (cliquelist[cnum].atomlist[i] == anum) return 0;
    }
    return 1;
}

static int find_atom_clique (atom_info *atoms, int anum, int cnum)
{
    int i;
    atom *atomlist = atoms->atomlist;

    for (i=0; i<atomlist[anum].cliquecount; i++) {
        if (atomlist[anum].cliquelist[i] == cnum) return 0;
    }
    return 1;
}

static int verify_subtour (atom_info *atoms)
{
    if (atoms->cliquecount == 1 && atoms->atomcount == 2 &&
        atoms->rhs == 2 && atoms->sense == 'G') {
        return 0;
    } else {
        return -1;
    }
}

static int atom_singleton (atom *a, vclique *cliquelist, int c, int nflipped)
{
    int i, j, ccnt = 0, fcnt = 0;

    if (cliquelist[c].flipped) {
        fcnt++;
        ccnt++;
    }
    
    for (i=0; i<a->cliquecount; i++) {
        j = a->cliquelist[i];
        if (cliquelist[j].flipped == 0) {
            if (j != c) return 0;
            ccnt++;
        } else {
            if (j == c) return 0;
            fcnt++;
        }
    }
    if (fcnt == nflipped && ccnt == 1) return 1;
    else return 0;
}

static int atom_doubleton (atom *a, vclique *cliquelist, int c1, int c2,
        int nflipped)
{
    int i, j, c1cnt = 0, c2cnt = 0, fcnt = 0;

    if (cliquelist[c1].flipped) {
        fcnt++;
        c1cnt++;
    }
    if (cliquelist[c2].flipped) {
        fcnt++;
        c2cnt++;
    }
    for (i=0; i<a->cliquecount; i++) {
        j = a->cliquelist[i];
        if (cliquelist[j].flipped == 0) {
            if      (j == c1) c1cnt++;
            else if (j == c2) c2cnt++;
            else              return 0;
        } else {
            if (j == c1 || j == c2) return 0;
            fcnt++;
        }
    }
    if (c1cnt == 1 && c2cnt == 1 && fcnt == nflipped) return 1;
    else return 0;
}

/* combs come in 4 varieties, depending on whether there is a handle cavity
 * (hc), and whether there is an outside cavity (oc).  Every comb has
 * (1) an even number of cliques (1 handle + odd teeth)
 * (2) at least 4 cliques
 * (3) rhs 3*(ncliques-1) + 1, sense 'G'
 * (4) 2*(ncliques-1) + hc + oc atoms
 * (5) a handle, which contains ncliques-1+hc atoms
 * (6) ncliques-1 teeth, which contain 2 atoms or atomcount-2 atoms
 * If we "invert" the teeth so that each tooth contains 2 atoms, then
 * (7) each tooth contains exactly 2 atoms, one of just the tooth, and one
 *     of just the tooth and the handle */

static int verify_chvatal_comb (atom_info *atoms)
{
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    int atomcount = atoms->atomcount;
    int i, handlelo, handlehi, toothlo, toothhi, handlecount, nflipped;
    int handle = -1;

    if (cliquecount % 2 != 0) return 1;
    if (cliquecount < 4) return 1;
    if (atoms->rhs != 3 * cliquecount - 2) return 1;
    if (atoms->sense != 'G') return 1;

    if (atomcount == 2*cliquecount-2) {
        handlelo = cliquecount-1;
        handlehi = cliquecount-1;
    } else if (atomcount == 2*cliquecount-1) {
        handlelo = cliquecount-1;
        handlehi = cliquecount;
    } else if (atomcount == 2*cliquecount) {
        handlelo = cliquecount;
        handlehi = cliquecount;
    } else {
        return 1;
    }

    toothlo = 2;
    toothhi = atomcount-2;

    handlecount = 0;
    nflipped = 0;
    for (i=0; i<cliquecount; i++) {
        if (cliquelist[i].atomcount >= handlelo &&
            cliquelist[i].atomcount <= handlehi) {
            cliquelist[i].flipped = 0;
            handlecount++;
            handle = i;
        } else if (cliquelist[i].atomcount == toothlo) {
            cliquelist[i].flipped = 0;
        } else if (cliquelist[i].atomcount == toothhi) {
            cliquelist[i].flipped = 1;
            nflipped++;
        } else {
            return 1;
        }
    }
    if (handlecount != 1) return 1;

    for (i=0; i<cliquecount; i++) {
        if (i != handle) {
            if (verify_clean_tooth (atoms, i, handle, nflipped)) {
                return 1;
            }
        }
    }

    return 0;
}

static int verify_clean_tooth (atom_info *atoms, int tooth, int handle,
			       int nflipped)
{
    int atomcount = atoms->atomcount;
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int cliqueatomcount = cliquelist[tooth].atomcount;
    int *cliqueatomlist = cliquelist[tooth].atomlist;
    int inside = 0, outside = 0, i, inmark;

    clear_atom_marks (atoms);

    for (i=0; i<cliqueatomcount; i++) {
        atomlist[cliqueatomlist[i]].mark = 1;
    }

    if (cliquelist[tooth].flipped == 0) {
        inmark = 1;
    } else {
        inmark = 0;
    }

    for (i=0; i<atomcount; i++) {
	if (atomlist[i].mark == inmark) {
	    if (atom_singleton (&atomlist[i], cliquelist, tooth, nflipped)) {
		if (outside == 1) return 1;
		outside = 1;
	    } else if (atom_doubleton (&atomlist[i], cliquelist, handle,
				       tooth, nflipped)) {
		if (inside == 1) return 1;
		inside = 1;
	    } else {
		return 1;
	    }
	}
    }

    if (inside == 1 && outside == 1) return 0;
    else return 1;
}

static int build_binesting (atom_info *atoms)
{
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    atom *a;
    int i, rval;

    clear_clique_marks (atoms);

    for (i=0; i<cliquecount; i++) {
        rval = build_cross_cliques (atoms, i);
        if (rval) return rval;
    }
    
    for (i=0; i<cliquecount; i++) {
        cliquelist[i].color = -1;
    }

    rval = binesting_dfs (atoms, 0, 0);
    if (rval) return rval;

    for (i=0; i<cliquecount; i++) {
        if (cliquelist[i].color != 0 && cliquelist[i].color != 1) {
            return 1;
        }
    }

    /* every clique is flipped so 0 is outside, guaranteeing that */
    /* the family nests.                                          */

    for (i=0; i<cliquecount; i++) {
        cliquelist[i].flipped = 0;
    }
    a = &atoms->atomlist[0];
    for (i=0; i<a->cliquecount; i++) {
        cliquelist[a->cliquelist[i]].flipped = 1;
    }
    
    rval = nest_atom_color (atoms, 0);
    if (rval) {
        fprintf (stderr, "nest_atom_color failed\n"); return rval;
    }
    rval = nest_atom_color (atoms, 1);
    if (rval) {
        fprintf (stderr, "nest_atom_color failed\n"); return rval;
    }
    return 0;
}

static int build_cross_cliques (atom_info *atoms, int cliq)
{
    vclique *cliquelist = atoms->cliquelist;
    vclique *c = &cliquelist[cliq];
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    int catomcount = c->atomcount;
    int *catomlist = c->atomlist;
    int crosscount = 0;
    int *crosslist = (int *) NULL;
    atom *a;
    int i, j, k, rval;

    crosslist = BB_SAFE_MALLOC (atoms->cliquecount, int);
    if (crosslist == (int *) NULL) {
        fprintf (stderr, "Out of memory in build_cross_cliques\n");
        rval = 1; goto CLEANUP;
    }
    
    for (i=0; i<catomcount; i++) {
        a = &atomlist[catomlist[i]];
        for (j=0; j<a->cliquecount; j++) {
            k = a->cliquelist[j];
            cliquelist[k].mark++;
        }
    }
    cliquelist[cliq].mark = 0;

    crosscount = 0;
    for (i=0; i<catomcount; i++) {
        a = &atomlist[catomlist[i]];
        for (j=0; j<a->cliquecount; j++) {
            k = a->cliquelist[j];
            if (cliquelist[k].mark > 0 &&
                cliquelist[k].mark < catomcount &&
                cliquelist[k].mark < cliquelist[k].atomcount &&
                catomcount + cliquelist[k].atomcount <
                cliquelist[k].mark + atomcount) {
                crosslist[crosscount++] = k;
            }
            cliquelist[k].mark = 0;
        }
    }

    if (crosscount == 0) {
        printf ("Warning: clique %d crosses no other cliques\n", cliq);
        c->crosslist = (int *) NULL;
    } else {
        c->crosslist = BB_SAFE_MALLOC (crosscount, int);
        if (c->crosslist == (int *) NULL) {
            fprintf (stderr, "Out of memory in build_cross_cliques\n");
            rval = 1; goto CLEANUP;
        }
        for (i=0; i<crosscount; i++) {
            c->crosslist[i] = crosslist[i];
        }
    }
    c->crosscount = crosscount;
    rval = 0;

CLEANUP:
    BB_IFFREE (crosslist, int);
    return rval;
}

static int binesting_dfs (atom_info *atoms, int cliq, int color)
{
    vclique *c = &atoms->cliquelist[cliq];
    int crosscount = c->crosscount;
    int *crosslist = c->crosslist;
    int i, rval;

    if (c->color == color) return 0;
    else if (c->color != -1) return 1;
    
    c->color = color;

    for (i=0; i<crosscount; i++) {
        rval = binesting_dfs (atoms, crosslist[i], 1-color);
        if (rval) return rval;
    }
    return 0;
}

static int nest_atom_color (atom_info *atoms, int color)
{
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    int atomcount = atoms->atomcount;
    vclique *family = &atoms->family[color];
    vclique *c;
    int i, colorcount, rval = 0;
    int *perm = (int *) NULL, *keys = (int *) NULL;

    perm = BB_SAFE_MALLOC (cliquecount, int);
    keys = BB_SAFE_MALLOC (cliquecount, int);
    if (perm == (int *) NULL ||
        keys == (int *) NULL) {
        fprintf (stderr, "Out of memory in nest_atom_color\n");
        rval = 1; goto CLEANUP;
    }

    family->parent = (vclique *) NULL;
    family->sibling = (vclique *) NULL;
    family->child = (vclique *) NULL;
    family->flipped = 0;
    family->mark = 0;
    family->atomcount = 0;
    family->atomlist = (int *) NULL;
    family->crosscount = 0;
    family->crosslist = (int *) NULL;
    family->color = color;
    
    colorcount = 0;
    for (i=0; i<cliquecount; i++) {
        c = &cliquelist[i];
        if (c->color == color) {
            if (c->flipped) keys[i] = atomcount - c->atomcount;
            else            keys[i] = c->atomcount;
            perm[colorcount++] = i;
        }
    }

    BButil_int_perm_quicksort (perm, keys, colorcount);

    clear_atom_marks (atoms);
    clear_clique_marks (atoms);
    
    for (i=0; i<colorcount; i++) {
        rval = nest_clique (atoms, perm[i], family);
        if (rval) goto CLEANUP;
    }

CLEANUP:
    BB_IFFREE (perm, int);
    BB_IFFREE (keys, int);
    return rval;
}

static int nest_clique (atom_info *atoms, int cliqnum, vclique *family)
{
    vclique *cliquelist = atoms->cliquelist;
    vclique *cliq = &cliquelist[cliqnum];
    atom *atomlist = atoms->atomlist;
    int catomcount, i, rval = 0;
    int *catomlist = (int *) NULL;
    vclique *c;
    atom *a;

    if (cliq->flipped) {
        rval = clique_complement (atoms, cliq, &catomcount, &catomlist);
        if (rval) {
            fprintf (stderr, "clique_complement failed\n");
            goto CLEANUP;
        }
    } else {
        catomcount = cliq->atomcount;
        catomlist = cliq->atomlist;
    }
    
    cliq->child = (vclique *) NULL;
    
    for (i=0; i<catomcount; i++) {
        a = &atomlist[catomlist[i]];
        if (a->mark != 0) {
            c = &cliquelist[a->mark-1];
            if (c->mark == 1) {
                CHILD_DEL(c);
                c->mark = 0;
                CHILD_ADD (cliq, c);
            }
        }
        a->mark = cliqnum+1;
    }
    cliq->mark = 1;
    CHILD_ADD (family,cliq);

CLEANUP:
    if (cliq->flipped) {
        BB_IFFREE (catomlist, int);
    }
    return rval;
}

static int clique_complement (atom_info *atoms, vclique *cliq, int *catomcount,
        int **catomlist)
{
    int atomcount = atoms->atomcount;
    int *cmark = (int *) NULL;
    int cnt, i, rval;

    *catomcount = atomcount - cliq->atomcount;
    *catomlist = BB_SAFE_MALLOC (*catomcount, int);
    cmark = BB_SAFE_MALLOC (atomcount, int);

    if (*catomlist == (int *) NULL ||
        cmark == (int *) NULL) {
        fprintf (stderr, "Out of memory in clique_complement\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<atomcount; i++) cmark[i] = 0;
    for (i=0; i<cliq->atomcount; i++) cmark[cliq->atomlist[i]] = 1;
    cnt = 0;
    for (i=0; i<atomcount; i++) {
        if (cmark[i] == 0) {
            (*catomlist)[cnt++] = i;
        }
    }
    if (cnt != *catomcount) {
        fprintf (stderr, "Lost some atoms\n");
        rval = 1; goto CLEANUP;
    }
    rval = 0;

CLEANUP:
    if (rval) BB_IFFREE (*catomlist, int);
    BB_IFFREE (cmark, int);
    return rval;
}

static void flip_marked_cliques (vclique *family)
{
    int cnt;
    vclique *cflip = (vclique *) NULL;
    vclique *c;

    for (;;) {
        cnt = 0;
        for (c = family->child; c; c = c->sibling) {
            if (c->mark != 0) {
                cnt++;
                cflip = c;
            }
        }
        if (cnt == 1) {
            cflip->mark = 0;
            flip_clique (family, cflip);
        } else {
            return;
        }
    }
}

static void flip_clique (vclique *family, vclique *cliq)
{
    vclique *c;

    CHILD_DEL(cliq);
    c = cliq->child;

    cliq->child = family->child;
    if (cliq->child) cliq->child->prevsib = &(cliq->child);

    cliq->sibling = c;
    if (c) c->prevsib = &(cliq->sibling);

    family->child = cliq;
    cliq->prevsib = &(family->child);

    for (c = family->child; c; c = c->sibling) {
        c->parent = family;
    }

    for (c = cliq->child; c; c = c->sibling) {
        c->parent = cliq;
    }

    cliq->flipped = ! cliq->flipped;
}

static int verify_star (atom_info *atoms, int handnum)
{
    vclique *handles = &atoms->family[handnum];
    vclique *teeth = &atoms->family[1-handnum];
    int tcnt, hcnt, tsum, tmult, acnt, i, rval, nflipped;
    vclique *c, *t;
    int *alist = (int *) NULL;

    clear_clique_marks (atoms);
    
    rval = flip_disjoint (teeth);
    if (rval) goto CLEANUP;

    rval = flip_concentric (handles);
    if (rval) goto CLEANUP;
    
    tcnt = child_count (teeth);
    if (tcnt < 3 || (tcnt % 2) != 1) {
        rval = 1; goto CLEANUP;
    }
    if (atoms->sense != 'G') {
        rval = 1; goto CLEANUP;
    }

    nflipped = 0;
    hcnt = family_count (handles, &nflipped) - 1;
    tsum = family_count (teeth, (int *) NULL) - 1;

    if (atoms->rhs != (tcnt+1) * hcnt + 2*tsum) {
        rval = 1; goto CLEANUP;
    }

    alist = BB_SAFE_MALLOC (atoms->atomcount, int);
    if (alist == (int *) NULL) {
        fprintf (stderr, "Out of memory in verify_star\n");
        rval = 1; goto CLEANUP;
    }
    
    clear_clique_marks (atoms);
    check_mark_concentric (handles->child, 1);

    for (c = teeth->child; c; c = c->sibling) {
        t = concentric_innermost (c);
        rval = collect_atom_mark_counts (atoms, t, alist, &acnt, nflipped);
        if (rval) goto CLEANUP;
        BButil_int_array_quicksort (alist, acnt);
        if (alist[0] != 0 || alist[acnt-1] != hcnt) {
            rval = 1; goto CLEANUP;
        }
        tmult = family_count (c, (int *) NULL);
        for (i=1; i<acnt; i++) {
            if (alist[i] - alist[i-1] > tmult) {
                rval = 1; goto CLEANUP;
            }
        }
    }

    rval = 0;

CLEANUP:
    BB_IFFREE (alist, int);
    return rval;
}

/* verify_bipartition rounds the required coefficient on degenerate teeth
 * touching at least three handles up to the next integer.  */

static int verify_bipartition (atom_info *atoms, int handnum)
{
    int teethnum = 1 - handnum;
    vclique *handles = &atoms->family[handnum];
    vclique *teeth = &atoms->family[teethnum];
    vclique *cliquelist = atoms->cliquelist;
    int rval, tsum, hsum, tmult, i, nflipped;
    vclique *c, *t;
    int *alist = (int *) NULL;

    clear_clique_marks (atoms);
    
    rval = flip_disjoint (teeth);
    if (rval) goto CLEANUP;

    rval = flip_disjoint_one (handles);
    if (rval) goto CLEANUP;

    if (atoms->sense != 'G') {
        rval = 1; goto CLEANUP;
    }

    clear_clique_marks (atoms);

    for (c = teeth->child; c; c = c->sibling) {
        t = concentric_innermost (c);
        for (i=0; i<t->crosscount; i++) {
            cliquelist[t->crosslist[i]].mark++;
        }
    }

    hsum = 0;
    for (c = handles->child; c; c = c->sibling) {
        tmult = family_count (c, (int *) NULL);
        if (tmult != 1 || (c->mark % 2) != 1) {
            rval = 1; goto CLEANUP;
        }
        hsum += c->mark + 1;
    }
    tsum = family_count (teeth, (int *) NULL) - 1;

    if (atoms->rhs != 2 * tsum + hsum) {
        rval = 1; goto CLEANUP;
    }

    clear_clique_marks (atoms);

    nflipped = 0;
    for (c = handles->child; c; c = c->sibling) {
        c->mark = 1;
        if (c->flipped) nflipped++;
    }

    for (c = teeth->child; c; c = c->sibling) {
        t = concentric_innermost (c);
        if (!tooth_has_cavity (atoms, t, nflipped)) {
            tmult = family_count (c, (int *) NULL);
            if (tmult < 2) {
                rval = 1; goto CLEANUP;
            }
        }
    }
    rval = 0;

CLEANUP:
    BB_IFFREE (alist, int);
    return rval;
}

static int flip_disjoint (vclique *family)
{
    int ccnt, rval;
    vclique *c;

    if (family == (vclique *) NULL || family->child == (vclique *) NULL) {
        return 1;
    }

    ccnt = child_count (family);
    if (ccnt == 1) {
        rval = check_mark_concentric (family->child, 1);
        if (rval) {
            flip_marked_cliques (family);
        }
    } else if (ccnt == 2) {
        rval = check_mark_concentric (family->child, 1);
        if (rval == 0) {
            check_mark_concentric (family->child, 0);
            rval = check_mark_concentric (family->child->sibling, 1);
        }
        if (rval) {
            flip_marked_cliques (family);
        }
    }

    for (c = family->child; c; c = c->sibling) {
        rval = check_mark_concentric (c, 0);
        if (rval) return 1;
    }

    return 0;
}

static int flip_disjoint_one (vclique *family)
{
    vclique *c;

    if (family == (vclique *) NULL || family->child == (vclique *) NULL) {
        return 1;
    }

    c = family->child;
    if (c->sibling == (vclique *) NULL) {
        if (c->child == (vclique *) NULL) {
            return 0;
        } else {
            flip_clique (family, c);
        }
    }

    for (c = family->child; c; c = c->sibling) {
        if (c->child != (vclique *) NULL) return 1;
    }

    return 0;
}

static int flip_concentric (vclique *family)
{
    int rval;

    if (family == (vclique *) NULL || family->child == (vclique *) NULL) {
        return 1;
    }

    if (child_count (family) == 2) {
        rval = check_mark_concentric (family->child, 1);
        if (rval == 0) {
            flip_marked_cliques (family);
        }
    }

    rval = check_mark_concentric (family, 0);
    return rval;
}

static int child_count (vclique *root)
{
    int n = 0;
    vclique *c;

    if (root == (vclique *) NULL) return -1;
    for (c = root->child; c; c = c->sibling) n++;

    return n;
}

/* family_count includes the root of the family.  Thus, when counting an
   entire family, subtract 1 for the dummy root */

static int family_count (vclique *root, int *nflipped)
{
    int n;
    vclique *c;

    if (root == (vclique *) NULL) return -1;

    n = 1;
    if (root->flipped && nflipped != (int *) NULL) (*nflipped) += 1;
    for (c = root->child; c; c = c->sibling) {
        n += family_count (c, nflipped);
    }
    return n;
}

static int check_mark_concentric (vclique *root, int mark)
{
    vclique *c;

    if (root == (vclique *) NULL) return 0;

    root->mark = mark;
    for (c = root->child; c; c = c->child) {
        if (c->sibling != (vclique *) NULL) return 1;
        c->mark = mark;
    }
    return 0;
}

static vclique *concentric_innermost (vclique *root)
{
    if (root == (vclique *) NULL) return (vclique *) NULL;
    while (root->child != (vclique *) NULL) {
        root = root->child;
    }
    return root;
}

static int collect_atom_mark_counts (atom_info *atoms, vclique *cliq,
        int *alist, int *p_acnt, int nflipped)
{
    int acnt = 0;
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int *catomlist = (int *) NULL;
    int *acliquelist;
    int catomcount, acliquecount, i, j, a, cnt, fcnt, rval;

    if (cliq->flipped) {
        rval = clique_complement (atoms, cliq, &catomcount, &catomlist);
        if (rval) {
            fprintf (stderr, "clique_complement failed\n");
            goto CLEANUP;
        }
    } else {
        catomcount = cliq->atomcount;
        catomlist = cliq->atomlist;
    }
        
    for (i=0; i<catomcount; i++) {
        a = catomlist[i];
        acliquecount = atomlist[a].cliquecount;
        acliquelist = atomlist[a].cliquelist;
        cnt = 0;
        fcnt = 0;
        for (j=0; j<acliquecount; j++) {
            if (cliquelist[acliquelist[j]].mark) {
                if (cliquelist[acliquelist[j]].flipped) fcnt++;
                else cnt++;
            }
        }
        alist[acnt++] = cnt + (nflipped - fcnt);
    }
    *p_acnt = acnt;

    rval = 0;

CLEANUP:
    if (cliq->flipped) {
        BB_IFFREE (catomlist, int);
    }
    return rval;
}

static int tooth_has_cavity (atom_info *atoms, vclique *cliq, int nflipped)
{
    atom *atomlist = atoms->atomlist;
    vclique *cliquelist = atoms->cliquelist;
    int *catomlist = (int *) NULL;
    int *acliquelist;
    int catomcount, acliquecount, i, j, a, cnt, fcnt, rval;

    if (cliq->flipped) {
        rval = clique_complement (atoms, cliq, &catomcount, &catomlist);
        if (rval) {
            fprintf (stderr, "clique_complement failed\n");
            goto CLEANUP;
        }
    } else {
        catomcount = cliq->atomcount;
        catomlist = cliq->atomlist;
    }
        
    for (i=0; i<catomcount; i++) {
        a = catomlist[i];
        acliquecount = atomlist[a].cliquecount;
        acliquelist = atomlist[a].cliquelist;
        cnt = 0;
        fcnt = 0;
        for (j=0; j<acliquecount; j++) {
            if (cliquelist[acliquelist[j]].mark) {
                if (cliquelist[acliquelist[j]].flipped) fcnt++;
                else cnt++;
            }
        }
        if (cnt + (nflipped - fcnt) == 0) {
            rval = 1; /* cavity found */
            goto CLEANUP;
        }
    }

    rval = 0;

CLEANUP:
    if (cliq->flipped) {
        BB_IFFREE (catomlist, int);
    }
    return rval;
}

static int verify_other (atom_info *atoms, BBdatagroup *getdat)
{
    int rval;
    graph *g = (graph *) NULL;

    rval = build_complete_graph (atoms->atomcount, &g);
    BBcheck_rval (rval, "build_complete_graph failed");
    compute_lhs (atoms, g);

    if (getdat) {
        rval = BButil_graph2dat_matrix (g->nodecount, g->edgecount,
                       g->elist, g->elen, 1000000, getdat);
        BBcheck_rval (rval, "BBtsp_graph2dat_matrix failed");
    } else {
        rval = verify_rhs (atoms, g);
        if (rval) goto CLEANUP;
    }

CLEANUP:
    if (g) free_graph (&g);
    return rval;
}

static int build_complete_graph (int nodecount, graph **p_g)
{
    graph *g = (graph *) NULL;
    int edgecount = (nodecount * (nodecount - 1)) / 2;
    int edgecount2;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int i, j; 

    g = BB_SAFE_MALLOC (1, graph);
    if (!g) {
        fprintf (stderr, "Out of memory in build_graph\n"); goto FAILURE;
    }
    g->nodecount = nodecount;
    g->edgecount = edgecount;
    g->elist = (int *) NULL;
    g->elen = (int *) NULL;

    g->elist = (int *) BB_SAFE_MALLOC (edgecount*2, int);
    if (!g->elist) {
        fprintf (stderr, "Out of memory in build_graph\n"); goto FAILURE;
    }

    g->elen = (int *) BB_SAFE_MALLOC (edgecount, int);
    if (!g->elen) {
        fprintf (stderr, "Out of memory in build_graph\n"); goto FAILURE;
    }

    edgecount2 = 0;
    elist = g->elist;
    elen = g->elen;
    for (i=0; i<nodecount; i++) {
        for (j=i+1; j<nodecount; j++) {
            elist[2*edgecount2] = i;
            elist[2*edgecount2+1] = j;
            elen[edgecount2] = 0;
            edgecount2++;
        }
    }
    if (edgecount2 != edgecount) {
        fprintf (stderr, "ERROR - found %d edges, expected %d\n",
                 edgecount2, edgecount);
        goto FAILURE;
    }

    *p_g = g;
    return 0;

FAILURE:
    if (g) {
        BB_IFFREE (g->elen, int);
        BB_IFFREE (g->elist, int);
        BB_FREE (g, graph);
    }
    *p_g = (graph *) NULL;
    return -1;
}

static void free_graph (graph **p_g)
{
    if (*p_g) {
        BB_IFFREE ((*p_g)->elist, int);
        BB_IFFREE ((*p_g)->elen, int);
        BB_FREE (*p_g, graph)
    }
}

static void compute_lhs (atom_info *atoms, graph *g)
{
    int edgecount = g->edgecount;
    int *elen = g->elen;
    int *elist = g->elist;
    int cliquecount = atoms->cliquecount;
    vclique *cliquelist = atoms->cliquelist;
    atom *atomlist = atoms->atomlist;
    vclique *clique;
    int i, j;

    for (i=0; i<edgecount; i++) elen[i] = 0;
    clear_atom_marks (atoms);
    for (i=0; i<cliquecount; i++) {
        clique = &cliquelist[i];
        for (j=0; j<clique->atomcount; j++) {
            atomlist[clique->atomlist[j]].mark = 1;
        }
        for (j=0; j<edgecount; j++) {
            if (atomlist[elist[2*j]].mark != atomlist[elist[2*j+1]].mark) {
                elen[j]++;
            }
        }
        for (j=0; j<clique->atomcount; j++) {
            atomlist[clique->atomlist[j]].mark = 0;
        }
    }
}

static int verify_rhs (atom_info *atoms, graph *g)
{
    int rval, foundtour;
    double upbound, optval;
    BBdatagroup vdat;

    BButil_init_datagroup (&vdat);

    if (atoms->sense != 'G') {
        fprintf (stderr, "Constraint is not a >=\n");
        rval = 1; goto CLEANUP;
    }

    upbound = (double) atoms->rhs;
    rval = BBheld_small_elist (g->nodecount, g->edgecount, g->elist,
                               g->elen, &upbound, &optval, &foundtour,
                               1, (int *) NULL, BB_HELDKARP_LIMIT, 1);
    if (rval == BB_HELDKARP_SEARCHLIMITEXCEEDED) {
        fprintf (stderr, "BBheld_small_elist search limit exceeded\n");
        rval = 1;
    } else if (rval) {
        fprintf (stderr, "BBheld_small_elist failed\n");
        rval = 1;
    } else if (foundtour) {
        fprintf (stderr, "BBheld_small_elist found better tour\n");
        rval = 1;
    } else {
        printf ("Held-Karp verified the cut\n"); fflush (stdout);
        rval = 0;
    }

CLEANUP:
    BButil_freedatagroup (&vdat);
    return rval;
}
      
static void cleanup_atom_info (BBtsp_lpcut_in *cut, atom_info *atoms)
{
    node *nodelist = atoms->nodelist;
    BBtsp_lpclique *cliques = cut->cliques;
    int cliquecount = cut->cliquecount;
    int i, k, j;

    for (i=0; i<cliquecount; i++) {
        BB_FOREACH_NODE_IN_CLIQUE (j, cliques[i], k) {
            nodelist[j].atom = -1;
        }
    }
}

static void clear_atom_marks (atom_info *atoms)
{
    atom *atomlist = atoms->atomlist;
    int atomcount = atoms->atomcount;
    int i;

    for (i=0; i<atomcount; i++) atomlist[i].mark = 0;
}

static void clear_clique_marks (atom_info *atoms)
{
    vclique *cliquelist = atoms->cliquelist;
    int cliquecount = atoms->cliquecount;
    int i;

    for (i=0; i<cliquecount; i++) cliquelist[i].mark = 0;
}

static int max_node (BBtsp_lpcut_in *c)
{
    int i, j, k;
    int maxn;
    int cliquecount = c->cliquecount;
    BBtsp_lpclique *cliques = c->cliques;
    int dominocount = c->dominocount;
    BBtsp_lpdomino *dominos = c->dominos;

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

#define BBtsp_DOMINORHS(c) (3*(c)->dominocount + 1)

static int test_dp_cut (int ncount, BBtsp_lpcut_in *c, int *yes_no)
{
    int rval = 0;
    int i, handlecount, ndomino;
    int *Acount = (int *) NULL;
    int *Bcount = (int *) NULL;
    int **A = (int **) NULL;
    int **B = (int **) NULL;
    int *handle = (int *) NULL;
    int chkvalid, chkcomb;

    if (yes_no) *yes_no = 0;

    if (c->dominocount <= 0 || c->cliquecount != 1) {
        printf ("DP-cut has no dominos or no handle\n"); fflush (stdout);
        if (yes_no) *yes_no = 0;
        goto CLEANUP;
    }

    if (c->rhs != BBtsp_DOMINORHS(c)) {
        printf ("DP-cut has wrong RHS\n"); fflush (stdout);
        if (yes_no) *yes_no = 0;
        goto CLEANUP;
    }

    rval = extract_dp_cut (c, &ndomino, &Acount, &A, &Bcount, &B,
                              &handlecount, &handle);
    BBcheck_rval (rval, "extract_dp_cut failed");

    rval = check_raw_domino (ncount, handlecount, handle, ndomino, 
                             Acount, A, Bcount, B, &chkvalid, &chkcomb);
    BBcheck_rval (rval, "check_raw_domino");

    if (chkvalid == 0) {
        printf ("Raw DP-cut is not valid\n"); fflush (stdout);
        if (yes_no) *yes_no = 0;
    } else {
        if (yes_no) *yes_no = 1;
    }

CLEANUP:
    BB_IFFREE (Acount, int);
    BB_IFFREE (Bcount, int);
    if (A) {
        for (i = 0; i < ndomino; i++) {
            BB_IFFREE (A[i], int);
        }
        BB_FREE (A, int *);
    }
    if (B) {
        for (i = 0; i < ndomino; i++) {
            BB_IFFREE (B[i], int);
        }
       BB_FREE (B, int *);
    }
    BB_IFFREE (handle, int);

    return rval;
}

static int extract_dp_cut (BBtsp_lpcut_in *c, int *pndomino, int **pAcount,
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

    rval = BButil_clique_to_array (&c->cliques[0], &handle, &handlecount);
    BBcheck_rval (rval, "BButil_clique_to_array failed");

    ndomino = c->dominocount;
    Acount = BB_SAFE_MALLOC (ndomino, int);
    BBcheck_NULL (Acount, "out of memory for Acount");
    Bcount = BB_SAFE_MALLOC (ndomino, int);
    BBcheck_NULL (Acount, "out of memory for Bcount");
    A = BB_SAFE_MALLOC (ndomino, int *);
    BBcheck_NULL (Acount, "out of memory for A");
    B = BB_SAFE_MALLOC (ndomino, int *);
    BBcheck_NULL (Acount, "out of memory for B");

    for (i = 0; i < ndomino; i++) {
        rval = BButil_clique_to_array (&c->dominos[i].sets[0],
                     &A[i], &Acount[i]);
        BBcheck_rval (rval, "BButil_clique_to_array failed");
        rval = BButil_clique_to_array (&c->dominos[i].sets[1],
                     &B[i], &Bcount[i]);
        BBcheck_rval (rval, "BButil_clique_to_array failed");
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
        BB_IFFREE (Acount, int);
        BB_IFFREE (Bcount, int);
        BB_IFFREE (A, int *);  /* Note:  Misses the sets */
        BB_IFFREE (B, int *);  /* Note:  Misses the sets */
        BB_IFFREE (handle, int);
    }

    return rval;
}

static int check_raw_domino (int ncount, int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B, int *valid, int *comb)
{
    int rval = 0;
    int *marks = (int *) NULL;
    int i, j, marker = 0;
    int ain, aout, bin, bout;

    *valid = 0;
    *comb = 0;

    marks = BB_SAFE_MALLOC (ncount, int);
    BBcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    if (dcount % 2 == 0) {
        fprintf (stderr, "Even number of dominos in a DP cut\n");
        goto CLEANUP;
    }

    if (hcount <= 0 || hcount >= ncount) goto CLEANUP;
    marker++;
    for (i = 0; i < hcount; i++) {
        if (hand[i] < 0 || hand[i] >= ncount) goto CLEANUP;
        if (marks[hand[i]] == marker) goto CLEANUP;
        marks[hand[i]] = marker;
    }
    for (i = 0; i < dcount; i++) {
        if (Acount[i] + Bcount[i] >= ncount) goto CLEANUP;
        if (Acount[i] <= 0 || Bcount[i] <= 0) goto CLEANUP;
        marker++;
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) goto CLEANUP;
            marks[A[i][j]] = marker;
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) goto CLEANUP;
            marks[B[i][j]] = marker;
        }
    }
    *valid = 1;

     /* Test for comb */

     /* First check in dominos intersect handle correctly */

    marker++;
    for (i = 0; i < hcount; i++) marks[hand[i]] = marker;
    for (i = 0; i < dcount; i++) {
        ain = aout = 0;
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) ain++;
            else                          aout++;
        }
        if (ain && aout) goto CLEANUP;

        bin = bout = 0;
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) bin++;
            else                          bout++;
        }
        if (bin && bout) goto CLEANUP;
        if ((aout && bout) || (ain && bin)) goto CLEANUP;
    }

    /* Now check that teeth are disjoint */

    marker++;
    for (i = 0; i < dcount; i++) {
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) goto CLEANUP;
            marks[A[i][j]] = marker;
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) goto CLEANUP;
            marks[B[i][j]] = marker;
        }
    }

    *comb = 1;

CLEANUP:
    BB_IFFREE (marks, int);
    return rval;
}

