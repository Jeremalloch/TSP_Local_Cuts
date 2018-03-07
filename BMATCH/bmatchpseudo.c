#include <stdio.h>
#include <stdlib.h>
#include "bmatch.h"

#define DEBUG 0

static int getcomblist (CCbmat_graph *G, CCbmat_node **comblist);
static int addcomblist (CCbmat_graph *G, CCbmat_node *p, int k,
    CCbmat_node **comblist);
static void shrinknode (CCbmat_node *px, CCbmat_node *p, CCbmat_edge *eold,
    CCbmat_node *pold);
static void shrinkroot (CCbmat_node *px, CCbmat_node *p, CCbmat_edge *eold1,
    CCbmat_edge *eold2, CCbmat_node *pold1, CCbmat_node *pold2);
static void update_surfs_shrink (CCbmat_node *px);
static void update_surfs_blowup (CCbmat_node *px);
static void update_surf_work (CCbmat_node *p, CCbmat_node *s);
static void grab_edges (CCbmat_heaps *H, CCbmat_node *px, CCbmat_node *p);
static void return_edges (CCbmat_node *px);
static void inandout (CCbmat_edge *e, CCbmat_node *p, CCbmat_node **in,
    CCbmat_node **out);
static void makecycle (CCbmat_node *proot, CCbmat_edge *blossom_edge);
static void reversecycle (CCbmat_node *p);
static CCbmat_edge *findedgeout(CCbmat_node *p);
static int toothend (CCbmat_edge *e, CCbmat_node *p);
static void fixteeth (CCbmat_node *p);

int CCbmat_freecombs (CCbmat_graph *G, CCbmat_node **pseudo_freelist)
{
    int i, n, count, rval = 0;
    CCbmat_node **comblist = (CCbmat_node **) NULL;

    comblist = CC_SAFE_MALLOC (G->ncount, CCbmat_node *);
    CCcheck_NULL (comblist, "out of memory for comblist");

    n = getcomblist (G, comblist);
    for (i=0, count=0; i<n; i++) {
        if (comblist[i]->ymat > 0.001)
            count++;
        CCbmat_node_free (comblist[i], pseudo_freelist);
    }
    printf ("USED %d COMBS (%d with nonzero multiplier)\n",n,count);

CLEANUP:
    CC_IFFREE (comblist, CCbmat_node *);
    return rval;
}

static int getcomblist (CCbmat_graph *G, CCbmat_node **comblist)
{
    int i,k = 0;
    CCbmat_node *p;

    for (i=0, p = G->nodelist; i < G->ncount; i++, p++) {
        k = addcomblist (G, p, k, comblist);
    }
    return k;
}

static int addcomblist (CCbmat_graph *G, CCbmat_node *p, int k,
        CCbmat_node **comblist)
{
    /* add combs whose smallest real vertex is p */
    /* k indexes the next free spot in comblist */

    CCbmat_node *p1;

    for (p1 = p->nest.parent; p1 && smallest(G, p1) == p;
                              p1 = p1->nest.parent) {
        comblist[k++] = p1;
    }
    return k;
}

static void shrinknode (CCbmat_node *px, CCbmat_node *p, CCbmat_edge *eold,
        CCbmat_node *pold)
{
    CCbmat_node *p1, *p2;
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;

    /* hook children of p to px except pold, and mark all matching edges
       out of p as teeth except eold and parentedge, and nest p in px */

    for (p1 = p->tree.child; p1; p1 = p2) {
        p2 = p1->tree.sibling;
        if (p1 != pold) {
            CCbmat_node_remove (p1);
            CCbmat_addon (p1, p1->parentedge, px);
        }
    }
    if (!p->status.pseudo) {
        for (pee = p->edgelist; pee; pee = pee->next) {
            pe = CCbmat_THISEDGE(pee);
            if (pe->x && pe != eold && pe != p->parentedge) {
                if (pe->ends[0] == p)
                    pe->istooth[1] = CCbmat_TRUE;
                else
                    pe->istooth[0] = CCbmat_TRUE;
            }
        }
    }
    p->nest.parent = px;
    p->nest.sibling = px->nest.child;
    px->nest.child = p;
    p->label = CCbmat_UNLABELED;
}

static void shrinkroot (CCbmat_node *px, CCbmat_node *p, CCbmat_edge *eold1,
        CCbmat_edge *eold2, CCbmat_node *pold1, CCbmat_node *pold2)
{
    CCbmat_node *p1, *p2;
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;

    /* safety measure */
    if (p->tree.parent == NULL)
        p->parentedge = NULL;

    for (p1 = p->tree.child; p1; p1 = p2) {
        p2 = p1->tree.sibling;
        if (p1 != pold1 && p1 != pold2) {
            CCbmat_node_remove (p1);
            CCbmat_addon (p1, p1->parentedge, px);
        }
    }
    if (!p->status.pseudo) {
        for (pee = p->edgelist; pee; pee = pee->next) {
            pe = CCbmat_THISEDGE(pee);
            if (pe->x && pe != eold1 && pe != eold2 && pe != p->parentedge) {
                if (pe->ends[0] == p)
                    pe->istooth[1] = CCbmat_TRUE;
                else
                    pe->istooth[0] = CCbmat_TRUE;
            }
        }
    }
    pe = p->parentedge;
    if (p->label == CCbmat_MINUSMAT) {
        assert (!p->status.pseudo);
        assert (surfcheck(pe->ends[0], p) == (pe->surfends[0] == p));
        if (pe->surfends[0] == p) {
            assert (!pe->istooth[1]);
            pe->istooth[1] = CCbmat_TRUE;
        }
        else {
            assert (!pe->istooth[0]);
            pe->istooth[0] = CCbmat_TRUE;
        }
    }
    p->nest.parent = px;
    p->nest.sibling = px->nest.child;
    px->nest.child = p;
    p1 = p->tree.parent;
    if (p1) {
        CCbmat_node_remove (p);
        CCbmat_addon (px, pe, p1);
    }
    p->label = CCbmat_UNLABELED;
}

static void update_surfs_shrink (CCbmat_node *px)
{
    update_surf_work (px, px);
}

static void update_surfs_blowup (CCbmat_node *px)
{
    CCbmat_node *p;

    for (p = px->nest.child; p; p = p->nest.sibling) {
        update_surf_work (p, p);
    }
}

static void update_surf_work (CCbmat_node *p, CCbmat_node *s)
{
    CCbmat_node *q1, *q2;
/*  This code really does the same thing as:

    for (q = p->nest.child; q; q = q->nest.sibling) {
        update_surf_work (q, s);
    }
*/

    p->surf = s;
    q1 = p;
    for (;;) {
        while ((q2 = q1->nest.child)) {
            q2->surf = s;
            q1 = q2;
        }
        if (q1 == p) return;
        while (!(q2 = q1->nest.sibling)) {
            q1 = q1->nest.parent;
            if (q1 == p) return;
        }
        q2->surf = s;
        q1 = q2;
    }
}

static void grab_edges (CCbmat_heaps *H, CCbmat_node *px, CCbmat_node *p)
{
    CCbmat_edgeptr *pe, *penext;
    CCbmat_edge *e;

    assert (p->nest.parent == px);

    for (pe = p->edgelist, p->edgelist = (CCbmat_edgeptr *) NULL;
         pe; pe = penext) {
        penext = pe->next;
        e = CCbmat_THISEDGE(pe);
        if (e->surfends[0] == p) {
            if (e->surfends[1]->nest.parent == px) {
                /* shrunken */
                CCbmat_heap_edgeshrunk (H, e);
                pe->next = p->edgelist;
                p->edgelist = pe;
            } else {
                /* not shrunken */
                pe->next = px->edgelist;
                px->edgelist = pe;
                e->surfends[0] = px;
            }
        } else {
            assert (e->surfends[1] == p);
            if (e->surfends[0]->nest.parent == px) {
                /* shrunken */
                CCbmat_heap_edgeshrunk (H, e);
                pe->next = p->edgelist;
                p->edgelist = pe;
            } else {
                /* not shrunken */
                pe->next = px->edgelist;
                px->edgelist = pe;
                e->surfends[1] = px;
            }
        }
    }
}

static void return_edges (CCbmat_node *px)
{
    CCbmat_edgeptr *pe, *penext;
    CCbmat_edge *e;
    CCbmat_node *p;

    /* return_edges assumes update_surfs_blowup has already been done */

    for (pe = px->edgelist; pe; pe = penext) {
        penext = pe->next;
        e = CCbmat_THISEDGE(pe);
        if (e->surfends[0] == px) {
            assert (e->surfends[1] != px);
            p = e->ends[0]->surf;
            e->surfends[0] = p;
        } else {
            assert (e->surfends[1] == px);
            p = e->ends[1]->surf;
            e->surfends[1] = p;
        }
        pe->next = p->edgelist;
        p->edgelist = pe;
    }
    px->edgelist = (CCbmat_edgeptr *) NULL;
}

void CCbmat_really_shrink_blossom (CCbmat_heaps *H, CCbmat_edge *e,
        CCbmat_node *proot, CCbmat_node **pseudo_freelist)
{
    CCbmat_node *p1 = e->surfends[0];
    CCbmat_node *p2 = e->surfends[1];
    CCbmat_node *p3, *px;
    CCbmat_node *pold1, *pold2;
    CCbmat_edge *eold1, *eold2;

#if    DEBUG
    printf ("really_shrink_blossom (edge %ld, node %ld)\n",
                                    e-edgelist,proot-G->nodelist);
    fflush (stdout);
#endif    /* DEBUG */
    /* e->slack == 0, e forms blossom */

    assert (CCbmat_surftest (e));

    px = CCbmat_node_alloc (pseudo_freelist);
    px->status.pseudo = CCbmat_TRUE;
    px->label = CCbmat_PLUSMAT;
    px->edgelist = (CCbmat_edgeptr *) NULL;
    for (p3 = p1, eold1 = e, pold1 = NULL; p3 != proot; pold1 = p3,
                  eold1 = p3->parentedge, p3 = p3->tree.parent) {
        shrinknode (px, p3, eold1, pold1);
    }
    for (p3 = p2, eold2 = e, pold2 = NULL; p3 != proot; pold2 = p3,
                  eold2 = p3->parentedge, p3 = p3->tree.parent) {
        shrinknode (px, p3, eold2, pold2);
    }
    shrinkroot (px, proot, eold1, eold2, pold1, pold2);
    proot->parentedge = e;

    update_surfs_shrink (px);
    /* edge grabbing needs to occur after all node shrinking, so that
       shrunken edges can be recognized */
    for (p3 = p1; p3 != proot; p3 = p3->tree.parent) {
        grab_edges (H, px, p3);
        CCbmat_heap_shrink_node (H, p3);
    }
    for (p3 = p2; p3 != proot; p3 = p3->tree.parent) {
        grab_edges (H, px, p3);
        CCbmat_heap_shrink_node (H, p3);
    }
    grab_edges (H, px, proot);
    CCbmat_heap_shrink_node (H, proot);

    CCbmat_heap_label_update (H, px);
}

static void inandout (CCbmat_edge *e, CCbmat_node *p, CCbmat_node **in, 
        CCbmat_node **out)
{
    if ((*out = e->surfends[0]) == p) {
        *in = surfsecond (e->ends[0]);
        *out = e->surfends[1];
    }
    else
        *in = surfsecond (e->ends[1]);
}

static void makecycle (CCbmat_node *proot, CCbmat_edge *blossom_edge)
{
    CCbmat_node *p1, *plast, *pnext;

    p1 = proot->tree.child;
    /* split into two paths */
    CCbmat_node_remove (p1);
    /* reverse one path */
    plast = proot;
    while (p1) {
        pnext = p1->tree.child;
        plast->tree.parent = p1;
        plast->parentedge = p1->parentedge;
        p1->tree.child = plast;
        plast = p1;
        p1 = pnext;
    }

    /* hook the paths together */
    for (p1 = proot; p1->tree.child; p1 = p1->tree.child)
        ;
    p1->tree.child = plast;
    plast->tree.parent = p1;
    plast->parentedge = blossom_edge;
}

static void reversecycle (CCbmat_node *p)
{
    CCbmat_node *p1, *p2, *p3;
    CCbmat_edge *esave = p->parentedge;

    p1 = p;
    p2 = p1->tree.child;
    do {
        p1->parentedge = p2->parentedge;
        p1->tree.parent = p2;
        p3 = p2->tree.child;
        p2->tree.child = p1;
        p1 = p2;
        p2 = p3;
    } while (p1 != p);
    p->tree.child->parentedge = esave;
}

static CCbmat_edge *findedgeout(CCbmat_node *p)
{
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;

    for (pee = p->edgelist; pee; pee = pee->next) {
        pe = CCbmat_THISEDGE(pee);
        assert (surfcheck (pe->ends[0], p) == (pe->surfends[0] == p));
        if (pe->surfends[0] == p) {
            assert (!surfcheck (pe->ends[1], p));
            if ((pe->x && !pe->istooth[1]) ||
                (!pe->x && pe->istooth[1]))
                return pe;
        } else {
            assert (surfcheck (pe->ends[1], p));
            if ((pe->x && !pe->istooth[0]) ||
                (!pe->x && pe->istooth[0]))
                return pe;
        }
    }
    return NULL;
}

static int toothend (CCbmat_edge *e, CCbmat_node *p)
{
    assert (surfcheck (e->ends[0], p) == (e->surfends[0] == p));
    if (e->surfends[0] == p)
        return e->istooth[1];
    else {
        assert (surfcheck (e->ends[1], p));
        return e->istooth[0];
    }
}

void CCbmat_pseudoexpand (CCbmat_heaps *H, CCbmat_node *p,
        CCbmat_node **pseudo_freelist)
{
    /* Strip the pseudonode out of the nest structure, clear up
       appropriate teeth, but watch out for pre-existing teeth */
      
    CCbmat_node *proot, *p1, *p2;
    CCbmat_edge *parent_edge, *child_edge, *blossom_edge;
    CCbmat_node *child_inside, *child_outside, *parent_inside, *parent_outside;
    int parity;
    int child_end, parent_end;

    assert(CCbmat_checklabel (p));
    proot = p->nest.child;
    parent_edge = p->parentedge;
    child_edge = findedgeout (p);
    assert (child_edge);
    blossom_edge = proot->parentedge;

    /* fix matching in polygon */
    inandout (child_edge, p, &child_inside, &child_outside);
    if (p->tree.child == NULL)
        child_outside = NULL;
    assert (p->tree.child == child_outside);
    inandout (parent_edge, p, &parent_inside, &parent_outside);
    /* check out the polygon case */
    if (!child_inside->status.pseudo && child_inside == parent_inside) {
        child_end = child_edge->ends[0] == child_inside ? 0 : 1;
        assert (child_edge->ends[child_end] == child_inside);
        parent_end = parent_edge->ends[0] == child_inside ? 0 : 1;
        assert (parent_edge->ends[parent_end] == child_inside);
        if (child_edge->istooth[1-child_end] !=
            parent_edge->istooth[1-parent_end]) {
            child_edge->istooth[1-child_end] = 
                                 !child_edge->istooth[1-child_end];
            parent_edge->istooth[1-parent_end] =
                                 !parent_edge->istooth[1-parent_end];
            p->label = CCbmat_PLUSMAT;
            assert(CCbmat_checklabel (parent_outside));
            CCbmat_heap_label_update (H, p);
            return;
        }
    }

    update_surfs_blowup (p);
    return_edges (p);

    fixteeth (proot);

    /* fix matching */
    CCbmat_flipinterval (child_inside, proot);
    makecycle (proot, blossom_edge);
    if ((child_edge->x ^ toothend(child_edge, child_inside)) ==
        (child_inside->parentedge->x ^
        toothend(child_inside->parentedge, child_inside))) {
        CCbmat_flipcycle2 (proot);
    }

    /* fix augmentation tree */
    if (child_outside)
        CCbmat_node_remove (child_outside);
    CCbmat_node_remove (p);
    if ((parent_edge->x ^ toothend(parent_edge, parent_inside)) ==
        (parent_inside->tree.child->parentedge->x ^
        toothend (parent_inside->tree.child->parentedge, parent_inside)))
        reversecycle (proot);

    for (p1 = child_inside->tree.child; p1 != parent_inside; p1 = p2) {
        p2 = p1->tree.child;
        p1->tree.parent = p1->tree.sibling = p1->tree.child = NULL;
    }
    parent_inside->tree.parent = NULL;
    child_inside->tree.child = NULL;
    CCbmat_addon (parent_inside, parent_edge, parent_outside);
    if (child_outside)
        CCbmat_addon (child_outside, child_edge, child_inside);

    /* fix labeling */
    parity = parent_outside->label;
    for (p1 = parent_inside; p1 != child_outside; p1 = p1->tree.child) {
        p1->label = parity = ((parity == CCbmat_PLUSMAT) ^ 
                      (p1->parentedge->istooth[0] == CCbmat_TRUE) ^ 
                      (p1->parentedge->istooth[1] == CCbmat_TRUE)) 
                   ? CCbmat_MINUSMAT : CCbmat_PLUSMAT;
    }

    /* fix pseudonode nest. Hopefully, we can successfully postpone this
       until near the end, so we can use it to also update the heaps */
    for (p1 = proot; p1; p1 = p2) {
        p1->nest.parent = NULL;
        p2 = p1->nest.sibling;
        p1->nest.sibling = NULL;
        CCbmat_heap_label_update (H, p1);
    }

    assert(CCbmat_checklabel (parent_outside));

    CCbmat_node_free (p, pseudo_freelist);
}

static void fixteeth (CCbmat_node *p)
{
    CCbmat_node *p1;
    CCbmat_edge *pe;
    CCbmat_edgeptr *pee;

    for (p1 = p->tree.child; p1; p1 = p1->tree.sibling)
        fixteeth (p1);
    if (!p->status.pseudo) {
        for (pee = p->edgelist; pee; pee = pee->next) {
            pe = CCbmat_THISEDGE(pee);
            if (pe->ends[0] == p)
                pe->istooth[1] = CCbmat_FALSE;
            else
                pe->istooth[0] = CCbmat_FALSE;
        }
    }
}

void CCbmat_fixmatching (CCbmat_node *p)
{
    CCbmat_edge *child_edge, *paredge;
    CCbmat_node *child_inside, *proot;
    int inside_end;
    int parity;

    if (!p->status.pseudo)
        return;
    child_edge = findedgeout (p);
    assert (surfcheck(child_edge->ends[0], p) ==
                     (child_edge->surfends[0] == p));
    inside_end = (child_edge->surfends[0] == p) ? 0 : 1;
    child_inside = surfbelow (child_edge->ends[inside_end], p);
    proot = p->nest.child;
    CCbmat_flipinterval (child_inside, proot);

    /* convert to cycle (using child and parent pointers) */
    makecycle (proot, proot->parentedge);
    paredge = child_inside->parentedge;
    parity = child_edge->x ^ paredge->x;
    if (child_inside->status.pseudo) {
        if (child_edge->istooth[1-inside_end])
            parity = !parity;
        assert (surfcheck (paredge->ends[0], child_inside) ==
                          (paredge->surfends[0] == child_inside));
        if (paredge->surfends[0] == child_inside)
            parity ^= paredge->istooth[1];
        else
            parity ^= paredge->istooth[0]; 
    }
    if (!parity)
        CCbmat_flipcycle2 (proot);
    update_surfs_blowup (p);
    return_edges (p);
}

void CCbmat_zeroblowup (CCbmat_node *p, CCbmat_node **pseudo_freelist)
{
    /* Strip the pseudonode out of the nest structure, clear up
       appropriate teeth, but watch out for pre-existing teeth */
      
    CCbmat_node *proot, *p1, *p2;
    CCbmat_edge *child_edge, *blossom_edge;
    CCbmat_node *child_inside, *child_outside;

    assert (p->status.pseudo);
    assert (p->ymat == 0);

    proot = p->nest.child;
    child_edge = findedgeout (p);
    assert (child_edge);
    blossom_edge = proot->parentedge;

    /* fix matching in polygon */
    inandout (child_edge, p, &child_inside, &child_outside);

    update_surfs_blowup(p);
    return_edges (p);
    fixteeth (proot);

    CCbmat_flipinterval (child_inside, proot);

    /* convert to cycle (using child and parent pointers) */
    makecycle (proot, blossom_edge);
    if ((child_edge->x ^ toothend(child_edge, child_inside)) ==
        (child_inside->parentedge->x ^
        toothend(child_inside->parentedge, child_inside)))
        CCbmat_flipcycle2 (proot);

    /* fix pseudonode nest */
    for (p1 = proot; p1; p1 = p2) {
        p1->nest.parent = NULL;
        p2 = p1->nest.sibling;
        p1->nest.sibling = NULL;
    }

    /* cleanup augmentation tree */
    p1 = proot;
    do {
        p2 = p1->tree.parent;
        p1->tree.parent = p1->tree.child = p1->tree.sibling
                                         = (CCbmat_node *) NULL;
        p1->parentedge = (CCbmat_edge *) NULL;
        p1->label = CCbmat_UNLABELED;
        p1 = p2;
    } while (p1 != proot);
    CCbmat_node_free (p, pseudo_freelist);
}

