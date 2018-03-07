#include <stdio.h>
#include "bmatch.h"

#define DEBUG 0

static int analyzedoubletooth (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
    CCbmat_node *p2, CCbmat_node **pseudo_freelist);
static int analyzetooth (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
    CCbmat_node *p2, CCbmat_node **pseudo_freelist);
static int analyzenoteeth (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
    CCbmat_node *p2, CCbmat_node **pseudo_freelist);
static int addnode (CCbmat_heaps *H, CCbmat_node *p1, CCbmat_edge *e,
    CCbmat_node *p2, int l, CCbmat_node **pseudo_freelist);
static CCbmat_node *get_blossom_root (CCbmat_node *p1, CCbmat_node *p2);
static CCbmat_node *finddeficient (CCbmat_node *p1, CCbmat_node *proot);
static int shrink_blossom (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
    CCbmat_node *p2, CCbmat_node **pseudo_freelist);
static void augmentpath2 (CCbmat_node *p);
static void flippath (CCbmat_node *p);
static void augmentblossom (CCbmat_node *p1, CCbmat_node *p2, CCbmat_edge *e,
   CCbmat_node *p3);

/* make tree rooted at new a descendent of old, connected by e */
void CCbmat_addon (CCbmat_node *new, CCbmat_edge *e, CCbmat_node *old)    
{
    new->tree.sibling = old->tree.child;
    old->tree.child = new;
    new->tree.parent = old;
    new->parentedge = e;
}

void CCbmat_node_remove (CCbmat_node *p)
{
    CCbmat_node *p2, *p3;

    if ((p2 = p->tree.parent)) {
        if ((p3 = p2->tree.child) == p)
            p2->tree.child = p->tree.sibling;
        else {
            while (p3->tree.sibling != p)
                p3 = p3->tree.sibling;
            p3->tree.sibling = p->tree.sibling;
        }
    }
    p->tree.sibling = p->tree.parent = NULL;
}

int CCbmat_analyze (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node **pseudo_freelist)
{
    CCbmat_node *p1, *p2;

    assert (e->slack == 0);
    p1 = e->surfends[0];
    p2 = e->surfends[1];
    assert (p1 != p2);
    if (e->istooth[0])
        if (e->istooth[1])
            return analyzedoubletooth (H, e, p1, p2, pseudo_freelist);
        else
            return analyzetooth (H, e, p2, p1, pseudo_freelist);
    else if (e->istooth[1])
            return analyzetooth (H, e, p1, p2, pseudo_freelist);
    else
        return analyzenoteeth (H, e, p1, p2, pseudo_freelist);
}

/* e is a tooth edge in both directions */
static int analyzedoubletooth (CCbmat_heaps *H, CCbmat_edge *e,
        CCbmat_node *p1, CCbmat_node *p2, CCbmat_node **pseudo_freelist)
{
    if (e->x) {
        if (p1->label == CCbmat_PLUSMAT)
            return addnode (H, p2, e, p1, CCbmat_MINUSMAT, pseudo_freelist);
        else if (p2->label == CCbmat_PLUSMAT)
            return addnode (H, p1, e, p2, CCbmat_MINUSMAT, pseudo_freelist);
        else {
            assert (0);
            return 0;
        }
    } else {
        if (p1->label == CCbmat_MINUSMAT)
            return addnode (H, p2, e, p1, CCbmat_PLUSMAT, pseudo_freelist);
        else if (p2->label == CCbmat_MINUSMAT)
            return addnode (H, p1, e, p2, CCbmat_PLUSMAT, pseudo_freelist);
        else {
            assert (0);
            return 0;
        }
    }
}

/* e is a tooth edge, p1 the handle, p2 the tooth */
static int analyzetooth (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
        CCbmat_node *p2, CCbmat_node **pseudo_freelist)
{
    if (e->x) {
        if (p1->label == CCbmat_PLUSMAT)
            return addnode (H, p2, e, p1, CCbmat_PLUSMAT, pseudo_freelist);
        else if (p2->label == CCbmat_MINUSMAT)
            return addnode (H, p1, e, p2, CCbmat_MINUSMAT, pseudo_freelist);
        else {
            assert (0);
            return 0;
        }
    } else {
        if (p1->label == CCbmat_MINUSMAT)
            return addnode (H, p2, e, p1, CCbmat_MINUSMAT, pseudo_freelist);
        else if (p2->label == CCbmat_PLUSMAT)
            return addnode (H, p1, e, p2, CCbmat_PLUSMAT, pseudo_freelist);
        else {
            assert (0);
            return 0;
        }
    }
}

/* e is not a tooth edge */
static int analyzenoteeth (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
        CCbmat_node *p2, CCbmat_node **pseudo_freelist)
{
    if (e->x) {
        if (p1->label == CCbmat_MINUSMAT)
            return addnode (H, p2, e, p1, CCbmat_PLUSMAT, pseudo_freelist);
        else if (p2->label == CCbmat_MINUSMAT)
            return addnode (H, p1, e, p2, CCbmat_PLUSMAT, pseudo_freelist);
        else {
            assert (0);
            return 0;
        }
    } else {
        if (p1->label == CCbmat_PLUSMAT)
            return addnode (H, p2, e, p1, CCbmat_MINUSMAT, pseudo_freelist);
        else if (p2->label == CCbmat_PLUSMAT)
            return addnode (H, p1, e, p2, CCbmat_MINUSMAT, pseudo_freelist);
        else {
            assert (0);
            return 0;
        }
    }
}

/* add p1 to the tree by using e to connect it to p2, and label it l */
static int addnode (CCbmat_heaps *H, CCbmat_node *p1, CCbmat_edge *e,
        CCbmat_node *p2, int l, CCbmat_node **pseudo_freelist)
{
#if    DEBUG
    printf ("addnode (node %d, edge %d, node %d, label %d\n",p1-nodelist,
        e-edgelist,p2-nodelist,l);
    fflush (stdout);
#endif    /* DEBUG */
    assert (l != CCbmat_UNLABELED);
    if (p1->label == CCbmat_UNLABELED) {
        CCbmat_addon (p1, e, p2);
        p1->label = l;
        CCbmat_heap_label_update (H, p1);
        if (l == CCbmat_MINUSMAT && p1->deficiency) {
#if    DEBUG
            printf ("    (deficient, augmented)\n");
            fflush (stdout);
#endif    /* DEBUG */
            augmentpath2 (p1);
            return 1;
        } else {
#if    DEBUG
            printf ("   (added to tree, expanding)\n");
            fflush (stdout);
#endif    /* DEBUG */
            return 0;
        }
    } else {
        assert (p1->label != l);
#if    DEBUG
        printf ("   (blossom forming, evaluating)\n");
        fflush (stdout);
#endif    /* DEBUG */
        return shrink_blossom (H, e, p1, p2, pseudo_freelist);
    }
}

static CCbmat_node *get_blossom_root (CCbmat_node *p1, CCbmat_node *p2)
{
    CCbmat_node *p;

    for (p = p1; p; p = p->tree.parent)
        p->status.inpath = CCbmat_TRUE;
    while (!p2->status.inpath)
        p2 = p2->tree.parent;
    for (p = p1; p; p = p->tree.parent)
        p->status.inpath = CCbmat_FALSE;
    return p2;
}

static CCbmat_node *finddeficient (CCbmat_node *p1, CCbmat_node *proot)
{
    while (p1 != proot) {
        if (p1->deficiency)
            return p1;
        p1 = p1->tree.parent;
    }
    return proot->deficiency ? proot : NULL;
}

static int shrink_blossom (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node *p1,
        CCbmat_node *p2, CCbmat_node **pseudo_freelist)
{
    CCbmat_node *proot, *p3;

    /* e forms a blossom from p1 to p2.  Check to see if we can augment,
       and if not, save the info for later use. */
    proot = get_blossom_root (p1, p2);
    if ((p3 = finddeficient (p1, proot))) {
        augmentblossom (p3, p1, e, p2);
        return 1;
    } else if ((p3 = finddeficient (p2, proot))) {
        augmentblossom (p3, p2, e, p1);
        return 1;
    } else {
        CCbmat_really_shrink_blossom (H, e, proot, pseudo_freelist);
        return 0;
    }
}

static void augmentpath2 (CCbmat_node *p)
{
    /* augment along the path from p to the root */
    p->deficiency--;
    flippath (p);
}

static void flippath (CCbmat_node *p)
{
    while (p->tree.parent) {
        p->parentedge->x = !p->parentedge->x;
        p = p->tree.parent;
    }
}

void CCbmat_flipinterval (CCbmat_node *pstart, CCbmat_node *pend)
{
    while (pstart != pend) {
        pstart->parentedge->x = !pstart->parentedge->x;
        pstart = pstart->tree.parent;
    }
}

void CCbmat_flipcycle2 (CCbmat_node *p)
{
    CCbmat_node *p1;

    p1 = p;
    do {
        p1->parentedge->x = !p1->parentedge->x;
        p1 = p1->tree.parent;
    } while (p1 != p);
}

static void augmentblossom (CCbmat_node *p1, CCbmat_node *p2, CCbmat_edge *e,
        CCbmat_node *p3)
{
    /* augment from p1 up to p2, across e to p3, and down to root */
    p1->deficiency--;
    CCbmat_flipinterval (p2, p1);
    e->x = !e->x;
    flippath (p3);
}

