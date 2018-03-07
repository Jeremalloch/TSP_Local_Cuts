#include <stdio.h>
#include "bmatch.h"

int CCbmat_checklabel (CCbmat_node *n)
{
    CCbmat_node *p, *p1, *p2;
    CCbmat_edge *e;

    assert (n->label != UNLABELED);
    if ((p = n->tree.parent)) {
        assert (p->label != UNLABELED);
        e = n->parentedge;
        assert (e);
        assert (CCbmat_surftest (e));
        p1 = e->surfends[0];
        p2 = e->surfends[1];
        assert ((p1 == n && p2 == p) || (p2 == n && p1 == p));
        if (e->istooth[0] == e->istooth[1]) {
            assert (n->label != p->label);
        } else {
            assert (n->label == p->label);
        }
    } else {
        assert (n->label == PLUSMAT);
    }
    for (p = n->tree.child; p; p=p->tree.sibling)
        CCbmat_checklabel (p);
    return 1;
}

#if 0
static void spit_out_combs (graph *G)
{
    int i, k;
    CCbmat_edge *e;
    CCbmat_node *n;

    printf ("Combs used:\n");
    for (i = 0, n = nodelist; i < nnodes; i++, n++)
        dumpacomb (n);
    printf ("Teeth used:\n");
    for (i=0, k=0, e=edgelist; i<nedges; i++,e++) {
        if (e->istooth[0] || e->istooth[1])
            printf ("%c%ld%s%ld",k++%8?' ':'\n',e->ends[0]-nodelist,
                e->istooth[0] ? (e->istooth[1] ? "<->" : "<--")
                          : "-->",e->ends[1]-nodelist);
    }
    for (e = G->newedges; e; e = e->next) {
        if (e->istooth[0] || e->istooth[1])
            printf ("%c%ld%s%ld",k++%8?' ':'\n',e->ends[0]-nodelist,
                e->istooth[0] ? (e->istooth[1] ? "<->" : "<--")
                          : "-->",e->ends[1]-nodelist);
    }
    printf ("\n");
    fflush (stdout);
}
#endif

#if 0
static long cutweight (CCbmat_node *n1, CCbmat_node *n2)
{
    long cw = 0;
    CCbmat_node *p;
    CCbmat_edge *e;

    if ((e = findedge (n1, n2))) {
        if (e->slack <= 0) {
            return e->weight;
        } else {
            return e->weight - e->slack;
        }
    }

    for (p = n1; p; p = p->nest.parent)
        p->status.inpath = TRUE;
    for (p = n2; p && !p->status.inpath; p = p->nest.parent)
        cw += p->ymat;
    for (; p; p = p->nest.parent)
        p->status.inpath = FALSE;
    for (p = n1; p && p->status.inpath; p = p->nest.parent) {
        cw += p->ymat;
        p->status.inpath = FALSE;
    }
    return cw;
}
#endif

#if 0
static void dumpacomb (CCbmat_node *p)
{
    /* dump combs whose smallest real vertex = p */
    CCbmat_node *p1;

    for (p1 = p->nest.parent; p1 && smallest(p1) == p; p1 = p1->nest.parent) {
        printcomb (p1);
        printf ("\n");
    }
    fflush (stdout);
}
#endif

#if 0
static void printcomb (CCbmat_node *p)
{
    CCbmat_node *p1;

    if (p->status.pseudo) {
        for (p1 = p->nest.child; p1; p1 = p1->nest.sibling)
            printcomb (p1);
    } else {
        printf ("%ld ",p-nodelist);
    }
    fflush (stdout);
}
#endif

void CCbmat_surftest (CCbmat_edge *e)
{
    assert (CCbmat_slowsurf (e->ends[0]) == e->ends[0]->surf);
    assert (e->ends[0]->surf == e->surfends[0]);
    assert (CCbmat_slowsurf (e->ends[1]) == e->ends[1]->surf);
    assert (e->ends[1]->surf == e->surfends[1]);
}

/* find the outermost pseudonode containing n */
CCbmat_node *CCbmat_slowsurf (CCbmat_node *n)
{
    CCbmat_node *n1;

    while ((n1 = n->nest.parent)) {
        n = n1;
    }
    return n;
}

