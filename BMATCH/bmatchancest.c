#include <stdio.h>
#include "macrorus.h"
#include "bmatch.h"

/* exported: CCbmat_ancest_init, CCbmat_ancest_checkout, CCbmat ancest_flush */
/* external: CCbmat_repair_world, findedge */

#define ANCEST_CHUNK ((1048576-16)/sizeof (CCbmat_anc))

static CCbmat_anc *find_anc (CCbmat_anc *x);
static CCbmat_anc *union_anc (CCbmat_anc *x, CCbmat_anc *y);
static CCbmat_anc *ancest_work (CCbmat_node *p, CCbmat_anc **repair_list,
    CCbmat_anc **free_list);
static int do_repairs (CCbmat_graph *G, CCbmat_heaps *H, CCbmat_anc *repair_list,
    CCbmat_anc **anc_freelist, CCbmat_edge **edge_freelist,
    CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply);
static int build_free_anc (CCbmat_anc **anc_freelist,
    CCbmat_anc_ptr **anc_supply);
static void number_tree (CCbmat_node *n, int *dfsnum);
static void get_surface (CCbmat_graph *G, CCbmat_node **surfacelist);
static void checkout_anc (CCbmat_anc *a, CCbmat_node *plca,
    CCbmat_anc **repair_list, CCbmat_anc **anc_freelist);

static CCbmat_anc *find_anc (CCbmat_anc *x)
{
    CCbmat_anc *p = x->parent;

    assert (x);
    assert (p);

    return (x == p) ? x : (x->parent = find_anc (p));
}

static CCbmat_anc *union_anc (CCbmat_anc *x, CCbmat_anc *y)
{
    CCbmat_anc *t;

    assert (x);
    assert (y);
    assert (x->parent == x);
    assert (y->parent == y);

    if (x->rank > y->rank) {
	CC_SWAP (x, y, t);
    } else if (x->rank == y->rank) {
	y->rank++;
    }
    x->parent = y;
    return y;
}

int CCbmat_ancest_init (CCbmat_graph *G, CCbmat_node **surfacelist,
        CCbmat_anc **anc_freelist, CCbmat_anc_ptr **anc_supply)
{
    int dfsnum = 0, rval = 0;
    CCbmat_node *n;

    *surfacelist = (CCbmat_node *) NULL;
    get_surface (G, surfacelist);
    CCbmat_get_sums (G);

    for (n = *surfacelist; n; n = n->next) {
	number_tree (n, &dfsnum);
    }
    if (!(*anc_freelist)) {
	rval = build_free_anc (anc_freelist, anc_supply);
        CCcheck_rval (rval, "build_free_anc failed");
    }

CLEANUP:
    return rval;
}

int CCbmat_ancest_checkout (CCbmat_graph *G, CCbmat_heaps *H, CCbmat_node *n1,
        CCbmat_node *n2, int w, CCbmat_node **surfacelist,
        CCbmat_anc **anc_freelist, CCbmat_edge **edge_freelist,
        CCbmat_node **pseudo_freelist, CCbmat_anc_ptr **anc_supply,
        CCbmat_edgeptr **edge_supply)
{
    CCbmat_anc *p;
    CCbmat_anc *p2;
    CCbmat_node *t;

    assert (*anc_freelist);
    p = *anc_freelist;
    *anc_freelist = p->next;

    if (n1->ancest.number > n2->ancest.number) {
	CC_SWAP(n1, n2, t);
    }
    assert (n1->ancest.number < n2->ancest.number);

    p->left = n1;
    p->right = n2;
    p->weight = w;
    p->next = n2->ancest.evallist;
    n2->ancest.evallist = p;
    if ((p2 = n1->ancest.addset)) {
	p2->rank = 1;
	p->rank = 0;
	p->parent = p2;
    } else {
	p->ancest = n1;
	p->rank = 0;
	n1->ancest.addset = p;
	p->parent = p;
    }

    if (!*anc_freelist) {
	return CCbmat_ancest_flush (G, H, surfacelist, anc_freelist,
                      edge_freelist, pseudo_freelist, anc_supply, edge_supply);
    } else {
	return 0;
    }
}

int CCbmat_ancest_flush (CCbmat_graph *G, CCbmat_heaps *H,
        CCbmat_node **surfacelist, CCbmat_anc **anc_freelist,
        CCbmat_edge **edge_freelist, CCbmat_node **pseudo_freelist,
        CCbmat_anc_ptr **anc_supply, CCbmat_edgeptr **edge_supply)
{
    CCbmat_node *p;
    CCbmat_anc *proot = (CCbmat_anc *) NULL;
    CCbmat_anc *s;
    int added;
    CCbmat_anc *repair_list = (CCbmat_anc *) NULL;

    putchar ('A'); fflush (stdout);
    repair_list = (CCbmat_anc *) NULL;

    for (p = *surfacelist; p; p = p->next) {
	s = ancest_work (p, &repair_list, anc_freelist);
	if (s) {
	    if (proot) {
		proot = union_anc (s, proot);
	    } else {
		proot = s;
	    }
	    proot->ancest = (CCbmat_node *) NULL;
	}
    }

    added = do_repairs (G, H, repair_list, anc_freelist, edge_freelist,
                        pseudo_freelist, edge_supply);
    CCbmat_ancest_init (G, surfacelist, anc_freelist, anc_supply);
    putchar ('Z'); fflush (stdout);
    return added;
}

static CCbmat_anc *ancest_work (CCbmat_node *p, CCbmat_anc **repair_list,
        CCbmat_anc **anc_freelist)
{
    CCbmat_node *q;
    CCbmat_anc *a, *anext;
    CCbmat_node *n;

    for (q = p->nest.child; q; q = q->nest.sibling) {
	a = ancest_work (q, repair_list, anc_freelist);
	if (a) {
	    if (p->ancest.addset) {
		p->ancest.addset = union_anc (a, p->ancest.addset);
	    } else {
		p->ancest.addset = a;
	    }
	    p->ancest.addset->ancest = p;
	}
    }

    for (a = p->ancest.evallist; a; a = anext) {
	n = find_anc(a)->ancest;
	anext = a->next;
	checkout_anc (a, n, repair_list, anc_freelist);
    }

    return p->ancest.addset;
}

static int do_repairs (CCbmat_graph *G, CCbmat_heaps *H, CCbmat_anc *repair_list,
        CCbmat_anc **anc_freelist, CCbmat_edge **edge_freelist,
        CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply)
{
    CCbmat_anc *p, *pnext;
    int added = 0;

    for (p = repair_list; p; p = pnext) {
	added += CCbmat_repair_world (G, H, p->left, p->right, p->weight,
                               edge_freelist, pseudo_freelist, edge_supply);
	pnext = p->next;
	p->next = *anc_freelist;
	*anc_freelist = p;
    }
    return added;
}

static int build_free_anc (CCbmat_anc **anc_freelist, CCbmat_anc_ptr **anc_supply)
{
    int rval = 0;
    CCbmat_anc *p, *pend;
    CCbmat_anc_ptr *a = (CCbmat_anc_ptr *) NULL;

    *anc_freelist = CC_SAFE_MALLOC (ANCEST_CHUNK, CCbmat_anc);
    CCcheck_NULL (*anc_freelist, "out of memory for anc_freelist");

    for (p = *anc_freelist, pend = p + ANCEST_CHUNK-1; p < pend; p++) {
	p->next = p+1;
    }
    pend->next = (CCbmat_anc *) NULL;

    a = CC_SAFE_MALLOC (1, CCbmat_anc_ptr);
    CCcheck_NULL (a, "out or memory for a");
    a->this = *anc_freelist;
    a->next = *anc_supply;
    *anc_supply = a;

CLEANUP:
    return rval;
}

static void number_tree (CCbmat_node *n, int *dfsnum)
{
    CCbmat_node *p;

    n->ancest.number = (*dfsnum)++;
    n->ancest.addset = (CCbmat_anc *) NULL;
    n->ancest.evallist = (CCbmat_anc *) NULL;
    for (p = n->nest.child; p; p = p->nest.sibling) {
	number_tree (p, dfsnum);
    }
}

static void get_surface (CCbmat_graph *G, CCbmat_node **surfacelist)
{
    CCbmat_node *n, *nend = G->nodelist + G->ncount;

    for (n = G->nodelist; n < nend; n++) {
	n->surf->status.inpath = CCbmat_FALSE;
    }
    *surfacelist = (CCbmat_node *) NULL;
    for (n = G->nodelist; n < nend; n++) {
	if (!n->surf->status.inpath) {
	    n->surf->next = *surfacelist;
	    *surfacelist = n->surf;
	    n->surf->status.inpath = CCbmat_TRUE;
	}
    }
    for (n = *surfacelist; n; n = n->next) {
	assert (n->status.inpath);
	n->status.inpath = CCbmat_FALSE;
    }
}

static void checkout_anc (CCbmat_anc *a, CCbmat_node *plca,
        CCbmat_anc **repair_list, CCbmat_anc **anc_freelist)
{
    int wbar;
    int lcaval;
    CCbmat_node *n1 = a->left;
    CCbmat_node *n2 = a->right;

    assert (plca == lca (n1, n2));
    lcaval = plca ? plca->sum_to_root : 0;
    wbar = a->weight - n1->sum_to_root - n2->sum_to_root + 2 * lcaval;
    if (wbar < 0 && !findedge (n1, n2)) {
	a->next = *repair_list;
	*repair_list = a;
    } else {
	a->next = *anc_freelist;
	*anc_freelist = a;
    }
}

