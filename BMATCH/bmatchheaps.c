#include <stdio.h>
#include "bmatch.h"

#define DEBUG 0

#define MINUS_NODE_HEAP 1
#define MINUS_PSEUDO_HEAP 2
#define PLUS_NODE_HEAP 3
#define MINUS2_MAT_HEAP 1
#define MINUS2_NOTMAT_HEAP 2
#define MINUS1_MAT_HEAP 3
#define MINUS1_NOTMAT_HEAP 4
#define PLUS1_MAT_HEAP 5
#define PLUS1_NOTMAT_HEAP 6
#define PLUS2_MAT_HEAP 7
#define PLUS2_NOTMAT_HEAP 8

static int heap_checkout (CCbmat_heaps *H, CCbmat_edge *e);
static void get_node_from_heap (CCbmat_heaps *H, CCbmat_node *n);
static int add_node_to_heap (CCbmat_heaps *H, CCbmat_node *n);
static void get_edge_from_heap (CCbmat_heaps *H, CCbmat_edge *n);
static int add_edge_to_heap (CCbmat_heaps *H, CCbmat_edge *n);
static void node_heap_init_struct (CCbmat_node_heap *h);
static int node_heap_init (CCbmat_node_heap *h);
static void node_heap_free (CCbmat_node_heap *h);
static int node_heap_expand (CCbmat_node_heap *h);
static int node_heap_insert (CCbmat_node *i, CCbmat_node_heap *h);
static void node_heap_delete (CCbmat_node *i, CCbmat_node_heap *h);
static void node_heap_siftup (CCbmat_node *i, int x, CCbmat_node_heap *h);
static void node_heap_siftdown (CCbmat_node *i, int x, CCbmat_node_heap *h);
static int node_heap_minchild (int x, CCbmat_node_heap *h);
static void node_heap_purge (CCbmat_node_heap *h, int delta);
static void edge_heap_init_struct (CCbmat_edge_heap *h);
static int edge_heap_init (CCbmat_edge_heap *h);
static void edge_heap_free (CCbmat_edge_heap *h);
static int edge_heap_expand (CCbmat_edge_heap *h);
static int edge_heap_insert (CCbmat_edge *i, CCbmat_edge_heap *h);
static void edge_heap_delete (CCbmat_edge *i, CCbmat_edge_heap *h);
static void edge_heap_siftup (CCbmat_edge *i, int x, CCbmat_edge_heap *h);
static void edge_heap_siftdown (CCbmat_edge *i, int x, CCbmat_edge_heap *h);
static int edge_heap_minchild (int x, CCbmat_edge_heap *h);
static void edge_heap_purge (CCbmat_edge_heap *h, int delta);
static int node_heap_verify_order (CCbmat_node_heap *h, int type);
static int edge_heap_verify_order (CCbmat_edge_heap *h, int type);

#if DEBUG
static int edge_correct_type (CCbmat_edge *e);
static int node_heap_contains (CCbmat_node *n, CCbmat_node_heap *h);
static int edge_heap_contains (CCbmat_edge *n, CCbmat_edge_heap *h);
#endif

int CCbmat_heap_verify (CCbmat_graph *G, CCbmat_heaps *H);

void CCbmat_init_heaps_struct (CCbmat_heaps *H)
{
    if (H) {
        node_heap_init_struct (&H->minus_node_heap);
        node_heap_init_struct (&H->minus_pseudo_heap);
        node_heap_init_struct (&H->plus_node_heap);

        edge_heap_init_struct (&H->minus2_mat_heap);
        edge_heap_init_struct (&H->minus1_mat_heap);
        edge_heap_init_struct (&H->plus1_mat_heap);
        edge_heap_init_struct (&H->plus2_mat_heap);
        edge_heap_init_struct (&H->minus2_notmat_heap);
        edge_heap_init_struct (&H->minus1_notmat_heap);
        edge_heap_init_struct (&H->plus1_notmat_heap);
        edge_heap_init_struct (&H->plus2_notmat_heap);
        H->delta = 0;
    }
}

int CCbmat_init_heaps (CCbmat_heaps *H)
{
    int rval = 0;

    rval = node_heap_init (&H->minus_node_heap);
    CCcheck_rval (rval, "node_heap_init failed");
    rval = node_heap_init (&H->minus_pseudo_heap);
    CCcheck_rval (rval, "node_heap_init failed");
    rval = node_heap_init (&H->plus_node_heap);
    CCcheck_rval (rval, "node_heap_init failed");

    rval = edge_heap_init (&H->minus2_mat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->minus1_mat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->plus1_mat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->plus2_mat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->minus2_notmat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->minus1_notmat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->plus1_notmat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");
    rval = edge_heap_init (&H->plus2_notmat_heap);
    CCcheck_rval (rval, "edge_heap_init failed");

    H->delta = 0;

CLEANUP:
    return rval;
}

void CCbmat_free_heaps (CCbmat_heaps *H)
{
    node_heap_free (&H->minus_node_heap);
    node_heap_free (&H->minus_pseudo_heap);
    node_heap_free (&H->plus_node_heap);
    edge_heap_free (&H->minus2_mat_heap);
    edge_heap_free (&H->minus1_mat_heap);
    edge_heap_free (&H->plus1_mat_heap);
    edge_heap_free (&H->plus2_mat_heap);
    edge_heap_free (&H->minus2_notmat_heap);
    edge_heap_free (&H->minus1_notmat_heap);
    edge_heap_free (&H->plus1_notmat_heap);
    edge_heap_free (&H->plus2_notmat_heap);
}


int CCbmat_heap_dispatch (CCbmat_heaps *H, int ybound, 
        CCbmat_node **pseudo_freelist)
{
    CCbmat_node *n, *n1;
    CCbmat_edge *e, *e1;
    int eval, minval = ybound - H->delta;

    assert (minval > 0);

    e = CCbmat_edge_heap_findmin (&H->plus1_mat_heap);
    if (e) {
        eval = - (e->slack + H->delta);
        if (eval < minval) minval = eval;
        if (eval == 0) {
            e1 = CCbmat_edge_heap_deletemin (&H->plus1_mat_heap);
            assert (e == e1);
            e1->heap.type = CCbmat_NOTHEAPED;
            e->slack += H->delta;
            assert (e->slack == 0);
            return CCbmat_analyze (H, e, pseudo_freelist);
        }
    }
    e = CCbmat_edge_heap_findmin (&H->minus1_notmat_heap);
    if (e) {
        eval = e->slack - H->delta;
        if (eval < minval) minval = eval;
        if (eval == 0) {
            e1 = CCbmat_edge_heap_deletemin (&H->minus1_notmat_heap);
            assert (e == e1);
            e1->heap.type = CCbmat_NOTHEAPED;
            e->slack -= H->delta;
            assert (e->slack == 0);
            return CCbmat_analyze (H, e, pseudo_freelist);
        }
    }
    e = CCbmat_edge_heap_findmin (&H->plus2_mat_heap);
    if (e) {
        eval = - (e->slack + 2*H->delta);
        if (eval/2 < minval) minval = eval/2;
        if (eval == 0) {
            e1 = CCbmat_edge_heap_deletemin (&H->plus2_mat_heap);
            assert (e == e1);
            e1->heap.type = CCbmat_NOTHEAPED;
            e->slack += 2*H->delta;
            assert (e->slack == 0);
            return CCbmat_analyze (H, e, pseudo_freelist);
        }
    }
    e = CCbmat_edge_heap_findmin (&H->minus2_notmat_heap);
    if (e) {
        eval = e->slack - 2*H->delta;
        if (eval/2 < minval) minval = eval/2;
        if (eval == 0) {
            e1 = CCbmat_edge_heap_deletemin (&H->minus2_notmat_heap);
            assert (e == e1);
            e1->heap.type = CCbmat_NOTHEAPED;
            e->slack -= 2*H->delta;
            assert (e->slack == 0);
            return CCbmat_analyze (H, e, pseudo_freelist);
        }
    }

    n = CCbmat_node_heap_findmin (&H->minus_pseudo_heap);
    if (n) {
        eval = n->ymat - H->delta;
        if (eval < minval) minval = eval;
        if (eval == 0) {
            n1 = CCbmat_node_heap_deletemin (&H->minus_pseudo_heap);
            assert (n == n1);
            n->heap.type = CCbmat_NOTHEAPED;
            n->ymat -= H->delta;
            assert (n->ymat == 0);
            CCbmat_pseudoexpand (H, n, pseudo_freelist);
            return 0;
        }
    }

    assert (minval > 0);
    H->delta += minval;

    if (H->delta == ybound) return 1;
    else return 0;
}

void CCbmat_clear_heaps (CCbmat_heaps *H)
{
    node_heap_purge (&H->minus_node_heap, -H->delta);
    node_heap_purge (&H->minus_pseudo_heap, -H->delta);
    node_heap_purge (&H->plus_node_heap, H->delta);

    edge_heap_purge (&H->minus2_mat_heap, -2*H->delta);
    edge_heap_purge (&H->minus1_mat_heap, -H->delta);
    edge_heap_purge (&H->plus1_mat_heap, H->delta);
    edge_heap_purge (&H->plus2_mat_heap, 2*H->delta);
    edge_heap_purge (&H->minus2_notmat_heap, -2*H->delta);
    edge_heap_purge (&H->minus1_notmat_heap, -H->delta);
    edge_heap_purge (&H->plus1_notmat_heap, H->delta);
    edge_heap_purge (&H->plus2_notmat_heap, 2*H->delta);

    H->delta = 0;
}

int CCbmat_heap_label_update (CCbmat_heaps *H, CCbmat_node *n)
{
    int rval = 0;
    CCbmat_edgeptr *pe;

    get_node_from_heap (H, n);

    switch (n->label) {
    case CCbmat_PLUSMAT:
        n->heap.type = PLUS_NODE_HEAP; 
        break;
    case CCbmat_MINUSMAT: 
        n->heap.type = n->status.pseudo ? MINUS_PSEUDO_HEAP : MINUS_NODE_HEAP;
        break;
    default:
        n->heap.type = CCbmat_NOTHEAPED;
        break;
    }

    rval = add_node_to_heap (H, n);
    CCcheck_rval (rval, "add_node_to_heap failed");

    for (pe = n->edgelist; pe; pe = pe->next) {
        rval = heap_checkout (H, CCbmat_THISEDGE(pe));
        CCcheck_rval (rval, "heap_checkout failed");
    }

CLEANUP:
    return rval;
}

void CCbmat_heap_edgeshrunk (CCbmat_heaps *H, CCbmat_edge *e)
{
    get_edge_from_heap (H, e);
    e->heap.type = CCbmat_NOTHEAPED;
}

void CCbmat_heap_shrink_node (CCbmat_heaps *H, CCbmat_node *n)
{
    get_node_from_heap (H, n);
    n->heap.type = CCbmat_NOTHEAPED;
}

static int heap_checkout (CCbmat_heaps *H, CCbmat_edge *e)
{
    int mult = 0, rval = 0;

    get_edge_from_heap (H, e);

    switch (e->surfends[0]->label) {
    case CCbmat_PLUSMAT: e->istooth[1] ? mult++ : mult--; break;
    case CCbmat_MINUSMAT: e->istooth[1] ? mult-- : mult++; break;
    }
    switch (e->surfends[1]->label) {
    case CCbmat_PLUSMAT: e->istooth[0] ? mult++ : mult--; break;
    case CCbmat_MINUSMAT: e->istooth[0] ? mult-- : mult++; break;
    }
    switch (mult) {
    case -2:
        e->heap.type = e->x ? MINUS2_MAT_HEAP : MINUS2_NOTMAT_HEAP;
        break;
    case -1:
        e->heap.type = e->x ? MINUS1_MAT_HEAP : MINUS1_NOTMAT_HEAP;
        break;
    case 0:
        e->heap.type = CCbmat_NOTHEAPED;
        break;
    case 1:
        e->heap.type = e->x ? PLUS1_MAT_HEAP : PLUS1_NOTMAT_HEAP;
        break;
    case 2:
        e->heap.type = e->x ? PLUS2_MAT_HEAP : PLUS2_NOTMAT_HEAP;
        break;
    }

    rval = add_edge_to_heap (H, e);
    CCcheck_rval (rval, "add_edge_to_heap failed");

CLEANUP:
    return rval;
}

static void get_node_from_heap (CCbmat_heaps *H, CCbmat_node *n)
{
    switch (n->heap.type) {
    case PLUS_NODE_HEAP:
        node_heap_delete (n, &H->plus_node_heap);
        n->ymat += H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        break;
    case MINUS_NODE_HEAP:
        node_heap_delete (n, &H->minus_node_heap);
        n->ymat -= H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        break;
    case MINUS_PSEUDO_HEAP:
        node_heap_delete (n, &H->minus_pseudo_heap);
        n->ymat -= H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->ymat >= 0);
        break;
    case CCbmat_NOTHEAPED:
        break;
    }
}

static int add_node_to_heap (CCbmat_heaps *H, CCbmat_node *n)
{
    int rval = 0;

    switch (n->heap.type) {
    case PLUS_NODE_HEAP:
        n->ymat -= H->delta;
        n->heap.key = 0;
        rval = node_heap_insert (n, &H->plus_node_heap);
        CCcheck_rval (rval, "node_heap_insert failed");
        break;
    case MINUS_NODE_HEAP:
        n->ymat += H->delta;
        n->heap.key = 0;
        rval = node_heap_insert (n, &H->minus_node_heap);
        CCcheck_rval (rval, "node_heap_insert failed");
        break;
    case MINUS_PSEUDO_HEAP:
        assert (n->ymat >= 0);
        n->ymat += H->delta;
        n->heap.key = n->ymat;
        rval = node_heap_insert (n, &H->minus_pseudo_heap);
        CCcheck_rval (rval, "node_heap_insert failed");
        break;
    case CCbmat_NOTHEAPED:
        break;
    }

CLEANUP:
    return rval;
}

static void get_edge_from_heap (CCbmat_heaps *H, CCbmat_edge *n)
{
    switch (n->heap.type) {
    case MINUS2_MAT_HEAP:
        edge_heap_delete (n, &H->minus2_mat_heap);
        n->slack -= 2*H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack <= 0);
        break;
    case MINUS2_NOTMAT_HEAP:
        edge_heap_delete (n, &H->minus2_notmat_heap);
        n->slack -= 2*H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack >= 0);
        break;
    case MINUS1_MAT_HEAP:
        edge_heap_delete (n, &H->minus1_mat_heap);
        n->slack -= H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack <= 0);
        break;
    case MINUS1_NOTMAT_HEAP:
        edge_heap_delete (n, &H->minus1_notmat_heap);
        n->slack -= H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack >= 0);
        break;
    case PLUS1_MAT_HEAP:
        edge_heap_delete (n, &H->plus1_mat_heap);
        n->slack += H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack <= 0);
        break;
    case PLUS1_NOTMAT_HEAP:
        edge_heap_delete (n, &H->plus1_notmat_heap);
        n->slack += H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack >= 0);
        break;
    case PLUS2_MAT_HEAP:
        edge_heap_delete (n, &H->plus2_mat_heap);
        n->slack += 2*H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack <= 0);
        break;
    case PLUS2_NOTMAT_HEAP:
        edge_heap_delete (n, &H->plus2_notmat_heap);
        n->slack += 2*H->delta;
        n->heap.type = CCbmat_NOTHEAPED;
        assert (n->slack >= 0);
        break;
    case CCbmat_NOTHEAPED:
        break;
    }
}

static int add_edge_to_heap (CCbmat_heaps *H, CCbmat_edge *n)
{
    int rval = 0;

    switch (n->heap.type) {
    case MINUS2_MAT_HEAP:
        assert (n->slack <= 0);
        n->slack += 2*H->delta;
        n->heap.key = 0;
        rval = edge_heap_insert (n, &H->minus2_mat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case MINUS2_NOTMAT_HEAP:
        assert (n->slack >= 0);
        n->slack += 2*H->delta;
        n->heap.key = n->slack;
        rval = edge_heap_insert (n, &H->minus2_notmat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case MINUS1_MAT_HEAP:
        assert (n->slack <= 0);
        n->slack += H->delta;
        n->heap.key = 0;
        rval = edge_heap_insert (n, &H->minus1_mat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case MINUS1_NOTMAT_HEAP:
        assert (n->slack >= 0);
        n->slack += H->delta;
        n->heap.key = n->slack;
        rval = edge_heap_insert (n, &H->minus1_notmat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case PLUS1_MAT_HEAP:
        assert (n->slack <= 0);
        n->slack -= H->delta;
        n->heap.key = - n->slack;
        rval = edge_heap_insert (n, &H->plus1_mat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case PLUS1_NOTMAT_HEAP:
        assert (n->slack >= 0);
        n->slack -= H->delta;
        n->heap.key = 0;
        rval = edge_heap_insert (n, &H->plus1_notmat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case PLUS2_MAT_HEAP:
        assert (n->slack <= 0);
        n->slack -= 2*H->delta;
        n->heap.key = - n->slack;
        rval = edge_heap_insert (n, &H->plus2_mat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case PLUS2_NOTMAT_HEAP:
        assert (n->slack >= 0);
        n->slack -= 2*H->delta;
        n->heap.key = 0;
        rval = edge_heap_insert (n, &H->plus2_notmat_heap);
        CCcheck_rval (rval, "edge_heap_insert failed");
        break;
    case CCbmat_NOTHEAPED:
        break;
    }

CLEANUP:
    return rval;
}

int CCbmat_heap_verify (CCbmat_graph *G, CCbmat_heaps *H)
{
    int i;
    CCbmat_node *p;
    CCbmat_edge *e;

    node_heap_verify_order(&H->plus_node_heap, PLUS_NODE_HEAP);
    node_heap_verify_order(&H->minus_node_heap, MINUS_NODE_HEAP);
    node_heap_verify_order(&H->minus_pseudo_heap, MINUS_PSEUDO_HEAP);

    edge_heap_verify_order(&H->minus2_notmat_heap, MINUS2_NOTMAT_HEAP);
    edge_heap_verify_order(&H->minus1_notmat_heap, MINUS1_NOTMAT_HEAP);
    edge_heap_verify_order(&H->plus1_notmat_heap, PLUS1_NOTMAT_HEAP);
    edge_heap_verify_order(&H->plus2_notmat_heap, PLUS2_NOTMAT_HEAP);
    edge_heap_verify_order(&H->minus2_mat_heap, MINUS2_MAT_HEAP);
    edge_heap_verify_order(&H->minus1_mat_heap, MINUS1_MAT_HEAP);
    edge_heap_verify_order(&H->plus1_mat_heap, PLUS1_MAT_HEAP);
    edge_heap_verify_order(&H->plus2_mat_heap, PLUS2_MAT_HEAP);

    for (i=0; i<G->ncount; i++) {
        for (p = &G->nodelist[i]; p; p = p->nest.parent) {
            if (p->nest.parent) {
                assert (p->heap.type == CCbmat_NOTHEAPED);
            } else if (p->label == CCbmat_PLUSMAT) {
                assert (p->heap.type == PLUS_NODE_HEAP);
                assert (node_heap_contains (p, &H->plus_node_heap));
            } else if (p->label == CCbmat_MINUSMAT) {
                if (p->status.pseudo) {
                    assert (p->heap.type == MINUS_PSEUDO_HEAP);
                    assert (node_heap_contains (p, &H->minus_pseudo_heap));
                } else {
                    assert (p->heap.type == MINUS_NODE_HEAP);
                    assert (node_heap_contains (p, &H->minus_node_heap));
                }
            } else {
                assert (p->heap.type == CCbmat_NOTHEAPED);
            }
        }
    }

    for (i=0, e = G->edgelist; i < G->ecount; i++, e++) {
        assert (edge_correct_type(e));
    }
    for (e = G->newedges; e; e = e->next) {
        assert (edge_correct_type(e));
    }
    return 1;
}

#if DEBUG
static int edge_correct_type (CCbmat_edge *e)
{
    int mult;

    if (e->surfends[0]->nest.parent) {
        assert (e->surfends[1]->nest.parent == e->surfends[0]->nest.parent);
        assert (e->heap.type == CCbmat_NOTHEAPED);
    } else {
        mult = 0;
        if (e->surfends[0]->label == CCbmat_PLUSMAT) {
            e->istooth[1] ? mult-- : mult++;
        } else if (e->surfends[0]->label == CCbmat_MINUSMAT) {
            e->istooth[1] ? mult++ : mult--;
        }
        if (e->surfends[1]->label == CCbmat_PLUSMAT) {
            e->istooth[0] ? mult-- : mult++;
        } else if (e->surfends[1]->label == CCbmat_MINUSMAT) {
            e->istooth[0] ? mult++ : mult--;
        }
        if (e->x) {
            switch (mult) {
            case -2:
                assert (e->heap.type == PLUS2_MAT_HEAP);
                assert (edge_heap_contains (e, &plus2_mat_heap));
                break;
            case -1:
                assert (e->heap.type == PLUS1_MAT_HEAP);
                assert (edge_heap_contains (e, &plus1_mat_heap));
                break;
            case 0:
                assert (e->heap.type == CCbmat_NOTHEAPED);
                break;
            case 1:
                assert (e->heap.type == MINUS1_MAT_HEAP);
                assert (edge_heap_contains (e, &minus1_mat_heap));
                break;
            case 2:
                assert (e->heap.type == MINUS2_MAT_HEAP);
                assert (edge_heap_contains (e, &minus2_mat_heap));
                break;
            default:
                assert (mult != mult);
                break;
            }
        } else {
            switch (mult) {
            case -2:
                assert (e->heap.type == PLUS2_NOTMAT_HEAP);
                assert (edge_heap_contains (e, &plus2_notmat_heap));
                break;
            case -1:
                assert (e->heap.type == PLUS1_NOTMAT_HEAP);
                assert (edge_heap_contains (e, &plus1_notmat_heap));
                break;
            case 0:
                assert (e->heap.type == CCbmat_NOTHEAPED);
                break;
            case 1:
                assert (e->heap.type == MINUS1_NOTMAT_HEAP);
                assert (edge_heap_contains (e, &minus1_notmat_heap));
                break;
            case 2:
                assert (e->heap.type == MINUS2_NOTMAT_HEAP);
                assert (edge_heap_contains (e, &minus2_notmat_heap));
                break;
            default:
                assert (mult != mult);
                break;
            }
        }
    }
    return 1;
}
#endif

/****  Old heaps.c file  ****/

#define INIT_HEAP_SIZE ((65536-16)/sizeof (CCbmat_node *))
#define HEAP_D 3
#define HEAP_UP(x) (((x)-1)/HEAP_D)
#define HEAP_DOWN(x) (((x)*HEAP_D)+1)

static void node_heap_init_struct (CCbmat_node_heap *h)
{
    h->entry = (CCbmat_node **) NULL;
    h->total_space = 0;
    h->size = 0;
}

static int node_heap_init (CCbmat_node_heap *h)
{
    int rval = 0;

    h->entry = CC_SAFE_MALLOC (INIT_HEAP_SIZE, CCbmat_node *);
    CCcheck_NULL (h->entry, "out of memory for heap entries");
    h->total_space = INIT_HEAP_SIZE;
    h->size = 0;

CLEANUP:
    return rval;
}

static void node_heap_free (CCbmat_node_heap *h)
{
   if (h) {
       CC_IFFREE (h->entry, CCbmat_node *);    
   }
}

static int node_heap_expand (CCbmat_node_heap *h)
{
    int rval = 0;

    h->total_space *= 2;
    h->entry = (CCbmat_node **) CCutil_reallocrus ((void *) h->entry,
                           h->total_space * sizeof (CCbmat_node *));
    CCcheck_NULL (h->entry, "out of memory for realloc of heap entries");

CLEANUP:
    return rval;
}

CCbmat_node *CCbmat_node_heap_findmin (CCbmat_node_heap *h)
{
    if (h->size == 0) return (CCbmat_node *) NULL;
    else return h->entry[0];
}

static int node_heap_insert (CCbmat_node *i, CCbmat_node_heap *h)
{
    int rval = 0;

    if (h->size >= h->total_space) {
        rval = node_heap_expand (h);
        CCcheck_rval (rval, "node_heap_expand failed");
    }
    h->size++;
    node_heap_siftup (i, h->size-1, h);

CLEANUP:
    return rval;
}

static void node_heap_delete (CCbmat_node *i, CCbmat_node_heap *h)
{
    CCbmat_node *j;

    h->size--;
    j = h->entry[h->size];
    h->entry[h->size] = (CCbmat_node *) NULL;

    if (j != i) {
        if (j->heap.key <= i->heap.key) {
            node_heap_siftup (j, i->heap.loc, h);
        } else {
            node_heap_siftdown (j, i->heap.loc, h);
        }
    }
}

CCbmat_node *CCbmat_node_heap_deletemin (CCbmat_node_heap *h)
{
    CCbmat_node *i;

    if (h->size == 0) return (CCbmat_node *) NULL;
    else {
        i = h->entry[0];
        node_heap_delete (i, h);
        return i;
    }
}

static void node_heap_siftup (CCbmat_node *i, int x, CCbmat_node_heap *h)
{
    int p;

    p = HEAP_UP(x);
    while (x && h->entry[p]->heap.key > i->heap.key) {
        h->entry[x] = h->entry[p];
        h->entry[p]->heap.loc = x;
        x = p;
        p = HEAP_UP(p);
    }
    h->entry[x] = i;
    i->heap.loc = x;
}

static void node_heap_siftdown (CCbmat_node *i, int x, CCbmat_node_heap *h)
{
    int c;

    c = node_heap_minchild (x, h);

    while (c >= 0 && h->entry[c]->heap.key < i->heap.key) {
        h->entry[x] = h->entry[c];
        h->entry[c]->heap.loc = x;
        x = c;
        c = node_heap_minchild (c, h);
    }
    h->entry[x] = i;
    i->heap.loc = x;
}

static int node_heap_minchild (int x, CCbmat_node_heap *h)
{
    int c = HEAP_DOWN(x);
    int cend;
    int minval;
    int minloc;

    if (c >= h->size) return -1;
    minval = h->entry[c]->heap.key;
    minloc = c;
    cend = c + HEAP_D;
    if (h->size < cend) cend = h->size;
    for (c++; c < cend; c++) {
        if (h->entry[c]->heap.key < minval) {
            minval = h->entry[c]->heap.key;
            minloc = c;
        }
    }
    return minloc;
}

static void node_heap_purge (CCbmat_node_heap *h, int delta)
{
    CCbmat_node **p, **pend;
    CCbmat_node *q;

    for (p = h->entry, pend = p + h->size; p < pend; p++) {
        q = *p;
        q->heap.type = CCbmat_NOTHEAPED;
        q->ymat += delta;
    }
    h->size = 0;
}

static void edge_heap_init_struct (CCbmat_edge_heap *h)
{
    if (h) {
        h->entry = (CCbmat_edge **) NULL;
        h->total_space = 0;
        h->size = 0;
    }
}

static int edge_heap_init (CCbmat_edge_heap *h)
{
    int rval = 0;

    h->entry = CC_SAFE_MALLOC (INIT_HEAP_SIZE, CCbmat_edge *);
    CCcheck_NULL (h->entry, "out of memory for edge hash entries");
    h->total_space = INIT_HEAP_SIZE;
    h->size = 0;

CLEANUP:
    return rval;
}

static void edge_heap_free (CCbmat_edge_heap *h)
{
   if (h) {
       CC_IFFREE (h->entry, CCbmat_edge *);    
   }
}

static int edge_heap_expand (CCbmat_edge_heap *h)
{
    int rval = 0;

    h->total_space *= 2;
    h->entry = (CCbmat_edge **) CCutil_reallocrus ((void *) h->entry,
                           h->total_space * sizeof (CCbmat_edge *));
    CCcheck_NULL (h->entry, "out of memory for realloc of edge heap entries");

CLEANUP:
    return rval;
}

CCbmat_edge *CCbmat_edge_heap_findmin (CCbmat_edge_heap *h)
{
    if (h->size == 0) return (CCbmat_edge *) NULL;
    else return h->entry[0];
}

static int edge_heap_insert (CCbmat_edge *i, CCbmat_edge_heap *h)
{
    int rval = 0;

    if (h->size >= h->total_space) {
        rval = edge_heap_expand (h);
        CCcheck_rval (rval, "edge_heap_expand failed");
    }
    h->size++;
    edge_heap_siftup (i, h->size-1, h);

CLEANUP:
    return rval;
}

static void edge_heap_delete (CCbmat_edge *i, CCbmat_edge_heap *h)
{
    CCbmat_edge *j;

    h->size--;
    j = h->entry[h->size];
    h->entry[h->size] = (CCbmat_edge *) NULL;

    if (j != i) {
        if (j->heap.key <= i->heap.key) {
            edge_heap_siftup (j, i->heap.loc, h);
        } else {
            edge_heap_siftdown (j, i->heap.loc, h);
        }
    }
}

CCbmat_edge *CCbmat_edge_heap_deletemin (CCbmat_edge_heap *h)
{
    CCbmat_edge *i;

    if (h->size == 0) return (CCbmat_edge *) NULL;
    else {
        i = h->entry[0];
        edge_heap_delete (i, h);
        return i;
    }
}

static void edge_heap_siftup (CCbmat_edge *i, int x, CCbmat_edge_heap *h)
{
    int p;

    p = HEAP_UP(x);
    while (x && h->entry[p]->heap.key > i->heap.key) {
        h->entry[x] = h->entry[p];
        h->entry[p]->heap.loc = x;
        x = p;
        p = HEAP_UP(p);
    }
    h->entry[x] = i;
    i->heap.loc = x;
}

static void edge_heap_siftdown (CCbmat_edge *i, int x, CCbmat_edge_heap *h)
{
    int c;

    c = edge_heap_minchild (x, h);

    while (c >= 0 && h->entry[c]->heap.key < i->heap.key) {
        h->entry[x] = h->entry[c];
        h->entry[c]->heap.loc = x;
        x = c;
        c = edge_heap_minchild (c, h);
    }
    h->entry[x] = i;
    i->heap.loc = x;
}

static int edge_heap_minchild (int x, CCbmat_edge_heap *h)
{
    int c = HEAP_DOWN(x);
    int cend;
    int minval;
    int minloc;
    
    if (c >= h->size) return -1;
    minval = h->entry[c]->heap.key;
    minloc = c;
    cend = c + HEAP_D;
    if (h->size < cend) cend = h->size;
    for (c++; c < cend; c++) {
        if (h->entry[c]->heap.key < minval) {
            minval = h->entry[c]->heap.key;
            minloc = c;
        }
    }
    return minloc;
}

static void edge_heap_purge (CCbmat_edge_heap *h, int delta)
{
    CCbmat_edge **p, **pend;
    CCbmat_edge *q;

    for (p = h->entry, pend = p + h->size; p < pend; p++) {
        q = *p;
        q->heap.type = CCbmat_NOTHEAPED;
        q->slack += delta;
    }
    h->size = 0;
}

static int node_heap_verify_order (CCbmat_node_heap *h, int type)
{
    int i;

    for (i=1; i<h->size; i++) {
        assert (h->entry[i]->heap.key >= h->entry[HEAP_UP(i)]->heap.key);
    }
    for (i=0; i<h->size; i++) {
        assert (h->entry[i]->heap.type == type);
        assert (h->entry[i]->heap.loc == i);
    }
    return 1;
}

static int edge_heap_verify_order (CCbmat_edge_heap *h, int type)
{
    int i;

    for (i=1; i<h->size; i++) {
        assert (h->entry[i]->heap.key >= h->entry[HEAP_UP(i)]->heap.key);
    }
    for (i=0; i<h->size; i++) {
        assert (h->entry[i]->heap.type == type);
        assert (h->entry[i]->heap.loc == i);
    }
    return 1;
}

#if DEBUG
static int node_heap_contains (CCbmat_node *n, CCbmat_node_heap *h)
{
    if (h->entry[n->heap.loc] == n) return 1;
    else return 0;
}
#endif

#if DEBUG
static int edge_heap_contains (CCbmat_edge *n, CCbmat_edge_heap *h)
{
    if (h->entry[n->heap.loc] == n) return 1;
    else return 0;
}
#endif

