/**************************************************************************/
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "pancake.h"

#define MAGICNODE 0
#define PANEPSILON 0.0001
#define PANFALSE 0
#define PANPOSITIVE 1
#define PANFEW 1
#define PANSHORT 1000
#define PANBIGNEG -10000000000.0
#define PANSWAP(x,y,temp) {temp = x; x = y; y = temp;}

static int panalloc (int ecount, int *elist);
static int buildpanadjlist (int ecount, int *elist);
static int pancakex (int ecount, int *elist, double *x);
static void panfree (void);

static int initpancake (int ecount);
static int buildfirsttree (int ecount);
static int decompositiontree (int ecount);
static int initdecompositiontree (int ecount);
static void drop (panedge *e, vaseknode *x);
static void throw (vaseknode *x, vaseknode *y, panedge *e);
static vaseknode *anc (vaseknode *v);
static void trickledown (int i);
static vaseknode *newcomp (vaseknode *v, double w);
static void hookup (vaseknode *parent, vaseknode *child);
static void distribute (void);
static void initdistribute (void);
static void split (vaseknode *a);
static void bruteforce (panedge *e);
static void update (panedge *e);
static void dealwith (panedge *e, vaseknode **pa);
static void attach (panedge *e);
static void magicrc (void);
static double min2 (panedge **elist);
static double findbound (void);
static double blnode (vaseknode *v, int *hit, int *cutarray, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param);
static void labeler (vaseknode *p, int *cutcount, int *cutarray);
static void freepancake (void);
static int pancakemain (int ecount);
static int blobsviolated (int *hit, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param);

static int vstep;
static vaseknode *vroot;
static vaseknode *vpannodes, *vnodehit;
static vaseknode *head, *tail;
static vaseknode *vnodestack;
static triomino *trisupply, *tristack;
static panedge *work;
static panedge **vheap;
static int vheapend;
static int vcomponentcount;
static panedge **panedgespace, *panedgelist;
static pannode *pannodelist;
static int nnodes;

int CCtest_blobcuts (int ncount, int ecount, int *elist, double *x, int *hit,
        int (*doit_fn) (double, int, int *, void *), void *pass_param)
{
    int rval = 0;

    /* printf ("CCtest_blobcuts ...\n"); fflush (stdout); */

    nnodes = ncount;

    pancakex (ecount, elist, x);

    rval = blobsviolated (hit, x, doit_fn, pass_param);
    if (rval) {
        fprintf (stderr, "blobsviolated failed\n");  goto CLEANUP;
    }
    freepancake ();
    panfree ();

CLEANUP:

    return rval;
}

static int pancakex (int ecount, int *elist, double *x)
{
    int i;
    int rval = 0;

    /* printf ("pancakex ...\n"); fflush (stdout); */

    rval = panalloc (ecount, elist);
    if (rval) {
        fprintf (stderr, "panalloc failed\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < ecount; i++) {
        panedgelist[i].panweight = -x[i];
    }

    pancakemain (ecount);

CLEANUP:

    return rval;
}

static int panalloc (int ecount, int *elist)
{
    int i;
    panedge *pa;
    int rval = 0;

    /* printf ("panalloc (%d)  ...\n", ecount); fflush (stdout); */

    pannodelist = (pannode *) malloc (nnodes * sizeof (pannode));
    panedgelist = (panedge *) malloc (ecount * sizeof (panedge));
    if (!pannodelist || !panedgelist) {
        fprintf (stderr, "out of memory in panalloc\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, pa = panedgelist; i < ecount; i++, pa++) {
        pa->ends[0] = pannodelist + elist[2*i];
        pa->ends[1] = pannodelist + elist[2*i+1];
    }

    rval = buildpanadjlist (ecount, elist);
    if (rval) {
        fprintf (stderr, "buildpanadjlist failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static void panfree (void)
{
    if (pannodelist) {
        free (pannodelist);
        pannodelist = NULL;
    }
    if (panedgelist) {
        free (panedgelist);
        panedgelist = NULL;
    }
    if (panedgespace) {
        free (panedgespace);
        panedgespace = NULL;
    }
}

static int buildpanadjlist (int ecount, int *elist)
{
    int i, *degrees = (int *) NULL, *pint;
    pannode *pv;
    panedge *pa, **edgespace;
    int rval = 0;

    /* printf ("buildpanadjlist ...\n"); fflush (stdout); */

    edgespace = (panedge **) malloc ((nnodes+(2*ecount)) * sizeof (panedge *));
    degrees = (int *) malloc (nnodes * sizeof (int));
    if (!edgespace || !degrees) {
        fprintf (stderr, "out of memory in buildpanadjlist\n");
        rval = 1;  goto CLEANUP;
    }

    panedgespace = edgespace;
    for (i = 0, pint = degrees; i < nnodes; i++) {
        *pint++ = 0;
    }

    for (i = 0; i < ecount; i++) {
        degrees[elist[2*i]]++;
        degrees[elist[2*i+1]]++;
    }

    for (i = 0, pint = degrees, pv = pannodelist; i < nnodes;
         i++, pint++, pv++) {
        if (*pint) {
            pv->edgelist = pv->goodedge = edgespace;
            edgespace += *pint + 1;
        } else
            pv->edgelist = pv->goodedge = (panedge **) NULL;
    }

    for (i = 0, pa = panedgelist; i < ecount; i++, pa++) {
        *(pa->ends[0]->goodedge++) = pa;
        *(pa->ends[1]->goodedge++) = pa;
    }

    for (i = 0, pv = pannodelist; i < nnodes; i++, pv++) {
        *(pv->goodedge) = (panedge *) NULL;
    }

CLEANUP:

    if (degrees) {
        free (degrees);
    }
    return rval;
}

/**************************************************************************/
/*                         Old dtree.c                                    */
/**************************************************************************/

static int pancakemain (int ecount)
{
    int rval = 0;

    rval = initpancake (ecount);
    if (rval) {
        fprintf (stderr, "initpancake failed\n"); goto CLEANUP;
    }
    rval = buildfirsttree (ecount);
    if (rval) {
        fprintf (stderr, "buildfirsttree failed\n"); goto CLEANUP;
    }

    findbound ();

CLEANUP:

    return rval;
}

static int initpancake (int ecount)
{
    int i;
    pannode *pn;
    vaseknode *pv;
    triomino *pt;
    int rval = 0;

    vpannodes = (vaseknode *) malloc ((2*nnodes-1) * sizeof (vaseknode));
    if (!vpannodes) {
        fprintf (stderr, "out of memory in initpancake\n");
        rval = 1; goto CLEANUP;
    } 

    for (i = nnodes, pn = pannodelist, pv = vpannodes; i; i--, pn++, pv++) {
        pv->parent = (vaseknode *) NULL;
        pv->child = (vaseknode *) NULL;
        pv->sibling = (vaseknode *) NULL;
        pv->adj = (triomino *) NULL;
        pv->n = 0;
        pv->b = 1;
        pv->anc = pv;
        pv->junk = (panedge *) NULL;
        pv->w = PANBIGNEG;
        pv->y = 0.0;
        pn->vptr = pv;
    }

    vnodestack = vpannodes + nnodes;
    for (i = nnodes - 2, pv = vnodestack; i; i--, pv++)
        pv->ptr = pv + 1;
    pv->ptr = (vaseknode *) NULL;

    trisupply = (triomino *) malloc ((2*ecount) * sizeof (triomino));
    if (!trisupply) {
        fprintf (stderr, "out of memory for trisupply\n");
        rval = 1;  goto CLEANUP;
    }
    tristack = trisupply;
    for (i = 2 * ecount - 1, pt = tristack; i; i--, pt++)
        pt->next = pt + 1;

    VNODEALLOC (head);

CLEANUP:

    return rval;
}

static void freepancake (void)
{
    if (vpannodes) {
        free (vpannodes);
        vpannodes = NULL;
    }
    if (trisupply) {
        free (trisupply);
        trisupply = NULL;
    }
    if (vheap) {
        free (vheap);
        vheap = NULL;
    }
}

static int buildfirsttree (int ecount)
{
    vaseknode *p;
    double parw;
    vaseknode *q;
    int rval = 0;

    rval = decompositiontree (ecount);
    if (rval) {
        fprintf (stderr, "decompositiontree failed\n"); goto CLEANUP;
    }
    distribute ();


    head->ptr = tail = vroot;
    p = head;
    do {
        p = p->ptr;
        parw = p->w;
        for (q = p->child; q != (vaseknode *) NULL; q = q->sibling) {
            q->mult = parw - q->w;
            if (q->child != (vaseknode *) NULL) {
                tail->ptr = q;
                tail = q;
            }
        }
    } while (p != tail);
    vroot->mult = -vroot->w;

    magicrc ();

CLEANUP:

    return rval;
}

static int decompositiontree (int ecount)
{
    panedge *e;
    int i;
    double w, ub;
    vaseknode *x, *y;
    int rval = 0;

    /* printf ("decompositiontree ...\n"); fflush (stdout); */

    rval = initdecompositiontree (ecount);
    if (rval) {
        fprintf (stderr, "initdecompositiontree failed\n"); goto CLEANUP;
    }

    for (vcomponentcount = nnodes - 2; vcomponentcount;) {
        for (;;) {
            e = *vheap;
            *vheap = vheap[vheapend--];
            trickledown (0);
            x = anc (e->ends[0]->vptr);
            y = anc (e->ends[1]->vptr);
            if (x != y)
                break;
            drop (e, x);
        }

        w = e->panweight;
        throw (x, y, e);

        ub = w + 0.01;

        while (vheapend >= 0 && (e = *vheap)->panweight < ub) {
            *vheap = vheap[vheapend--];
            trickledown (0);
            x = anc (e->ends[0]->vptr);
            y = anc (e->ends[1]->vptr);
            if (x != y)
                throw (x, y, e);
            else
                drop (e, x);
        }

        for (; vnodehit != (vaseknode *) NULL; vnodehit = vnodehit->ptr) {
            if (vnodehit->n)
                vroot = newcomp (vnodehit, w);
        }
    }
    i = vheapend + 1;
    while (i) {
	drop (vheap[--i], vroot);
    }

CLEANUP:

    return rval;
}

static int initdecompositiontree (int ecount)
{
    int i;
    panedge *pe, **ph;
    pannode *pm;
    int rval = 0;

    vnodehit = (vaseknode *) NULL;
    work = (panedge *) NULL;

    for (i = ecount, pe = panedgelist; i; i--, pe++)
        pe->tag = PANFALSE;

    vheap = (panedge **) malloc (ecount * sizeof (panedge *));
    if (!vheap) {
        fprintf (stderr, "out of memory for vheap\n");
        rval = 1;  goto CLEANUP;
    }

    pm = pannodelist + MAGICNODE;
    for (i = ecount, pe = panedgelist, ph = vheap; i; i--, pe++)
        if (pe->ends[0] != pm && pe->ends[1] != pm)
            *(ph++) = pe;
        else
            pe->rc = pe->panweight;

    vheapend = (ph - vheap) - 1;
    for (i = vheapend / 2; i >= 0; i--) {
        trickledown (i);
    }

CLEANUP:

    return rval;
}

static void drop (panedge *e, vaseknode *x)
{
    e->a[0] = e->ends[0]->vptr;
    e->a[1] = e->ends[1]->vptr;
    e->roof = x;
    e->rc = PANPOSITIVE;
    e->next = work;
    work = e;
}

static void throw (vaseknode *x, vaseknode *y, panedge *e)
{
    e->a[0] = x;
    e->a[1] = y;

    e->rc = 0.0;
    attach (e);

    if (!(x->n)) {
        x->n = 1;
        x->ptr = vnodehit;
        vnodehit = x;
    }
    if (!(y->n)) {
        y->n = 1;
        y->ptr = vnodehit;
        vnodehit = y;
    }
    e->next = work;
    work = e;
}

static vaseknode *anc (vaseknode *v)
{
    vaseknode *hand, *va;

    hand = v;
    while (hand != hand->anc)
	hand = hand->anc;

    va = hand;
    for (hand = v; hand != va; hand = hand->anc)
        hand->anc = va;

    return va;
}

static void trickledown (int i)
{
    panedge *memo;
    int k, minchild;

    memo = vheap[i];

    while ((k = (2 * i) + 2) <= vheapend) {
        minchild = (vheap[k - 1]->panweight <= vheap[k]->panweight ? k - 1 : k);

        if (memo->panweight > vheap[minchild]->panweight) {
            vheap[i] = vheap[minchild];
            i = minchild;
        } else {
            vheap[i] = memo;
            return;
        }
    }
    if (k - 1 == vheapend && memo->panweight > vheap[vheapend]->panweight) {
        vheap[i] = vheap[vheapend];
        i = vheapend;
    }
    vheap[i] = memo;
}

static vaseknode *newcomp (vaseknode *v, double w)
{
    vaseknode *new, *stack;
    triomino *t;

    VNODEALLOC (new);
    new->parent = (vaseknode *) NULL;
    new->child = (vaseknode *) NULL;
    new->sibling = (vaseknode *) NULL;
    new->anc = new;
    new->n = 0;
    new->b = 0;
    new->w = w;
    new->adj = (triomino *) NULL;
    new->junk = (panedge *) NULL;
    new->tag = PANFALSE;

    hookup (new, v);
    v->qtr = (vaseknode *) NULL;

    do {
        stack = v->qtr;
        for (t = v->adj; t != (triomino *) NULL; t = t->next) {
            t->edge->roof = new;
            v = t->end;
            if (v->n) {
                v->qtr = stack;
                stack = v;
                hookup (new, v);
                vcomponentcount--;
            }
        }
    } while ((v = stack) != (vaseknode *) NULL);

    return new;
}

static void hookup (vaseknode *parent, vaseknode *child)
{
    child->n = 0;
    child->parent = parent;
    child->anc = parent;
    child->sibling = parent->child;
    parent->child = child;
}

static void distribute (void)
{
    vaseknode *active;
    panedge *e, *f;

    /* printf ("distribute ...\n"); fflush (stdout); */

    initdistribute ();
    active = (vaseknode *) NULL;
    vroot->n = 0;

    for (e = work, work = (panedge *) NULL; e != (panedge *) NULL; e = f) {
        f = e->next;
        e->top = vroot;
        dealwith (e, &active);
    }

    while (work) {
        for (; active != (vaseknode *) NULL; active = active->ptr)
            if (active->n < PANFEW)
                active->n = -1;
            else {
                active->n = 0;
                split (active);
            }
        for (e = work, work = (panedge *) NULL; e != (panedge *) NULL;
                                                                    e = f) {
            f = e->next;
            if (e->top->n >= 0) {
                update (e);
                dealwith (e, &active);
            } else
                bruteforce (e);
        }
        vstep = vstep / 2;
    }
}

static void initdistribute (void)
{
    vaseknode *stack, *finger, *x;
    int maxd, twice, d;

    maxd = 0;

    vroot->d = 0;
    vroot->ptr = (vaseknode *) NULL;

    for (finger = vroot; finger != (vaseknode *) NULL; finger = stack) {
        stack = finger->ptr;
        finger->anc = vroot;
        if ((x = finger->child) != (vaseknode *) NULL) {
            d = finger->d + 1;
            do {
                x->d = d;
                x->ptr = stack;
                stack = x;
                x = x->sibling;
            } while (x != (vaseknode *) NULL);
            if (d > maxd)
                maxd = d;
        }
    }
    vstep = 1;
    twice = 2;
    while (twice < maxd) {
	vstep = twice;
	twice = vstep + vstep;
    }
}

static void split (vaseknode *a)
{
    int mid, bot;
    vaseknode *stack, *hand, *foot, *memo, *x;

    mid = vstep + a->d;
    bot = vstep + mid;

    a->qtr = (vaseknode *) NULL;
    for (hand = a; hand != (vaseknode *) NULL; hand = stack) {
        stack = hand->qtr;
        if (hand->d == mid) {
            memo = hand->qtr;
            hand->qtr = (vaseknode *) NULL;
            for (foot = hand; foot != (vaseknode *) NULL; foot = stack) {
                stack = foot->qtr;
                foot->anc = hand;
                if (foot->d != bot) {
                    for (x = foot->child; x != (vaseknode *) NULL;
                         x = x->sibling) {
                        x->qtr = stack;
                        stack = x;
                    }
                }
            }
            hand->qtr = memo;
        } else
            for (x = hand->child; x != (vaseknode *) NULL; x = x->sibling) {
                x->qtr = stack;
                stack = x;
            }
    }
}

static void bruteforce (panedge *e)
{
    vaseknode *x, *y, *nx, *ny;
    int dx, dy;

    x = e->a[0];
    y = e->a[1];

    if (x == y) {
        printf ("Tough luck Pal 1.\n");
        exit (1);
    }
    dx = x->d;
    dy = y->d;

    while (dx > dy) {
        x = x->parent;
        dx--;
    }
    if (x == y) {
        printf ("Tough luck Pal 2.\n");
        exit (1);
    }
    while (dy > dx) {
        y = y->parent;
        dy--;
    }
    if (x == y) {
        printf ("Tough luck Pal 3.\n");
        exit (1);
    }
    nx = x->parent;
    ny = y->parent;
    while (nx != ny) {
	x = nx;
	y = ny;
	nx = x->parent;
	ny = y->parent;
    }

    e->a[0] = x;
    e->a[1] = y;

    e->roof = nx;
    /* if (e->rc > 0.0) { e->next = nx->junk; nx->junk = e; } else { e->tag =
     * PANFALSE; attach (e); } */
    e->next = nx->junk;
    nx->junk = e;
}

static void update (panedge *e)
{
    vaseknode *x, *y, *v;

    x = e->a[0]->anc;
    y = e->a[1]->anc;
    v = e->top;

    if (x == v) {
        if (y != v)
            e->a[1] = y;
    } else if (y == v)
        e->a[0] = x;
    else if (x != y) {
        e->a[0] = x;
        e->a[1] = y;
    } else {
        e->top = x;
        if (x->d > e->roof->d)
            e->roof = x;
    }
}

static void dealwith (panedge *e, vaseknode **pa)
{
    if ((e->roof->d) - (e->a[0]->d) < PANSHORT &&
        (e->roof->d) - (e->a[1]->d) < PANSHORT) {
        bruteforce (e);
    } else {
        e->next = work;
        work = e;
        if (!e->top->n) {
            e->top->ptr = *pa;
            *pa = e->top;
        }
        (e->top->n)++;
    }
}

static void attach (panedge *e)
{
    triomino *cell;

    TRIALLOC (cell);
    cell->edge = e;
    cell->end = e->a[1];
    cell->next = e->a[0]->adj;
    e->a[0]->adj = cell;

    TRIALLOC (cell);
    cell->edge = e;
    cell->end = e->a[0];
    cell->next = e->a[1]->adj;
    e->a[1]->adj = cell;
}

static void magicrc (void)
{
    double a;
    panedge **pee, *pe;

    a = min2 ((pannodelist + MAGICNODE)->edgelist);

    for (pee = (pannodelist + MAGICNODE)->edgelist;
                (pe = *pee) != (panedge *) NULL; pee++)
        pe->rc -= a;

    (pannodelist + MAGICNODE)->vptr->y += a;
}

static double min2 (panedge **elist)
{
    double minweight, minweight2, td;
    panedge *e;

    if (elist == (panedge **) NULL ||
     elist[0] == (panedge *) NULL || elist[1] == (panedge *) NULL) {
        fprintf (stderr, "Vertex has degree < two\n");
        exit (1);
    }
    minweight = elist[0]->rc;
    minweight2 = elist[1]->rc;
    if (minweight > minweight2) {
        PANSWAP (minweight, minweight2, td);
    }
    for (elist += 2; (e = *elist) != (panedge *) NULL; elist++) {
        if (e->rc < minweight2) {
            minweight2 = e->rc;
            if (minweight > minweight2) {
                PANSWAP (minweight, minweight2, td);
            }
        }
    }
    return minweight2;
}

static double findbound (void)
{
    vaseknode *p, *q, *stack;
    triomino *tri;
    panedge **pee, *pe;
    double tree_bound = 0.0;
    double star_bound = 0.0;
    double edge_bound = 0.0;

    vroot->ptr = (vaseknode *) NULL;
    vroot->n = 1;
    vroot->b = 0;


    for (p = stack = vroot; p; p = stack) {
        if (p->n) {
            p->n = 0;
            q = p->child;
            if (q)
                for (; q; q = q->sibling) {
                    q->ptr = stack;
                    stack = q;
                    q->n = 1;
                    q->b = 0;
                }
            else {
                stack = p->ptr;
                (p->parent->b)++;
                star_bound += p->y;
            }
            for (tri = p->adj; tri; tri = tri->next)
		edge_bound += tri->edge->rc;
        } else {
            stack = p->ptr;
            if (stack)
                (p->parent->b) += p->b;
            tree_bound -= (p->mult) * ((p->b) - 1);
        }
    }
    star_bound *= 2.0;
    edge_bound /= 2.0;

    for (pee = ((pannodelist + MAGICNODE)->edgelist);
                 (pe = *pee) != (panedge *) NULL; pee++)
        if (pe->rc < 0.0)
            edge_bound += pe->rc;
    star_bound += (pannodelist + MAGICNODE)->vptr->y * 2;

    return tree_bound + star_bound + edge_bound;
}

static int blobsviolated (int *hit, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param)
{
    int *cutarray = (int *) NULL;
    int rval = 0;

    /* printf ("blobsviolated ....\n"); fflush (stdout); */

    if (hit) *hit = 0;

    cutarray = (int *) malloc (nnodes * sizeof (int));
    if (!cutarray) {
        fprintf (stderr, "out of memory for cutarray\n");
        rval = 1; goto CLEANUP;
    }

    blnode (vroot, hit, cutarray, x, doit_fn, pass_param);

CLEANUP:
  
    if (cutarray) free (cutarray);
    return rval;
}

#define CUTTOLERANCE 0.01

static double blnode (vaseknode *v, int *hit, int *cutarray, double *x,
        int (*doit_fn) (double, int, int *, void *), void *pass_param)
{
    double w = 0.0;
    double t;
    panedge *e;
    vaseknode *c;
    int rval = 0;
    int cutcount;

    if (!v->child)
        return 0.0;
    else {
        for (e = v->junk; e; e = e->next)
            w += (x[e - panedgelist]);
        for (c = v->child; c; c = c->sibling)
            w += blnode (c, hit, cutarray, x, doit_fn, pass_param);
        t = v->b;
        if (w > t - 1.0 + CUTTOLERANCE) {
            cutcount = 0;
            labeler (v, &cutcount, cutarray);
            if (doit_fn) {
                rval = doit_fn (1.9, cutcount, cutarray, pass_param);
                if (rval) {
                    fprintf (stderr, "doit_fn failed\n");
                    exit (1);
                }
            }
            (*hit)++;
        }
        return w;
    }
}

static void labeler (vaseknode *p, int *cutcount, int *cutarray)
{
    vaseknode *c;

    if (!p->child) {
        cutarray[*cutcount] = p - vpannodes;
        (*cutcount)++;
    } else {
        for (c = p->child; c; c = c->sibling) {
            labeler (c, cutcount, cutarray);
        }
    }
}

