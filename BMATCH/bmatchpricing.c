#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bmatch.h"

#define DEBUG 0

static CCbmat_edge *build_edge (CCbmat_node *p1, CCbmat_node *p2,
    CCbmat_edge **edge_freelist, CCbmat_edgeptr **edge_supply);
static void contract_edge (CCbmat_edge *old, CCbmat_edge *new, CCbmat_node *pdel,
    CCbmat_node *p);
static int sum_to_root (CCbmat_node *p);
static int get_sum_work (CCbmat_node *n);
static void clear_sum_work (CCbmat_node *n);
static void initlist (CCbmat_graph *G, CCbmat_node *head, CCbmat_node *tail,
    CCbmat_node *head2, CCbmat_node *tail2);
static int price_scan (CCbmat_graph *G, CCdatagroup *D, CCbmat_heaps *H,
    CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply);

#if DEBUG
static int sum_to_level (CCbmat_node *p, CCbmat_node *ptop);
#endif

static CCbmat_edge *build_edge (CCbmat_node *p1, CCbmat_node *p2,
        CCbmat_edge **edge_freelist, CCbmat_edgeptr **edge_supply)
{
    CCbmat_edge *newe;

    assert (!p1->status.pseudo && !p2->status.pseudo);
    assert (p1->surf != p2->surf);

    newe = CCbmat_edge_alloc (edge_freelist, edge_supply);
    if (!newe) {
        fprintf (stderr, "CCbmat_edge_alloc failed\n");
        return (CCbmat_edge *) NULL;
    }
    newe->ends[0] = p1;
    newe->ends[1] = p2;
    newe->surfends[0] = p1->surf;
    newe->surfends[1] = p2->surf;
    newe->x = 0;
    newe->count = 0;
    newe->status.needschange = 0;
    newe->istooth[0] = newe->istooth[1] = 0;
    newe->slack = 0;
    newe->heap.type = CCbmat_NOTHEAPED;

    newe->ptrs[0].next = p1->surf->edgelist;
    p1->surf->edgelist = &(newe->ptrs[0]);
    newe->ptrs[1].next = p2->surf->edgelist;
    p2->surf->edgelist = &(newe->ptrs[1]);

    return newe;
}

static void contract_edge (CCbmat_edge *old, CCbmat_edge *new, CCbmat_node *pdel,
        CCbmat_node *p)
{
    CCbmat_edgeptr *pe, **ppe;
    int oldend, newend;

    /* contract old.  old connects pdel to p.  pdel has degree 2.
       new is the other edge at pdel */

    assert (p != pdel);

    oldend = (old->surfends[0] == p) ? 0 : 1;
    assert (old->surfends[oldend] == p);
    assert (old->surfends[1-oldend] == pdel);

    newend = (new->surfends[0] == pdel) ? 0 : 1;
    assert (new->surfends[newend] == pdel);
    assert (new->surfends[1-newend] != pdel);

    assert (newend == 1);
    assert (oldend == 0);

    for (ppe = &(p->edgelist), pe = *ppe; pe; ppe = &(pe->next), pe = *ppe) {
        if (pe == &(old->ptrs[oldend])) {
            new->ptrs[newend].next = pe->next;
            *ppe = &(new->ptrs[newend]);
            new->ends[newend] = old->ends[oldend];
            new->surfends[newend] = p;
            new->istooth[1-newend] = old->istooth[1-oldend];
            return;
        }
    }
    fprintf (stderr, "Help!!! I shouldn't be here\n");
    abort ();
}

CCbmat_node *lca (CCbmat_node *p1, CCbmat_node *p2)
{
    CCbmat_node *p;

    if (p1->surf != p2->surf) return (CCbmat_node *) NULL;

    for (p = p1; p; p = p->nest.parent) {
        p->status.inpath = CCbmat_TRUE;
    }
    while (p2 && !p2->status.inpath) {
        p2 = p2->nest.parent;
    }
    for (p = p1; p; p = p->nest.parent) {
        p->status.inpath = CCbmat_FALSE;
    }
    return p2;
}

#if DEBUG
static int sum_to_level (CCbmat_node *p, CCbmat_node *ptop)
{
    int val = 0;

    while (p && p != ptop) {
        val += p->ymat;
        p = p->nest.parent;
    }
    return val;
}
#endif

static int sum_to_root (CCbmat_node *p)
{
    if (p) {
        p->sum_to_root = p->ymat + sum_to_root (p->nest.parent);
        return p->sum_to_root;
    } else {
        return 0;
    }
}

int CCbmat_repair_world (CCbmat_graph *G, CCbmat_heaps *H, CCbmat_node *p1,
       CCbmat_node *p2, int weight, CCbmat_edge **edge_freelist,
       CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply)
{
    CCbmat_edge *newe1, *newe2;
    CCbmat_node *newp;
    int rval;
    int fix = 0;
    CCbmat_node *plca;
    int lcaval;
    int slack;

    plca = lca (p1, p2);
    slack = weight - sum_to_root (p1) - sum_to_root (p2);
    lcaval = plca ? plca->sum_to_root : 0;
    /* the sum_to_root's on the previous line will have made */
    /* plca->sum_to_root accurate */
    assert (lcaval == sum_to_root (plca));
    slack += 2 * lcaval;
    if (slack >= 0) {
        putchar ('-'); fflush (stdout);
        return 0;
    }
    putchar ('+'); fflush (stdout);

    newp = CCbmat_node_alloc (pseudo_freelist);
    newp->deficiency = 2;

    newe1 = build_edge(p1, newp, edge_freelist, edge_supply);
    newe2 = build_edge(p2, newp, edge_freelist, edge_supply);
    /* Note: newe1 or newe2 could be NULL if out of memory */
    if (slack & 0x1) {
        newe2->slack = 1;
        fix = 1;
        assert (CCbmat_slowsurf(p1) != CCbmat_slowsurf(p2));
    }

    CCbmat_augment2 (H, newp, lcaval + (((-slack)+fix)/2), pseudo_freelist);
    rval = (newp->ymat < lcaval + (((-slack)+fix)/2));

    while (p1->surf == p2->surf && p1->surf->ymat == 0) {
        putchar ('!'); fflush (stdout);
        CCbmat_zeroblowup (p1->surf, pseudo_freelist);
    }

    if (p1->surf == p2->surf) {
        printf ("!!! Bill was right !!!\n");
        printf ("p1 and p2 are both in same pseudonode\n");
        printf ("and it has weight %d\n",p1->surf->ymat);
    }
    assert (p1->surf != p2->surf && p1->surf == CCbmat_slowsurf(p1) &&
                                    p2->surf == CCbmat_slowsurf(p2));    
    if (newp->nest.parent) {
        printf ("!!! Bill was right !!!\n");
        printf ("newp is in pseudonode\n");
        printf ("and it has weight %d\n",newp->surf->ymat);
    }
    assert (!newp->nest.parent);
    assert (newp->surf == newp);
    assert (newe1->slack == 0);
    assert (newe2->slack == 0);

    assert (newe1->ends[0] == p1);
    assert (newe1->ends[1] == newp);
    assert (newe2->ends[0] == p2);
    assert (newe2->ends[1] == newp);

    contract_edge (newe2, newe1, newp->surf, p2->surf);

    newe1->slack = slack + 2 * (newp->ymat - lcaval) - fix;
    newe1->weight = weight;

    newe1->next = G->newedges;
    G->newedges = newe1;

    assert ((rval && newe1->x && newe2->x && newe1->slack < 0) ||
        (!rval&&!newe1->x &&!newe2->x && newe1->slack == 0));
    assert (newe1->slack == newe1->weight - p1->ymat - p2->ymat
     - (newe1->istooth[1] ? -1 : 1)*(sum_to_level(p1->nest.parent,lca(p1, p2)))
     - (newe1->istooth[0] ? -1 : 1)*(sum_to_level(p2->nest.parent,lca(p1,p2))));

    CCbmat_node_free (newp, pseudo_freelist);
    CCbmat_edge_free (newe2, edge_freelist);

    return 1;
}

int CCbmat_priceedges (char *f, CCbmat_graph *G, CCbmat_heaps *H,
        CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply)
{
    FILE *fin;
    int n, m, i, x, y, w;
    int nadded, rval = 0;
    CCbmat_node *surfacelist = (CCbmat_node *) NULL;
    CCbmat_anc *anc_freelist = (CCbmat_anc *) NULL;
    CCbmat_anc_ptr *anc_supply = (CCbmat_anc_ptr *) NULL, *a, *anext;
    CCbmat_edge *edge_freelist = (CCbmat_edge *) NULL;

    printf ("\nStarting to scan edge file for negative edges\n");
    fflush (stdout);

    if ((fin = fopen (f, "r")) == NULL) {
        perror (f);
        fprintf (stderr, "Unable to open pricing file, aborting pricing.\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCbmat_ancest_init (G, &surfacelist, &anc_freelist, &anc_supply);
    CCcheck_rval (rval, "CCbmat_ancest_init failed");

top:
    CCbmat_get_sums (G);

    fscanf (fin, "%d %d",&n, &m);
    if (n != G->ncount) {
        fprintf (stderr, "Pricing file contains wrong number of nodes\n");
        rval = 1; goto CLEANUP;
    }

    nadded = 0;
    for (i=0; i<m; i++) {
        fscanf (fin, "%d %d %d",&x,&y,&w);
        if (!findedge(G->nodelist+x, G->nodelist+y)) {
            w = w +w;
            if (w-G->nodelist[x].sum_to_root-G->nodelist[y].sum_to_root < 0) {
                nadded += CCbmat_ancest_checkout(G, H, G->nodelist+x,
                   G->nodelist+y, w, &surfacelist, &anc_freelist,
                   &edge_freelist, pseudo_freelist, &anc_supply, edge_supply);
            }
        }
        if (i % 1000 == 99) {putchar ('%'); fflush (stdout); }
    }
    nadded += CCbmat_ancest_flush (G, H, &surfacelist, &anc_freelist,
                  &edge_freelist, pseudo_freelist, &anc_supply, edge_supply);
    printf ("pass complete, %d edges added\n",nadded);
    fflush (stdout);
    if (nadded) {
        rewind (fin);
        goto top;
    }

CLEANUP:
    /* free the CCbmat_anc structs in the anc_supply list and free */
    /* the anc_ptr's in the list itself                           */
    a = anc_supply;
    while (a) {
        anext = a->next;
        CC_IFFREE (a->this, CCbmat_anc);
        CC_FREE (a, CCbmat_anc_ptr);
        a = anext;
    }
    return rval;
}

int CCbmat_pricedats (CCdatagroup *dat, CCbmat_graph *G, CCbmat_heaps *H, 
        CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply)
{
    int i, rval = 0, inorm;

    printf ("\n"); fflush (stdout);
    CCutil_dat_getnorm (dat, &inorm);
    if ((inorm & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        fprintf (stderr, "dat file does not support kd-tree norms\n");
        rval = 1; goto CLEANUP;
    }
    if ((inorm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "code only set up for 2-dimensional norms\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < G->ncount; i++) {
/*
        int ix, iy;
        ix = (int) dat->x[i];
        iy = (int) dat->y[i];
        if (dat->x[i] != (double) ix || dat->y[i] != (double) iy) {
            fprintf (stderr, "code only set up for integer x,y coordinates\n");
            rval = 1; goto CLEANUP;
        }
*/
        G->nodelist[i].co[0] = dat->x[i]; 
        G->nodelist[i].co[1] = dat->y[i];
    }
      
    CCbmat_get_sums (G);

    printf ("Starting to scan dat file for negative edges\n");

    for (i = 0; i < G->ncount; i++) {
        G->nodelist[i].status.marked = 1;
    }

    rval = price_scan (G, dat, H, pseudo_freelist, edge_supply);
    CCcheck_rval (rval, "price_scan failed");

CLEANUP:
    return rval;
}

void CCbmat_get_sums (CCbmat_graph *G)
{
    CCbmat_node *n, *nend = G->nodelist + G->ncount;

    for (n = G->nodelist; n < nend; n++)
        get_sum_work (n);
    for (n = G->nodelist; n < nend; n++)
        clear_sum_work (n);
}

static int get_sum_work (CCbmat_node *n)
{
    if (!n) return 0;
    else if (n->status.inpath) return n->sum_to_root;
    else {
        n->status.inpath = CCbmat_TRUE;
        n->sum_to_root = n->ymat + get_sum_work (n->nest.parent);
        return n->sum_to_root;
    }
}

static void clear_sum_work (CCbmat_node *n)
{
    if (n && n->status.inpath) {
        n->status.inpath = CCbmat_FALSE;
        clear_sum_work (n->nest.parent);
    }
}

static void initlist (CCbmat_graph *G, CCbmat_node *head, CCbmat_node *tail,
        CCbmat_node *head2, CCbmat_node *tail2)
{
    CCbmat_node *n;
    CCbmat_node *nodeend = G->nodelist + G->ncount;
    CCbmat_node *p;
    int bound;

    head->sort.order = -CCbmat_MAXIWEIGHT;
    tail->sort.order = CCbmat_MAXIWEIGHT;
    head->sort.next = tail;
    head->sort.prev = 0;
    tail->sort.next = 0;
    tail->sort.prev = &(head->sort.next);
    head2->sort.order = CCbmat_MAXIWEIGHT;
    tail2->sort.order = -CCbmat_MAXIWEIGHT;
    head2->sort.next = tail2;
    head2->sort.prev = 0;
    tail2->sort.next = 0;
    tail2->sort.prev = &(head2->sort.next);

    for (n = nodeend-1; n >= G->nodelist; n--) {
        bound = 2 * n->co[0] - n->sum_to_root;
        for (p = head->sort.next; p->sort.order < bound; p = p->sort.next)
        ;
        n->sort.order = bound;
        n->sort.next = p;
        n->sort.prev = p->sort.prev;
        *(n->sort.prev) = n;
        p->sort.prev = &n->sort.next;
    }
}

static int price_scan (CCbmat_graph *G, CCdatagroup *D, CCbmat_heaps *H,
        CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply)
{
    int w, rval = 0;
    double bound;
    int added, nodeschecked, edgeschecked, totaladded = 0;
    CCbmat_node *n1, *n2;
    CCbmat_node *nodeend = G->nodelist + G->ncount;
    CCbmat_node high_fakehead, high_faketail, low_fakehead, low_faketail;
    CCbmat_node *surfacelist = (CCbmat_node *) NULL;
    CCbmat_anc *anc_freelist = (CCbmat_anc *) NULL;
    CCbmat_edge *edge_freelist = (CCbmat_edge *) NULL;
    CCbmat_anc_ptr *anc_supply = (CCbmat_anc_ptr *) NULL, *a, *anext;

    initlist (G, &high_fakehead, &high_faketail, &low_fakehead, &low_faketail);
    rval = CCbmat_ancest_init (G, &surfacelist, &anc_freelist, &anc_supply);
    CCcheck_rval (rval, "CCbmat_ancest_init failed");

    do {
        added = 0;
        nodeschecked = 0;
        edgeschecked = 0;
        for (n1 = G->nodelist; n1 != nodeend; n1++) {
            *(n1->sort.prev) = n1->sort.next;
            n1->sort.next->sort.prev = n1->sort.prev;
            if (n1->status.marked) {
                nodeschecked++;
                n1->status.marked = 0;
                n1->sum_to_root = sum_to_root (n1);
                bound = 2*n1->co[0] + n1->sum_to_root;
                for (n2 = high_fakehead.sort.next; n2->sort.order <= bound &&
                          !n1->status.marked; n2 = n2->sort.next) {
                    edgeschecked++;
                    if (!n2->status.marked) {
                        w = 2*CCutil_dat_edgelen (n1-G->nodelist,
                                                  n2-G->nodelist, D);
                        if (w - n1->sum_to_root - n2->sum_to_root < 0) {
                            added += CCbmat_ancest_checkout (G, H, n1, n2, w,
                                &surfacelist, &anc_freelist, &edge_freelist,
                                pseudo_freelist, &anc_supply, edge_supply);
                        }
                    }
                }
                bound = 2*n1->co[0] - n1->sum_to_root;
                for (n2 = low_fakehead.sort.next; n2->sort.order >= bound &&
                         !n1->status.marked; n2 = n2->sort.next) {
                    edgeschecked++;
                    if (!n2->status.marked) {
                         w = 2*CCutil_dat_edgelen (n1-G->nodelist,
                                                  n2-G->nodelist, D);
                        if (w - n1->sum_to_root - n2->sum_to_root < 0) {
                            added += CCbmat_ancest_checkout (G, H, n1, n2, w,
                               &surfacelist, &anc_freelist, &edge_freelist,
                               pseudo_freelist, &anc_supply, edge_supply);
                        }
                    }
                }
            }
            bound = 2*n1->co[0] + n1->sum_to_root;
            for (n2 = low_fakehead.sort.next; n2->sort.order > bound;
                                              n2 = n2->sort.next)
            ;
            n1->sort.order = bound;
            n1->sort.next = n2;
            n1->sort.prev = n2->sort.prev;
            *(n1->sort.prev) = n1;
            n2->sort.prev = &n1->sort.next;
            if ((n1 - G->nodelist) % 100 == 99) {
                putchar ('%'); fflush (stdout);
            }
        }
        added += CCbmat_ancest_flush (G, H, &surfacelist, &anc_freelist,
              &edge_freelist, pseudo_freelist, &anc_supply, edge_supply);
        totaladded += added;
        printf ("\nForward pass complete, %d nodes checked, %d edges checked\n",
                       nodeschecked,edgeschecked);
        printf ("    %d edges added, total %d edges added\n",added,totaladded);
        fflush (stdout);
        if (added == 0) break;
        added = 0;
        nodeschecked = 0;
        edgeschecked = 0;
        for (n1 = nodeend-1; n1 >= G->nodelist; n1--) {
            *(n1->sort.prev) = n1->sort.next;
            n1->sort.next->sort.prev = n1->sort.prev;
            if (n1->status.marked) {
                nodeschecked++;
                n1->status.marked = 0;
                n1->sum_to_root = sum_to_root (n1);
                bound = 2*n1->co[0] + n1->sum_to_root;
                for (n2 = high_fakehead.sort.next; n2->sort.order <= bound
                         && !n1->status.marked; n2 = n2->sort.next) {
                    edgeschecked++;
                    if (!n2->status.marked) {
                        w = 2*CCutil_dat_edgelen (n1-G->nodelist,
                                                  n2-G->nodelist, D);
                        if (w - n1->sum_to_root - n2->sum_to_root < 0) {
                            added += CCbmat_ancest_checkout (G,H,n1,n2,w,
                             &surfacelist, &anc_freelist, &edge_freelist,
                             pseudo_freelist, &anc_supply, edge_supply);
                        }
                    }
                }
                bound = 2*n1->co[0] - n1->sum_to_root;
                for (n2 = low_fakehead.sort.next; n2->sort.order >= bound
                         && !n1->status.marked; n2 = n2->sort.next) {
                    edgeschecked++;
                    if (!n2->status.marked) {
                        w = 2*CCutil_dat_edgelen (n1-G->nodelist,
                                                  n2-G->nodelist, D);
                        if (w - n1->sum_to_root - n2->sum_to_root < 0) {
                            added += CCbmat_ancest_checkout (G,H,n1,n2,w,
                             &surfacelist, &anc_freelist, &edge_freelist,
                             pseudo_freelist, &anc_supply, edge_supply);
                        }
                    }
                }
            }
            bound = 2*n1->co[0] - n1->sum_to_root;
            for (n2 = high_fakehead.sort.next; n2->sort.order < bound;
                                          n2 = n2->sort.next)
            ;
            n1->sort.order = bound;
            n1->sort.next = n2;
            n1->sort.prev = n2->sort.prev;
            *(n1->sort.prev) = n1;
            n2->sort.prev = &n1->sort.next;
            if ((nodeend - n1) % 100 == 0) {
                putchar ('%'); fflush (stdout);
            }
        }
        added += CCbmat_ancest_flush (G, H, &surfacelist, &anc_freelist,
                    &edge_freelist, pseudo_freelist, &anc_supply, edge_supply);
        totaladded += added;
        printf ("\nBackward pass complete, %d nodes checked, %d edges checked\n",
                                    nodeschecked,edgeschecked);
        printf ("    %d edges added, total %d edges added\n",added,totaladded);
        fflush (stdout);
    } while (added);

CLEANUP:
    a = anc_supply;  /* free anc_supply */
    while (a) {
        anext = a->next;
        CC_IFFREE (a->this, CCbmat_anc);
        CC_FREE (a, CCbmat_anc_ptr);
        a = anext;
    }
    return rval;
}

