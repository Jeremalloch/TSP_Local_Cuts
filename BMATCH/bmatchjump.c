#include <stdio.h>
#include <stdlib.h>
#include "bmatch.h"

static void initmat (CCbmat_graph *G, int bdegree, int *degrees);
static int chkmat (CCbmat_graph *G, int bdegree, int *degrees);
static int bfraction (CCbmat_graph *G);
static int augment (CCbmat_graph *G, CCbmat_node *n, int PLUS, int MINUS);
static int expand (CCbmat_graph *G, CCbmat_node *n, int PLUS, int MINUS);
static CCbmat_node *dualchange (CCbmat_node *n, int PLUS);
static void minalpha (CCbmat_node *n, CCbmat_node **new, int *alpha, int PLUS);
static void subalpha (CCbmat_node *n, int alpha, int PLUS);
static void augmentpath (CCbmat_node *n);
static void augmentpath1 (CCbmat_node *n);
static void augmentplushalf (CCbmat_node *n, CCbmat_edge *e);
static void augmentminushalf (CCbmat_node *n, CCbmat_edge *e);
static void flipcycle (CCbmat_edge *e, int v);
static void augmentpair (CCbmat_node *n, CCbmat_node *n1, CCbmat_edge *e);
static void setflag (CCbmat_node *n, int v);
static CCbmat_node *findflag (CCbmat_node *n);
static CCbmat_node *findhole (CCbmat_node *n, CCbmat_node *n2);
static void flipseg (CCbmat_node *n, CCbmat_node *n2);
static void stringup (CCbmat_node *n, CCbmat_node *n1, CCbmat_node *n2,
    CCbmat_edge *e);

int CCbmat_jump (CCbmat_graph *G, int bdegree, int *degrees)
{
    int v, rval = 0;

    initmat (G, bdegree, degrees);
    rval = bfraction (G);
    CCcheck_rval (rval, "bfraction failed");
    v = chkmat (G, bdegree, degrees);
    printf ("Fractional Matching: %d\n", v / 2); fflush (stdout);

CLEANUP:
    return rval;
}

static void initmat (CCbmat_graph *G, int bdegree, int *degrees)
{
    int i;

    for (i=0; i<G->ncount; i++) {
        G->nodelist[i].ymat = 0;   
        if (degrees) {
            G->nodelist[i].deficiency = degrees[i];
        } else {
            G->nodelist[i].deficiency = bdegree;
        }
        G->nodelist[i].parentedge = NULL;
        G->nodelist[i].label = 0;
    }
    
    for (i=0; i<G->ecount; i++) {
        G->edgelist[i].zmat = 0;
        G->edgelist[i].x = CCbmat_ZERO;
        G->edgelist[i].next = NULL;
    }
}

static int chkmat (CCbmat_graph *G, int bdegree, int *degrees)
{
    int i, k, n;
    int val = 0, dualval = 0;
    CCbmat_edgeptr *pe;

    if (degrees) {
        printf ("Checking fractional matching with mixed degrees!\n");
        fflush (stdout);
    }

    for (i=0; i<G->ncount; i++) {
        k = (degrees ? degrees[i]: bdegree);
        n = 0;
        for (pe = G->nodelist[i].edgelist; pe; pe = pe->next) {
            n += CCbmat_THISEDGE(pe)->x;
        }
        if (n != 2*k ) {
            printf ("Not a matching.  Node %d has 2 %d-degree %d.\n", i, k, n);
            exit (1);
        }
        dualval += k * G->nodelist[i].ymat;
    }
    for (i=0; i<G->ecount; i++) {
        switch (G->edgelist[i].x) {
        case CCbmat_TWO:
            val += G->edgelist[i].weight;
            break;
        case CCbmat_ONE:
            val += (G->edgelist[i].weight / 2);
            break;
        }
        dualval -= G->edgelist[i].zmat;
    }

    if (val <= dualval-1 || val >= dualval+1) {
        printf ("The primal and dual objective values differ.\n");
        printf ("primal = %d  dual = %d\n", val, dualval);
        exit (1);
    }
    return val;
}

static int bfraction (CCbmat_graph *G)
{
    int i, rval = 0;
    int PLUS = 1, MINUS = 2;

    for (i=0; i<G->ncount; i++) {
        while (G->nodelist[i].deficiency) {
            rval = augment (G, &G->nodelist[i], PLUS, MINUS);
            CCcheck_rval (rval, "augment failed");
            PLUS += 2;
            MINUS += 2;
        }
    }

CLEANUP:
    return rval;
}

static int augment (CCbmat_graph *G, CCbmat_node *n, int PLUS, int MINUS)
{
    CCbmat_node *auglist;

    n->label = PLUS;
    n->parentedge = NULL;
    n->deficiency--;
    if (expand (G, n, PLUS, MINUS)) return 0;
    while ((auglist = dualchange (n, PLUS)) != NULL) {
        while (auglist) {
            if (expand (G, auglist, PLUS, MINUS)) return 0;
            auglist = auglist -> next;
        }
    }
    return 1;
}

static int expand (CCbmat_graph *G, CCbmat_node *n, int PLUS, int MINUS)
{
    CCbmat_node *n1;
    CCbmat_edge *e;
    CCbmat_edgeptr *pe;

    if (n->label < PLUS) {
        fprintf (stderr,"Error, expanding an unlabeled node\n");
        exit (1);
    } else if (n->label == PLUS) {
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            if (e->x == CCbmat_ONE) {
                augmentplushalf (n, e);
                return 1;
            }
        }
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            switch (e->x) {
            case CCbmat_ZERO:
                n1 = e->ends[0];
                if (n1 == n) n1 = e->ends[1];
                if (n1 -> ymat + n -> ymat - e->zmat  != e->weight) break;
                if (n1->label < PLUS) {       /* n1 has no label */
                    n1->label = MINUS;
                    n1->parentedge = e;
                    if (n1->deficiency) {
                        augmentpath (n1);
                        return 1;
                    }
                    if (expand (G, n1, PLUS, MINUS)) return 1;
                } else if (n1->label == PLUS) {
                    augmentpair (n, n1, e);
                    return 1;
                }
                break;
            case CCbmat_TWO:
                break;
            }
        }
    } else {        /* MINUS */
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            if (e->x  == CCbmat_ONE) {
                augmentminushalf (n, e);
                return 1;
            }
        }
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            switch (e->x) {
            case CCbmat_ZERO: break;
            case CCbmat_TWO:
                if (e->zmat != 0) break;
                n1 = e->ends[0];
                if (n1 == n) n1 = e->ends[1];
                if (n1 -> label < PLUS) {          /* n1 has no label */
                    n1 -> label = PLUS;
                    n1 -> parentedge = e;
                    if (expand (G, n1, PLUS, MINUS)) return 1;
                } else if (n1->label == MINUS) {
                    augmentpair (n, n1, e);
                    return 1;
                }
                break;
            }
        }
    }
    return 0;
}

static CCbmat_node *dualchange (CCbmat_node *n, int PLUS)
{
    CCbmat_node *new = NULL;
    int alpha = CCbmat_MAXIWEIGHT;

    minalpha (n, &new, &alpha, PLUS);
    if (alpha == CCbmat_MAXIWEIGHT) {
        printf ("Whoops, dual change required, but no candidate edges\n");
        return NULL;
    }
    alpha /= 2;
    subalpha (n, alpha, PLUS);
    return new;
}

static void minalpha (CCbmat_node *n, CCbmat_node **new, int *alpha, int PLUS)
{
    int minv = CCbmat_MAXIWEIGHT;
    int thisv;
    CCbmat_node *n1;
    CCbmat_edge *e;
    CCbmat_edgeptr *pe;

    if (n->label < PLUS) {
        fprintf (stderr, "minalpha at unlabeled node\n");
        exit (1);
    } else if (n->label == PLUS) {
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            if (e->x == CCbmat_ZERO) {
                n1 = e->ends[0];
                if (n1 == n)
                n1 = e->ends[1];
                if (n1->label < PLUS) {      /* n1 is unlabeled */
                    thisv = e->weight - n->ymat - n1->ymat;
                    thisv += thisv;
                    if (thisv < minv) minv = thisv;
                } else if (n1 -> label == PLUS) {
                    thisv = e->weight - n->ymat - n1->ymat;
                    if (thisv < minv) minv = thisv;
                } else {                     /* n1 has a minus label */
                    if (n1 -> parentedge == e) minalpha (n1, new, alpha, PLUS);
                }
            }
        }
    } else {                       /* MINUS case */
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            if (e->x == CCbmat_TWO) {
                n1 = e->ends[0];
                if (n1 == n) n1 = e->ends[1];
                if (n1->label < PLUS) {
                    thisv = e->zmat + e->zmat;
                    if (thisv < minv) minv = thisv;
                } else if (n1->label == PLUS) {
                    if (n1 -> parentedge == e) minalpha (n1, new, alpha, PLUS);
                } else {           /* n1 has a MINUS label */
                    thisv = e->zmat;
                    if (thisv < minv) minv = thisv;
                }
            }
        }
    }

    if (minv < *alpha) {
        *new = n;
        n->next = NULL;
        *alpha = minv;
    } else if (minv == *alpha) {
        n->next = *new;
        *new = n;
    }
}

static void subalpha (CCbmat_node *n, int alpha, int PLUS)
{
    CCbmat_edge *e;
    CCbmat_edgeptr *pe;
    CCbmat_node *n1;

    if (n->label < PLUS) {
        fprintf (stderr, "subalpha at unlabeled node\n");
        exit (1);
    } else if (n->label == PLUS) {
        n->ymat += alpha;
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            if (e->x == CCbmat_TWO) {
                e->zmat += alpha;
            } else if (e->x == CCbmat_ZERO) {
                n1 = e->ends[0];
                if (n1 == n) n1 = e->ends[1];
                if (n1->parentedge == e && n1->label >= PLUS)
                subalpha (n1, alpha, PLUS);
            }
        }
    } else {        /* MINUS */
        n->ymat -= alpha;
        for (pe = n->edgelist; pe; pe = pe->next) {
            e = CCbmat_THISEDGE(pe);
            if (e->x == CCbmat_TWO) {
                e->zmat -= alpha;
                n1 = e -> ends[0];
                if (n1 == n) n1 = e->ends[1];
                if (n1 -> parentedge == e && n1->label >= PLUS) {
                    subalpha (n1, alpha, PLUS);
                }
            }
        }
    }
}

static void augmentpath (CCbmat_node *n)
{
    n->deficiency--;
    augmentpath1 (n);
}

static void augmentpath1 (CCbmat_node *n)
{
    CCbmat_node *n1;
    CCbmat_edge *e;

    while ((e = n->parentedge) != NULL) {
        e->x = 2 - e->x;
        n1 = e->ends[0];
        if (n1 == n) n1 = e->ends[1];
        n = n1;
    }
}

static void augmentplushalf (CCbmat_node *n, CCbmat_edge *e)
{
    CCbmat_edge *e1 = e->next;

    if (e1->ends[0] == n || e1->ends[1] == n) e = e1;
    flipcycle (e, CCbmat_TWO);
    augmentpath1 (n);
}

static void augmentminushalf (CCbmat_node *n, CCbmat_edge *e)
{
    CCbmat_edge *e1 = e->next;

    if (e1->ends[0] == n || e1->ends[1] == n) e = e1;
    flipcycle (e, CCbmat_ZERO);
    augmentpath1 (n);
}

static void flipcycle (CCbmat_edge *e, int v)
{
    CCbmat_edge *e1 = e;
    CCbmat_edge *e2;

    do {
        e1->x = v;
        v = 2 - v;
        e2 = e1 -> next;
        e1 -> next = NULL;
        e1 = e2;
    } while (e1 != e);
}

static void augmentpair (CCbmat_node *n, CCbmat_node *n1, CCbmat_edge *e)
{
    CCbmat_node *n2, *n3;

    setflag (n, CCbmat_FALSE);
    setflag (n1, CCbmat_TRUE);
    n2 = findflag (n);
    if ((n3 = findhole (n, n2)) != NULL) {
        n3 -> deficiency--;
        flipseg (n, n3);
        e->x = 2 - e->x;
        augmentpath1 (n1);
        return;
    }
    if ((n3 = findhole (n1, n2)) != NULL) {
        n3 -> deficiency--;
        flipseg (n1, n3);
        e->x = 2 - e->x;
        augmentpath1 (n);
        return;
    }
    stringup (n, n1, n2, e);
    augmentpath1 (n2);
}

static void setflag (CCbmat_node *n, int v)
{
    CCbmat_edge *e;
    CCbmat_node *n1;

    n -> status.inpath = v;
    while ((e = n->parentedge) != NULL) {
        n1 = e->ends[0];
        if (n1 == n)
        n1 = e->ends[1];
        n = n1;
        n->status.inpath = v;
    }
}

static CCbmat_node *findflag (CCbmat_node *n)
{
    CCbmat_edge *e;
    CCbmat_node *n1;

    while (!n->status.inpath) {
        if ((e = n->parentedge) == NULL) {
            fprintf (stderr,"Find, ran out of parents\n");
            exit (1);
        }
        n1 = e -> ends[0];
        if (n1 == n) n1 = e->ends[1];
        n = n1;
    }
    return n;
}

static CCbmat_node *findhole (CCbmat_node *n, CCbmat_node *n2)
{
    CCbmat_edge *e;
    CCbmat_node *n1;

    while (n != n2) {
        if (n->deficiency) return n;
        e = n->parentedge;
        n1 = e->ends[0];
        if (n1 == n) n1 = e->ends[1];
        n = n1;
    }
    return n2->deficiency ? n2 : NULL;
}

static void flipseg (CCbmat_node *n, CCbmat_node *n2)
{
    CCbmat_edge *e;
    CCbmat_node *n1;

    while (n != n2) {
        e = n->parentedge;
        e->x = 2 - e->x;
        n1 = e -> ends[0];
        if (n1 == n) n1 = e->ends[1];
        n = n1;
    }
}

static void stringup (CCbmat_node *n, CCbmat_node *n1, CCbmat_node *n2,
        CCbmat_edge *e)
{
    CCbmat_node *n3;
    CCbmat_edge *preve, *thise, *savee;

    preve = e;
    while (n != n2) {
        thise = n -> parentedge;
        thise -> x = CCbmat_ONE;
        preve -> next = thise;
        preve = thise;
        n3 = thise -> ends[0];
        if (n3 == n) n3 = thise -> ends[1];
        n = n3;
    }
    savee = preve;
    preve = e;
    while (n1 != n2) {
        thise = n1->parentedge;
        thise -> x = CCbmat_ONE;
        thise -> next = preve;
        preve = thise;
        n3 = thise -> ends[0];
        if (n3 == n1) n3 = thise -> ends[1];
        n1 = n3;
    }
    savee->next = preve;
    e->x = CCbmat_ONE;
}
