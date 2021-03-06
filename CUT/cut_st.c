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
/*         MINIMUM S-T CUTS IN DIRECTED AND UNDIRECTED GRAPHS               */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: July 25, 1996                                                     */
/*                                                                          */
/*  An implementation of the Push-Relabel Flow Algorithm described in       */
/*   A. Goldberg and R. Tarjan, "A new approach to the maximum-flow         */
/*   problem",  Journal of the ACM 35 (1988) 921-940.                       */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCcut_mincut_st (int ncount, int ecount, int *elist, double *ecap,  */
/*      int s, int t, double *value, int **cut, int *cutcount)              */
/*    COMPUTES the min st-cut in a directed or undirected graph.            */
/*      -ncount is the number of nodes in the graph.                        */
/*      -ecount is the number of directed (undirected) edges.               */
/*      -elist gives the edges in node node format (interpreted as          */
/*           tail head when compiled for directed graphs).                  */
/*      -ecap gives the capacities of the edges.                            */
/*      -s is the name of the source node.                                  */
/*      -t is the name of the sink node.                                    */
/*      -value returns the capacity of the minimum cut.                     */
/*      -cut (if not NULL) returns a list of nodes in a a minimum cut (it   */
/*       returns the side that contains t); it will be allocated to an      */
/*       array of the appropriate size.                                     */
/*      -cutcount returns the number of nodes in the listed cut, if cut     */
/*       is not NULL (if cut is NULL, then cutcount can be NULL).           */
/*                                                                          */
/*    NOTES:                                                                */
/*      Returns 0 if it worked and 1 otherwise (for example, when one       */
/*      of the mallocs failed). The nodes in the graph should be named      */
/*      0 through #nodes - 1.                                               */
/*                                                                          */
/*      Define UNDIRECTED_GRAPH to compile the code for undirected          */
/*      graphs. (This appears to be the way to go for tsp instances.)       */
/*                                                                          */
/*      Two node selection rules are implemented: queue and highest         */
/*      label. One of QUEUE_PRF and HIGHEST_LABEL_PRF must be defined       */
/*      (but not both).                                                     */
/*                                                                          */
/*      The code can carry out global relabelings via a backwards           */
/*      breadth-first-search from the sink. The frequency of the            */
/*      relabelings is controlled by the defined constant                   */
/*      GLOBAL_RELABEL_FREQUENCY. A relabling will occur after each         */
/*      #nodes * GLOBAL_RELABEL_FREQUENCY nodes have been processed.        */
/*      A resonable choice for the constant is 1.                           */
/*                                                                          */
/*      Defining USE_GAP turns on the gap heuristic of Derigs and           */
/*      Meyer for determing nodes that can be labeled to ncount.            */
/*      This can be used with either the queue or highest label             */
/*      variants.                                                           */
/*                                                                          */
/*      To use this code for maxflows, allow nodes with labels up to        */
/*      2*ncount to become active, or implement an algorithm to             */
/*      decompose the preflow to create a flow.                             */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"

#define UNDIRECTED_GRAPH

#define HIGHEST_LABEL_PRF
#undef  QUEUE_PRF

#define USE_GAP
#define GLOBAL_RELABEL_FREQUENCY 1.0

#define CUT_INFINITY (1<<30)
#define PRF_EPSILON 0.000000001
#define GOING_IN  0
#define GOING_OUT 1

#ifdef QUEUE_PRF
#define ADD_TO_ACTIVE(n) {                                                \
    if (!(n)->active) {                                                   \
        (n)->qnext = (node *) NULL;                                       \
        if (qtail) {                                                      \
            qtail->qnext = (n);                                           \
        } else {                                                          \
            qhead = (n);                                                  \
        }                                                                 \
        qtail = (n);                                                      \
        (n)->active = 1;                                                  \
    }                                                                     \
}
#endif

#ifdef HIGHEST_LABEL_PRF
#define ADD_TO_ACTIVE(n) {                                                \
    if (!(n)->active) {                                                   \
        (n)->highnext = high[(n)->flowlabel];                             \
        high[(n)->flowlabel] = (n);                                       \
        if (G->highest < (n)->flowlabel)                                  \
            G->highest = (n)->flowlabel;                                  \
        (n)->active = 1;                                                  \
    }                                                                     \
}
#endif

#ifdef UNDIRECTED_GRAPH
#define RELABEL_BODY(n)                                                   \
    for (rele = (n)->out; rele; rele = rele->outnext) {                   \
        if (rele->cap - rele->flow > PRF_EPSILON &&                       \
                (relt = rele->ends[1]->flowlabel) < relm)                 \
            relm = relt;                                                  \
    }                                                                     \
    for (rele = (n)->in; rele; rele = rele->innext) {                     \
        if (rele->cap + rele->flow > PRF_EPSILON &&                       \
                 (relt = rele->ends[0]->flowlabel) < relm)                \
            relm = relt;                                                  \
    }                                                                     \
    (n)->flowlabel = ++relm;
#else
#define RELABEL_BODY(n)                                                   \
    for (rele = (n)->out; rele; rele = rele->outnext) {                   \
        if (rele->cap - rele->flow > PRF_EPSILON &&                       \
                (relt = rele->ends[1]->flowlabel) < relm)                 \
            relm = relt;                                                  \
    }                                                                     \
    for (rele = (n)->in; rele; rele = rele->innext) {                     \
        if (rele->flow > PRF_EPSILON &&                                   \
                 (relt = rele->ends[0]->flowlabel) < relm)                \
            relm = relt;                                                  \
    }                                                                     \
    (n)->flowlabel = ++relm;
#endif

#ifdef USE_GAP
#define RELABEL(n) {                                                      \
    int relm = CUT_INFINITY;                                              \
    edge *rele;                                                           \
    int relt, relold = (n)->flowlabel;                                    \
                                                                          \
    RELABEL_BODY(n)                                                       \
                                                                          \
    if ((n)->levelprev) {                                                 \
        (n)->levelprev->levelnext = (n)->levelnext;                       \
    } else {                                                              \
        level[relold] = (n)->levelnext;                                   \
    }                                                                     \
    if ((n)->levelnext)                                                   \
        (n)->levelnext->levelprev = (n)->levelprev;                       \
                                                                          \
    if (relm < ncount) {                                                  \
        if (level[relm]) {                                                \
            level[relm]->levelprev = (n);                                 \
            (n)->levelnext = level[relm];                                 \
            (n)->levelprev = (node *) NULL;                               \
            level[relm] = (n);                                            \
        } else {                                                          \
            (n)->levelprev = (node *) NULL;                               \
            (n)->levelnext = (node *) NULL;                               \
            level[relm] = (n);                                            \
        }                                                                 \
        if (!level[relold]) {                                             \
            relold++;                                                     \
            while (level[relold]) {                                       \
                node *relno;                                              \
                for (relno = level[relold]; relno;                        \
                                            relno = relno->levelnext) {   \
                    relno->flowlabel = ncount;                            \
                }                                                         \
                level[relold] = (node *) NULL;                            \
                relold++;                                                 \
            }                                                             \
        }                                                                 \
    }                                                                     \
}
#else
#define RELABEL(n) {                                                      \
    int relm = CUT_INFINITY;                                              \
    edge *rele;                                                           \
    int relt;                                                             \
                                                                          \
    RELABEL_BODY(n)                                                       \
}
#endif /* USE_GAP */

typedef struct edge {
    struct node     *ends[2];
    struct edge     *innext;
    struct edge     *outnext;
    double          cap;
    double          flow;
} edge;

typedef struct node {
    struct  node    *qnext;
    struct  node    *tnext;
#ifdef USE_GAP
    struct  node    *levelnext;
    struct  node    *levelprev;
#endif
#ifdef HIGHEST_LABEL_PRF
    struct  node    *highnext;
#endif
    struct  edge    *in;
    struct  edge    *out;
    struct  edge    *incurrent;
    struct  edge    *outcurrent;
    double           excess;
    int              magiclabel;
    int              flowlabel;
    char             inout;
    char             active;
} node;

typedef struct graph {
    struct  node    *nodelist;
    struct  edge    *edgelist;
#ifdef USE_GAP
    struct  node   **level;
#endif
#ifdef HIGHEST_LABEL_PRF
    struct  node   **high;
    int              highest;
#endif
    int              nnodes;
    int              nedges;
    int              magicnum;
} graph;



static void
    setlabels (graph *G, node *s, node *t),
    backwards_bfs (node *s, int K, graph *G),
    init_graph (graph *G),
    free_graph (graph *G);
static int
    grab_the_cut (graph *G, node *n, int **cut, int *cutcount),
    buildgraph (graph *G, int ncount, int ecount, int *elist, double *gap);
static double
    flow (graph *G, node *s, node *t);


int CCcut_mincut_st (int ncount, int ecount, int *elist, double *ecap,
            int s, int t, double *value, int **cut, int *cutcount)
{
    int rval = 0;
    graph G;

    init_graph (&G);

    if (cut) {
        *cut = (int *) NULL;
        if (cutcount) {
            *cutcount = 0;
        } else {
            fprintf (stderr, "cut is specified but not cutcount\n");
            rval = 1; goto CLEANUP;
        }
    }

    rval = buildgraph (&G, ncount, ecount, elist, ecap);
    if (rval) {
        fprintf (stderr, "Buildgraph failed\n"); goto CLEANUP;
    }
    *value = flow (&G, G.nodelist + s, G.nodelist + t);
    if (cut) {
/*
        {
            int i;
            double tval = 0.0;

            for (i = 0; i < ecount; i++) {
                if (elist[2*i] == t || elist[2*i+1] == t) tval += ecap[i];
            }
            if (tval <= *value + PRF_EPSILON &&
                tval >= *value - PRF_EPSILON) {
                *cut = CC_SAFE_MALLOC (1, int);
                CCcheck_NULL (*cut, "out of memory for cut");
                (*cut)[0] = t;
                *cutcount = 1;
                printf ("BING!\n"); fflush (stdout);
            }
        }
*/

        rval = grab_the_cut (&G, G.nodelist + t, cut, cutcount);
        CCcheck_rval (rval, "grab_the_cut failed");

        {
            int i, joke[100000];
            double z = 0.0;

            for (i = 0; i < ncount; i++) joke[i] = 0;
            for (i = 0; i < *cutcount; i++) joke[(*cut)[i]] = 1;
            for (i = 0; i < ecount; i++) {
                if (joke[elist[2*i]] != joke[elist[2*i+1]]) {
                    z += ecap[i];
                }
            }
            if (z - *value > 0.001 || *value - z > 0.001) {
                printf ("mincut_st = %.4f, actually = %.4f\n", *value, z);
                fflush (stdout);
                printf ("Cut:\n");
                for (i = 0; i < *cutcount; i++) {
                    printf (" %2d", (*cut)[i]);
                    if (i % 10 == 9) printf ("\n");
                }
                printf ("\n");
                exit (1);
            }
        }
    }

CLEANUP:

    free_graph (&G);
    return rval;
}

static double flow (graph *G, node *s, node *t)
{
#ifdef QUEUE_PRF
    node *qhead = (node *) NULL;
    node *qtail = (node *) NULL;
#endif
    node *n;
    edge *e;
    int count, myrnd;
    int i;
    int ncount = G->nnodes;
    edge *edgelist = G->edgelist;
    node *nodelist = G->nodelist;
#ifdef USE_GAP
    node **level = G->level;
#endif
#ifdef HIGHEST_LABEL_PRF
    node **high = G->high;
#endif

/*
    printf ("Find cut separating %d and %d ...\n", s - nodelist, t - nodelist);
    fflush (stdout);
*/

    for (i = 0; i < ncount; i++) {
        nodelist[i].excess = 0.0;
        nodelist[i].active = 0;
#ifdef HIGHEST_LABEL_PRF
        high[i] = (node *) NULL;
#endif
    }
#ifdef HIGHEST_LABEL_PRF
    G->highest = 0;
#endif

    for (i = G->nedges - 1; i >= 0; i--)
        edgelist[i].flow = 0.0;

    t->active = 1;              /* a lie, which keeps s and t off the */
    s->active = 1;              /* active int                         */

    for (e = s->out; e; e = e->outnext) {
        if (e->cap > 0.0) {
            e->flow = e->cap;
            e->ends[1]->excess += e->cap;
#ifdef QUEUE_PRF
            ADD_TO_ACTIVE(e->ends[1]);
#endif
#ifdef HIGHEST_LABEL_PRF
            e->ends[1]->active = 1;
#endif
        }
    }
#ifdef UNDIRECTED_GRAPH
    for (e = s->in; e; e = e->innext) {
        if (e->cap > 0.0) {
            e->flow = -e->cap;
            e->ends[0]->excess += e->cap;
#ifdef QUEUE_PRF
            ADD_TO_ACTIVE(e->ends[0]);
#endif
#ifdef HIGHEST_LABEL_PRF
            e->ends[0]->active = 1;
#endif
        }
    }
#endif


    setlabels (G, s, t);
    count = 0;
    myrnd = (int) (GLOBAL_RELABEL_FREQUENCY * ncount);

#ifdef QUEUE_PRF
    while (qhead) {
        n = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = (node *) NULL;
        n->active = 0;
        if (n->flowlabel >= ncount)
            continue;
#endif

#ifdef HIGHEST_LABEL_PRF
    while (G->highest) {
        n = high[G->highest];
        n->active = 0;
        high[G->highest] = high[G->highest]->highnext;
        if (!high[G->highest]) {
            G->highest--;
            while (G->highest && (high[G->highest] == (node *) NULL))
                G->highest--;
        }
#endif

        if (count == myrnd) {
            setlabels (G, s, t);
            if (n->flowlabel >= ncount)
                continue;
            count = 0;
        } else
            count++;

        if (n->inout == GOING_IN)
            goto DO_ME_IN;

        if (n->outcurrent) {
            while (n->excess > 0.0) {
                e = n->outcurrent;
                { /* PUSH OUT */
                    double rf = e->cap - e->flow;
                    node *n1 = e->ends[1];
                    if (n->flowlabel == n1->flowlabel + 1 && rf > 0.0) {
                        if (n->excess <= rf) {
                            e->flow += n->excess;
                            n1->excess += n->excess;
                            n->excess = 0.0;
                            ADD_TO_ACTIVE(n1);
                        } else {
                            e->flow += rf;
                            n1->excess += rf;
                            n->excess -= rf;
                            ADD_TO_ACTIVE(n1);
                            n->outcurrent = e->outnext;
                            if (!n->outcurrent) {
                                n->outcurrent = n->out;
                                n->inout = GOING_IN;
                                goto DO_ME_IN;
                            }
                        }
                    } else {
                        n->outcurrent = e->outnext;
                        if (!n->outcurrent) {
                            n->outcurrent = n->out;
                            n->inout = GOING_IN;
                            goto DO_ME_IN;
                        }
                    }
                }
            }
        }
DO_ME_IN:
        if (n->incurrent) {
            while (n->excess > 0.0) {
                e = n->incurrent;
                { /* PUSH IN */
#ifdef UNDIRECTED_GRAPH
                    double rf = e->cap + e->flow;
#else
                    double rf = e->flow;
#endif
                    node *n1 = e->ends[0];
                    if (n->flowlabel == n1->flowlabel + 1 && rf > 0.0) {
                        if (n->excess <= rf) {
                            e->flow -= n->excess;
                            n1->excess += n->excess;
                            n->excess = 0.0;
                            ADD_TO_ACTIVE(n1);
                        } else {
                            e->flow -= rf;
                            n1->excess += rf;
                            n->excess -= rf;
                            ADD_TO_ACTIVE(n1);
                            n->incurrent = e->innext;
                            if (!n->incurrent) {
                                n->incurrent = n->in;
                                n->inout = GOING_OUT;
                                RELABEL(n);
                                break;
                            }
                        }
                    } else {
                        n->incurrent = e->innext;
                        if (!n->incurrent) {
                            n->incurrent = n->in;
                            n->inout = GOING_OUT;
                            RELABEL(n);
                            break;
                        }
                    }
                }
            }
        } else {
            /* n->in is NULL */
            n->inout = GOING_OUT;
            RELABEL(n);
        }
        if (n->excess > 0.0 && n->flowlabel < ncount) {
            ADD_TO_ACTIVE(n);
        }
    }

    return t->excess;
}

static void setlabels (graph *G, node *s, node *t)
{
    node *n;
    int ncount = G->nnodes;
    int num = ++(G->magicnum);
    int i;
    /* static int duke = 0; */

    t->flowlabel = 0;
    backwards_bfs (t, num, G);
    if (s->magiclabel == num) {
        printf ("Help - s should not get a label\n");
        s->flowlabel = ncount;
    }

    for (i = G->nnodes, n = G->nodelist; i;  i--, n++) {
        n->outcurrent = n->out;
        n->incurrent = n->in;
        n->inout = GOING_OUT;
        if (n->magiclabel != num) {
            n->flowlabel = ncount;
        }
    }
}

static void backwards_bfs (node *s, int K, graph *G)
{
    node *this, *next, *tail;
    edge *e;
    int  dist;
#ifdef USE_GAP
    node dummy;
    node **level = G->level;
#endif
#ifdef HIGHEST_LABEL_PRF
    node **high = G->high;
#endif

    s->magiclabel = K;
    next = s;
    s->tnext = (node *) NULL;
    dist = s->flowlabel;

#ifdef USE_GAP
    {
        int i;
        for (i = 0; level[i]; i++)
            level[i] = (node *) NULL;
        level[dist] = s;
        s->levelnext = (node *) NULL;
    }
#endif
#ifdef HIGHEST_LABEL_PRF
    {
        int i;
        for (i = 0; i <= G->highest; i++)
            high[i] = (node *) NULL;
        G->highest = 0;
    }
#endif

    do {
#ifdef USE_GAP
        level[dist]->levelprev = (node *) NULL;
        level[dist + 1] = &dummy;
        dummy.levelprev = (node *) NULL;
#endif
        dist++;
        for (this = next, next = (node *) NULL; this; this = this->tnext) {
            for (e = this->out; e; e = e->outnext) {
                tail = e->ends[1];
#ifdef UNDIRECTED_GRAPH
                if (tail->magiclabel != K && e->cap + e->flow > PRF_EPSILON) {
#else
                if (tail->magiclabel != K && e->flow > PRF_EPSILON) {
#endif
                    tail->flowlabel = dist;
                    tail->tnext = next;
                    next = tail;
                    tail->magiclabel = K;
#ifdef USE_GAP
                    tail->levelnext = level[dist];
                    level[dist]->levelprev = tail;
                    level[dist] = tail;
#endif
#ifdef HIGHEST_LABEL_PRF
                    if (tail->active) {
                        tail->highnext = high[dist];
                        high[dist] = tail;
                    }
#endif
                }
            }
            for (e = this->in; e; e = e->innext) {
                tail = e->ends[0];
                if (tail->magiclabel != K && e->cap - e->flow > PRF_EPSILON) {
                    tail->flowlabel = dist;
                    tail->tnext = next;
                    next = tail;
                    tail->magiclabel = K;
#ifdef USE_GAP
                    tail->levelnext = level[dist];
                    level[dist]->levelprev = tail;
                    level[dist] = tail;
#endif
#ifdef HIGHEST_LABEL_PRF
                    if (tail->active) {
                        tail->highnext = high[dist];
                        high[dist] = tail;
                    }
#endif
                }
            }
        }
#ifdef USE_GAP
        if (dummy.levelprev) {
            dummy.levelprev->levelnext = (node *) NULL;
            level[dist]->levelprev = (node *) NULL;
        } else {
            level[dist] = (node *) NULL;
        }
#endif
#ifdef HIGHEST_LABEL_PRF
        if (high[dist])
            G->highest = dist;
#endif
    } while (next);
}

static int grab_the_cut (graph *G, node *n, int **cut, int *cutcount)
{
    int rval = 0;
    edge *e;
    node *q, *top;
    int count = 0;
    int i, num;
    node *nodelist = G->nodelist;
    int *tcut = (int *) NULL;

    *cut = (int *) NULL;
    *cutcount = 0;

    tcut = CC_SAFE_MALLOC (G->nnodes, int);
    CCcheck_NULL (tcut, "out of memory for tcut");

    G->magicnum++;
    num = G->magicnum;
    tcut[count++] = (int) (n - nodelist);
    q = n;
    q->magiclabel = num;
    q->tnext = (node *) NULL;

    while (q) {
        top = q;
        q = q->tnext;
        for (e = top->out; e; e = e->outnext) {
#ifdef UNDIRECTED_GRAPH
            if (e->cap + e->flow > PRF_EPSILON &&
                e->ends[1]->magiclabel != num) {
#else
            if (e->flow > PRF_EPSILON && e->ends[1]->magiclabel != num) {
#endif
                tcut[count++] = (int) (e->ends[1] - nodelist);
                e->ends[1]->magiclabel = num;
                e->ends[1]->tnext = q;
                q = e->ends[1];
            }
        }
        for (e = top->in; e; e = e->innext) {
            if (e->cap - e->flow > PRF_EPSILON &&
                e->ends[0]->magiclabel != num) {
                tcut[count++] = (int) (e->ends[0] - nodelist);
                e->ends[0]->magiclabel = num;
                e->ends[0]->tnext = q;
                q = e->ends[0];
            }
        }
    }

    *cut = CC_SAFE_MALLOC (count, int);
    CCcheck_NULL (*cut, "out of memory for cut");
    for (i = 0; i < count; i++) {
        (*cut)[i] = tcut[i];
    }
    *cutcount = count;

CLEANUP:

    if (rval) {
        CC_IFFREE (*cut, int);
    }
    CC_IFFREE (tcut, int);
    return rval;
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist,
                       double *ecap)
{
    int i;
    edge *edgelist;
    node *nodelist;

    G->nodelist = (node *) NULL;
    G->edgelist = (edge *) NULL;
#ifdef USE_GAP
    G->level = (node **) NULL;
#endif
#ifdef HIGHEST_LABEL_PRF
    G->high = (node **) NULL;
#endif

    G->magicnum = 0;
    G->nnodes = ncount;
    G->nedges = ecount;
    G->nodelist = CC_SAFE_MALLOC (ncount, node);
    G->edgelist = CC_SAFE_MALLOC (ecount, edge);
    if (!G->nodelist || !G->edgelist) {
        fprintf (stderr, "Out of memory in buildgraph\n");
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->edgelist, edge);
        return 1;
    }
#ifdef USE_GAP
    G->level = CC_SAFE_MALLOC (ncount + 1, node *);
    if (!G->level) {
        fprintf (stderr, "Out of memory in buildgraph\n");
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->edgelist, edge);
        return 1;
    }
    for (i = 0; i < ncount; i++)
        G->level[i] = (node *) NULL;
    G->level[ncount] = (node *) NULL;  /* A guard dog for a while loop */
#endif

#ifdef HIGHEST_LABEL_PRF
    G->high = CC_SAFE_MALLOC (ncount, node *);
    if (!G->high) {
        fprintf (stderr, "Out of memory in buildgraph\n");
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->edgelist, edge);
        return 1;
    }
#endif

    nodelist = G->nodelist;
    edgelist = G->edgelist;

    for (i = 0; i < ncount; i++) {
        nodelist[i].in = (edge *) NULL;
        nodelist[i].out = (edge *) NULL;
        nodelist[i].magiclabel = 0;
    }

    for (i = 0; i < ecount; i++) {
        int tail = elist[2 * i];
        int head = elist[(2 * i) + 1];
        if (tail < 0 || tail >= ncount ||
            head < 0 || head >= ncount) {
            fprintf (stderr, "Edge list in wrong format: Edge %d = [%d, %d]\n",
                         i, tail, head);
            return 1;
        }
        edgelist[i].ends[0] = nodelist + tail;
        edgelist[i].ends[1] = nodelist + head;
        edgelist[i].cap = ecap[i];
        edgelist[i].outnext = nodelist[tail].out;
        nodelist[tail].out = &(edgelist[i]);
        edgelist[i].innext = nodelist[head].in;
        nodelist[head].in = &(edgelist[i]);
    }
    return 0;
}

static void init_graph (graph *G)
{
    if (G) {
        G->nodelist = (node *) NULL;
        G->edgelist = (edge *) NULL;
#ifdef USE_GAP
        G->level = (node **) NULL;
#endif
#ifdef HIGHEST_LABEL_PRF
        G->high = (node **) NULL;
#endif
    }
}

static void free_graph (graph *G)
{
    CC_IFFREE (G->nodelist, node);
    CC_IFFREE (G->edgelist, edge);
#ifdef USE_GAP
    CC_IFFREE (G->level, node *);
#endif
#ifdef HIGHEST_LABEL_PRF
    CC_IFFREE (G->high, node *);
#endif
}
