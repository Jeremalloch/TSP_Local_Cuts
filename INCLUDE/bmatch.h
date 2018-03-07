/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE's bmatch code                             */
/*                                                                          */
/*  (c) Copyright 2012 by David Applegate and William Cook                 */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __BMATCH_H
#define __BMATCH_H

#define NDEBUG
#include <assert.h>
#include "util.h"

#define CCbmat_MAXIWEIGHT 2147483647

#define CCbmat_ZERO 0
#define CCbmat_ONE 1
#define CCbmat_TWO 2

#define CCbmat_FALSE 0
#define CCbmat_TRUE 1

#define CCbmat_UNLABELED 0
#define CCbmat_PLUSMAT 1
#define CCbmat_MINUSMAT 2

typedef struct CCbmat_node {
    struct CCbmat_edgeptr *edgelist;
    struct {
        struct CCbmat_node *parent;   /* jump does without this structure */
        struct CCbmat_node *child;
        struct CCbmat_node *sibling;
    } tree;
    struct {
        struct CCbmat_node *parent;   /* jump doesn't have a nest */
        struct CCbmat_node *child;
        struct CCbmat_node *sibling;
    } nest;
    int ymat;                  /* used by both */
    int label;                 /* used by both */
    struct CCbmat_edge *parentedge;   /* used by both */
    struct {
        unsigned char pseudo;
        unsigned char inpath;  /* used by both */
        unsigned char marked;
    } status;
    char deficiency;
    struct CCbmat_node *next; /* used in jump, and for pseudonode free list */
                             /* and in ancest, to form surface node list   */
    double co[2];  /* coords for pricing */
    struct {
        struct CCbmat_node *next;
        struct CCbmat_node **prev;
        double order;
    } sort;
    int sum_to_root;
    struct CCbmat_node *surf;
    struct {
        int number;
        struct CCbmat_anc *addset;
        struct CCbmat_anc *evallist;
    } ancest;
    struct {
        int type;
        int loc;
        int key;
    } heap;
} CCbmat_node;

typedef struct CCbmat_anc {
    struct CCbmat_anc *parent;
    int rank;
    struct CCbmat_node *ancest;
    struct CCbmat_node *left;
    struct CCbmat_node *right;
    int weight;
    struct CCbmat_anc *next;
} CCbmat_anc;

typedef struct CCbmat_anc_ptr {
    struct CCbmat_anc *this;
    struct CCbmat_anc_ptr *next;
} CCbmat_anc_ptr;

typedef struct CCbmat_edgeptr {
    struct CCbmat_edge *this;
    struct CCbmat_edgeptr *next;
} CCbmat_edgeptr;

#define CCbmat_THISEDGE(x) ((x)->this)
#define CCbmat_SETTHIS(x)  ((x)->ptrs[0].this = (x)->ptrs[1].this = (x))

typedef struct CCbmat_edge {
    struct CCbmat_edgeptr ptrs[2];
    CCbmat_node *ends[2];
    CCbmat_node *surfends[2];
    int x;            /* edge value (often meaningless) */
    char count;
    struct {
        unsigned char needschange;
    } status;
    unsigned char istooth[2];
    int zmat;                   /* z weight for fractional matching */
    int slack;                  /* slack for integral matching */
    int weight;                 /* weight for integral matching */
    struct CCbmat_edge *next;    /* for fractional matching and keeping track */
                                /* of added edges and edge free list */
    struct {
        int type;
        int loc;
        int key;
    } heap;
    int filler;
} CCbmat_edge;

typedef struct CCbmat_graph {
    int ncount, ecount;
    struct CCbmat_node *nodelist;
    struct CCbmat_edge *edgelist;
    struct CCbmat_edge *newedges;
} CCbmat_graph;

typedef struct CCbmat_node_heap {
    CCbmat_node **entry;
    int total_space;
    int size;
} CCbmat_node_heap;

typedef struct CCbmat_edge_heap {
    CCbmat_edge **entry;
    int total_space;
    int size;
} CCbmat_edge_heap;

typedef struct CCbmat_heaps {
    struct CCbmat_node_heap minus_node_heap;
    struct CCbmat_node_heap minus_pseudo_heap;
    struct CCbmat_node_heap plus_node_heap;
    struct CCbmat_edge_heap minus2_mat_heap;
    struct CCbmat_edge_heap minus1_mat_heap;
    struct CCbmat_edge_heap plus1_mat_heap;
    struct CCbmat_edge_heap plus2_mat_heap;
    struct CCbmat_edge_heap minus2_notmat_heap;
    struct CCbmat_edge_heap minus1_notmat_heap;
    struct CCbmat_edge_heap plus1_notmat_heap;
    struct CCbmat_edge_heap plus2_notmat_heap;
    int delta;
} CCbmat_heaps;

#define CCbmat_NOTHEAPED 0

/* bmatchancest.c */
int CCbmat_ancest_init (CCbmat_graph *G, CCbmat_node **surfacelist,
    CCbmat_anc **anc_freelist, CCbmat_anc_ptr **anc_supply);
int CCbmat_ancest_checkout (CCbmat_graph *G, CCbmat_heaps *H, CCbmat_node *n1,
    CCbmat_node *n2, int w, CCbmat_node **surfacelist, CCbmat_anc **anc_freelist,
    CCbmat_edge **edge_freelist, CCbmat_node **pseudo_freelist,
    CCbmat_anc_ptr **anc_supply, CCbmat_edgeptr **edge_supply);
int CCbmat_ancest_flush (CCbmat_graph *G, CCbmat_heaps *H,
    CCbmat_node **surfacelist, CCbmat_anc **anc_freelist,
    CCbmat_edge **edge_freelist, CCbmat_node **pseudo_freelist,
    CCbmat_anc_ptr **ancest_supply, CCbmat_edgeptr **edge_supply);

/* bmatchheaps.c */
void CCbmat_init_heaps_struct (CCbmat_heaps *H);
void CCbmat_clear_heaps (CCbmat_heaps *H);
void CCbmat_heap_edgeshrunk (CCbmat_heaps *H, CCbmat_edge *e);
void CCbmat_heap_shrink_node (CCbmat_heaps *H, CCbmat_node *n);
int CCbmat_init_heaps (CCbmat_heaps *H);
void CCbmat_free_heaps (CCbmat_heaps *H);
int CCbmat_heap_dispatch (CCbmat_heaps *H, int ybound,
    CCbmat_node **pseudo_freelist);
int CCbmat_heap_label_update (CCbmat_heaps *H, CCbmat_node *n);
CCbmat_node *CCbmat_node_heap_findmin (CCbmat_node_heap *h);
CCbmat_edge *CCbmat_edge_heap_findmin (CCbmat_edge_heap *h);
CCbmat_node *CCbmat_node_heap_deletemin(CCbmat_node_heap *h);
CCbmat_edge *CCbmat_edge_heap_deletemin(CCbmat_edge_heap *h);

/* bmatchjump.c */
int CCbmat_jump (CCbmat_graph *G, int bdegree, int *degrees);

/* bmatch.c */
int CCbmat_onematch (int ncount, int ecount, int *elist, int *elen,
    long *val, CCdatagroup *dat, char *pedgefname, int *mcount, int **melist,
    int **melen);
int CCbmat_twomatch (int ncount, int ecount, int *elist, int *elen,
    long *val, CCdatagroup *dat, char *pedgefname, int *mcount, int **melist,
    int **melen);
int CCbmat_bmatch (int ncount, int ecount, int *elist, int *elen,
    int *degrees, long *val, CCdatagroup *dat, char *pedgefname, int *mcount,
    int **melist, int **melen);
int CCbmat_onematch_geom (int ncount, CCdatagroup *dat, long *val, int *mcount,
    int **melist, int **melen, CCrandstate *rstate);
int CCbmat_twomatch_geom (int ncount, CCdatagroup *dat, long *val, int *mcount,
    int **melist, int **melen, CCrandstate *rstate);
int CCbmat_bmatch_geom (int ncount, CCdatagroup *dat, int *degrees, long *val,
    int *mcount, int **melist, int **melen, CCrandstate *rstate);
void CCbmat_augment2 (CCbmat_heaps *H, CCbmat_node *n, int ybound,
    CCbmat_node **pseudo_freelist);

/* bmatchchk.c */
void CCbmat_surftest (CCbmat_edge *e);
int CCbmat_checklabel (CCbmat_node *n);
CCbmat_node *CCbmat_slowsurf (CCbmat_node *n);

/* bmatchexpand.c */
void CCbmat_addon (CCbmat_node *new, CCbmat_edge *e, CCbmat_node *old);
void CCbmat_node_remove (CCbmat_node *p);
void CCbmat_flipinterval (CCbmat_node *pstart, CCbmat_node *pend);
void CCbmat_flipcycle2 (CCbmat_node *p);
int CCbmat_analyze (CCbmat_heaps *H, CCbmat_edge *e, CCbmat_node **pseudo_freelist);

/* bmatchpseudo.c */
void CCbmat_zeroblowup (CCbmat_node *p, CCbmat_node **pseudo_freelist);
void CCbmat_fixmatching (CCbmat_node *p);
void CCbmat_really_shrink_blossom (CCbmat_heaps *H, CCbmat_edge *e,
    CCbmat_node *proot, CCbmat_node **pseudo_freelist);
void CCbmat_pseudoexpand (CCbmat_heaps *H, CCbmat_node *p,
    CCbmat_node **pseudo_freelist);
int CCbmat_freecombs (CCbmat_graph *G, CCbmat_node **pseudo_freelist);

/* bmatchutil.c */
int CCbmat_matchingedges (CCbmat_graph *G, long *pval, int *pmcount,
    int **pmelist, int **pmelen);
int CCbmat_dumpmatch (CCbmat_graph *G, long *val);
int CCbmat_surfcheck (CCbmat_node *n, CCbmat_node *p);
CCbmat_edge *findedge (CCbmat_node *n1, CCbmat_node *n2);
CCbmat_node *smallest (CCbmat_graph *G, CCbmat_node *p);
CCbmat_node *surfsecond (CCbmat_node *n);
CCbmat_node *surfbelow (CCbmat_node *n, CCbmat_node *p);

/* bmatchalloc.c */
void CCbmat_edge_free(CCbmat_edge *p, CCbmat_edge **edge_freelist);
void CCbmat_node_free (CCbmat_node *p, CCbmat_node **pseudo_freelist);
CCbmat_edge *CCbmat_edge_alloc (CCbmat_edge **edge_freelist,
    CCbmat_edgeptr **edge_supply);
CCbmat_node *CCbmat_node_alloc (CCbmat_node **pseudo_freelist);

/* bmatchpricing.c */
void CCbmat_get_sums (CCbmat_graph *G);
int CCbmat_repair_world (CCbmat_graph *G, CCbmat_heaps *H, CCbmat_node *p1,
    CCbmat_node *p2, int weight, CCbmat_edge **edge_freelist,
    CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply);
int CCbmat_priceedges (char *f, CCbmat_graph *G, CCbmat_heaps *H,
    CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply);
int CCbmat_pricedats (CCdatagroup *dat, CCbmat_graph *G, CCbmat_heaps *H,
    CCbmat_node **pseudo_freelist, CCbmat_edgeptr **edge_supply);
CCbmat_node *lca (CCbmat_node *p1, CCbmat_node *p2);

#define CCbmat_MAX 1
#define CCbmat_EUCLIDEAN 2

#endif  /* __BMATCH_H */
