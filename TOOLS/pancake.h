#ifndef __PANCAKE_H
#define __PANCAKE_H

int CCtest_blobcuts (int ncount, int ecount, int *elist, double *x, int *hit,
        int (*doit_fn) (double, int, int *, void *), void *pass_param);

typedef struct pannode {
        struct panedge **edgelist;
        struct panedge **goodedge;
        int degree;
        struct vaseknode *vptr;
} pannode;

typedef struct panedge {
        double panweight;
        pannode *ends[2];
        int elim;
        int tag;
        struct vaseknode *a[2];
        struct vaseknode *roof;
        struct vaseknode *top;
        struct panedge *next;
        double rc;
        double realrc;
} panedge;

typedef struct vaseknode {
        struct vaseknode *parent;
        struct vaseknode *child;
        struct vaseknode *sibling;
        struct vaseknode *anc;
        struct vaseknode *ptr;
        struct vaseknode *qtr;
        struct vaseknode *listpointer;
        int b;
        int d;
        int n;
        int tag;
        int fringe;
        double y;
        double w;
        double mult;
        struct panedge *tree;
        struct panedge *junk;
        struct triomino *adj;
        struct triomino *scan;
} vaseknode;

typedef struct triomino {
        panedge *edge;
        vaseknode *end;
        struct triomino *next;
} triomino;

#define VNODEALLOC(vrequest) {                                          \
        if (vnodestack == NULL) {                                       \
                printf ("Ran out of vnode supply\n");                   \
                exit (1);                                               \
        }                                                               \
        vrequest = vnodestack;                                          \
        vnodestack = vnodestack->ptr;                                   \
    }

#define VNODEFREE(vreturn) {                                            \
        vreturn->ptr = vnodestack;                                      \
        vnodestack = vreturn;                                           \
    }

#define TRIALLOC(trequest) {                                            \
        trequest = tristack;                                            \
        tristack = tristack->next;                                      \
    }

#define TRIFREE(treturn) {                                              \
        treturn->next = tristack;                                       \
        tristack = treturn;                                             \
    }
#endif
