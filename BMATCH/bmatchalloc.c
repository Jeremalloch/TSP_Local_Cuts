#include <stdio.h>
#include <stdlib.h>
#include "bmatch.h"

#define EDGE_ALLOC_CHUNK ((65536 - 16) / sizeof (CCbmat_edge))

CCbmat_edge *CCbmat_edge_alloc (CCbmat_edge **edge_freelist,
        CCbmat_edgeptr **edge_supply)
{
    CCbmat_edge *p = (CCbmat_edge *) NULL;

    if (!(*edge_freelist)) {
        CCbmat_edge *pend;
        CCbmat_edgeptr *e;

        *edge_freelist = CC_SAFE_MALLOC (EDGE_ALLOC_CHUNK, CCbmat_edge);
        if (*edge_freelist == (CCbmat_edge *) NULL) {
            fprintf (stderr, "out of memory for edge_freelist\n");
            goto CLEANUP;
        }
        for (p = *edge_freelist, pend = p + EDGE_ALLOC_CHUNK-1; p < pend; p++) {
            p->next = p+1;
        }
        pend->next = (CCbmat_edge *) NULL;
  
        e = CC_SAFE_MALLOC (1, CCbmat_edgeptr);
        if (e == (CCbmat_edgeptr *) NULL) {
            fprintf (stderr, "out of memory for e\n");
            goto CLEANUP;
        }
        e->this = *edge_freelist;
        e->next = *edge_supply;
        *edge_supply = e;
    }

    p = *edge_freelist;
    *edge_freelist = p->next;
    CCbmat_SETTHIS(p);

CLEANUP:
    return p;
}

void CCbmat_edge_free(CCbmat_edge *p, CCbmat_edge **edge_freelist)
{
    p->next = *edge_freelist;
    *edge_freelist = p;
}

CCbmat_node *CCbmat_node_alloc (CCbmat_node **pseudo_freelist)
{
    CCbmat_node *px;

    if (*pseudo_freelist) {
        px = *pseudo_freelist;
        *pseudo_freelist = (*pseudo_freelist)->next;
    } else {
        px = CC_SAFE_MALLOC (1, CCbmat_node);
        if (!px) {
            fprintf (stderr, "out of memory for px\n");
            goto CLEANUP;
        }
    }
    px->edgelist = (CCbmat_edgeptr *) NULL;
    px->parentedge = (CCbmat_edge *) NULL;
    px->tree.parent = px->tree.sibling = px->tree.child = (CCbmat_node *) NULL;
    px->nest.parent = px->nest.sibling = px->nest.child = (CCbmat_node *) NULL;
    px->ymat = 0;
    px->status.pseudo = CCbmat_FALSE;
    px->status.inpath = CCbmat_FALSE;
    px->deficiency = 0;
    px->label = CCbmat_UNLABELED;
    px->surf = px;
    px->heap.type = CCbmat_NOTHEAPED;
CLEANUP:
    return px;
}

void CCbmat_node_free (CCbmat_node *p, CCbmat_node **pseudo_freelist)
{
    p->next = *pseudo_freelist;
    *pseudo_freelist = p;
}

