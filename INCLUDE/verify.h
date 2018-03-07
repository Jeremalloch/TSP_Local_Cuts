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

#ifndef __VERIFY_H
#define __VERIFY_H

#include "tsp.h"

#define CC_TYPE_SUBTOUR 1
#define CC_TYPE_COMB 2
#define CC_TYPE_STAR 4
#define CC_TYPE_BIPARTITION 8
#define CC_TYPE_OTHER 16
#define CC_TYPE_DIRTY_COMB 32
#define CC_TYPE_DOMINO 64
#define CC_TYPE_2P 128
#define CC_TYPE_ALL    (CC_TYPE_SUBTOUR     | CC_TYPE_COMB  | CC_TYPE_STAR | \
                        CC_TYPE_BIPARTITION | CC_TYPE_OTHER | \
                        CC_TYPE_DIRTY_COMB)
#define CC_TYPE_SIMPLE (CC_TYPE_SUBTOUR     | CC_TYPE_COMB  | CC_TYPE_STAR | \
                        CC_TYPE_BIPARTITION | CC_TYPE_DIRTY_COMB)

typedef struct CCverify_cutclass {
    int type;
    int nhandles;
    int nfamilies;
    int *cliques;
    int *inverted;
    int *family_start;
} CCverify_cutclass;


int
    CCverify_cut (CCtsp_lpcut_in *cut, int ncount, int check_types, int *type,
        int use_tsp, int *did_use_tsp, char *cname),
    CCverify_edgeclone (CCtsp_lpcut_in *cut, int ncount, int *reduced,
        CCtsp_lpcut_in *new),
    CCverify_backbone_cut (CCtsp_lpcut_in *cut, CCtsp_lpcut_in *new,
        int *no_outside, int ncount),
    CCverify_classify (CCtsp_lpcut_in *cut, int check_types,
        CCverify_cutclass *cclass),
    CCverify_dump_tsp (CCtsp_lpcut_in *cut, char *cname);

void
    CCverify_initcutclass (CCverify_cutclass *cclass),
    CCverify_freecutclass (CCverify_cutclass *cclass);

#define CCtsp_VERIFY_WORK        'A'
#define CCtsp_VERIFY_NO          'N'
#define CCtsp_VERIFY_RECEIVE     'R'
#define CCtsp_VERIFY_SEND        'S'
#define CCtsp_VERIFY_SAVE        'T'
#define CCtsp_VERIFY_WAIT        'W'
#define CCtsp_VERIFY_YES         'Y'
#define CCtsp_VERIFY_EXIT        'X'

#endif  /* __VERIFY_H */
