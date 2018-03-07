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
/*                   Find the Union of a Set of Edge Set                    */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by: Applegate, Bixby, Chvatal, and Cook                         */
/*  Date: May 3, 1999                                                       */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edge_file_union (int ncount, int nfiles, char **flist,       */
/*      int *ecount, int **elist, int **elen, int *foundtour,               */
/*      int *besttourlen)                                                   */
/*    MERGES a list of edge sets.                                           */
/*     -ncount is the number of nodes.                                      */
/*     -nfiles is the number of files to be read.                           */
/*     -flist is the list of the files.                                     */
/*     -ecount, elist, elen returns the merged edge set.                    */
/*     -foundtour will return a 1 if at least one of the files is the       */
/*      edgeset of a tour (it can be NULL).                                 */
/*     -besttour will return the length of the best tour amongst the        */
/*      edge sets (it can be NULL)                                          */
/*     Returns a nonzero value if there was and error.                      */
/*                                                                          */
/*  int CCutil_edge_file_intersect (int ncount, char *fname1, char *fname2, */
/*        int *ecount, int **elist, int **elen)                             */
/*    INTERSECTS the two edge sets in files fname1 and fname2               */
/*                                                                          */
/*  int CCutil_edge_file_difference (int ncount, char *fname1, char *fname2,*/
/*      int *ecount, int **elist, int **elen, int justfirst)                */
/*    SYMMETRIC DIFFERENCE of edge sets in files fname1 and fname2          */
/*     -justfirst if nonzero then only output A - B                         */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

int CCutil_edge_file_union (int ncount, int nfiles, char **flist, int *ecount,
        int **elist, int **elen, int *foundtour, double *besttourlen) 
{
    int i, j, rval = 0;
    CCutil_edgehash h;
    int tcount;
    int *tlist = (int *) NULL;
    int *telen = (int *) NULL;

    *ecount = 0;
    *elist = (int *) NULL;
    *elen  = (int *) NULL;
    if (foundtour) *foundtour = 0;
    if (besttourlen) *besttourlen = CCutil_MAXDOUBLE;

    rval = CCutil_edgehash_init (&h, 2*ncount);
    if (rval) {
        fprintf (stderr, "CCutil_edgehash_init failed\n"); goto CLEANUP;
    }

    for (i = 0; i < nfiles; i++) {
        rval = CCutil_getedgelist (ncount, flist[i], &tcount, &tlist,
                                   &telen, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getedgelist failed\n"); goto CLEANUP;
        }

        for (j = 0; j < tcount; j++) {
            rval = CCutil_edgehash_set (&h, tlist[2*j], tlist[2*j+1],
                                        telen[j]);
            if (rval) {
                fprintf (stderr, "CCutil_edgehash_set failed\n"); 
                goto CLEANUP;
            }
        }

        if (foundtour && (tcount == ncount)) {
            int yesno;
            rval = CCutil_edge_to_cycle (ncount, tlist, &yesno, (int *) NULL);
            if (rval) {
                fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
            }
            if (yesno) {
                *foundtour = 1;
                if (besttourlen) {
                    double len = 0.0;
                    for (j = 0; j < tcount; j++) {
                        len += (double) telen[j];
                    }
                    if (len < *besttourlen) {
                        *besttourlen = len;
                    }
                }
            }

        }

        CC_IFFREE (tlist, int);
        CC_IFFREE (telen, int);
    }

    rval = CCutil_edgehash_getall (&h, ecount, elist, elen);
    if (rval) {
        fprintf (stderr, "CCutil_edgehash_getall failed\n");
        goto CLEANUP;
    }


CLEANUP:

    if (rval) {
        *ecount = 0;
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
        if (foundtour) *foundtour = 0;
        if (besttourlen) *besttourlen = CCutil_MAXDOUBLE;
    }
    CC_IFFREE (tlist, int);
    CC_IFFREE (telen, int);
    CCutil_edgehash_free (&h);

    return rval;
}

int CCutil_edge_file_intersect (int ncount, char *fname1, char *fname2,
        int *ecount, int **elist, int **elen) 
{
    int j, test, len, rval = 0;
    CCutil_edgehash h;
    int tcount, scount;
    int *tlist = (int *) NULL;
    int *telen = (int *) NULL;
    int *slist = (int *) NULL;
    int *selen = (int *) NULL;

    *ecount = 0;
    *elist = (int *) NULL;
    *elen  = (int *) NULL;

    rval = CCutil_edgehash_init (&h, 10*ncount);
    CCcheck_rval (rval, "CCutil_edgehash_init failed");

    rval = CCutil_getedgelist (ncount, fname1, &tcount, &tlist, &telen, 0);
    CCcheck_rval (rval, "CCutil_getedgelist failed");
   
    for (j = 0; j < tcount; j++) {
        rval = CCutil_edgehash_set (&h, tlist[2*j], tlist[2*j+1], telen[j]);
        CCcheck_rval (rval, "CCutil_edgehash_set failed");
    }

    CC_IFFREE (tlist, int);
    CC_IFFREE (telen, int);

    rval = CCutil_getedgelist (ncount, fname2, &tcount, &tlist, &telen, 0);
    CCcheck_rval (rval, "CCutil_getedgelist failed");

    slist = CC_SAFE_MALLOC (2 * tcount, int);
    CCcheck_NULL (tlist, "out of memory for slist");
    selen = CC_SAFE_MALLOC (tcount, int);
    CCcheck_NULL (telen, "out of memory for selen");
    scount = 0;

    for (j = 0; j < tcount; j++) {
        test = CCutil_edgehash_find (&h, tlist[2*j], tlist[2*j+1], &len);
        if (test != -1) {
            if (len != telen[j]) {
                fprintf (stderr, "edge with different lengths in two files\n");
                fprintf (stderr, "%d %d %d vs %d\n", tlist[2*j], tlist[2*j+1],
                                                     len, telen[j]);
                rval = 1;  goto CLEANUP;
            }
            slist[2*scount] = tlist[2*j];
            slist[2*scount+1] = tlist[2*j+1];
            selen[scount] = telen[j];
            scount++;
        }
    }

    if (scount) {
        *elist = CC_SAFE_MALLOC (2*scount, int);
        CCcheck_NULL (*elist, "out of memory for elist");
        *elen = CC_SAFE_MALLOC (scount, int);
        CCcheck_NULL (*elen, "out of memory for elen");
        for (j = 0; j < scount; j++) {
            (*elist)[2*j] = slist[2*j];
            (*elist)[2*j+1] = slist[2*j+1];
            (*elen)[j] = selen[j];
        }
        *ecount = scount;
    }
         
CLEANUP:

    if (rval) {
        *ecount = 0;
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    CC_IFFREE (tlist, int);
    CC_IFFREE (telen, int);
    CC_IFFREE (slist, int);
    CC_IFFREE (selen, int);
    CCutil_edgehash_free (&h);

    return rval;
}

int CCutil_edge_file_difference (int ncount, char *fname1, char *fname2,
        int *ecount, int **elist, int **elen, int justfirst) 
{
    int j, test, len, rval = 0, tcount1, tcount2, scount;
    CCutil_edgehash h;
    int *tlist1 = (int *) NULL, *telen1 = (int *) NULL;
    int *tlist2 = (int *) NULL, *telen2 = (int *) NULL;
    int *slist  = (int *) NULL,  *selen = (int *) NULL;

    *ecount = 0;
    *elist = (int *) NULL;
    *elen  = (int *) NULL;

    rval = CCutil_getedgelist (ncount, fname1, &tcount1, &tlist1, &telen1, 0);
    CCcheck_rval (rval, "CCutil_getedgelist failed");
    rval = CCutil_getedgelist (ncount, fname2, &tcount2, &tlist2, &telen2, 0);
    CCcheck_rval (rval, "CCutil_getedgelist failed");

    slist = CC_SAFE_MALLOC (2 * (tcount1 + tcount2), int);
    CCcheck_NULL (slist, "out of memory for slist");
    selen = CC_SAFE_MALLOC (tcount1 + tcount2, int);
    CCcheck_NULL (selen, "out of memory for selen");
    scount = 0;

    /* if justfirst is 1, then only output A - B */

    if (justfirst == 0) {
        rval = CCutil_edgehash_init (&h, 10*ncount);
        CCcheck_rval (rval, "CCutil_edgehash_init failed");
        for (j = 0; j < tcount1; j++) {
            rval = CCutil_edgehash_set (&h, tlist1[2*j], tlist1[2*j+1],
                                            telen1[j]);
            CCcheck_rval (rval, "CCutil_edgehash_set failed");
        }

        for (j = 0; j < tcount2; j++) {
            test = CCutil_edgehash_find (&h, tlist2[2*j], tlist2[2*j+1], &len);
            if (test == -1) {
                slist[2*scount] = tlist2[2*j];
                slist[2*scount+1] = tlist2[2*j+1];
                selen[scount] = telen2[j];
                scount++;
            }
        }
        CCutil_edgehash_free (&h);
    }

    rval = CCutil_edgehash_init (&h, 10*ncount);
    CCcheck_rval (rval, "CCutil_edgehash_init failed");
    for (j = 0; j < tcount2; j++) {
        rval = CCutil_edgehash_set (&h, tlist2[2*j], tlist2[2*j+1], telen2[j]);
        CCcheck_rval (rval, "CCutil_edgehash_set failed");
    }

    for (j = 0; j < tcount1; j++) {
        test = CCutil_edgehash_find (&h, tlist1[2*j], tlist1[2*j+1], &len);
        if (test == -1) {
            slist[2*scount] = tlist1[2*j];
            slist[2*scount+1] = tlist1[2*j+1];
            selen[scount] = telen1[j];
            scount++;
        }
    }

    if (scount) {
        *elist = CC_SAFE_MALLOC (2*scount, int);
        CCcheck_NULL (*elist, "out of memory for elist");
        *elen = CC_SAFE_MALLOC (scount, int);
        CCcheck_NULL (*elen, "out of memory for elen");
        for (j = 0; j < scount; j++) {
            (*elist)[2*j] = slist[2*j];
            (*elist)[2*j+1] = slist[2*j+1];
            (*elen)[j] = selen[j];
        }
        *ecount = scount;
    }
         
CLEANUP:
    if (rval) {
        *ecount = 0;
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    CC_IFFREE (tlist1, int);
    CC_IFFREE (telen1, int);
    CC_IFFREE (tlist2, int);
    CC_IFFREE (telen2, int);
    CC_IFFREE (slist, int);
    CC_IFFREE (selen, int);
    CCutil_edgehash_free (&h);
    return rval;
}
