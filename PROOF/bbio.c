/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--2007 by David Applegate, Robert Bixby,              */
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
/*  Functions for reading/writing/allocating proofs                         */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date:  June 1, 2007                                                     */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int BBio_read_short_proof ()                                            */
/*  int BBio_write_lpcut_in ()                                              */
/*  int BBio_write_lpclique ()                                              */
/*  int BBio_write_lpdomino ()                                              */
/*  int BBio_read_lpcut_in ()                                               */
/*  int BBio_read_lpclique ()                                               */
/*  int BBio_read_lpdomino ()                                               */
/*  int BBio_read_skeleton ()                                               */
/*  int BBio_write_skeleton ()                                              */
/*  void BBio_init_proofnode ()                                             */
/*  void BBio_free_proofnode ()                                             */
/*  void BBio_init_cutproof ()                                              */
/*  void BBio_free_cutproof ()                                              */
/*                                                                          */
/****************************************************************************/

#include "bbproof.h"

static int
    read_short_proof_work (BB_SFILE *f, BBcutproof *p, int silent),
    read_cuts (BB_SFILE *f, int ncount, BBtsp_lpcuts *cuts);

static void
    free_lpcuts (BBtsp_lpcuts *cuts);

void BBio_init_proofnode (BBproofnode *p) 
{
    if (p) {
        p->number = -1;
        p->parent = -1;
        p->child0 = -1;
        p->child1 = -1;
        p->side   = -1;
        p->branch = (BBtsp_branchobj *) NULL;
        p->exact_dual = (BBtsp_bigdual *) NULL;
    }
}

void BBio_free_proofnode (BBproofnode *p) 
{
    if (p) {
        if (p->branch) {
            BButil_free_branchobj (p->branch);
            BB_FREE (p->branch, BBtsp_branchobj);
        }
        if (p->exact_dual) {
            BButil_free_bigdual (&p->exact_dual);
        }
        BBio_init_proofnode (p);
    }
}

int BBio_read_short_proof (BBcutproof **pp, BBcuttype **pcutlist,
        char *fname)
{
    int rval = 0;
    BB_SFILE *f = (BB_SFILE *) NULL;
    int i, ind, cutcount, proofcount;
    BBcutproof *p;
    BBcuttype *clist = (BBcuttype *) NULL;

    *pcutlist = (BBcuttype *) NULL;
    *pp = (BBcutproof *) NULL;

    p = BB_SAFE_MALLOC (1, BBcutproof);
    BBcheck_NULL (p, "out of memory for p");
    BBio_init_cutproof (p);

    f = BBsafe_sopen (fname, "r");
    if (!f) {
        fprintf (stderr, "could not safe open %s for reading\n", fname);
        rval =1; goto CLEANUP;
    }  

    rval = read_short_proof_work (f, p, 0);
    BBcheck_rval (rval, "read_short_proof_work failed");

    /* Get cutlist if present */

    rval = BBsafe_sread_int (f, &cutcount);
    BBcheck_rval (rval, "BBsafe_sread_int failed");

    if (cutcount) {
        if (cutcount != p->cuts->cutcount) {
            fprintf (stderr, "cutlist does not match the proof\n");
            rval = 1; goto CLEANUP;
        }

        clist = BB_SAFE_MALLOC (cutcount, BBcuttype);
        BBcheck_NULL (clist, "out of memory for clist");
        for (i = 0; i < cutcount; i++) {
            rval = BBsafe_sread_int (f, &clist[i].class);
            BBcheck_rval (rval, "BBsafe_sread_int failed");
            rval = BBsafe_sread_int (f, &clist[i].isomorph);
            BBcheck_rval (rval, "BBsafe_sread_int failed");
            rval = BBsafe_sread_int (f, &clist[i].hk);
            BBcheck_rval (rval, "BBsafe_sread_int failed");
            clist[i].proof = (BBcutproof *) NULL;
        }
        *pcutlist = clist;

        rval = BBsafe_sread_int (f, &proofcount);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        printf ("File has %d cut proofs\n", proofcount); fflush (stdout);
        for (i = 0; i < proofcount; i++) {
            rval = BBsafe_sread_int (f, &ind);
            BBcheck_rval (rval, "BBsafe_sread_int failed");
            clist[ind].proof = BB_SAFE_MALLOC (1, BBcutproof);
            BBcheck_NULL (clist[ind].proof, "out of memory for proof");
            printf ("Read proof for cut %d\n", ind); fflush (stdout);
            rval = read_short_proof_work (f, clist[ind].proof, 1);
            BBcheck_rval (rval, "read_short_proof_work failed for cut");
        }
    }
    *pp = p;

CLEANUP:
    if (f) BBsafe_sclose (f);
    return rval;
}

static int read_short_proof_work (BB_SFILE *f, BBcutproof *p, int silent)
{
    int rval = 0;
    int i, j, gotbranch, segcount, cutcount;
    BBtsp_bigdual *d = (BBtsp_bigdual *) NULL;
    BBtsp_branchobj *b = (BBtsp_branchobj *) NULL;
    BBproofnode *t;

    p->probname = BB_SAFE_MALLOC (BBtsp_PROB_FILE_NAME_LEN, char);
    BBcheck_NULL (p->probname, "out of memory for probname");

    for (i = 0; i < BBtsp_PROB_FILE_NAME_LEN; i++) {
        rval = BBsafe_sread_char (f, &p->probname[i]);
        BBcheck_rval (rval, "BBsafe_sread_char failed");
    }

    rval = BBsafe_sread_int (f, &p->ncount);
    BBcheck_rval (rval, "BBsafe_sread_int failed");

    p->cuts = BB_SAFE_MALLOC (1, BBtsp_lpcuts);
    BBcheck_NULL (p->cuts, "out of memory for cuts");
    BButil_init_tsp_lpcuts_struct (p->cuts);

    rval = read_cuts (f, p->ncount, p->cuts);
    BBcheck_rval (rval, "read_cuts failed");
 
    if (!silent) {
        printf ("Number of cuts: %d\n", p->cuts->cutcount); fflush (stdout);
    }

    p->tour = BB_SAFE_MALLOC (p->ncount, int);
    BBcheck_NULL (p->tour, "out of memory for tour");

    for (i = 0; i < p->ncount; i++) {
        rval = BBsafe_sread_int (f, &p->tour[i]);
        BBcheck_rval (rval, "BBsafe_swrite_int failed");
    }

    rval = BBbigguy_sread (f, &p->lpbound);
    BBcheck_rval (rval, "BBbigguy_sread failed");

    if (!silent) {
        printf ("LP bound: %12f\n", BBbigguy_bigguytod (p->lpbound));
        fflush (stdout);
    }

    rval = BBsafe_sread_int (f, &p->bbcount);
    BBcheck_rval (rval, "BBsafe_sread_int failed");

    t = BB_SAFE_MALLOC (p->bbcount, BBproofnode);
    BBcheck_NULL (t, "out of memory for lproof");
    p->lproof = t;

    for (i = 0; i < p->bbcount; i++) {
        rval = BBsafe_sread_int (f, &t[i].number);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &t[i].parent);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &t[i].side);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &t[i].child0);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &t[i].child1);
        BBcheck_rval (rval, "BBsafe_sread_int failed");

        rval = BBsafe_sread_int (f, &gotbranch);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        if (gotbranch) {
            b = BB_SAFE_MALLOC (1, BBtsp_branchobj);
            BBcheck_NULL (b, "out of memory for branchobj");
            BButil_init_branchobj (b);
            rval = BBsafe_sread_int (f, &b->ends[0]);
            BBcheck_rval (rval, "BBsafe_sread_int failed");
            rval = BBsafe_sread_int (f, &b->ends[1]);
            BBcheck_rval (rval, "BBsafe_sread_int failed");

            rval = BBsafe_sread_int (f, &segcount);
            BBcheck_rval (rval, "BBsafe_sread_int failed");
            if (segcount) {
                b->clique = BB_SAFE_MALLOC (1, BBtsp_lpclique);
                BBcheck_NULL (b->clique, "out or memory for clique");
                BButil_init_lpclique (b->clique);
                b->clique->segcount = segcount;
                b->clique->nodes = BB_SAFE_MALLOC (segcount, BBtsp_segment);
                BBcheck_NULL (b->clique->nodes, "out of memory for nodes");
                for (j = 0; j < segcount; j++) {
                    rval = BBsafe_sread_int (f, &b->clique->nodes[j].lo);
                    BBcheck_rval (rval, "BBsafe_sread_int failed");
                    rval = BBsafe_sread_int (f, &b->clique->nodes[j].hi);
                    BBcheck_rval (rval, "BBsafe_sread_int failed");
                }
            } else {
                b->clique = (BBtsp_lpclique *) NULL;
            }
            t[i].branch = b;
        } else {
            t[i].branch = (BBtsp_branchobj *) NULL;
        }
        rval = BBsafe_sread_int (f, &cutcount);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        if (cutcount) {
            d = BB_SAFE_MALLOC (1, BBtsp_bigdual);
            BBcheck_NULL (d, "out of memory for d");
            d->cutcount = cutcount;
            d->node_pi = BB_SAFE_MALLOC (p->ncount, BBbigguy);
            BBcheck_NULL (d->node_pi, "out of memory for node_pi");
            d->cut_pi = BB_SAFE_MALLOC (cutcount, BBbigguy);
            BBcheck_NULL (d->cut_pi, "out of memory for cut_pi");
            for (j = 0; j < p->ncount; j++) {
                rval = BBbigguy_sread (f, &d->node_pi[j]);
                BBcheck_rval (rval, "BBbigguy_sread failed");
            }
            for (j = 0; j < d->cutcount; j++) {
                rval = BBbigguy_sread (f, &d->cut_pi[j]);
                BBcheck_rval (rval, "BBbigguy_sread failed");
            }

            t[i].exact_dual = d;
        } else {
            t[i].exact_dual = (BBtsp_bigdual *) NULL;
        }
    }
    if (!silent) {
        printf ("BBcount: %d\n", p->bbcount); fflush (stdout);
    }

CLEANUP:
    return rval;
}

#define BB_EXTRA_CLIQUES 1000

static int read_cuts (BB_SFILE *f, int ncount, BBtsp_lpcuts *c)
{
    int rval = 0;
    int i, j, nbits, cbits, dbits, cliqcount, dominocount, cutcount;
    int twodom_cliquecount, twodom0, twodom1;
    BBtsp_lpcut *u;
    char version;

    rval = BBsafe_sread_char (f, &version);
    BBcheck_rval (rval, "BBsafe_sread_char failed");
    if (version != 3) {
        fprintf (stderr, "Possibly old cut version %d\n", (unsigned) version);
        rval = 1;  goto CLEANUP;
    }

    rval = BBsafe_sread_int (f, &i);
    BBcheck_rval (rval, "BBsafe_sread_int failed");
    if (i != ncount) {
        fprintf (stderr, "Cuts do not match the proof\n");
        rval = 1;  goto CLEANUP;
    }
    nbits = BBsafe_sbits (ncount);

    rval = BBsafe_sread_int (f, &cliqcount);
    BBcheck_rval (rval, "BBsafe_sread_int failed");

    c->cliques = BB_SAFE_MALLOC (cliqcount+BB_EXTRA_CLIQUES, BBtsp_lpclique);
    BBcheck_NULL (c->cliques, "out of memory for c->cliques");
    c->cliquespace = cliqcount+BB_EXTRA_CLIQUES;
    for (i = 0; i < cliqcount; i++) {
        rval = BBio_read_lpclique (f, &c->cliques[i], ncount);
        BBcheck_rval (rval, "BBio_read_lpclique failed");
    }
    c->cliqueend = cliqcount;

    rval = BBsafe_sread_int (f, &dominocount);
    BBcheck_rval (rval, "BBsafe_sread_int failed");
    if (dominocount > 0) {
        c->dominos  = BB_SAFE_MALLOC (dominocount, BBtsp_lpdomino);
        for (i = 0; i < dominocount; i++) {
            rval = BBio_read_lpdomino (f, &c->dominos[i], ncount);
            BBcheck_rval (rval, "BBio_read_lpdomino failed");
        }
    }
    c->dominoend = dominocount;

    cbits = BBsafe_sbits (cliqcount);
    dbits = BBsafe_sbits (dominocount);
    rval = BBsafe_sread_int (f, &cutcount);
    BBcheck_rval (rval, "BBsafe_sread_int failed");
    c->cuts = BB_SAFE_MALLOC (cutcount+BB_EXTRA_CLIQUES, BBtsp_lpcut);
    BBcheck_NULL (c->cuts, "out of memory for c->cuts");
    c->cutspace = cutcount+BB_EXTRA_CLIQUES;
    for (i = 0; i < cutcount; i++) {
        u = &c->cuts[i];
        rval = BBsafe_sread_int (f, &u->cliquecount);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &twodom_cliquecount);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &u->dominocount);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &twodom0);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_int (f, &twodom1);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        if (twodom_cliquecount > 0 || twodom0 > 0 || twodom1 > 0 ) {
            fprintf (stderr, "Not set up for 2-dom cuts\n");
            rval = 1; goto CLEANUP;
        }
        rval = BBsafe_sread_int (f, &u->rhs);
        BBcheck_rval (rval, "BBsafe_sread_int failed");
        rval = BBsafe_sread_char (f, &u->sense);
        u->branch = 0;
        u->cliques = BB_SAFE_MALLOC (u->cliquecount, int);
        BBcheck_NULL (u->cliques, "out of memory for u->cliques");
        for (j = 0; j < u->cliquecount; j++) {
            rval = BBsafe_sread_bits (f, &u->cliques[j], cbits);
            BBcheck_rval (rval, "BBsafe_sread_bits failed");
        }
        if (u->dominocount > 0) {
            u->dominos = BB_SAFE_MALLOC (u->dominocount, int);
            BBcheck_NULL (u->dominos, "out of memory for u->dominos");
            for (j = 0; j < u->dominocount; j++) {
                rval = BBsafe_sread_bits (f, &u->dominos[j], dbits);
                BBcheck_rval (rval, "BBsafe_sread_bits failed");
            }
        } else {
            u->dominos = (int *) NULL;
        }
        rval = BBio_read_skeleton (f, &u->skel, ncount);
        BBcheck_rval (rval, "BBio_read_skeleton failed");
    }
    c->cutcount = cutcount;

CLEANUP:
    return rval;
}

int BBio_write_lpcut_in (BB_SFILE *f, BBtsp_lpcut_in *c, int ncount)
{
    int i, rval = 0;

    rval = BBsafe_swrite_int (f, c->cliquecount);
    BBcheck_rval (rval, "BBsafe_swrite_int failed");

    rval = BBsafe_swrite_int (f, c->dominocount);
    BBcheck_rval (rval, "BBsafe_swrite_int failed");

    for (i = 0; i < c->cliquecount; i++) {
        rval = BBio_write_lpclique (f, &c->cliques[i], ncount);
        BBcheck_rval (rval, "BBio_write_lpclique failed");
    }

    for (i = 0; i < c->dominocount; i++) {
        rval = BBio_write_lpdomino (f, &c->dominos[i], ncount);
        BBcheck_rval (rval, "BBio_write_lpdomino failed");
    }

    rval = BBsafe_swrite_int (f, c->rhs);
    BBcheck_rval (rval, "BBsafe_swrite_int failed");

    rval = BBsafe_swrite_char (f, c->sense);
    BBcheck_rval (rval, "BBsafe_swrite_char failed");

    rval = BBio_write_skeleton (f, &c->skel, ncount);
    BBcheck_rval (rval, "BBio_write_skeleton failed");
    
CLEANUP:

    return rval;
}

int BBio_write_lpclique (BB_SFILE *f, BBtsp_lpclique *c, int ncount)
{
    int rval = 0;
    int i;
    int nbits = BBsafe_sbits (ncount);

    rval = BBsafe_swrite_bits (f, c->segcount, nbits);
    BBcheck_rval (rval, "BBsafe_swrite_bits failed");

    for (i=0; i < c->segcount; i++) {
        rval = BBsafe_swrite_bits (f, c->nodes[i].lo, nbits);
        BBcheck_rval (rval, "BBsafe_swrite_bits failed");

        rval = BBsafe_swrite_bits (f, c->nodes[i].hi, nbits);
        BBcheck_rval (rval, "BBsafe_swrite_bits failed");
    }

CLEANUP:

    return rval;
}

int BBio_write_lpdomino (BB_SFILE *f, BBtsp_lpdomino *c, int ncount)
{
    int rval = 0;
    int k;

    for (k = 0; k < 2; k++) {
        rval = BBio_write_lpclique (f, &(c->sets[k]), ncount);
        BBcheck_rval (rval, "BBio_write_lpcplique failed");
    }

CLEANUP:

    return rval;
}

int BBio_read_lpcut_in (BB_SFILE *f, BBtsp_lpcut_in *c, int ncount)
{
    int i, rhs, ncliques, ndominos, rval = 0;
    char sense;

    BButil_init_lpcut_in (c);
    
    rval = BBsafe_sread_int (f, &ncliques);
    BBcheck_rval (rval, "BBsafe_sread_int failed");

    rval = BBsafe_sread_int (f, &ndominos);
    BBcheck_rval (rval, "BBsafe_sread_int failed");

    c->cliquecount = ncliques;
    if (ncliques > 0) {
        c->cliques = BB_SAFE_MALLOC (ncliques, BBtsp_lpclique);
        BBcheck_NULL (c->cliques, "out of memory for cliques");
        for (i=0; i<ncliques; i++) {
            BButil_init_lpclique (&c->cliques[i]);
        }
        for (i = 0; i < ncliques; i++) {
            rval = BBio_read_lpclique (f, &c->cliques[i], ncount);
            BBcheck_rval (rval, "BBio_read_lpclique failed");
        }
    } else {
        c->cliques = (BBtsp_lpclique *) NULL;
    }

    c->dominocount = ndominos;
    if (ndominos > 0) {
        c->dominos = BB_SAFE_MALLOC (ndominos, BBtsp_lpdomino);
        BBcheck_NULL (c->dominos, "out of memory for dominos");
        for (i=0; i<ndominos; i++) {
            BButil_init_lpdomino (&c->dominos[i]);
        }
        for (i = 0; i < ndominos; i++) {
            rval = BBio_read_lpdomino (f, &c->dominos[i], ncount);
            BBcheck_rval (rval, "BBio_read_lpdomino failed");
        }
    } else {
        c->dominos = (BBtsp_lpdomino *) NULL;
    }

    rval = BBsafe_sread_int (f, &rhs);
    BBcheck_rval (rval, "BBsafe_sread_int failed");
    c->rhs = rhs;

    rval = BBsafe_sread_char (f, &sense);
    BBcheck_rval (rval, "BBsafe_sread_int failed");
    c->sense = sense;

    c->branch = 0;

    rval = BBio_read_skeleton (f, &c->skel, ncount);
    BBcheck_rval (rval, "BBio_read_skeleton failed");
                  
    rval = 0;

CLEANUP:

    if (rval) {
        BButil_free_lpcut_in (c);
    }
    return rval;
}

int BBio_read_lpclique (BB_SFILE *f, BBtsp_lpclique *c, int ncount)
{
    int nbits = BBsafe_sbits (ncount);
    int i, lo, hi, rval = 0;

    rval = BBsafe_sread_bits (f, &c->segcount, nbits);
    BBcheck_rval (rval, "BBsafe_sread_bits failed");
    c->nodes = BB_SAFE_MALLOC (c->segcount, BBtsp_segment);
    BBcheck_NULL (c->nodes, "out of memory for c->nodes");

    for (i = 0; i < c->segcount; i++) {
        rval = BBsafe_sread_bits (f, &lo, nbits);
        BBcheck_rval (rval, "BBsafe_sread_bits failed");
        rval = BBsafe_sread_bits (f, &hi, nbits);
        BBcheck_rval (rval, "BBsafe_sread_bits failed");
        c->nodes[i].lo = lo;
        c->nodes[i].hi = hi;
    }

CLEANUP:
    return rval;
}

int BBio_read_lpdomino (BB_SFILE *f, BBtsp_lpdomino *d, int ncount)
{
    int k, rval = 0;

    for (k = 0; k < 2; k++) {
        rval = BBio_read_lpclique (f, &(d->sets[k]), ncount);
        BBcheck_rval (rval, "BBio_read_lpclique failed");
    }

CLEANUP:
    return rval;
}

static void free_lpcuts (BBtsp_lpcuts *cuts)
{
    int i;

    if (cuts->cuts) {
        for (i=0; i < cuts->cutcount; i++) {
            BB_IFFREE (cuts->cuts[i].cliques, int);
            BB_IFFREE (cuts->cuts[i].dominos, int);
            BButil_free_skeleton (&cuts->cuts[i].skel);
        }
        BB_FREE (cuts->cuts, BBtsp_lpcut);
    }

    if (cuts->cliques) {
        for (i=0; i < cuts->cliqueend; i++) {
            BB_IFFREE (cuts->cliques[i].nodes, BBtsp_segment);
        }
        BB_FREE (cuts->cliques, BBtsp_lpclique);
    }
    if (cuts->dominos) {
        for (i=0; i < cuts->dominoend; i++) {
            BButil_free_lpdomino (&cuts->dominos[i]);
        }
        BB_FREE (cuts->dominos, BBtsp_lpdomino);
    }
}

void BBio_init_cutproof (BBcutproof *p)
{
    if (p) {
        p->ncount = 0;
        p->bbcount = 0;
        p->tour = (int *) NULL;
        p->probname = (char *) NULL;
        p->lproof = (BBproofnode *) NULL;
        p->cuts = (BBtsp_lpcuts *) NULL;
        p->lpbound = BBbigguy_ZERO ;
    }
}

void BBio_free_cutproof (BBcutproof *p)
{
    int i;

    if (p) {
        BB_IFFREE (p->tour, int);
        BB_IFFREE (p->probname, char);
        if (p->lproof) {
            for (i = 0; i < p->bbcount; i++) {
                BBio_free_proofnode (&p->lproof[i]);
            }
            BB_FREE (p->lproof, BBproofnode);
        }
        if (p->cuts) free_lpcuts (p->cuts);
        BB_IFFREE (p->cuts, BBtsp_lpcuts);
        BBio_init_cutproof (p);
    }
}

#define SKELETON_WILD 0

int BBio_read_skeleton (BB_SFILE *f, BBtsp_skeleton *skel, int ncount)
{
    int rval;
    int atomcount;
    int i;
    int nbits = BBsafe_sbits (ncount);
    char type;
    
    BButil_init_skeleton (skel);

    rval = BBsafe_sread_char (f, &type);
    BBcheck_rval (rval, "BBsafe_sread_char failed");

    switch (type) {
    case SKELETON_WILD:
        rval = BBsafe_sread_bits (f, &atomcount, nbits);
        BBcheck_rval (rval, "BBsafe_sread_bits failed");
        skel->atomcount = atomcount;
        
        if (atomcount == 0) {
            skel->atoms = (int *) NULL;
            rval = 0;
            goto CLEANUP;
        }
        
        skel->atoms = BB_SAFE_MALLOC (atomcount, int);
        BBcheck_NULL (skel->atoms, "out of memory for skel->atoms");

        for (i=0; i<atomcount; i++) {
            rval = BBsafe_sread_bits (f, &skel->atoms[i], nbits);
            BBcheck_rval (rval, "BBsafe_sread_bits failed");
        }
        break;
    default:
        fprintf (stderr, "Unknown skeleton type %ud\n", (unsigned) type);
        rval = 1; goto CLEANUP;
    }
    rval = 0;

CLEANUP:
    if (rval) {
        BButil_free_skeleton (skel);
    }
    return rval;
}

int BBio_write_skeleton (BB_SFILE *f, BBtsp_skeleton *skel, int ncount)
{
    int rval;
    int atomcount = skel->atomcount;
    int i;
    int nbits = BBsafe_sbits (ncount);

    rval = BBsafe_swrite_char (f, SKELETON_WILD);
    BBcheck_rval (rval, "BBsafe_swrite_char failed");
   
    rval = BBsafe_swrite_bits (f, atomcount, nbits);
    BBcheck_rval (rval, "BBsafe_swrite_bits failed");
   
    for (i=0; i<atomcount; i++) {
        rval = BBsafe_swrite_bits (f, skel->atoms[i], nbits);
        BBcheck_rval (rval, "BBsafe_swrite_bits failed");
    }

    rval = 0;

CLEANUP:
    return rval;
}

