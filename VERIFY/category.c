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
/*  Old code from Dave, looking for ismorphic cuts.                         */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "tsp.h"
#include "util.h"
#include "verify.h"

static int textin = 0;
static int probin = 0;
static int poolin = 0;
static char *cutname = (char *) NULL;
static int in_ncount = 0;
static int run_silently = 0;

typedef struct category {
    int type;
    int count;
    double weight;
    int cut[1];
} category;

typedef struct category_pile {
    CCgenhash table;
    int table_init;
    int *mark;
} category_pile;
    
static int
    categorize_text_pool (category_pile *pile, char *poolname, int nodecount),
    categorize_binary_pool (category_pile *pile, char *poolname),
    categorize_problem (category_pile *pile, char *probname, int silent),
    categorize_cut (category_pile *pile, CCtsp_lpcut_in *cut, double wgt),
    pile_create (category_pile *pile, int nodecount),
    parseargs (int ac, char **av);

static void
    usage (char *fname),
    pile_init (category_pile *pile),
    pile_free (category_pile *pile),
    category_freehash (void *k, void *d, void *u_data);


int main (int ac, char **av)
{
    int rval;
    CCutil_timer z;
    category_pile pile;

    CCutil_init_timer (&z, "Categorization");
    CCutil_start_timer (&z);

    pile_init (&pile);

    rval = parseargs (ac, av);
    if (rval) return -1;

    if (textin) {
        rval = categorize_text_pool (&pile, cutname, in_ncount);
        if (rval) {
            fprintf (stderr, "categorize_text_pool failed\n");
            goto CLEANUP;
        }
    } else if (probin) {
        rval = categorize_problem (&pile, cutname, run_silently);
        if (rval) {
            fprintf (stderr, "categorize_problem failed\n");
            goto CLEANUP;
        }
    } else if (poolin) {
        rval = categorize_binary_pool (&pile, cutname);
        if (rval) {
            fprintf (stderr, "categorize_binary_pool failed\n");
            goto CLEANUP;
        }
    }

  CLEANUP:
    pile_free (&pile);
    printf ("Categorization completed in %.2f seconds\n",
            CCutil_stop_timer (&z, 1));
    fflush (stdout);
    return rval;
}

static int categorize_text_pool (category_pile *pile, char *poolname,
                                 int nodecount)
{
    int rval;
    int *tour = (int *) NULL;
    CCtsp_lpcut_in *cuts = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *cnext;
    int cutcount;
    int i;

    rval = pile_create (pile, nodecount);
    if (rval) {
        fprintf (stderr, "pile_create failed\n");
        goto CLEANUP;
    }
    
    tour = CC_SAFE_MALLOC (nodecount, int);
    if (tour == (int *) NULL) {
        fprintf (stderr, "Out of memory in categorize_text_pool\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<nodecount; i++) tour[i] = i;
    
    rval = CCtsp_file_cuts (poolname, &cuts, &cutcount, nodecount, tour);
    if (rval) {
        fprintf (stderr, "CCtsp_file_cuts failed\n");
        goto CLEANUP;
    }

    while (cuts) {
        cnext = cuts->next;
        rval = categorize_cut (pile, cuts, 1.0);
        if (rval) {
            fprintf (stderr, "categorize_cut failed\n");
            goto CLEANUP;
        }
        CCtsp_free_lpcut_in (cuts);
        CC_FREE (cuts, CCtsp_lpcut_in);
        cuts = cnext;
    }

  CLEANUP:
    CC_IFFREE (tour, int);
    while (cuts) {
        cnext = cuts->next;
        CCtsp_free_lpcut_in (cuts);
        CC_FREE (cuts, CCtsp_lpcut_in);
        cuts = cnext;
    }
    return rval;
}

static int categorize_binary_pool (category_pile *pile, char *poolname)
{
    int rval;
    int ncount = 0;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    int i;
    CCtsp_lpcut_in cut;

    rval = CCtsp_init_cutpool (&ncount, poolname, &pool);
    if (rval) {
        fprintf (stderr, "CCtsp_init_cutpool failed\n");
        goto CLEANUP;
    }

    rval = pile_create (pile, ncount);
    if (rval) {
        fprintf (stderr, "pile_create failed\n");
        goto CLEANUP;
    }
    
    for (i=0; i<pool->cutcount; i++) {
        rval = CCtsp_lpcut_to_lpcut_in (pool, &(pool->cuts[i]), &cut);
        if (rval) {
            fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
            goto CLEANUP;
        }
        rval = categorize_cut (pile, &cut, 1.0);
        if (rval) {
            fprintf (stderr, "categorize_cut failed\n");
            goto CLEANUP;
        }
        CCtsp_free_lpcut_in (&cut);
    }

  CLEANUP:
    if (pool) {
        CCtsp_free_cutpool (&pool);
    }
    return rval;
}

static int categorize_problem (category_pile *pile, char *probname, int silent)
{
    int i;
    int rval;
    CCtsp_lp *lp = (CCtsp_lp *) NULL;
    CCtsp_lpcuts *cuts;
    CCtsp_lpcut_in cut;
    int ncount = 0;

    lp = CC_SAFE_MALLOC (1, CCtsp_lp);
    if (lp == (CCtsp_lp *) NULL) {
        fprintf (stderr, "Out of memory in categorize_problem\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_init_tsp_lp_struct (lp);
    rval = CCtsp_read_probfile (lp, probname, (char *) NULL, &ncount, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_read_probfile failed\n");
        goto CLEANUP;
    }

    rval = pile_create (pile, ncount);
    if (rval) {
        fprintf (stderr, "pile_create failed\n");
        goto CLEANUP;
    }
    
    rval = CCtsp_build_lpadj (&lp->graph, 0, lp->graph.ecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n");
        goto CLEANUP;
    }
    
    rval = CCtsp_add_branchhistory_to_lp (lp);
    if (rval) {
        fprintf (stderr, "CCtsp_add_branchhistory_to_lp failed\n");
        goto CLEANUP;
    }

    cuts = &lp->cuts;
    for (i=0; i<cuts->cutcount; i++) {
        if (cuts->cuts[i].branch == 0) {
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &cut);
            if (rval) {
                fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
                goto CLEANUP;
            }
            rval = categorize_cut (pile, &cut, 1.0);
            if (rval) {
                fprintf (stderr, "categorize_cut failed\n");
                goto CLEANUP;
            }
            CCtsp_free_lpcut_in (&cut);
        }
    }

  CLEANUP:
    CCtsp_free_tsp_lp_struct (&lp);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "NPt:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'N': probin = 1; break;
        case 'P': poolin = 1; break;
        case 't': textin = 1; in_ncount = atoi(boptarg); break;
        default: usage (av[0]); return 1;
        }
    }
    if (textin == 0 && probin == 0 && poolin == 0) {
        usage (av[0]);
        return 1;
    }
    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    cutname = av[boptind++];
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [NPt:] cutfile\n", fname);
    fprintf (stderr, "   -N        cutfile is problem file\n");
    fprintf (stderr, "   -P        cutfile is a cut pool\n");
    fprintf (stderr, "   -t n      cutfile is text file, n is the nodecount\n");
    fprintf (stderr, "   One of -N, -P, or -t must be specified\n");
}

typedef struct cat_clique {
    int cat;
    int mult;
    int nnodes;
    int *nodes;
} cat_clique;

typedef struct cat_node {
    int cat;
    int mult;
    int ncliques;
    int *cliques;
} cat_node;

typedef struct cat_work {
    int ncliques;
    int ncliquecats;
    cat_clique *c;
    int nnodes;
    int nnodecats;
    cat_node *n;
    category *cbest;
    category *cwork;
} cat_work;

static int cat_clique_cmp (cat_node *n, cat_clique *a, cat_clique *b)
{
    int i;
    int nnodes;
    
    if (a->cat != b->cat) return a->cat - b->cat;
    if (a->nnodes != b->nnodes) return a->nnodes - b->nnodes;
    if (a->mult != b->mult) return a->mult - b->mult;
    nnodes = a->nnodes;
    for (i=0; i<nnodes; i++) {
        if (n[a->nodes[i]].cat != n[b->nodes[i]].cat) {
            return n[a->nodes[i]].cat - n[b->nodes[i]].cat;
        }
    }
    return 0;
}

static int cat_node_cmp (cat_clique *c, cat_node *a, cat_node *b)
{
    int i;
    int ncliques;
    
    if (a->cat != b->cat) return a->cat - b->cat;
    if (a->ncliques != b->ncliques) return a->ncliques - b->ncliques;
    if (a->mult != b->mult) return a->mult - b->mult;
    ncliques = a->ncliques;
    for (i=0; i<ncliques; i++) {
        if (c[a->cliques[i]].cat != c[b->cliques[i]].cat) {
            return c[a->cliques[i]].cat - c[b->cliques[i]].cat;
        }
    }
    return 0;
}

static void collapse_mults (cat_work *cw)
{
    int ncliques = cw->ncliques;
    int nnodes = cw->nnodes;
    cat_clique *c = cw->c;
    cat_node *n = cw->n;
    int i;
    int j;
    int k;
    cat_clique ctmp;
    cat_node ntmp;
    int itmp;
    int cval;

    for (i=0; i<ncliques; i++) {
        c[i].cat = 0;
        c[i].mult = 1;
    }
    for (i=0; i<nnodes; i++) {
        n[i].cat = i;
    }
    for (i=0; i<ncliques; i++) {
        ctmp = c[i];
        for (j=0; j<ctmp.nnodes; j++) {
            itmp = ctmp.nodes[j];
            for (k=j; k>0 && itmp < ctmp.nodes[k-1]; k--) {
                ctmp.nodes[k] = ctmp.nodes[k-1];
            }
            if (k != j) ctmp.nodes[k] = itmp;
        }
        for (j=i; j>0 && cat_clique_cmp (n, &ctmp, &c[j-1]) < 0; j--) {
            c[j] = c[j-1];
        }
        if (j != i) c[j] = ctmp;
    }

    for (i=1, j=0, k=1; i<ncliques; i++) {
        if (cat_clique_cmp (n, c[j], c[i])) {
            c[j].mult = k;
            j++;
            k = 1;
            if (j != i) c[j] = c[i];
        } else {
            k++;
        }
    }
    c[j].mult = k;
    j++;
    ncliques = j;
    cw->ncliques = j;

    for (i=0; i<ncliques; i++) {
        c[i].cat = i;
    }
    for (i=0; i<nnodes; i++) {
        n[i].cat = 0;
        n[i].mult = 1;
    }
    for (i=0; i<nnodes; i++) {
        ntmp = n[i];
        for (j=0; j<ntmp.ncliques; j++) {
            itmp = ntmp.cliques[j];
            for (k=j; k>0 && itmp < ntmp.cliques[k-1]; k--) {
                ntmp.cliques[k] = ntmp.cliques[k-1];
            }
            if (k != j) ntmp.cliques[k] = itmp;
        }
        for (j=i; j>0 && cat_node_cmp (&ntmp, &n[j-1]) < 0; j--) {
            n[j] = n[j-1];
        }
        if (j != i) n[j] = ntmp;
    }
    
    for (i=1, j=0, k=1; i<nnnodes; i++) {
        if (cat_nnode_cmp (c, n[j], n[i])) {
            n[j].mult = k;
            j++;
            k = 1;
            if (j != i) n[j] = n[i];
        } else {
            k++;
        }
    }
    n[j].mult = k;
    j++;
    nnnodes = j;
    cw->nnnodes = j;
}

static void refine_step (cat_work *cw)
{
    int ncliques = cw->ncliques;
    int nnodes = cw->nnodes;
    cat_clique *c = cw->c;
    cat_clique *n = cw->n;
    int i;
    int j;
    int k;
    cat_clique ctmp;
    cat_node ntmp;
    int itmp;
    int cval;

    for (i=0; i<ncliques; i++) {
        ctmp = c[i];
        for (j=0; j<ctmp.nnodes; j++) {
            itmp = ctmp.nodes[j];
            for (k=j; k>0 && itmp < ctmp.nodes[k-1]; k--) {
                ctmp.nodes[k] = ctmp.nodes[k-1];
            }
            if (k != j) ctmp.nodes[k] = itmp;
        }
        for (j=i; j>0 && cat_clique_cmp (n, &ctmp, &c[j-1]) < 0; j--) {
            c[j] = c[j-1];
        }
        if (j != i) c[j] = ctmp;
    }

    j = 0;
    for (i=1; i<ncliques; i++) {
        /* this strange route is because cat_clique_cmp looks at cat */
        k = cat_clique_cmp (&c[i-1], &c[i]);
        c[i-1].cat = j;
        if (k) {
            j++;
        }
    }
    c[ncliques-1].cat = j;
    cw->ncliquecats = j+1;
    
    for (i=0; i<nnodes; i++) {
        ntmp = n[i];
        for (j=0; j<ntmp.ncliques; j++) {
            itmp = ntmp.cliques[j];
            for (k=j; k>0 && itmp < ntmp.cliques[k-1]; k--) {
                ntmp.cliques[k] = ntmp.cliques[k-1];
            }
            if (k != j) ntmp.cliques[k] = itmp;
        }
        for (j=i; j>0 && cat_node_cmp (&ntmp, &n[j-1]) < 0; j--) {
            n[j] = n[j-1];
        }
        if (j != i) n[j] = ntmp;
    }

    j = 0;
    for (i=1; i<nnodes; i++) {
        /* this strange route is because cat_node_cmp looks at cat */
        k = cat_node_cmp (&n[i-1], &n[i]);
        n[i-1].cat = j;
        if (k) {
            j++;
        }
    }
    n[nnodes-1].cat = j;
    cw->nnodecats = j+1;
}

static void refine_category (cat_work *cw)
{
    int oldnc;
    int oldnn;

    if (cw->ncliquecats == cw->ncliques && cw->nnodecats == cw->nnodes) {
        return;
    }
    do {
        oldnc = cw->ncliquecats;
        oldnn = cw->nnodecats;
        refine_step (cw);
    } while (oldnc < cw->ncliquecats || oldnn < cw->nnodecats);
    if (cw->ncliquecats == cw->ncliques && cw->nnodecats == cw->nnodes) {
        return;
    }

}

static category *build_category (CCtsp_lpcut_in *cut, int *mark)
{
    cat_clique *c = (cat_clique *) NULL;
    cat_node *n = (cat_node *) NULL;
    int *space = (int *) NULL;
    category *ret = (category *) NULL;
    int ncliques = cut->cliquecount;
    int nnodes = cut->skel.atomcount;
    CCtsp_lpclique *cl;
    int *p;
    int i;
    int j;
    int k;
    int n;
    int catsize;
    cat_work cw;
    
    c = CC_SAFE_MALLOC (ncliques, cat_clique);
    n = CC_SAFE_MALLOC (nnodes, cat_node);
    space = CC_SAFE_MALLOC (int, 2 * ncliques * nnodes);
    if (!c || !n || !space) {
        fprintf (stderr, "Out of memory in build_category\n");
        goto CLEANUP;
    }

    p = space;
    for (i=0; i<ncliques; i++) {
        c[i].nnodes = 0;
        c[i].cat = 0;
        c[i].mult = 1;
        c[i].nodes = p;
        p += nnodes;
    }
    for (i=0; i<nnodes; i++) {
        n[i].ncliques = 0;
        n[i].cat = 0;
        n[i].mult = 1;
        n[i].cliques = p;
        p += nnodes;
    }
    for (i=0; i<nnodes; i++) {
        mark[cut->skel.atoms[i]] = i;
    }

    catsize = 0;
    for (i=0; i<cut->cliquecount; i++) {
        cl = &cut->cliques[i];
        CC_FOREACH_NODE_IN_CLIQUE(j,cl,k) {
            n = mark[j];
            if (n >= 0) {
                c[i].nodes[c[i].nnodes++] = n;
                n[n].cliques[n[n].ncliques++] = i;
                catsize++;
            }
        }
        catsize++;
    }
    catsize++;

    for (i=0; i<ncliques; i++) {
        c[i].cat = 0;
    }
    for (i=0; i<nnodes; i++) {
        n[i].cat = 0;
    }

    cw.ncliques = ncliques;
    cw.ncliquecats = 1;
    cw.c = c;
    cw.nnodes = nnodes;
    cw.nnodecats = 1;
    cw.n = n;
    cw.cbest = category_alloc (catsize);
    cw.cwork = category_alloc (catsize);
    if (!cw.cbest || !cw.cwork) {
        goto CLEANUP;
    }

    cw.cbest->type = -1;
    cw.cwork->type = -1;

    refine_category (&cw);

    ret = cw.cbest;
    cw.cbest = (category *) NULL;
    
CLEANUP:
    for (i=0; i<nnodes; i++) {
        mark[cut->skel.atoms[i]] = -1;
    }
    CC_IFFREE (c, cat_clique);
    CC_IFFREE (n, cat_node);
    CC_IFFREE (space, int);
    CC_IFFREE (cw.cwork, category);
    CC_IFFREE (cw.cbest, category);
    return ret;
}

static int categorize_cut (category_pile *pile, CCtsp_lpcut_in *cut,
                           double wgt)
{
    category *c;
    category *c2;
    int rval;
    int type;

    c = build_category (cut, pile->mark);
    if (!c) {
        fprintf (stderr, "build_category failed\n");
        return -1;
    }

    c2 = (category *) CCutil_genhash_lookup (&pile->table, c);
    if (c2) {
        c2->weight += wgt;
        c2->count++;
        category_freehash (c, (void *) NULL, (void *) NULL);
        return 0;
    }
    c->weight = wgt;
    c->count = 1;

    rval = CCverify_cut (cut, CC_TYPE_SIMPLE, &type);
    if (rval) {
        c->type = 'o';
    } else {
        switch (type) {
        case CC_TYPE_SUBTOUR: c->type = 's'; break;
        case CC_TYPE_COMB: c->type = 'c'; break;
        case CC_TYPE_STAR: c->type = 't'; break;
        case CC_TYPE_BIPARTITION: c->type = 'b'; break;
        default: c->type = '?'; break;
        }
    }
    rval = CCutil_genhash_insert (&pile->table, c, c);
    return rval;
}

static void pile_init (category_pile *pile)
{
    pile->table_init = 0;
}

static int category_cmp (void *vc1, void *vc2, void *u_data)
{
    int *c1 = ((category *) vc1)->cut;
    int *c2 = ((category *) vc2)->cut;
    int i;

    for (i=0; c1[i] != -2 && c2[i] != -2; i++) {
        if (c1[i] != c2[i]) return c1[i] - c2[i];
    }
    if (c1[i] != c2[i]) return c1[i] - c2[i];
    return 0;
}

static unsigned int category_hash (void *vc, void *u_data)
{
    int *c = ((category *) vc)->cut;
    int i;
    unsigned int v = 0;

    for (i=0; c[i] != -2; i++) {
        v = v * 37 + c[i];
    }
    return v;
}

static void category_freehash (void *k, void *d, void *u_data)
{
    category *c = (category *) k;
    free (c);
}

static int pile_create (category_pile *pile, int nodecount)
{
    int rval;
    int i;

    pile->mark = CC_SAFE_MALLOC (nodecount, int);
    if (!pile->mark) {
        fprintf (stderr, "Out of memory in pile_create\n");
        return -1;
    }
    for (i=0; i<nodecount; i++) {
        pile->mark[i] = -1;
    }
    
    rval = CCutil_genhash_init (&pile->table, 100000, category_cmp,
                                category_hash, (void *) NULL, 0.8, 0.4);
    if (rval) {
        fprintf (stderr, "CCutil_genhash_init failed\n");
        free (pile->mark);
        return rval;
    }
    pile->table_init = 1;
    return 0;
}

static void pile_free (category_pile *pile)
{
    if (pile->table_init) {
        CCutil_genhash_free (&pile->table, category_freehash);
        pile->table_init = 0;
    }
}

