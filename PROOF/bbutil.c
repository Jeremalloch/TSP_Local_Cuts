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
/*  Various functions modfied from Concorde for TSP proof                   */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 2, 2007                                                       */
/*                                                                          */
/****************************************************************************/

#include "bbproof.h"

void *BButil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }
    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void BButil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }
    free (p);
}

double BButil_zeit (void)
{
    struct rusage ru;
    getrusage (RUSAGE_SELF, &ru);
    return ((double) ru.ru_utime.tv_sec) +
            ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

void BButil_int_array_quicksort (int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) return;
    BB_SWAP (len[0], len[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[0];

    while (1) {
        do i++; while (i < n && len[i] < t);
        do j--; while (len[j] > t);
        if (j < i) break;
        BB_SWAP (len[i], len[j], temp);
    }
    BB_SWAP (len[0], len[j], temp);

    BButil_int_array_quicksort (len, j);
    BButil_int_array_quicksort (len + i, n - i);
}

void BButil_int_perm_quicksort (int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) return;
    BB_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        BB_SWAP (perm[i], perm[j], temp);
    }
    BB_SWAP (perm[0], perm[j], temp);

    BButil_int_perm_quicksort (perm, len, j);
    BButil_int_perm_quicksort (perm + i, len, n - i);
}


static double
    dtrunc (double);

static int
    edgelen_nonorm (int i, int j, BBdatagroup *dat),
    max_edgelen (int i, int j, BBdatagroup *dat),
    man_edgelen (int i, int j, BBdatagroup *dat),
    euclid_edgelen (int i, int j, BBdatagroup *dat),
    euclid_ceiling_edgelen (int i, int j, BBdatagroup *dat),
    euclid3d_edgelen (int i, int j, BBdatagroup *dat),
    geographic_edgelen (int i, int j, BBdatagroup *dat),
    geom_edgelen (int i, int j, BBdatagroup *dat),
    att_edgelen (int i, int j, BBdatagroup *dat),
    matrix_edgelen (int i, int j, BBdatagroup *dat),
    sparse_edgelen (int i, int j, BBdatagroup *dat);


#ifndef M_PI
#define M_PI 3.14159265358979323846264
#endif

int BButil_dat_edgelen (int i, int j, BBdatagroup *dat)
{
    return (dat->edgelen)(i, j, dat);
}

int BButil_dat_setnorm (BBdatagroup *dat, int norm)
{
    switch (norm) {
    case BB_EUCLIDEAN_CEIL:
        dat->edgelen = euclid_ceiling_edgelen;
        break;
    case BB_EUCLIDEAN:
        dat->edgelen = euclid_edgelen;
        break;
    case BB_MAXNORM:
        dat->edgelen = max_edgelen;
        break;
    case BB_MANNORM:
        dat->edgelen = man_edgelen;
        break;
    case BB_EUCLIDEAN_3D:
        dat->edgelen = euclid3d_edgelen;
        break;
    case BB_GEOGRAPHIC:
        dat->edgelen = geographic_edgelen;
        break;
    case BB_GEOM:
        dat->edgelen = geom_edgelen;
        break;
    case BB_ATT:
        dat->edgelen = att_edgelen;
        break;
    case BB_MATRIXNORM:
        dat->edgelen = matrix_edgelen;
        break;
    case BB_SPARSE:
        dat->edgelen = sparse_edgelen;
        break;
    default:
        fprintf (stderr, "ERROR:  Unknown NORM %d.\n", norm);
        return 1;
    }
    dat->norm = norm;

    return 0;
}

void BButil_dat_getnorm (BBdatagroup *dat, int *norm)
{
    (*norm) = dat->norm;
}

static int edgelen_nonorm (int i, int j, BBdatagroup *dat)
{
    fprintf (stderr, "BButil_dat_edgelen has been called with no norm set\n");
    fprintf (stderr, "This is a FATAL ERROR\n");
    if (i != 0 || j != 0 || dat != (BBdatagroup *) NULL) {
        /* so the compiler won't complain about unused variables */
        fprintf (stderr, "This is a FATAL ERROR\n");
        exit (1);
    }
    return -1;
}

/* Several variables that would normally be called y1 and y2 are called
   yy1 and yy2 to avoid conflict with the bessel functions */

static int max_edgelen (int i, int j, BBdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;
    t1 += 0.5;
    t2 += 0.5;

    return (int) (t1 < t2 ? t2 : t1);
}

static int man_edgelen (int i, int j, BBdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;

    return (int) (t1 + t2 + 0.5);
}


static int euclid_edgelen (int i, int j, BBdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
    return temp;
}

static int euclid3d_edgelen (int i, int j, BBdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    double t3 = dat->z[i] - dat->z[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2 + t3 * t3) + 0.5);
    return temp;
}

static int euclid_ceiling_edgelen (int i, int j, BBdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
/*
    int rd;
    double max;

    max = sqrt (t1 * t1 + t2 * t2);
    rd = (int) max;
    return (((max - rd) > .000000001) ? rd + 1 : rd);
*/
    return (int) (ceil (sqrt (t1 * t1 + t2 * t2)));
}

#define GH_PI (3.141592)

static int geographic_edgelen (int i, int j, BBdatagroup *dat)
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;
    double x1 = dat->x[i], x2 = dat->x[j], yy1 = dat->y[i], yy2 = dat->y[j];

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (yy1);
    min = yy1 - deg;
    longi = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (yy2);
    min = yy2 - deg;
    longj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                + 1.0);
    return dd;
}

static int geom_edgelen (int i, int j, BBdatagroup *dat)
{
    double lati, latj, longi, longj;
    double q1, q2, q3, q4, q5;

    lati = M_PI * dat->x[i] / 180.0;
    latj = M_PI * dat->x[j] / 180.0;

    longi = M_PI * dat->y[i] / 180.0;
    longj = M_PI * dat->y[j] / 180.0;

    q1 = cos (latj) * sin(longi - longj);
    q3 = sin((longi - longj)/2.0);
    q4 = cos((longi - longj)/2.0);
    q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
    q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
    return (int) (6378388.0 * atan2(sqrt(q1*q1 + q2*q2), q5) + 1.0);
}

#if 0
static int geom_edgelen (int i, int j, BBdatagroup *dat)
{
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;

    lati = M_PI * (dat->x[i] / 180.0);
    latj = M_PI * (dat->x[j] / 180.0);

    longi = M_PI * (dat->y[i] / 180.0);
    longj = M_PI * (dat->y[j] / 180.0);

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378388.0 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                + 1.0);
    return dd;
}
#endif

static int att_edgelen (int i, int j, BBdatagroup *dat)
{
    double xd = dat->x[i] - dat->x[j];
    double yd = dat->y[i] - dat->y[j];
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

static int matrix_edgelen (int i, int j, BBdatagroup *dat)
{
    if (i > j)
        return (dat->adj[i])[j];
    else
        return (dat->adj[j])[i];
}

static int sparse_edgelen (int i, int j, BBdatagroup *dat)
{
    int *adj;
    int k, deg;

    if (i > j) {
        BB_SWAP (i, j, k);
    }
    adj = dat->adj[i];
    deg = dat->degree[i];

    for (k = 0; k < deg; k++) {
        if (adj[k] == j) {
            return dat->len[i][k];
        }
    }
    return dat->default_len;
}

void BButil_init_datagroup (BBdatagroup *dat)
{
    dat->x = (double *) NULL;
    dat->y = (double *) NULL;
    dat->z = (double *) NULL;
    dat->adj = (int **) NULL;
    dat->adjspace = (int *) NULL;
    dat->len      = (int **) NULL;
    dat->lenspace = (int *) NULL;
    dat->degree   = (int *) NULL;
    dat->norm = 0;
    dat->default_len = 100000;
    dat->sparse_ecount = 0;
    dat->edgelen = edgelen_nonorm;
}

void BButil_freedatagroup (BBdatagroup *dat)
{
    BB_IFFREE (dat->x, double);
    BB_IFFREE (dat->y, double);
    BB_IFFREE (dat->z, double);
    BB_IFFREE (dat->adj, int *);
    BB_IFFREE (dat->adjspace, int);
    BB_IFFREE (dat->len, int *);
    BB_IFFREE (dat->lenspace, int);
    BB_IFFREE (dat->degree, int);
}

void BButil_init_lpcut (BBtsp_lpcut *c)
{
    if (c) {
        c->cliquecount = 0;
        c->dominocount = 0;
        c->rhs = 0;
        c->sense = 'X';
        c->branch = 0;
        c->cliques = (int *) NULL;
        c->dominos = (int *) NULL;
        BButil_init_skeleton (&c->skel);
    }
}

void BButil_init_lpcut_in (BBtsp_lpcut_in *c)
{
    if (c) {
        c->cliquecount = 0;
        c->dominocount = 0;
        c->rhs            = 0;
        c->sense          = 'X';
        c->branch         = 0;
        c->cliques        = (BBtsp_lpclique *) NULL;
        c->dominos        = (BBtsp_lpdomino *) NULL;
        BButil_init_skeleton (&c->skel);
    }
}

void BButil_free_lpcut_in (BBtsp_lpcut_in *c)
{
    int i;

    if (c != (BBtsp_lpcut_in *) NULL) {
        if (c->cliques != (BBtsp_lpclique *) NULL) {
            for (i = 0; i < c->cliquecount; i++) {
                BButil_free_lpclique (&c->cliques[i]);
            }
            BB_IFFREE (c->cliques, BBtsp_lpclique);
        }
        if (c->dominos != (BBtsp_lpdomino *) NULL) {
            for (i = 0; i < c->dominocount; i++) {
                BButil_free_lpdomino (&c->dominos[i]);
            }
            BB_IFFREE (c->dominos, BBtsp_lpdomino);
        }
        BButil_free_skeleton (&c->skel);
    }
}

int BButil_lpcut_to_lpcut_in (BBtsp_lpcuts *cuts, BBtsp_lpcut *c,
        BBtsp_lpcut_in *new)
{
    int i;
    BBtsp_lpclique *cl;
    BBtsp_lpdomino *dom;
    int rval = 0;

    BButil_init_lpcut_in (new);
    
    new->rhs = c->rhs;
    new->sense = c->sense;
    new->branch = c->branch;

    new->cliques = BB_SAFE_MALLOC (c->cliquecount, BBtsp_lpclique);
    BBcheck_NULL (new->cliques, "out of memory for cliques");
    if (c->dominocount > 0) {
        new->dominos = BB_SAFE_MALLOC (c->dominocount, BBtsp_lpdomino);
        BBcheck_NULL (new->dominos, "out of memory for dominos");
    }

    for (i = 0; i < c->cliquecount; i++) {
        cl = &(cuts->cliques[c->cliques[i]]);
        rval = BButil_copy_lpclique (cl, &new->cliques[i]);
        BBcheck_rval (rval, "BButil_copy_lpclique failed");
        new->cliquecount++;
    }

    if (c->dominocount > 0) {
        for (i = 0; i < c->dominocount; i++) {
            dom = &(cuts->dominos[c->dominos[i]]);
            rval = BButil_copy_lpdomino (dom, &new->dominos[i]);
            BBcheck_rval (rval , "BButil_copy_lpdomino failed");
            new->dominocount++;
        }
    }

    rval = BButil_copy_skeleton (&c->skel, &new->skel);
    BBcheck_rval (rval, "BButil_copy_skeleton failed");

    if (new->cliquecount != c->cliquecount ||
        new->dominocount != c->dominocount) {
        fprintf (stderr, "error in counts\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    if (rval) BButil_free_lpcut_in (new);
    return rval;
}

void BButil_init_lpclique (BBtsp_lpclique *c)
{
    if (c) {
        c->segcount = 0;
        c->nodes = (BBtsp_segment *) NULL;
    }
}

void BButil_free_lpclique (BBtsp_lpclique *c)
{
    if (c) {
        BB_IFFREE (c->nodes, BBtsp_segment);
        c->segcount = 0;
    }
}

void BButil_init_lpdomino (BBtsp_lpdomino *c)
{
    if (c) {
        BButil_init_lpclique (&(c->sets[0]));
        BButil_init_lpclique (&(c->sets[1]));
    }
}

void BButil_free_lpdomino (BBtsp_lpdomino *c)
{
    if (c) {
        BButil_free_lpclique (&(c->sets[0]));
        BButil_free_lpclique (&(c->sets[1]));
    }
}

void BButil_free_bigdual (BBtsp_bigdual **d) 
{
    if (d) {                 
        if (*d) {                                                           
            BB_IFFREE ((*d)->node_pi, BBbigguy);                            
            BB_IFFREE ((*d)->cut_pi,  BBbigguy);
            BB_IFFREE ((*d), BBtsp_bigdual);                                
        }                                                                   
    }                                                                       
}

void BButil_init_branchobj (BBtsp_branchobj *b)
{
    if (b) {
        b->depth     = 0;
        b->rhs       = 0;
        b->ends[0]   = -1;
        b->ends[1]   = -1;
        b->sense     = 'X';
        b->clique    = (BBtsp_lpclique *) NULL;
    }
}

void BButil_free_branchobj (BBtsp_branchobj *b)
{
    if (b) {
        b->depth     = 0;
        b->rhs       = 0;
        b->ends[0]   = -1;
        b->ends[1]   = -1;
        b->sense     = 'X';
        if (b->clique) {
            BButil_free_lpclique (b->clique);
            BB_FREE (b->clique, BBtsp_lpclique);
        }
    }
}

void BButil_init_tsp_lpcuts_struct (BBtsp_lpcuts *c)
{
    if (c) {
        c->cutcount        = 0;
        c->savecount       = 0;
        c->cliqueend       = 0;
        c->cutspace        = 0;
        c->cliquespace     = 0;
        c->cuts            = (BBtsp_lpcut *) NULL;
        c->cliques         = (BBtsp_lpclique *) NULL;
        c->dominoend       = 0;
        c->dominospace     = 0;
        c->dominos         = (BBtsp_lpdomino *) NULL;
    }
}

void BButil_init_skeleton (BBtsp_skeleton *skel)
{
    skel->atomcount = 0;
    skel->atoms = (int *) NULL;
}

void BButil_free_skeleton (BBtsp_skeleton *skel)
{
    if (skel != (BBtsp_skeleton *) NULL) {
        skel->atomcount = 0;
        BB_IFFREE (skel->atoms, int);
    }
}

int BButil_copy_skeleton (BBtsp_skeleton *old, BBtsp_skeleton *new)
{
    int i;

    BButil_init_skeleton (new);

    if (old->atomcount == 0) return 0;
    new->atoms = BB_SAFE_MALLOC (old->atomcount, int);
    if (new->atoms == (int *) NULL) {
        fprintf (stderr, "Out of memory in BButil_copy_skeleton\n");
        return 1;
    }
    for (i=0; i<old->atomcount; i++) {
        new->atoms[i] = old->atoms[i];
    }
    new->atomcount = old->atomcount;

    return 0;
}

int BButil_array_to_lpclique (int *ar, int acount, BBtsp_lpclique *cliq)
{
    int i, nseg, rval = 0;

    /* Function will alter the order on the array */

    BButil_int_array_quicksort (ar, acount);
    nseg = 0;
    i = 0;
    while (i < acount) {
        while (i < (acount - 1) && ar[i + 1] == (ar[i] + 1))
            i++;
        i++;
        nseg++;
    }

    cliq->nodes = BB_SAFE_MALLOC (nseg, BBtsp_segment);
    BBcheck_NULL (cliq->nodes, "out of memory for cliq->nodes");
    cliq->segcount = nseg;

    nseg = 0;
    i = 0;
    while (i < acount) {
        cliq->nodes[nseg].lo = ar[i];
        while (i < (acount - 1) && ar[i + 1] == (ar[i] + 1))
            i++;
        cliq->nodes[nseg].hi = ar[i++];
        nseg++;
    }

CLEANUP:
    return rval;
}

int BButil_clique_to_array (BBtsp_lpclique *c, int **ar, int *count)
{
    int rval = 0;
    int i, j, tmp;
    int k = 0;

    *ar = (int *) NULL;

    *count = 0;
    for (i = 0; i < c->segcount; i++) {
        (*count) += (c->nodes[i].hi - c->nodes[i].lo + 1);
    }

    if (*count) {
        *ar = BB_SAFE_MALLOC (*count, int);
        if (!(*ar)) {
            fprintf (stderr, "out of memory in BButil_clique_to_array\n");
            rval = 1; goto CLEANUP;
        }
        BB_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
            (*ar)[k++] = j;
        }
    }

CLEANUP:

    return rval;
}

int BButil_copy_lpclique (BBtsp_lpclique *c, BBtsp_lpclique *new)
{
    int k;
    BBtsp_segment *s = (BBtsp_segment *) NULL;

    BButil_init_lpclique (new);
    if (c->segcount) {
        s = BB_SAFE_MALLOC (c->segcount, BBtsp_segment);
        if (!s) {
            fprintf (stderr, "out of memory in copy_lpclique\n");
            return 1;
        }
        for (k = 0; k < c->segcount; k++) {
            s[k].lo = c->nodes[k].lo;
            s[k].hi = c->nodes[k].hi;
        }
    }
    new->segcount = c->segcount;
    new->nodes = s;
    return 0;
}

void BButil_lpclique_compare (BBtsp_lpclique *a, BBtsp_lpclique *b, int *diff)
{
    int i;

    if (a->segcount != b->segcount) {
        *diff = 1; return;
    } else {
        for (i = 0; i < a->segcount; i++) {
            if (a->nodes[i].lo != b->nodes[i].lo) {
                *diff = 1; return;
            }
            if (a->nodes[i].hi != b->nodes[i].hi) {
                *diff = 1; return;
            }
        }
    }
    *diff = 0; return;
}

int BButil_copy_lpdomino (BBtsp_lpdomino *c, BBtsp_lpdomino *new)
{
    int k;
    int rval = 0;
     
    BButil_init_lpdomino (new);
    for (k = 0; k < 2; k++) {
        rval = BButil_copy_lpclique (&(c->sets[k]), &(new->sets[k]));
        if (rval) {
            fprintf (stderr, "BButil_copy_lpclique failed\n");
            BButil_free_lpdomino (new);
            goto CLEANUP;
        }   
    }   
        
CLEANUP:    
    return rval;
}  

#define MATRIX_LOWER_DIAG_ROW  0
#define MATRIX_UPPER_ROW       1
#define MATRIX_UPPER_DIAG_ROW  2
#define MATRIX_FULL_MATRIX     3

int BButil_gettsplib (char *datname, int *ncount, BBdatagroup *dat)
{
    char buf[256], key[256], field[256];
    char *p;
    FILE *in;
    int matrixform = MATRIX_LOWER_DIAG_ROW;
    int norm = -1;

    BButil_init_datagroup (dat);
    *ncount = -1;

    if ((in = fopen (datname, "r")) == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        return 1;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        while (*p != '\0') {
            if (*p == ':')
                *p = ' ';
            p++;
        }
        p = buf;
        if (sscanf (p, "%s", key) != EOF) {
            p += strlen (key);
            while (*p == ' ')
                p++;
            if (!strcmp (key, "NAME")) {
                printf ("Problem Name: %s", p);
            } else if (!strcmp (key, "TYPE")) {
                printf ("Problem Type: %s", p);
                if (sscanf (p, "%s", field) == EOF || strcmp (field, "TSP")) {
                    fprintf (stderr, "Not a TSP problem\n");
                    return 1;
                }
            } else if (!strcmp (key, "COMMENT")) {
                printf ("%s", p);
            } else if (!strcmp (key, "DIMENSION")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in DIMENSION line\n");
                    return 1;
                }
                *ncount = atoi (field);
                printf ("Number of Nodes: %d\n", *ncount);
            } else if (!strcmp (key, "EDGE_WEIGHT_TYPE")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in EDGE_WEIGHT_TYPE line\n");
                    return 1;
                }
                if (!strcmp (field, "EXPLICIT")) {
                    norm = BB_MATRIXNORM;
                    printf ("Explicit Lengths (BB_MATRIXNORM)\n");
                } else if (!strcmp (field, "EUC_2D")) {
                    norm = BB_EUCLIDEAN;
                    printf ("Rounded Euclidean Norm (BB_EUCLIDEAN)\n");
                } else if (!strcmp (field, "EUC_3D")) {
                    norm = BB_EUCLIDEAN_3D;
                    printf ("Rounded Euclidean 3D Norm (BB_EUCLIDEAN_3D)\n");
                } else if (!strcmp (field, "MAX_2D")) {
                    norm = BB_MAXNORM;
                    printf ("Max Norm (BB_MAXNORM)\n");
                } else if (!strcmp (field, "MAN_2D")) {
                    norm = BB_MANNORM;
                    printf ("L1 Norm (BB_MANNORM)\n");
                } else if (!strcmp (field, "GEO")) {
                    norm = BB_GEOGRAPHIC;
                    printf ("Geographical Norm (BB_GEOGRAPHIC)\n");
                } else if (!strcmp (field, "GEOM")) {
                    norm = BB_GEOM;
                    printf ("Geographical Norm in Meters (BB_GEOM)\n");
                } else if (!strcmp (field, "ATT")) {
                    norm = BB_ATT;
                    printf ("ATT Norm (BB_ATT)\n");
                } else if (!strcmp (field, "CEIL_2D")) {
                    norm = BB_EUCLIDEAN_CEIL;
                    printf ("Rounded Up Euclidean Norm (BB_EUCLIDEAN_CEIL)\n");
                } else {
                    fprintf (stderr, "ERROR: Not set up for norm %s\n", field);
                    return 1;
                }
                if (BButil_dat_setnorm (dat, norm)) {
                    fprintf (stderr, "ERROR: Couldn't set norm %d\n", norm);
                    return 1;
                }
            } else if (!strcmp (key, "EDGE_WEIGHT_FORMAT")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in EDGE_WEIGHT_FORMAT line\n");
                    return 1;
                }
                if (!strcmp (field, "LOWER_DIAG_ROW")) {
                    matrixform = MATRIX_LOWER_DIAG_ROW;
                } else if (!strcmp (field, "UPPER_ROW")) {
                    matrixform = MATRIX_UPPER_ROW;
                } else if (!strcmp (field, "UPPER_DIAG_ROW")) {
                    matrixform = MATRIX_UPPER_DIAG_ROW;
                } else if (!strcmp (field, "FULL_MATRIX")) {
                    matrixform = MATRIX_FULL_MATRIX;
                } else if (strcmp (field, "FUNCTION")) {
                    fprintf (stderr, "Cannot handle format: %s\n", field);
                    return 1;
                }
            } else if (!strcmp (key, "NODE_COORD_SECTION")) {
                int i;
                if (*ncount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (dat->x != (double *) NULL) {
                    fprintf (stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    BButil_freedatagroup (dat);
                    return 1;
                }
                if ((norm & BB_NORM_SIZE_BITS) == BB_D2_NORM_SIZE) {
                    dat->x = BB_SAFE_MALLOC (*ncount, double);
                    if (!dat->x) {
                        BButil_freedatagroup (dat);
                        return 1;
                    }
                    dat->y = BB_SAFE_MALLOC (*ncount, double);
                    if (!dat->y) {
                        BButil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0; i < *ncount; i++) {
                        fscanf (in, "%*d %lf %lf", &(dat->x[i]), &(dat->y[i]));
                    }
                } else if ((norm & BB_NORM_SIZE_BITS) == BB_D3_NORM_SIZE) {
                    dat->x = BB_SAFE_MALLOC (*ncount, double);
                    if (!dat->x) {
                        BButil_freedatagroup (dat);
                        return 1;
                    }
                    dat->y = BB_SAFE_MALLOC (*ncount, double);
                    if (!dat->y) {
                        BButil_freedatagroup (dat);
                        return 1;
                    }
                    dat->z = BB_SAFE_MALLOC (*ncount, double);
                    if (!dat->z) {
                        BButil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0; i < *ncount; i++) {
                        fscanf (in, "%*d %lf %lf %lf",
                               &(dat->x[i]), &(dat->y[i]), &(dat->z[i]));
                    }
                } else {
                    fprintf (stderr, "ERROR: Node coordinates with norm %d?\n",
                                 norm);
                    return 1;
                }
            } else if (!strcmp (key, "EDGE_WEIGHT_SECTION")) {
                int i, j;
                if (*ncount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (dat->adj != (int **) NULL) {
                    fprintf (stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    BButil_freedatagroup (dat);
                    return 1;
                }
                if ((norm & BB_NORM_SIZE_BITS) == BB_MATRIX_NORM_SIZE) {
                    dat->adj = BB_SAFE_MALLOC (*ncount, int *);
                    dat->adjspace = BB_SAFE_MALLOC ((*ncount)*(*ncount+1)/2,
                                                    int);
                    if (dat->adj == (int **) NULL ||
                        dat->adjspace == (int *) NULL) {
                        BButil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0, j = 0; i < *ncount; i++) {
                        dat->adj[i] = dat->adjspace + j;
                        j += (i+1);
                    }
                    if (matrixform == MATRIX_LOWER_DIAG_ROW) {
                        for (i = 0; i < *ncount; i++) {
                            for (j = 0; j <= i; j++)
                                fscanf (in, "%d", &(dat->adj[i][j]));
                        }
                    } else if (matrixform == MATRIX_UPPER_ROW ||
                               matrixform == MATRIX_UPPER_DIAG_ROW ||
                               matrixform == MATRIX_FULL_MATRIX) {
                        int **tempadj = (int **) NULL;
                        int *tempadjspace = (int *) NULL;
                        tempadj = BB_SAFE_MALLOC (*ncount, int *);
                        tempadjspace = BB_SAFE_MALLOC ((*ncount) * (*ncount),
                                                       int);
                        if (tempadj == (int **) NULL ||
                            tempadjspace == (int *) NULL) {
                            BB_IFFREE (tempadj, int *);
                            BB_IFFREE (tempadjspace, int);
                            BButil_freedatagroup (dat);
                            return 1;
                        }
                        for (i = 0; i < *ncount; i++) {
                            tempadj[i] = tempadjspace + i * (*ncount);
                            if (matrixform == MATRIX_UPPER_ROW) {
                                tempadj[i][i] = 0;
                                for (j = i + 1; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            } else if (matrixform == MATRIX_UPPER_DIAG_ROW) {
                                for (j = i; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            } else {
                                for (j = 0; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            }
                        }
                        for (i = 0; i < *ncount; i++) {
                            for (j = 0; j <= i; j++)
                                dat->adj[i][j] = tempadj[j][i];
                        }
                        BB_FREE (tempadjspace, int);
                        BB_FREE (tempadj, int *);
                    }
                } else {
                    fprintf (stderr, "ERROR: Matrix with norm %d?\n",
                             norm);
                    return 1;
                }
            } else if (!strcmp (key, "FIXED_EDGES_SECTION")) {
                fprintf (stderr, "ERROR: Not set up for fixed edges\n");
                return 1;
            }
        }
    }
    fclose (in);

    if (dat->x == (double *) NULL && dat->adj == (int **) NULL) {
        fprintf (stderr, "ERROR: Didn't find the data\n");
        return 1;
    } else {
        return 0;
    }
}

int BButil_datagroup_perm (int ncount, BBdatagroup *dat, int *perm)
{
    int i, j;

    if (dat->norm == BB_SPARSE) {
        fprintf (stderr, "perm not set up for BB_SPARSE\n");
        return 1;
    }
    
    if (dat->x != (double *) NULL) {
        double *tempx;

        tempx = BB_SAFE_MALLOC (ncount, double);
        if (!tempx)
            return 1;
        for (i = 0; i < ncount; i++) {
            tempx[i] = dat->x[perm[i]];
        }
        BB_FREE (dat->x, double);
        dat->x = tempx;
    }
    if (dat->y != (double *) NULL) {
        double *tempy;

        tempy = BB_SAFE_MALLOC (ncount, double);
        if (!tempy)
            return 1;
        for (i = 0; i < ncount; i++) {
            tempy[i] = dat->y[perm[i]];
        }
        BB_FREE (dat->y, double);
        dat->y = tempy;
    }
    if (dat->z != (double *) NULL) {
        double *tempz;

        tempz = BB_SAFE_MALLOC (ncount, double);
        if (!tempz)
            return 1;
        for (i = 0; i < ncount; i++) {
            tempz[i] = dat->z[perm[i]];
        }
        BB_FREE (dat->z, double);
        dat->z = tempz;
    }
    if (dat->adj != (int **) NULL) {
        int **tempadj = (int **) NULL;
        int *tempadjspace = (int *) NULL;

        tempadj = BB_SAFE_MALLOC (ncount, int *);
        tempadjspace = BB_SAFE_MALLOC (ncount * (ncount+1) / 2, int);
        
        if (tempadj == (int **) NULL ||
            tempadjspace == (int *) NULL) {
            BB_IFFREE (tempadj, int *);
            BB_IFFREE (tempadjspace, int);
            return 1;
        }

        for (i = 0, j = 0; i < ncount; i++) {
            tempadj[i] = tempadjspace + j;
            j += (i+1);
        }
        for (i = 0, j = 0; i < ncount; i++) {
            for (j = 0; j <= i; j++) {
                if (perm[i] <  perm[j])
                    tempadj[i][j] = dat->adj[perm[j]][perm[i]];
                else
                    tempadj[i][j] = dat->adj[perm[i]][perm[j]];
            }
        }
        BB_FREE (dat->adj, int *);
        BB_FREE (dat->adjspace, int);
        dat->adj = tempadj;
        dat->adjspace = tempadjspace;
    }
    return 0;
}

int BButil_graph2dat_matrix (int ncount, int ecount, int *elist, int *elen,
        int defaultlen, BBdatagroup *dat)
{
    int i;
    int j;
    int k;
    int rval;
    
    BButil_init_datagroup (dat);
    dat->adj = BB_SAFE_MALLOC (ncount, int *);
    dat->adjspace = BB_SAFE_MALLOC (ncount * (ncount+1) / 2, int);
    if (dat->adj == (int **) NULL ||
        dat->adjspace == (int *) NULL) {
        fprintf (stderr, "Our of memory in BButil_graph2dat\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i=0, j=0; i<ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }
    for (i=0, k=0; i<ncount; i++) {
        for (j=0; j<i; j++) {
            dat->adj[i][j] = defaultlen;
        }
        dat->adj[i][i] = 0;
    }
    for (i=0; i<ecount; i++) {
        j = elist[2*i];
        k = elist[2*i+1];
        if (j < k) dat->adj[k][j] = elen[i];
        else dat->adj[j][k] = elen[i];
    }
    if (BButil_dat_setnorm (dat, BB_MATRIXNORM)) {
        fprintf (stderr, "BButil_dat_setnorm failed\n");
        rval = 1;
        goto CLEANUP;
    }
    rval = 0;

CLEANUP:
    if (rval) {
        BButil_freedatagroup (dat);
    }
    return rval;
}

int BButil_print_command (int ac, char **av)
{
    int rval = 0;
    int i, cmdlen = 0;
    char *cmdout = (char *) NULL;

    for (i=0; i<ac; i++) {
        cmdlen += strlen(av[i]) + 1;
    }
    cmdout = BB_SAFE_MALLOC (cmdlen, char);
    BBcheck_NULL (cmdout, "out of memory in print_command");

    cmdlen = 0;
    for (i=0; i<ac; i++) {
        strcpy (cmdout + cmdlen, av[i]);
        cmdlen += strlen(av[i]);
        cmdout[cmdlen] = ' ';
        cmdlen++;
    }
    cmdout[cmdlen-1] = '\0';
    printf ("%s\n", cmdout); fflush (stdout);

CLEANUP:

    BB_IFFREE (cmdout, char);
    return rval;
}

int BButil_bix_getopt (int ac, char **av, const char *def, int *p_optind,
        char **p_optarg)
{
    int c;
    char *sp = av[*p_optind];
    char bwarn[2];

    if (*p_optind < 1 || *p_optind >= ac) {
        *p_optind = ac;
        return (EOF);
    }
    if ((int) *sp != (int) '-')
        return (EOF);
    if ((int) *(sp + 1) == (int) '-') {
        (*p_optind)++;
        return (EOF);
    }
    (av[*p_optind])++;
    sp++;
    while ((int) *sp != (int) *def && (int) *def != (int) '\0')
            def++;
    if ((int) *def == (int) '\0') {
        *p_optind = ac;
        bwarn[0] = *sp;                          /* Bico: February 8, 1995 */
        bwarn[1] = '\0';
        printf ("Illegal option: -%s\n", bwarn);
        return BB_BIX_GETOPT_UNKNOWN;
    }
    if ((int) *(def + 1) != (int) ':') {
        c = *sp;
        if ((int) *(sp + 1) != (int) '\0')
            *sp = '-';
        else
            (*p_optind)++;
        return (c);
    } else {
        if ((int) *(sp + 1) != (int) '\0') {
            *p_optarg = sp + 1;
            c = *sp;
            (*p_optind)++;
            return (c);
        } else if (*p_optind >= ac - 1) {
            *p_optind = ac;
            return (EOF);
        } else {
            *p_optarg = av[*p_optind + 1];
            c = *sp;
            *p_optind += 2;
            return (c);
        }
    }
}

void BButil_printlabel (void)
{
    char buf[1024];

    gethostname (buf, 1024);
    printf ("Host: %s  Current process id: %d\n", buf, (int) getpid());
    fflush (stdout);
}

