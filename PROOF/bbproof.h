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

#ifndef  __BBPROOF_H
#define  __BBPROOF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
extern int gethostname (char *, int);

#define BB_SAFE_MALLOC(nnum,type)                                          \
    (type *) BButil_allocrus (((size_t) (nnum)) * sizeof (type))

#define BB_FREE(object,type) {                                             \
    BButil_freerus ((void *) (object));                                    \
    object = (type *) NULL;                                                \
}

#define BB_IFFREE(object,type) {                                           \
    if ((object)) BB_FREE ((object),type);                                 \
}

#define BBcheck_rval(rval,msg) {                                           \
    if ((rval)) {                                                          \
        fprintf (stderr, "%s\n", (msg));                                   \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define BBcheck_NULL(item,msg) {                                           \
    if ((!item)) {                                                         \
        fprintf (stderr, "%s\n", (msg));                                   \
        rval = 1;                                                          \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define BB_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define BB_SBUFFER_SIZE (4000)
#define BB_SFNAME_SIZE (32)
#define BB_BIX_GETOPT_UNKNOWN -3038
#define BBtsp_PROB_FILE_NAME_LEN 128

#define BB_KD_NORM_TYPE    128            /* Kdtrees work      */
#define BB_X_NORM_TYPE     256            /* Old nearest works */
#define BB_JUNK_NORM_TYPE  512            /* Nothing works     */

#define BB_D2_NORM_SIZE      1024         /* x,y coordinates   */
#define BB_D3_NORM_SIZE      2048         /* x,y,z coordinates */
#define BB_MATRIX_NORM_SIZE  4096         /* adj matrix        */

#define BB_NORM_BITS      (BB_KD_NORM_TYPE | BB_X_NORM_TYPE | \
                           BB_JUNK_NORM_TYPE)
#define BB_NORM_SIZE_BITS (BB_D2_NORM_SIZE | BB_D3_NORM_SIZE | \
                           BB_MATRIX_NORM_SIZE)

#define BB_MAXNORM        (0 |   BB_KD_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_EUCLIDEAN_CEIL (1 |   BB_KD_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_EUCLIDEAN      (2 |   BB_KD_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_EUCLIDEAN_3D   (3 |    BB_X_NORM_TYPE |     BB_D3_NORM_SIZE)
#define BB_USER           (4 | BB_JUNK_NORM_TYPE |                   0)
#define BB_ATT            (5 |    BB_X_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_GEOGRAPHIC     (6 |    BB_X_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_MATRIXNORM     (7 | BB_JUNK_NORM_TYPE | BB_MATRIX_NORM_SIZE)
#define BB_DSJRANDNORM    (8 | BB_JUNK_NORM_TYPE |                   0)
#define BB_CRYSTAL        (9 |    BB_X_NORM_TYPE |     BB_D3_NORM_SIZE)
#define BB_SPARSE        (10 | BB_JUNK_NORM_TYPE |                   0)
#define BB_EUCTOROIDAL   (16 | BB_JUNK_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_GEOM          (17 |    BB_X_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_MANNORM       (18 |   BB_KD_NORM_TYPE |     BB_D2_NORM_SIZE)
#define BB_SUBDIVISION   (99 | BB_JUNK_NORM_TYPE |                   0)

#define BB_GEOGRAPHIC_SCALE (6378.388 * 3.14 / 180.0)    /*  see edgelen.c  */
#define BB_GEOM_SCALE (6378388.0 * 3.14 / 180.0)         /*  see edgelen.c  */
#define BB_ATT_SCALE (.31622)                            /*  sqrt(1/10)     */

/* For X-NORMS, scales are such that |x[i] - x[j]| * scale <= edgelen(i,j). */
/* Geographic is slightly off, since the fractional part of x[i] is really  */
/* really minutes, not fractional degrees.                                  */

typedef struct BB_SFILE {
    int           status;
    int           desc;
    int           type;
    int           chars_in_buffer;
    int           current_buffer_char;     /* only used for reading */
    int           bits_in_last_char;       /* writing: number of empty bits in
                                            * buffer[chars_in_buffer];
                                            * reading: number of full bits in
                                            * buffer[?] */
    int           pos;
    char          fname[BB_SFNAME_SIZE];
    char          hname[BB_SFNAME_SIZE];
    unsigned char buffer[BB_SBUFFER_SIZE];
} BB_SFILE;

typedef struct BBdatagroup {
    int    (*edgelen) (int i, int j, struct BBdatagroup *dat);
    double  *x;
    double  *y;
    double  *z;
    int    **adj;
    int     *adjspace;
    int    **len;
    int     *lenspace;
    int     *degree;
    int      norm;
    int      default_len;     /* for edges not in sparse graph   */
    int      sparse_ecount;   /* number of edges in sparse graph */
} BBdatagroup;

typedef struct BB_SPORT {
    unsigned short port;
    int t;
} BB_SPORT;

typedef struct BBbigguy {
    unsigned short ihi;
    unsigned short ilo;
    unsigned short fhi;
    unsigned short flo;
} BBbigguy;

typedef struct BBtsp_bigdual {
    int           cutcount;
    BBbigguy     *node_pi;
    BBbigguy     *cut_pi;
} BBtsp_bigdual;

typedef struct BBtsp_lpnode {
    int                 deg;
    int                 mark;
    struct BBtsp_lpadj *adj;
} BBtsp_lpnode;

typedef struct BBtsp_lpedge {
    int       ends[2];   /* ends[0] should always be < ends[1] */
    int       fixed;
    int       branch;    /* < 0 means set to 0 and > 0 means set to 1 */
    int       len;
    int       coef;      /* should be maintained at zero */
    int       coefnext;  /* should be maintained at -2 */
} BBtsp_lpedge;

typedef struct BBtsp_lpadj {
    int       to;
    int       edge;
} BBtsp_lpadj;

typedef struct BBtsp_lpgraph {
    int              ncount;
    int              espace;
    int              ecount;
    int              nodemarker;
    BBtsp_lpnode    *nodes;
    BBtsp_lpedge    *edges;
    BBtsp_lpadj     *adjspace;
    int              adjstart;
    int              adjend;
} BBtsp_lpgraph;

typedef struct BBtsp_predge {
    int        ends[2];
    int        len;
    double     rc;
} BBtsp_predge;

typedef struct BBtsp_segment {
    int lo;
    int hi;
} BBtsp_segment;

typedef struct BBtsp_lpclique {
    int                   segcount;
    struct BBtsp_segment *nodes;
} BBtsp_lpclique;

typedef struct BBtsp_lpdomino {
    BBtsp_lpclique        sets[2];
} BBtsp_lpdomino;

#define BB_FOREACH_NODE_IN_CLIQUE(i,c,tmp) \
    for(tmp=0;tmp<(c).segcount;tmp++) \
        for(i=(c).nodes[tmp].lo;i<=(c).nodes[tmp].hi;i++)

typedef struct BBtsp_skeleton {
    int  atomcount;
    int *atoms;
} BBtsp_skeleton;

typedef struct BBtsp_lpcut {
    int                   cliquecount;
    int                   dominocount;
    int                   age;
    int                   rhs;
    char                  sense;
    char                  branch;
    int                  *cliques;
    int                  *dominos;
    BBtsp_skeleton        skel;
} BBtsp_lpcut;

typedef struct BBtsp_lpcut_in {
    int                    cliquecount;
    int                    dominocount;
    int                    rhs;
    char                   sense;
    char                   branch;
    BBtsp_lpclique        *cliques;
    BBtsp_lpdomino        *dominos;
    BBtsp_skeleton         skel;
} BBtsp_lpcut_in;

typedef struct BBtsp_lpcuts {
    int             cutcount;
    int             savecount;
    int             cliqueend;
    int             cutspace;
    int             cliquespace;
    BBtsp_lpcut    *cuts;
    BBtsp_lpclique *cliques;
    int             dominoend;
    int             dominospace;
    BBtsp_lpdomino *dominos;
} BBtsp_lpcuts;

typedef struct BBtsp_branchobj {
    int             depth;
    int             rhs;
    int             ends[2];
    char            sense;
    BBtsp_lpclique *clique;
} BBtsp_branchobj;

typedef struct BBproofnode {
    int number;
    int parent;
    int side;
    int child0;
    int child1;
    BBtsp_branchobj *branch;
    BBtsp_bigdual *exact_dual;
} BBproofnode;

typedef struct BBcutproof  {
    int ncount;
    int bbcount;
    int *tour;
    char *probname;
    BBproofnode *lproof;
    BBtsp_lpcuts *cuts;
    BBbigguy lpbound;
} BBcutproof;

typedef struct BBcuttype {
    int class;
    int isomorph;
    int hk;
    BBcutproof *proof;
} BBcuttype;


/******************** bbbsafe.c ********************/

BB_SFILE
    *BBsafe_sopen (const char *f, const char *s),
    *BBsafe_snet_open (const char *hname, unsigned short p),
    *BBsafe_snet_receive (BB_SPORT *s);

BB_SPORT
    *BBsafe_snet_listen (unsigned short p);

int
    BBsafe_swrite (BB_SFILE *f, char *buf, int size),
    BBsafe_swrite_bits (BB_SFILE *f, int x, int xbits),
    BBsafe_swrite_ubits (BB_SFILE *f, unsigned int x, int xbits),
    BBsafe_swrite_char (BB_SFILE *f, char x),
    BBsafe_swrite_string (BB_SFILE *f, const char *x),
    BBsafe_swrite_short (BB_SFILE *f, short x),
    BBsafe_swrite_ushort (BB_SFILE *f, unsigned short x),
    BBsafe_swrite_int (BB_SFILE *f, int x),
    BBsafe_swrite_uint (BB_SFILE *f, unsigned int x),
    BBsafe_swrite_double (BB_SFILE *f, double x),
    BBsafe_sread (BB_SFILE *f, char *buf, int size),
    BBsafe_sread_bits (BB_SFILE *f, int *x, int xbits),
    BBsafe_sread_ubits (BB_SFILE *f, unsigned int *x, int xbits),
    BBsafe_sread_char (BB_SFILE *f, char *x),
    BBsafe_sread_string (BB_SFILE *f, char *x, int maxlen),
    BBsafe_sread_short (BB_SFILE *f, short *x),
    BBsafe_sread_ushort (BB_SFILE *f, unsigned short *x),
    BBsafe_sread_int (BB_SFILE *f, int *x),
    BBsafe_sread_uint (BB_SFILE *f, unsigned int *x),
    BBsafe_sread_double (BB_SFILE *f, double *x),
    BBsafe_sflush (BB_SFILE *f),
    BBsafe_sclose (BB_SFILE *f),
    BBsafe_sbits (unsigned int x);

void
    BBsafe_snet_unlisten (BB_SPORT *s);


/******************** bbbigguy.c ********************/

extern const BBbigguy BBbigguy_MINBIGGUY;
extern const BBbigguy BBbigguy_MAXBIGGUY;
extern const BBbigguy BBbigguy_ZERO;
extern const BBbigguy BBbigguy_ONE;

void
    BBbigguy_addmult (BBbigguy *x, BBbigguy y, int m);

int
    BBbigguy_cmp (BBbigguy x, BBbigguy y),
    BBbigguy_swrite (BB_SFILE *f, BBbigguy x),
    BBbigguy_sread (BB_SFILE *f, BBbigguy *x);

double
    BBbigguy_bigguytod (BBbigguy x);

BBbigguy
    BBbigguy_itobigguy (int d),
    BBbigguy_dtobigguy (double d),
    BBbigguy_ceil (BBbigguy x);

#define BBbigguy_add(x,y) (BBbigguy_addmult(x,y,1))
#define BBbigguy_sub(x,y) (BBbigguy_addmult(x,y,-1))


/******************** bbio.c ********************/

void
    BBio_init_proofnode (BBproofnode *p),
    BBio_free_proofnode (BBproofnode *p),
    BBio_init_cutproof (BBcutproof *p),
    BBio_free_cutproof (BBcutproof *p);

int
    BBio_read_short_proof (BBcutproof **pp, BBcuttype **pcutlist, char *fname),
    BBio_read_lpcut_in (BB_SFILE *f, BBtsp_lpcut_in *c, int ncount),
    BBio_read_lpclique (BB_SFILE *f, BBtsp_lpclique *c, int ncount),
    BBio_read_lpdomino (BB_SFILE *f, BBtsp_lpdomino *d, int ncount),
    BBio_read_skeleton (BB_SFILE *f, BBtsp_skeleton *skel, int ncount),
    BBio_write_lpcut_in (BB_SFILE *f, BBtsp_lpcut_in *c, int ncount),
    BBio_write_lpclique (BB_SFILE *f, BBtsp_lpclique *c, int ncount),
    BBio_write_lpdomino (BB_SFILE *f, BBtsp_lpdomino *c, int ncount),
    BBio_write_skeleton (BB_SFILE *f, BBtsp_skeleton *skel, int ncount);


/******************** bbutil.c ********************/

void
   *BButil_allocrus (size_t size),
    BButil_freerus (void *p),
    BButil_int_array_quicksort (int *len, int n),
    BButil_int_perm_quicksort (int *perm, int *len, int n),
    BButil_printlabel (void),
    BButil_dat_getnorm (BBdatagroup *dat, int *norm),
    BButil_init_datagroup (BBdatagroup *dat),
    BButil_freedatagroup (BBdatagroup *dat),
    BButil_init_tsp_lpcuts_struct (BBtsp_lpcuts *c),
    BButil_init_lpgraph_struct (BBtsp_lpgraph *g),
    BButil_free_lpgraph (BBtsp_lpgraph *g),
    BButil_lpclique_compare (BBtsp_lpclique *a, BBtsp_lpclique *b, int *diff),
    BButil_init_lpcut (BBtsp_lpcut *c),
    BButil_init_lpcut_in (BBtsp_lpcut_in *c),
    BButil_init_lpclique (BBtsp_lpclique *c),
    BButil_free_lpclique (BBtsp_lpclique *c),
    BButil_init_branchobj (BBtsp_branchobj *b),
    BButil_free_branchobj (BBtsp_branchobj *b),
    BButil_free_lpcut_in (BBtsp_lpcut_in *c),
    BButil_init_lpdomino (BBtsp_lpdomino *c),
    BButil_free_lpdomino (BBtsp_lpdomino *c),
    BButil_init_skeleton (BBtsp_skeleton *skel),
    BButil_free_skeleton (BBtsp_skeleton *skel),
    BButil_free_bigdual (BBtsp_bigdual **d);

int
    BButil_gettsplib (char *datname, int *ncount, BBdatagroup *dat),
    BButil_datagroup_perm (int ncount, BBdatagroup *dat, int *perm),
    BButil_graph2dat_matrix (int ncount, int ecount, int *elist, int *elen,
        int defaultlen, BBdatagroup *dat),
    BButil_dat_edgelen (int i, int j, BBdatagroup *dat),
    BButil_dat_setnorm (BBdatagroup *dat, int norm),
    BButil_array_to_lpclique (int *ar, int acount, BBtsp_lpclique *cliq),
    BButil_clique_to_array (BBtsp_lpclique *c, int **ar, int *count),
    BButil_copy_lpclique (BBtsp_lpclique *c, BBtsp_lpclique *new),
    BButil_copy_lpdomino (BBtsp_lpdomino *c, BBtsp_lpdomino *new),
    BButil_lpcut_to_lpcut_in (BBtsp_lpcuts *cuts, BBtsp_lpcut *c,
        BBtsp_lpcut_in *new),
    BButil_construct_skeleton (BBtsp_lpcut_in *c, int nodecount),
    BButil_copy_skeleton (BBtsp_skeleton *old, BBtsp_skeleton *new),
    BButil_print_command (int ac, char **av),
    BButil_bix_getopt (int argc, char **argv, const char *def, int *p_optind,
        char **p_optarg);

double
    BButil_zeit (void);


/******************** bbprice.c ********************/

int
    BBprice_price (BBtsp_lpcuts *cuts, BBtsp_bigdual *exact_dual,
        BBdatagroup *dat, int ncount, int *elist, int ecount, int *fixlist,
        int fcount, BBbigguy *bound, int silent),
    BBprice_elim (BBtsp_lpcuts *cuts, BBtsp_bigdual *exact_dual,
        BBdatagroup *dat, double upperbound, BBbigguy exact_lowerbound,
        int ncount, int **elist, int *ecount, int **fixlist, int *fcount,
        int silent);


/******************** bbcuts.c ********************/

#define BB_NOCLASS 99999
#define BB_TYPE_SUBTOUR 1
#define BB_TYPE_COMB 2
#define BB_TYPE_STAR 4
#define BB_TYPE_BIPARTITION 8
#define BB_TYPE_OTHER 16
#define BB_TYPE_DOMINO 64 

#define BBtsp_VERIFY_WORK        'A'
#define BBtsp_VERIFY_NO          'N'
#define BBtsp_VERIFY_RECEIVE     'R'
#define BBtsp_VERIFY_SEND        'S'
#define BBtsp_VERIFY_SAVE        'T'
#define BBtsp_VERIFY_WAIT        'W'
#define BBtsp_VERIFY_YES         'Y'
#define BBtsp_VERIFY_EXIT        'X'
#define BBtsp_VERIFY_PORT ((unsigned short) 24870)

int
    BBcuts_verify (BBtsp_lpcut_in *cut, int ncount, int check_types, int *type,
        BBdatagroup *getdat),
    BBcuts_backbone (BBtsp_lpcut_in *cut, BBtsp_lpcut_in *new, int *no_outside);


/******************** bbheld.c ********************/

#define BB_HELDKARP_ERROR               -1
#define BB_HELDKARP_SEARCHLIMITEXCEEDED  1

int
    BBheld_small_elist (int ncount, int ecount, int *elist, int *elen,
        double *upbound, double *optval, int *foundtour, int anytour,
        int *tour_elist, int nodelimit, int silent);

#endif /* __BBPROOF_H */
