#include <stdio.h>
#include "verify.h"
#include "tsp.h"
#include "util.h"

#define VERIFY_STAT_OPEN     0
#define VERIFY_STAT_WORK     1
#define VERIFY_STAT_DONE     2
#define VERIFY_STAT_DONE_TSP 3
#define VERIFY_STAT_FAILED   4

static int save_remaining (CCtsp_lpcuts *pool, int *status, int ncount,
    char *probname);
static int receive_cut (CC_SFILE *s, int *status, CCtsp_lpcuts *pool,
    double *cumtime, int *result, const char *grunt);
static int send_cut (CC_SFILE *s, int *current, int *status, CCtsp_lpcuts *pool,
    int ncount);

int main (int ac, char **av)
{
    int rval = 0;
    int i, ncuts = 0, result;
    CC_SPORT *lport = (CC_SPORT *) NULL;
    CC_SFILE *s;
    int *status = (int *) NULL;
    char request;
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) NULL;
    char *probname = (char *) NULL;
    char grunt[1024];
    int ncount = 0;
    int donecount = 0;
    int current = 0;
    int hkcount = 0;
    int bccount = 0;
    int xxcount = 0;
    double cumtime = 0.0;

    if (ac != 2) {
        fprintf (stderr, "Usage: %s pool_file (nonsimple)\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_print_command (ac, av);
    CCcheck_rval (rval, "CCutil_print_command failed");

    CCutil_printlabel ();
    probname = CCtsp_problabel (av[1]);
    CCcheck_NULL (probname, "CCtsp_problabel failed");

    rval = CCtsp_init_cutpool (&ncount, av[1], &pool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    printf ("Total %d cuts\n", pool->cutcount);
    fflush (stdout);

    ncuts = pool->cutcount;
    status = CC_SAFE_MALLOC (ncuts, int);
    CCcheck_NULL (status, "out of memory for status array");

    for (i = 0; i < ncuts; i++) status[i] = VERIFY_STAT_OPEN;

/*
    int possible_count = 0;
    for (i = 0; i < ncuts; i++) {
        if (pool->cuts[i].rhs == CCtsp_COMBRHS (&pool->cuts[i])) {
            possible_count++;
        }
    }
    printf ("%d cuts have comb-rhs\n", possible_count); fflush (stdout);
    exit (1);
*/

/*
    {
        CCtsp_lpcut_in c;
        for (i = 0; i < ncuts; i++) {
            printf ("Cut %d\n\n", i);
            CCtsp_init_lpcut_in (&c);
            rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[i], &c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed"); 
            CCtsp_print_lpcut_in (&c);
            CCtsp_free_lpcut_in (&c);
            printf ("\n--------------------------------------------------\n\n");
        }
        exit (1);
    }
*/

/*
    {
        int total = 0;
        int j, k, hit, bad;
        char *hits = (char *) NULL;
        int *ccounts = (int *) NULL;
        int *tcounts = (int *) NULL;

        hits = CC_SAFE_MALLOC (ncuts, char);
        CCcheck_NULL (hits, "out of memory for hits");
        for (i = 0; i < ncuts; i++) hits[i] = 0;

        ccounts = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (ccounts, "out of memory for ccounts");
        tcounts = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (tcounts, "out of memory for tcounts");

        for (i = 0; i < ncuts; i++) {
            if (!hits[i]) {
                hit = 0;
                for (k = 0; k < pool->cuts[i].cliquecount; k++) {
                    CCtsp_clique_count (
                       &pool->cliques[pool->cuts[i].cliques[k]], &ccounts[k]);
                }
                CCutil_int_array_quicksort (ccounts, pool->cuts[i].cliquecount);
                for (j = i+1; j < ncuts; j++) {
                    if ((pool->cuts[i].skel.atomcount == 
                         pool->cuts[j].skel.atomcount) &&
                        (pool->cuts[i].rhs ==
                         pool->cuts[j].rhs) &&
                        (pool->cuts[i].cliquecount ==
                         pool->cuts[j].cliquecount)) {
                        for (k = 0; k < pool->cuts[j].cliquecount; k++) {
                            CCtsp_clique_count (
                                &pool->cliques[pool->cuts[j].cliques[k]],
                                &tcounts[k]);
                        }
                        CCutil_int_array_quicksort (tcounts,
                                                    pool->cuts[j].cliquecount);
                        bad = 0;
                        for (k = 0; k < pool->cuts[j].cliquecount; k++) {
                            if (ccounts[k] != tcounts[k]) {
                                bad = 1; break;
                            }
                        }
                        
                        if (!bad) {
                            hit++;
                            hits[j] = 1;
                            printf ("Hit-%d ", j);
                        }
                    }
                }
                if (hit) {
                    printf ("\n");
                    CCtsp_lpcut_in c;
                    total += hit;
                    printf ("Cut[%d] = %d, total = %d (%d atoms, %d cliques, %d rhs)\n",
                        i, hit, total, pool->cuts[i].skel.atomcount,
                        pool->cuts[i].cliquecount, pool->cuts[i].rhs);
                    fflush (stdout);

                    CCtsp_init_lpcut_in (&c);
                    rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[i], &c);
                    CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed"); 
                    CCtsp_print_lpcut_in (&c);
                    CCtsp_free_lpcut_in (&c);
                }
            }
        }
        printf ("TOTAL HITS: %d\n", total);
        exit (1);
    }
*/

/*
    {
        int total = 0;
        int j, hit, diff;
        char *hits = (char *) NULL;

        hits = CC_SAFE_MALLOC (ncuts, char);
        CCcheck_NULL (hits, "out of memory for hits");

        for (i = 0; i < ncuts; i++) hits[i] = 0;

        for (i = 0; i < ncuts; i++) {
            if (!hits[i]) {
                hit = 0;
                for (j = i+1; j < ncuts; j++) {
                    CCtsp_compare_skeletons (&pool->cuts[i].skel,
                                             &pool->cuts[j].skel, &diff);
                    if (!diff) {
                        hit++;
                        hits[j] = 1;
                    }
                }
                if (hit) {
                    total += hit;
                    printf ("Cut[%d] = %d, total = %d\n", i, hit, total);
                    fflush (stdout);
                }
            }
        }
        printf ("TOTAL HITS: %d\n", total);
        exit (1);
    }
*/

/*
    {
        int buk[1000];

        for (i = 0; i < 1000; i++) buk[i] = 0;
        for (i = 0; i < ncuts; i++) {
            buk[pool->cuts[i].skel.atomcount]++;;
        }
        for (i = 0; i < 1000; i++) {
            if (buk[i]) {
                printf ("%d: %d\n", i, buk[i]); fflush (stdout);
            }
        }
        exit (1);
    }
*/

/*
#define BUK_CNT 100000
    {
        int buk[BUK_CNT];
        int j, cnt;
        CCtsp_lpcut_in c;

        for (i = 0; i < BUK_CNT; i++) buk[i] = 0;

        for (i = 0; i < ncuts; i++) {
            CCtsp_init_lpcut_in (&c);
            rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[i], &c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed"); 

            for (j = 0; j < c.cliquecount; j++) {
                cnt = 0;
                CCtsp_clique_count (&c.cliques[j], &cnt);
                if (cnt < 0 || cnt >= BUK_CNT) {
                    printf ("WHAT: %d\n", cnt); fflush (stdout);
                    exit (1);
                }
                buk[cnt]++;;
            }

            CCtsp_free_lpcut_in (&c);
        }
        for (i = 0; i < 20; i++) {
            if (buk[i]) {
                printf ("Clique-size %d: %d\n", i, buk[i]); fflush (stdout);
            }
        }
        exit (1);
    }
*/

/*
#define BUK_CNT 100000
    {
        int buk[BUK_CNT];
        int j, k, cnt;
        CCtsp_lpcut_in c;

        for (i = 0; i < BUK_CNT; i++) buk[i] = 0;

        for (i = 0; i < ncuts; i++) {
            CCtsp_init_lpcut_in (&c);
            rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[i], &c);
            CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed"); 

            j = CCverify_cut (&c, CC_TYPE_OTHER, &k, 0, (int *) NULL,
                            (char *) NULL);
            buk[j]++;

            CCtsp_free_lpcut_in (&c);
        }
        for (i = 0; i < 100; i++) {
            if (buk[i]) {
                printf ("Clique-size %d: %d\n", i, buk[i]); fflush (stdout);
            }
        }
        exit (1);
    }
*/


    printf ("\nBEGINNING CUT-VERIFICATION NET PROCESSING\n\n");
    fflush (stdout);

    lport = CCutil_snet_listen (CCtsp_VERIFY_PORT);
    if (lport == (CC_SPORT *) NULL) {                                           
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }

    do {
        s = CCutil_snet_receive (lport);
        if (!s) {
            fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
            continue;
        }

        if (CCutil_sread_string (s, grunt, 1023)) {
            fprintf (stderr, "CCutil_sread_char string, abort con\n");
            CCutil_sclose (s);
            continue;
        }

        if (CCutil_sread_char (s, &request)) {
            fprintf (stderr, "CCutil_sread_char failed, abort con\n");
            CCutil_sclose (s);
            continue;
        }

        switch (request) {
        case CCtsp_VERIFY_SEND:
            if (current < ncuts) {
                rval = send_cut (s, &current, status, pool, ncount);
                if (rval) {
                    fprintf (stderr, "send_cut failed - abort con\n");
                }
                printf ("Sent cut %d to %s (remaining %d)\n",
                         current, grunt, ncuts - donecount);
                fflush (stdout);
            } else {
                rval = CCutil_swrite_char (s, CCtsp_VERIFY_NO);
                if (rval) {
                    fprintf (stderr, "VERIFY_NO write failed - abort con\n");
                }
            }
            CCutil_sclose (s);
            break;
        case CCtsp_VERIFY_RECEIVE:
            rval = receive_cut (s, status, pool, &cumtime, &result, grunt);
            if (rval) {
                fprintf (stderr, "receive_cut failed - abort connection\n");
                CCutil_sclose (s);
                break;
            } else {
                CCutil_sclose (s);
                donecount++;
                switch (result) {
                case VERIFY_STAT_DONE:     hkcount++; break;
                case VERIFY_STAT_DONE_TSP: bccount++; break;
                case VERIFY_STAT_FAILED:   xxcount++; break;
                default:  
                    fprintf (stderr, "Unknown return status\n");
                    rval = 1;  goto CLEANUP;
                }
                printf ("  %d with HK, %d with BC, %d FAILURES, %d remaining\n",
                     hkcount, bccount, xxcount, ncuts - donecount);
                break;
            }
        case CCtsp_VERIFY_SAVE:
            CCutil_sclose (s);
            rval = save_remaining (pool, status, ncount, probname);
            CCcheck_rval (rval, "save_remaining failed");
            break;
        case CCtsp_VERIFY_EXIT:
            printf ("Shutting down the verification boss\n"); fflush (stdout);
            CCutil_sclose (s);
            goto CLEANUP;
        default:
            fprintf (stderr, "Invalid request %c\n", request);
        }
    } while (donecount != ncuts);

    printf ("Completed %d cuts in %.2f seconds\n", ncuts, cumtime);
    fflush (stdout);
    printf ("%d cuts verified with HK\n", hkcount);
    printf ("%d cuts verified with with BFS-tsp\n", bccount);
    printf ("%d cuts failed verification\n", xxcount);
    fflush (stdout);

CLEANUP:

    if (donecount > 0 && (donecount < ncuts || xxcount > 0)) {
        rval = save_remaining (pool, status, ncount, probname);
        CCcheck_rval (rval, "save_remaining failed");
    }

    CCutil_snet_unlisten (lport);
    CC_IFFREE (status, int);
    CC_IFFREE (probname, char);
    if (pool) CCtsp_free_cutpool (&pool);

    return rval;
}

static int save_remaining (CCtsp_lpcuts *pool, int *status, int ncount,
    char *probname)
{
    int rval = 0;
    int i;
    CCtsp_lpcuts *newpool = (CCtsp_lpcuts *) NULL;
    char buf[1024];

    rval = CCtsp_init_cutpool (&ncount, (char *) NULL, &newpool);
    CCcheck_rval (rval, "CCtsp_init_cutpool failed");

    for (i = 0; i < pool->cutcount; i++) {
        if (status[i] != VERIFY_STAT_DONE && 
            status[i] != VERIFY_STAT_DONE_TSP) {
            rval = CCtsp_add_to_cutpool (newpool, pool, &pool->cuts[i]);
            CCcheck_rval (rval, "CCtsp_add_to_cutpool failed");
        }
    }

    sprintf (buf, "%s.remain", probname);
    printf ("Saving the remaining %d cuts to %s\n", newpool->cutcount, buf);
    fflush (stdout);

    rval = CCtsp_write_cutpool (ncount, buf, newpool);
    CCcheck_rval (rval, "CCtsp_write_cutpool failed");

CLEANUP:

    if (newpool) CCtsp_free_cutpool (&newpool);
    return rval;
}

static int receive_cut (CC_SFILE *s, int *status, CCtsp_lpcuts *pool,
        double *cumtime, int *result, const char *grunt)
{
    int rval = 0;
    int id, check;
    double rtime;
    char buf[1024];

    *result = -1;

    rval = CCutil_sread_int (s, &id);
    CCcheck_rval (rval, "CCutil_sread_int failed (id)");

    rval = CCutil_sread_double (s, &rtime);
    CCcheck_rval (rval, "CCutil_sread_double failed (rtime)");

    (*cumtime) += rtime;

    if (status[id] != VERIFY_STAT_WORK) {
        fprintf (stderr, "ERROR: Received cut that is not at work\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_sread_int (s, &check);
    CCcheck_rval (rval, "CCutil_sread_int failed (check)");

    if (check == 0) {
        *result = VERIFY_STAT_DONE;
        status[id] = *result; 
        printf ("Verified Cut %d (HK on %s), %.2f seconds, total %.2f\n",
                  id, grunt, rtime, *cumtime);
        fflush (stdout);
    } else if (check == 1) {
        *result = VERIFY_STAT_DONE_TSP;
        status[id] = *result;
        printf ("Verified Cut %d (BC on %s), %.2f seconds, total %.2f\n",
                  id, grunt, rtime, *cumtime);
        fflush (stdout);
    } else {
        CCtsp_lpcut_in c;

        *result = VERIFY_STAT_FAILED;
        status[id] = *result;

        fprintf (stderr, "FAILURE ON %s\n", grunt);
        fprintf (stderr, "Verification of cut %d failed, %.2f seconds\n",
                          id, rtime);
        sprintf (buf, "v%d.tsp", id);

        CCtsp_init_lpcut_in (&c);
        rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[id], &c);
        CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed"); 

        rval = CCverify_dump_tsp (&c, buf);
        CCcheck_rval (rval, "CCverify_dump_tsp failed");
        CCtsp_free_lpcut_in (&c);
    }

CLEANUP:

    return rval;
}

static int send_cut (CC_SFILE *s, int *current, int *status,
       CCtsp_lpcuts *pool, int ncount)
{
    int rval = 0;
    CCtsp_lpcut_in c;

    CCtsp_init_lpcut_in (&c);

    rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[*current], &c);
    CCcheck_rval (rval, "CCtsp_lpcut_to_lpcut_in failed");

    rval = CCutil_swrite_char (s, CCtsp_VERIFY_YES);
    CCcheck_rval (rval, "CCutil_swrite_int failed (YES)");

    rval = CCutil_swrite_int (s, *current);
    CCcheck_rval (rval, "CCutil_swrite_int failed (cut id)");

    rval = CCutil_swrite_int (s, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed (ncount)");

    rval = CCtsp_write_lpcut_in (s, &c, ncount);
    CCcheck_rval (rval, "CCtsp_write_lpcut_in failed\n");

    status[*current] = VERIFY_STAT_WORK;

    (*current)++;

CLEANUP:

    CCtsp_free_lpcut_in (&c);
    return rval;
}
