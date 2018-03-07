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
/*  Boss for bbgrunt.  Reads proof files and verifies those cuts that need  */
/*  to be checked with the Held-Karp solver.                                */
/****************************************************************************/

#include "bbproof.h"

#define VERIFY_STAT_OPEN     0
#define VERIFY_STAT_WORK     1
#define VERIFY_STAT_DONE     2
#define VERIFY_STAT_FAILED   4

static int
    receive_cut (BB_SFILE *s, int *status, double *cumtime, int *result,
        const char *grunt),
    send_cut (BB_SFILE *s, int current, int *status, BBtsp_lpcuts *pool,
        int ncount, int *other);

int main (int ac, char **av)
{
    int i, result, nother, rval = 0;
    int donecount = 0, current = 0, hkcount = 0, xxcount = 0;
    int *status = (int *) NULL, *other = (int *) NULL;
    char request, grunt[1024];
    double cumtime = 0.0;
    BB_SPORT *lport = (BB_SPORT *) NULL;
    BB_SFILE *s;
    BBcutproof *p = (BBcutproof *) NULL;
    BBcuttype *cutlist = (BBcuttype *) NULL;

    if (ac != 2) {
        fprintf (stderr, "Usage: %s proof_file\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    rval = BButil_print_command (ac, av);
    BBcheck_rval (rval, "BButil_print_command failed");
    BButil_printlabel ();

    rval = BBio_read_short_proof (&p, &cutlist, av[1]);
    BBcheck_rval (rval, "read_short_proof failed");

    printf ("Total Cuts: %d\n", p->cuts->cutcount); fflush (stdout);

    nother = 0;
    for (i = 0; i < p->cuts->cutcount; i++) {
        if (cutlist[i].class == BB_NOCLASS && cutlist[i].isomorph == -1 &&
                cutlist[i].hk == 1) {
            nother++;
        }
    }
    if (nother == 0) {
        printf ("No cuts need to be verified\n"); fflush (stdout);
        goto CLEANUP;
    }

    printf ("Cuts to be verified with Held-Karp: %d\n", nother);
    fflush (stdout);
    other = BB_SAFE_MALLOC (nother, int);
    BBcheck_NULL (other, "out of memory for other");

    nother = 0;
    for (i = 0; i < p->cuts->cutcount; i++) {
        if (cutlist[i].class == BB_NOCLASS && cutlist[i].isomorph == -1 &&
                cutlist[i].hk == 1) {
            other[nother++] = i;
        }
    }
    
    status = BB_SAFE_MALLOC (nother, int);
    BBcheck_NULL (status, "out of memory for status array");
    for (i = 0; i < nother; i++) status[i] = VERIFY_STAT_OPEN;

    printf ("\nBEGINNING CUT-VERIFICATION NET PROCESSING\n\n");
    fflush (stdout);

    lport = BBsafe_snet_listen (BBtsp_VERIFY_PORT);
    if (lport == (BB_SPORT *) NULL) {                                           
        fprintf (stderr, "BBsafe_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }

    do {
        s = BBsafe_snet_receive (lport);
        if (!s) {
            fprintf (stderr, "BBsafe_snet_receive failed, ignoring\n");
            continue;
        }

        if (BBsafe_sread_string (s, grunt, 1023)) {
            fprintf (stderr, "BBsafe_sread_char string, abort con\n");
            BBsafe_sclose (s);
            continue;
        }

        if (BBsafe_sread_char (s, &request)) {
            fprintf (stderr, "BBsafe_sread_char failed, abort con\n");
            BBsafe_sclose (s);
            continue;
        }

        switch (request) {
        case BBtsp_VERIFY_SEND:
            if (current < nother) {
                rval = send_cut (s, current, status, p->cuts, p->ncount, other);
                if (rval) {
                    fprintf (stderr, "send_cut failed - abort con\n");
                } else {
                    printf ("Sent cut %d (id %d) to %s (remaining %d)\n",
                         current, other[current], grunt, nother - donecount);
                    fflush (stdout);
                    current++;
                }
            } else {
                rval = BBsafe_swrite_char (s, BBtsp_VERIFY_NO);
                if (rval) {
                    fprintf (stderr, "VERIFY_NO write failed - abort con\n");
                }
            }
            BBsafe_sclose (s);
            break;
        case BBtsp_VERIFY_RECEIVE:
            rval = receive_cut (s, status, &cumtime, &result, grunt);
            if (rval) {
                fprintf (stderr, "receive_cut failed - abort connection\n");
                BBsafe_sclose (s);
                break;
            } else {
                BBsafe_sclose (s);
                donecount++;
                switch (result) {
                case VERIFY_STAT_DONE:     hkcount++; break;
                case VERIFY_STAT_FAILED:   xxcount++; break;
                default:  
                    fprintf (stderr, "Unknown return status\n");
                    rval = 1;  goto CLEANUP;
                }
                printf ("  %d with HK, %d FAILURES, %d remaining\n",
                     hkcount, xxcount, nother - donecount);
                break;
            }
        case BBtsp_VERIFY_EXIT:
            printf ("Shutting down the verification boss\n"); fflush (stdout);
            BBsafe_sclose (s);
            goto CLEANUP;
        default:
            fprintf (stderr, "Invalid request %c\n", request);
        }
    } while (donecount != nother);

    printf ("Completed %d cuts in %.2f seconds\n", nother, cumtime);
    printf ("%d cuts verified with HK\n", hkcount);
    printf ("%d cuts failed verification\n", xxcount);
    fflush (stdout);

CLEANUP:
    BBsafe_snet_unlisten (lport);
    BB_IFFREE (status, int);
    if (cutlist) {
        for (i = 0; i < p->cuts->cutcount; i++) {
            if (cutlist[i].proof) {
                BBio_free_cutproof (cutlist[i].proof);
                BB_FREE (cutlist[i].proof, BBcutproof);
            }
        }
        BB_FREE (cutlist, BBcuttype);
    }
    if (p) BBio_free_cutproof (p);

    return rval;
}

static int receive_cut (BB_SFILE *s, int *status, double *cumtime,
        int *result, const char *grunt)
{
    int id, check, rval = 0;
    double rtime;

    *result = -1;

    rval = BBsafe_sread_int (s, &id);
    BBcheck_rval (rval, "BBsafe_sread_int failed (id)");

    rval = BBsafe_sread_double (s, &rtime);
    BBcheck_rval (rval, "BBsafe_sread_double failed (rtime)");

    (*cumtime) += rtime;

    if (status[id] != VERIFY_STAT_WORK) {
        fprintf (stderr, "ERROR: Received cut that is not at work\n");
        rval = 1; goto CLEANUP;
    }

    rval = BBsafe_sread_int (s, &check);
    BBcheck_rval (rval, "BBsafe_sread_int failed (check)");

    if (check == 0) {
        *result = VERIFY_STAT_DONE;
        status[id] = *result; 
        printf ("Verified Cut %d (HK on %s), %.2f seconds, total %.2f\n",
                  id, grunt, rtime, *cumtime);
        fflush (stdout);
    } else {
        *result = VERIFY_STAT_FAILED;
        status[id] = *result;

        fprintf (stderr, "FAILURE ON %s\n", grunt);
        fprintf (stderr, "Verification of cut %d failed, %.2f seconds\n",
                          id, rtime);
        rval = 1;  goto CLEANUP;
    }

CLEANUP:
    return rval;
}

static int send_cut (BB_SFILE *s, int current, int *status,
       BBtsp_lpcuts *pool, int ncount, int *other)
{
    int rval = 0;
    BBtsp_lpcut_in c;

    BButil_init_lpcut_in (&c);

    rval = BButil_lpcut_to_lpcut_in (pool, &pool->cuts[other[current]], &c);
    BBcheck_rval (rval, "BButil_lpcut_to_lpcut_in failed");

    rval = BBsafe_swrite_char (s, BBtsp_VERIFY_YES);
    BBcheck_rval (rval, "BBsafe_swrite_int failed (YES)");

    rval = BBsafe_swrite_int (s, current);
    BBcheck_rval (rval, "BBsafe_swrite_int failed (cut id)");

    rval = BBsafe_swrite_int (s, ncount);
    BBcheck_rval (rval, "BBsafe_swrite_int failed (ncount)");

    rval = BBio_write_lpcut_in (s, &c, ncount);
    BBcheck_rval (rval, "BBio_write_lpcut_in failed\n");

    status[current] = VERIFY_STAT_WORK;

CLEANUP:
    BButil_free_lpcut_in (&c);
    return rval;
}
