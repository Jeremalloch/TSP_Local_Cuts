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
/*  Grunt (worker) for bbboss.  Verifies cuts with Held-Karp solver.        */
/****************************************************************************/

#include "bbproof.h"

static int
    receive_cut (BB_SFILE *s, int *ncount, int *id, BBtsp_lpcut_in *c),
    send_cut (BB_SFILE *s, int id, double rtime, int ccheck);

int main (int ac, char **av)
{
    char task, myname[128], *bosshost = (char *) NULL;
    double szeit, rtime = 0.0;
    int id = -1, k, ccheck, ctype, ncount = 0, rval = 0;
    BB_SFILE *s = (BB_SFILE *) NULL;
    BBtsp_lpcut_in c;

    if (ac != 2) {
        fprintf (stderr, "Usage: %s boss\n", av[0]);
        rval = 1; goto CLEANUP;
    }
    BButil_printlabel ();
    bosshost = av[1];

    rval = gethostname (myname, 127);
    BBcheck_rval (rval, "gethostname failed");
    printf ("Machine Name: %s\n", myname); fflush (stdout);

    while (1) {
        k = 0;
        do {
            s = BBsafe_snet_open (bosshost, BBtsp_VERIFY_PORT);
            if (!s) {
                fprintf (stderr, "BBsafe_snet_open failed\n");
                sleep (100);
            }
            k++;
        } while (!s && k < 5);

        if (!s) {
            fprintf (stderr, "Could not connect in %d trys.\n", k);
            goto CLEANUP;
        }

        rval = BBsafe_swrite_string (s, myname);
        BBcheck_rval (rval, "BBsafe_swrite_string failed (NAME)");

        rval = BBsafe_swrite_char (s, BBtsp_VERIFY_SEND);
        BBcheck_rval (rval, "BBsafe_swrite_char failed (SEND)");

        rval = BBsafe_sread_char (s, &task);
        BBcheck_rval (rval, "BBsafe_sread_char failed (task)");

        if (task == BBtsp_VERIFY_YES) {
            rval = receive_cut (s, &ncount, &id, &c);
            BBcheck_rval (rval, "receive_cut failed");
        } else {
            printf ("No more work - shut down\n");
            goto CLEANUP;
        }

        BBsafe_sclose (s);

        szeit = BButil_zeit ();
        rval = BBcuts_verify (&c, ncount, BB_TYPE_OTHER, &ctype,
                                         (BBdatagroup *) NULL);
        if (rval) {
            fprintf (stderr, "BBcuts_verify failed\n");
            ccheck = -1;
        } else {
            printf ("Cut %d is verified\n", id);  fflush (stdout);
            if (ctype != BB_TYPE_OTHER) {
                printf ("Interesting -- simple form of cut in list\n");
                fflush (stdout);
            }
            ccheck = 0;
        }

        rtime = BButil_zeit () - szeit;

        k = 0;
        do {
            s = BBsafe_snet_open (bosshost, BBtsp_VERIFY_PORT);
            if (!s) {
                fprintf (stderr, "BBsafe_snet_open failed\n");
                sleep (100);
            }
            k++;
        } while (!s && k < 5);

        if (!s) {
            fprintf (stderr, "Could not connect in %d trys.\n", k);
            goto CLEANUP;
        }

        rval = BBsafe_swrite_string (s, myname);
        BBcheck_rval (rval, "BBsafe_swrite_string failed (NAME)");

        rval = BBsafe_swrite_char (s, BBtsp_VERIFY_RECEIVE);
        BBcheck_rval (rval, "BBsafe_swrite_char failed (RECEIVE)");

        rval = send_cut (s, id, rtime, ccheck);
        BBcheck_rval (rval, "send_cut failed");

        BButil_free_lpcut_in (&c);
        BBsafe_sclose (s);
    }

CLEANUP:
    if (s != (BB_SFILE *) NULL) {
        BBsafe_sclose (s);
    }
    return rval;
}

static int receive_cut (BB_SFILE *s, int *ncount, int *id, BBtsp_lpcut_in *c)
{
    int n = 0, rval = 0;

    BButil_init_lpcut_in (c);

    rval = BBsafe_sread_int (s, id);
    BBcheck_rval (rval, "BBsafe_sread_int failed (id)");

    printf ("Cut id: %d\n", *id); fflush (stdout);

    rval = BBsafe_sread_int (s, &n);
    BBcheck_rval (rval, "BBsafe_sread_int failed (id)");

    rval = BBio_read_lpcut_in (s, c, n);
    BBcheck_rval (rval, "BBio_read_lpcut_in failed");

    if (ncount) *ncount = n;

CLEANUP:
    return rval;
}

static int send_cut (BB_SFILE *s, int id, double rtime, int ccheck)
{
    int rval = 0;

    rval = BBsafe_swrite_int (s, id);
    BBcheck_rval (rval, "BBsafe_swrite_int failed (id)");

    rval = BBsafe_swrite_double (s, rtime);
    BBcheck_rval (rval, "BBsafe_swrite_double failed (rtime)");

    rval = BBsafe_swrite_int (s, ccheck);
    BBcheck_rval (rval, "BBsafe_swrite_int failed (ccheck)");
    
CLEANUP:
    return rval;
}
