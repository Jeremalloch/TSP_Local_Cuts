#include <stdio.h>
#include "verify.h"
#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static int receive_cut (CC_SFILE *s, int *ncount, int *id, CCtsp_lpcut_in *c);
static int send_cut (CC_SFILE *s, int id, double rtime, int ccheck);

int main (int ac, char **av)
{
    char *bosshost = (char *) NULL;
    double rtime = 0.0;
    int id = -1;
    int rval = 0;
    CC_SFILE *s = (CC_SFILE *) NULL;
    double szeit;
    char task;
    int k, ccheck, ctype, ncount = 0;
    CCtsp_lpcut_in c;
    char buf[1024];
    char myname[128];

    CCtsp_init_lpcut_in (&c);

    if (ac != 2) {
        fprintf (stderr, "Usage: %s boss\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    CCutil_printlabel ();
    bosshost = av[1];

    rval = gethostname (myname, 127);
    CCcheck_rval (rval, "gethostname failed");
    printf ("Machine Name: %s\n", myname); fflush (stdout);

    while (1) {
        k = 0;
        do {
            s = CCutil_snet_open (bosshost, CCtsp_VERIFY_PORT);
            if (!s) {
                fprintf (stderr, "CCutil_snet_open failed\n");
                sleep (100);
            }
            k++;
        } while (!s && k < 5);

        if (!s) {
            fprintf (stderr, "Could not connect in %d trys.\n", k);
            goto CLEANUP;
        }

        rval = CCutil_swrite_string (s, myname);
        CCcheck_rval (rval, "CCutil_swrite_string failed (NAME)");

        rval = CCutil_swrite_char (s, CCtsp_VERIFY_SEND);
        CCcheck_rval (rval, "CCutil_swrite_char failed (SEND)");

        rval = CCutil_sread_char (s, &task);
        CCcheck_rval (rval, "CCutil_sread_char failed (task)");

        if (task == CCtsp_VERIFY_YES) {
            rval = receive_cut (s, &ncount, &id, &c);
            CCcheck_rval (rval, "receive_cut failed");
        } else {
            printf ("No more work - shut down\n");
            goto CLEANUP;
        }

        CCutil_sclose (s);

        szeit = CCutil_zeit ();

        /* Handle the Verification */

        CCtsp_print_lpcut_in (&c);

        sprintf (buf, "v%d", id);

        rval = CCverify_cut (&c, ncount, CC_TYPE_ALL , &ctype, 0, &ccheck, buf);
        if (rval) {
            fprintf (stderr, "CCverify_cut failed\n");
            ccheck = -1;
        } else {
            printf ("Cut %d is verified (%d)\n", id, ccheck);  fflush (stdout);
            if (ctype != CC_TYPE_OTHER) {
                printf ("Interesting -- simple form of cut in list\n");
                fflush (stdout);
            }
        }

        rtime = CCutil_zeit () - szeit;

        k = 0;
        do {
            s = CCutil_snet_open (bosshost, CCtsp_VERIFY_PORT);
            if (!s) {
                fprintf (stderr, "CCutil_snet_open failed\n");
                sleep (100);
            }
            k++;
        } while (!s && k < 5);

        if (!s) {
            fprintf (stderr, "Could not connect in %d trys.\n", k);
            goto CLEANUP;
        }

        rval = CCutil_swrite_string (s, myname);
        CCcheck_rval (rval, "CCutil_swrite_string failed (NAME)");

        rval = CCutil_swrite_char (s, CCtsp_VERIFY_RECEIVE);
        CCcheck_rval (rval, "CCutil_swrite_char failed (RECEIVE)");

        rval = send_cut (s, id, rtime, ccheck);
        CCcheck_rval (rval, "send_cut failed");

        CCtsp_free_lpcut_in (&c);
        CCtsp_init_lpcut_in (&c);
        CCutil_sclose (s);
    }

CLEANUP:

    if (s != (CC_SFILE *) NULL) {
        CCutil_sclose (s);
    }
    return rval;
}

static int receive_cut (CC_SFILE *s, int *ncount, int *id, CCtsp_lpcut_in *c)
{
    int rval = 0;
    int n = 0;

    CCtsp_init_lpcut_in (c);

    rval = CCutil_sread_int (s, id);
    CCcheck_rval (rval, "CCutil_sread_int failed (id)");

    printf ("Cut id: %d\n", *id); fflush (stdout);

    rval = CCutil_sread_int (s, &n);
    CCcheck_rval (rval, "CCutil_sread_int failed (id)");

    rval = CCtsp_read_lpcut_in (s, c, n);
    CCcheck_rval (rval, "CCutil_read_lpcut_in failed");

    if (ncount) *ncount = n;

CLEANUP:
    return rval;
}

static int send_cut (CC_SFILE *s, int id, double rtime, int ccheck)
{
    int rval = 0;

    rval = CCutil_swrite_int (s, id);
    CCcheck_rval (rval, "CCutil_swrite_int failed (id)");

    rval = CCutil_swrite_double (s, rtime);
    CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");

    rval = CCutil_swrite_int (s, ccheck);
    CCcheck_rval (rval, "CCutil_swrite_int failed (ccheck)");
    
CLEANUP:
    return rval;
}
