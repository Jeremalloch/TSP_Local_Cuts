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
/*  Binary Input/Output routines from UTIL/safe_io.h                        */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 13, 1995                                                 */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  BB_SFILE *BBsafe_sopen (char *f, char *s)                               */
/*      Opens a file for buffered binary I/O.  The buffered binary I/O      */
/*      routines using BB_SFILE's attempt to be machine independent,        */
/*      and only compatible with themselves.  Comparable to the stdio       */
/*      routine fopen().  If the file already exists and is being           */
/*      opened for output, the old file is renamed by prepending an O       */
/*      to is name.                                                         */
/*    f - the filename to open.  "stdin" means descriptor 0, "stdout"       */
/*        descriptor 1, and "stderr" descriptor 2.  "-" means               */
/*        descriptor 0 or 1, depending on wither the file is opened         */
/*        for reading or writing.                                           */
/*    s - the mode to open, either "r" for input, or "w" for output.        */
/*    returns a pointer to the opened file, or NULL if there is an          */
/*        error.                                                            */
/*                                                                          */
/*  int BBsafe_swrite (BB_SFILE *f, char *buf, int size)                    */
/*      writes to a buffered binary I/O file.                               */
/*    f - the BB_SFILE to write to                                          */
/*    buf - the data to write                                               */
/*    size - the number of bytes to write.                                  */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_bits (BB_SFILE *f, int x, int xbits)                  */
/*      writes bits to a buffered binary I/O file.                          */
/*    f - the BB_SFILE to write to                                          */
/*    x - an int containing the data to write                               */
/*    xbits - the number of bits to write.  The lowest order xbits          */
/*            bits of x will be written.                                    */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_ubits (BB_SFILE *f, unsigned int x, int xbits)        */
/*      writes bits to a buffered binary I/O file.                          */
/*    f - the BB_SFILE to write to                                          */
/*    x - an unsigned int int containing the data to write                  */
/*    xbits - the number of bits to write.  The lowest order xbits          */
/*            bits of x will be written.                                    */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_char (BB_SFILE *f, char x)                            */
/*      writes a char to a buffered binary I/O file.                        */
/*    f - the BB_SFILE to write to                                          */
/*    x - the char to write                                                 */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_string (BB_SFILE *f, const char *s)                   */
/*      writes a string to a buffered binary I/O file.                      */
/*    f - the BB_SFILE to write to                                          */
/*    s - the string to write.  The array of characters in s up to and      */
/*        including the first NULL are written.                             */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_short (BB_SFILE *f, short x)                          */
/*      writes a short to a buffered binary I/O file.                       */
/*    f - the BB_SFILE to write to                                          */
/*    x - the short to write                                                */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_ushort (BB_SFILE *f, unsigned short x)                */
/*      writes an unsigned short to a buffered binary I/O file.             */
/*    f - the BB_SFILE to write to                                          */
/*    x - the unsigned short to write                                       */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_int (BB_SFILE *f, int x)                              */
/*      writes an int to a buffered binary I/O file.                        */
/*    f - the BB_SFILE to write to                                          */
/*    x - the int to write                                                  */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_uint (BB_SFILE *f, unsigned int x)                    */
/*      writes an unsigned int to a buffered binary I/O file.               */
/*    f - the BB_SFILE to write to                                          */
/*    x - the unsigned int to write                                         */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_swrite_double (BB_SFILE *f, double x)                        */
/*      writes a double to a buffered binary I/O file.                      */
/*    f - the BB_SFILE to write to                                          */
/*    x - the double to write                                               */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread (BB_SFILE *f, char *buf, int size)                     */
/*      reads from a buffered binary I/O file.                              */
/*    f - the BB_SFILE to read from.                                        */
/*    buf - a buffer in which to store the data read.  buf should have      */
/*          space for size characters.                                      */
/*    size - the number of bytes to read.                                   */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_bits (BB_SFILE *f, int *x, int xbits)                  */
/*      reads bits from a buffered binary I/O file.                         */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the bits read (in the low-order           */
/*        xbits bits).                                                      */
/*    xbits - the number of bits read.                                      */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_ubits (BB_SFILE *f, unsigned int *x, int xbits)        */
/*      reads bits from a buffered binary I/O file.                         */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the bits read (in the low-order           */
/*        xbits bits).                                                      */
/*    xbits - the number of bits read.                                      */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_char (BB_SFILE *f, char *x)                            */
/*      reads a char from a buffered binary I/O file.                       */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the char read                             */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_string (BB_SFILE *f, char *x, int maxlen)              */
/*      reads a string from a buffered binary I/O file.                     */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the string read.                          */
/*    maxlen - the maximum number of characters to read.                    */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_short (BB_SFILE *f, short *x)                          */
/*      reads a short from a buffered binary I/O file.                      */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the short read                            */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_ushort (BB_SFILE *f, unsigned short *x)                */
/*      reads an unsigned short from a buffered binary I/O file.            */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the unsigned short read                   */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_int (BB_SFILE *f, int *x)                              */
/*      reads an int from a buffered binary I/O file.                       */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the int read                              */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_uint (BB_SFILE *f, unsigned int *x)                    */
/*      reads an unsigned int from a buffered binary I/O file.              */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the unsigned int read                     */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sread_double (BB_SFILE *f, double *x)                        */
/*      reads a double from a buffered binary I/O file.                     */
/*    f - the BB_SFILE to read from.                                        */
/*    x - on return, will contain the double read                           */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sflush (BB_SFILE *f)                                         */
/*      flushes the buffer of a buffered binary I/O file.                   */
/*    f - the BB_SFILE to flush                                             */
/*    returns 0 if succesful, nonzero if failure.                           */
/*                                                                          */
/*  int BBsafe_sclose (BB_SFILE *f)                                         */
/*      closes a BB_SFILE.                                                  */
/*    f - the BB_SFILE to close                                             */
/*    returns 0 for success, nonzero for failure.                           */
/*                                                                          */
/*  int BBsafe_sbits (unsigned int x)                                       */
/*      computes the number of bits necessary to represent all              */
/*      nonnegative integers <= x                                           */
/*    x - a number                                                          */
/*    returns the number of bits necessary to represent x.                  */
/*                                                                          */
/*  BB_SFILE *BBsafe_snet_open (char *h, unsigned short p)                  */
/*      Opens a network connection to a port on a remote host               */
/*    h - the name of the host to connect to                                */
/*    p - the port on the host to connect to                                */
/*    returns a BB_SFILE (opened for input and output) for buffered         */
/*            binary I/O to the specified port on the remote host,          */
/*            or NULL if there is a failure.                                */
/*    Only exists if BB_NETREADY is defined                                 */
/*                                                                          */
/*  BB_SFILE *BBsafe_snet_receive (BB_SPORT *s)                             */
/*      Accepts a network connection on a port.                             */
/*    s - the BB_SPORT to accept a connection from.  Must be the            */
/*        returned result of a successfull BBsafe_snet_listen call.         */
/*    returns a BB_SFILE (opened for input and output) for buffered         */
/*        binary I/O on the specified port, or NULL if there is a           */
/*        failure.                                                          */
/*    Only exists if BB_NETREADY is defined                                 */
/*                                                                          */
/*  BB_SPORT *BBsafe_snet_listen (unsigned short p)                         */
/*      Prepares to accept network connections on a port.                   */
/*    p - the port on which to accept connections.                          */
/*    returns a BB_SPORT for accepting connections on the specified         */
/*        port.  This return value is passed to BBsafe_snet_receive to      */
/*        accept a connection.  Returns NULL if there is a failure.         */
/*    Only exists if BB_NETREADY is defined                                 */
/*                                                                          */
/*  void BBsafe_snet_unlisten (BB_SPORT *s)                                 */
/*      Ceases accepting network connections from an BB_SPORT.              */
/*    s - the BB_SPORT to close.                                            */
/*    Only exists if BB_NETREADY is defined                                 */
/*                                                                          */
/****************************************************************************/

#include "bbproof.h"

static BB_SFILE
    *sopen_write (const char *f),
    *sopen_read (const char *f),
    *sdopen (int t),
    *sdopen_write (int t),
    *sdopen_read (int t),
    *sdopen_readwrite (int t);

static int
    swrite_buffer (BB_SFILE *f),
    sread_buffer (BB_SFILE *f),
    prepare_write (BB_SFILE *f),
    prepare_read (BB_SFILE *f);

static void
    sinit (BB_SFILE *s);


/* BBsafe_sopen interprets filenames "stdin" as descriptor 0, "stdout" as
 * descriptor 1, and "stderr" as descriptor 2.  "-" is interpreted as
 * 0 or 1, depending on whether the file is opened for reading or writing.
 *
 * BBsafe_sclose doesn't close descriptors 0, 1, and 2.
 */

/* When writing, written data extends from buffer[0] bit 7 through
 * buffer[chars_in_buffer-1] bit bits_in_last_char.  Empty space extends
 * from buffer[chars_in_buffer-1] bit bits_in_last_char-1 through
 * buffer[BB_SBUFFER_SIZE-1] bit 0.
 *
 * When reading, read data extends from buffer[0] bit 7 through
 * buffer[current_buffer_char] bit bits_in_last_char.  unread data
 * extends from buffer[current_buffer_char] bit bits_in_last_char-1
 * through buffer[chars_in_buffer-1] bit 0.  Empty space extends from
 * buffer[chars_in_buffer] bit 7 through buffer[BB_SBUFFER_SIZE-1] bit 0.
 */

/* If the routines detect an error, they return -1.
 */

#define SREAD 1
#define SWRITE 2
#define SRW_EMPTY 3
#define SRW_READ 4
#define SRW_WRITE 5

#define TFILE 1
#define TDESC 2
#define TNET 3

#define NBITMASK(n) ((1<<(n))-1)
#define BITRANGE(x,start,length) (((x) >> (start)) & NBITMASK(length))
#define BITS_PER_CHAR (8)

#ifndef O_BINARY
#define O_BINARY 0
#endif
#ifndef O_EXCL
#define O_EXCL 0
#endif

BB_SFILE *BBsafe_sopen (const char *f, const char *s)
{
    if (strcmp (s, "r") == 0) {
        return sopen_read (f);
    } else if (strcmp (s, "w") == 0) {
        return sopen_write (f);
    } else {
        fprintf (stderr, "Need to specify read/write in BBsafe_sopen\n");
        return (BB_SFILE *) NULL;
    }
}

static BB_SFILE *sopen_write (const char *f)
{
    BB_SFILE *s = (BB_SFILE *) NULL;
    int t;
    char fbuf[BB_SFNAME_SIZE];
    char fbuf_N[BB_SFNAME_SIZE + 32];
    char fbuf_Nx[BB_SFNAME_SIZE + 64];

    strncpy (fbuf, f, sizeof (fbuf) - 12);
    fbuf[sizeof (fbuf) - 12] = '\0';
    sprintf (fbuf_N, "N%s", fbuf);
    sprintf (fbuf_Nx, "N%s~", fbuf);

    if (strcmp (f, "stdout") == 0 || strcmp (f, "-") == 0) {
        s = sdopen_write (1);
    } else if (strcmp (f, "stderr") == 0) {
        s = sdopen_write (2);
    } else {
        t = open (fbuf_N, O_WRONLY | O_CREAT | O_BINARY | O_EXCL, 0644);
        if (t == -1 && errno == EEXIST) {
            fprintf (stderr, "%s already exists, renaming to %s\n",
                          fbuf_N, fbuf_Nx);
            if (rename (fbuf_N, fbuf_Nx)) {
                perror (fbuf_Nx);
                fprintf (stderr, "Couldn't rename %s to %s\n", fbuf_N,
                         fbuf_Nx);
                return (BB_SFILE *) NULL;
            }
            t = open (fbuf_N, O_WRONLY | O_CREAT | O_BINARY | O_EXCL, 0644);
        }
        if (t == -1) {
            perror (fbuf_N);
            fprintf (stderr, "Couldn't open %s for output\n", fbuf_N);
            return (BB_SFILE *) NULL;
        }
        s = sdopen_write (t);
        if (!s) {
            close (t);
        } else {
            s->type = TFILE;
        }
    }
    if (s) {
        strncpy (s->fname, fbuf, sizeof (s->fname));
        s->fname[sizeof (s->fname)-1] = '\0';
    }
    return s;
}

static BB_SFILE *sopen_read (const char *f)
{
    BB_SFILE *s = (BB_SFILE *) NULL;
    int t;

    if (strcmp (f, "stdin") == 0 || strcmp (f, "-") == 0) {
        s = sdopen_read (0);
    } else {
        t = open (f, O_RDONLY | O_BINARY, 0644);
        if (t == -1) {
            perror (f);
            fprintf (stderr, "Couldn't open for input\n");
            s = (BB_SFILE *) NULL;
        }
        s = sdopen_read (t);
        if (!s) {
            close (t);
        } else {
            s->type = TFILE;
        }
    }
    if (s) {
        strncpy (s->fname, f, sizeof (s->fname));
        s->fname[sizeof (s->fname)-1] = '\0';
    }
    return s;
}

static BB_SFILE *sdopen (int t)
{
    BB_SFILE *s = (BB_SFILE *) NULL;

    if (t < 0) {
        fprintf (stderr, "Invalid descriptor %d\n", t);
        return (BB_SFILE *) NULL;
    }

    s = BB_SAFE_MALLOC (1, BB_SFILE);
    if (s == (BB_SFILE *) NULL) {
        return (BB_SFILE *) NULL;
    }
    sinit (s);

    s->desc = t;
    s->type = TDESC;
    sprintf (s->fname, "descriptor %d", t);
    return s;
}

static BB_SFILE *sdopen_write (int t)
{
    BB_SFILE *s = (BB_SFILE *) NULL;

    s = sdopen (t);
    if (s) {
        s->status = SWRITE;
    }

    return s;
}

static BB_SFILE *sdopen_read (int t)
{
    BB_SFILE *s = (BB_SFILE *) NULL;

    s = sdopen (t);
    if (s) {
        s->status = SREAD;
    }

    return s;
}

static BB_SFILE *sdopen_readwrite (int t)
{
    BB_SFILE *s = (BB_SFILE *) NULL;

    s = sdopen (t);
    if (s) {
        s->status = SRW_EMPTY;
    }

    return s;
}

int BBsafe_swrite (BB_SFILE *f, char *buf, int size)
{
    int i;

    for (i=0; i<size; i++) {
        if (BBsafe_swrite_char (f, buf[i])) return -1;
    }
    return 0;
}

int BBsafe_swrite_bits (BB_SFILE *f, int x, int xbits)
{
    if (x < 0) {
        fprintf (stderr, "BBsafe_swrite_bits cannot write negative numbers\n");
        return -1;
    }
    return BBsafe_swrite_ubits (f, (unsigned int) x, xbits);
}

int BBsafe_swrite_ubits (BB_SFILE *f, unsigned int x, int xbits)
{
    int getbits;
    unsigned int v;

    if (prepare_write (f)) return -1;

    while (xbits) {
        if (f->bits_in_last_char == 0) {
            if (f->chars_in_buffer == BB_SBUFFER_SIZE) {
                if (swrite_buffer (f)) return -1;
            }
            f->buffer[f->chars_in_buffer++] = 0;
            f->bits_in_last_char = BITS_PER_CHAR;
        }
        getbits = f->bits_in_last_char;
        if (getbits > xbits)
            getbits = xbits;
        xbits -= getbits;
        f->bits_in_last_char -= getbits;
        v = BITRANGE (x, xbits, getbits);
        f->buffer[f->chars_in_buffer - 1] =
            (unsigned int) f->buffer[f->chars_in_buffer - 1] |
            (unsigned int) (v << f->bits_in_last_char);
    }
    return 0;
}

int BBsafe_swrite_char (BB_SFILE *f, char x)
{
    unsigned char ux = (unsigned char) x;
    
    if (prepare_write (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->chars_in_buffer + 1 > BB_SBUFFER_SIZE) {
        if (swrite_buffer (f)) return -1;
    }
    f->buffer[f->chars_in_buffer++] = ((unsigned int) ux) & 0xff;
    return 0;
}

int BBsafe_swrite_string (BB_SFILE *f, const char *s)
{
    int rval;

    while (*s) {
        rval = BBsafe_swrite_char (f, *s);
        if (rval)
            return rval;
        s++;
    }
    BBsafe_swrite_char (f, (unsigned char) 0);
    return 0;
}

int BBsafe_swrite_short (BB_SFILE *f, short x)
{
    return BBsafe_swrite_ushort (f, (unsigned short) x);
}

int BBsafe_swrite_ushort (BB_SFILE *f, unsigned short x)
{
    if (prepare_write (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->chars_in_buffer + 2 > BB_SBUFFER_SIZE) {
        if (swrite_buffer (f)) return -1;
    }

    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 8) & 0xff;
    f->buffer[f->chars_in_buffer++] = ((unsigned int) x) & 0xff;
    return 0;
}

int BBsafe_swrite_int (BB_SFILE *f, int x)
{
    return BBsafe_swrite_uint (f, (unsigned int) x);
}

int BBsafe_swrite_uint (BB_SFILE *f, unsigned int x)
{
    if (prepare_write (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->chars_in_buffer + 4 > BB_SBUFFER_SIZE) {
        if (swrite_buffer (f)) return -1;
    }

    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 24) & 0xff;
    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 16) & 0xff;
    f->buffer[f->chars_in_buffer++] = (((unsigned int) x) >> 8) & 0xff;
    f->buffer[f->chars_in_buffer++] = ((unsigned int) x) & 0xff;
    return 0;
}

int BBsafe_swrite_double (BB_SFILE *f, double x)
{
    unsigned short e;
    unsigned int m1;
    unsigned int m2;

    e = 128;

    if (x < 0) {
        e = (unsigned int) e + 256;
        x = -x;
    }

    if (x >= 1.0) {
#define MUNCH_HI_EXP(x,e,v,lv) if (x >= v) {e = (unsigned int) e + lv; x *= 1/v;}
        MUNCH_HI_EXP(x,e,18446744073709551616.0,64);
        MUNCH_HI_EXP(x,e,4294967296.0,32);
        MUNCH_HI_EXP(x,e,65536.0, 16);
        MUNCH_HI_EXP(x,e,256.0, 8);
        MUNCH_HI_EXP(x,e,16.0, 4);
        MUNCH_HI_EXP(x,e,4.0, 2);
        MUNCH_HI_EXP(x,e,2.0, 1);
#undef MUNCH_HI_EXP
        x /= 2;
        e = (unsigned int) e + 1;
    } else if (x < 0.5) {
#define MUNCH_LO_EXP(x,e,v,lv) if (x < 1/v) {e = (unsigned int) e - lv; x *= v;}
        MUNCH_LO_EXP(x,e,18446744073709551616.0,64);
        MUNCH_LO_EXP(x,e,4294967296.0,32);
        MUNCH_LO_EXP(x,e,65536.0, 16);
        MUNCH_LO_EXP(x,e,256.0, 8);
        MUNCH_LO_EXP(x,e,16.0, 4);
        MUNCH_LO_EXP(x,e,4.0, 2);
        MUNCH_LO_EXP(x,e,2.0, 1);
#undef MUNCH_LP_EXP
    }
    x *= 4294967296.0;
    m1 = (unsigned int) x;
    m2 = (unsigned int) ((x - m1) * 4294967296.0);
    if (BBsafe_swrite_ushort (f, e)) return -1;
    if (BBsafe_swrite_uint (f, m1)) return -1;
    if (BBsafe_swrite_uint (f, m2)) return -1;
    return 0;
}

int BBsafe_sread (BB_SFILE *f, char *buf, int size)
{
    int i;

    for (i=0; i<size; i++) {
        if (BBsafe_sread_char (f, &buf[i])) return -1;
    }
    return 0;
}

int BBsafe_sread_bits (BB_SFILE *f, int *x, int xbits)
{
    unsigned int ux = 0;
    int rval;
    
    rval = BBsafe_sread_ubits (f, &ux, xbits);
    *x = (int) ux;
    return rval;
}

int BBsafe_sread_ubits (BB_SFILE *f, unsigned int *x, int xbits)
{
    int getbits;
    unsigned int v;

    if (prepare_read (f)) return -1;

    *x = 0;
    while (xbits) {
        if (f->bits_in_last_char == 0) {
            if (f->current_buffer_char + 1 == f->chars_in_buffer) {
                if (sread_buffer (f)) return -1;
            }
            f->current_buffer_char++;
            f->bits_in_last_char = BITS_PER_CHAR;
        }
        getbits = f->bits_in_last_char;
        if (getbits > xbits)
            getbits = xbits;
        f->bits_in_last_char -= getbits;
        xbits -= getbits;
        v = BITRANGE ((unsigned int) f->buffer[f->current_buffer_char],
                      f->bits_in_last_char, getbits);
        *x |= v << xbits;
    }
    return 0;
}

int BBsafe_sread_char (BB_SFILE *f, char *x)
{
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = (char) (f->buffer[++f->current_buffer_char]);
    return 0;
}

int BBsafe_sread_string (BB_SFILE *f, char *x, int maxlen)
{
    int i, rval;

    maxlen--;
    for (i = 0; i < maxlen;  i++, x++) {
        rval = BBsafe_sread_char (f, x);
        if (rval)
            return rval;
        if (*x == 0)
            return 0;
    }
    *x = 0;
    return 0;
}

int BBsafe_sread_short (BB_SFILE *f, short *x)
{
    unsigned short ux = 0;
    int rval;

    rval = BBsafe_sread_ushort (f, &ux);
    *x = (short) ux;
    return rval;
}

int BBsafe_sread_ushort (BB_SFILE *f, unsigned short *x)
{
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = ((unsigned int) f->buffer[++f->current_buffer_char]) << 8;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = (unsigned int) *x | ((unsigned int) f->buffer[++f->current_buffer_char]);
    return 0;
}

int BBsafe_sread_int (BB_SFILE *f, int *x)
{
    unsigned int ux = 0;
    int rval;

    rval = BBsafe_sread_uint (f, &ux);
    *x = (int) ux;
    return rval;
}

int BBsafe_sread_uint (BB_SFILE *f, unsigned int *x)
{
    if (prepare_read (f)) return -1;

    f->bits_in_last_char = 0;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x = ((unsigned int) f->buffer[++f->current_buffer_char]) << 24;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 16;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x |= ((unsigned int) f->buffer[++f->current_buffer_char]) << 8;
    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        if (sread_buffer (f)) return -1;
    }
    *x |= ((unsigned int) f->buffer[++f->current_buffer_char]);
    return 0;
}

int BBsafe_sread_double (BB_SFILE *f, double *x)
{
    unsigned short e;
    unsigned int m1;
    unsigned int m2;

    if (BBsafe_sread_ushort (f, &e)) return -1;
    if (BBsafe_sread_uint (f, &m1)) return -1;
    if (BBsafe_sread_uint (f, &m2)) return -1;

    *x = ((m2 / 4294967296.0) + m1) / 4294967296.0;

    if ((unsigned int) e >= 256) {
        *x = -*x;
        e = (unsigned int) e - 256;
    }

    if ((unsigned int) e > 128) {
#define UNMUNCH_HI_EXP(x,e,v,lv) if ((unsigned int) e >= (unsigned int) (128 + lv)) \
                                     {e = (unsigned int) e - lv; x *= v;}
        UNMUNCH_HI_EXP(*x,e,18446744073709551616.0,64);
        UNMUNCH_HI_EXP(*x,e,4294967296.0,32);
        UNMUNCH_HI_EXP(*x,e,65536.0, 16);
        UNMUNCH_HI_EXP(*x,e,256.0, 8);
        UNMUNCH_HI_EXP(*x,e,16.0, 4);
        UNMUNCH_HI_EXP(*x,e,4.0, 2);
        UNMUNCH_HI_EXP(*x,e,2.0, 1);
#undef UNMUNCH_HI_EXP
    } else if ((unsigned int) e < 128) {
#define UNMUNCH_LO_EXP(x,e,v,lv) if ((unsigned int) e <= (unsigned int) (128 - lv)) \
                                     {e = (unsigned int) e + lv; x *= 1/v;}
        UNMUNCH_LO_EXP(*x,e,18446744073709551616.0,64);
        UNMUNCH_LO_EXP(*x,e,4294967296.0,32);
        UNMUNCH_LO_EXP(*x,e,65536.0, 16);
        UNMUNCH_LO_EXP(*x,e,256.0, 8);
        UNMUNCH_LO_EXP(*x,e,16.0, 4);
        UNMUNCH_LO_EXP(*x,e,4.0, 2);
        UNMUNCH_LO_EXP(*x,e,2.0, 1);
#undef UNMUNCH_LO_EXP
    }

    return 0;
}

int BBsafe_sflush (BB_SFILE *f)
{
    int rval;
    
    if (f == (BB_SFILE *) NULL) {
        rval = -1;
    } else if (f->status == SREAD || f->status == SRW_READ) {
        f->bits_in_last_char = 0;
        rval = 0;
    } else if (f->status == SWRITE || f->status == SRW_WRITE) {
        rval = swrite_buffer (f);
    } else if (f->status == SRW_EMPTY) {
        rval = 0;
    } else {
        fprintf (stderr, "Buffer %s has invalid status %d\n", f->fname,
                 f->status);
        rval = -1;
    }
        
    return rval;
}

int BBsafe_sclose (BB_SFILE *f)
{
    int retval = 0;
    char fbuf_O[BB_SFNAME_SIZE + 32];
    char fbuf_N[BB_SFNAME_SIZE + 32];

    if (!f) return -1;

    if ((f->status == SWRITE || f->status == SRW_WRITE) &&
        f->chars_in_buffer) {
        if (swrite_buffer (f)) retval = -1;
    }

    if (f->desc >= 3) {
        if (close (f->desc)) {
            perror ("close");
            fprintf (stderr, "Unable to close swrite file %s\n", f->fname);
            retval = -1;
        }
        if (f->status == SWRITE && f->type == TFILE) {
            sprintf (fbuf_N, "N%s", f->fname);
            sprintf (fbuf_O, "O%s", f->fname);
            rename (f->fname, fbuf_O);
            if (rename (fbuf_N, f->fname)) {
                perror (f->fname);
                fprintf (stderr, "Couldn't rename %s to %s\n",
                                               fbuf_N, f->fname);
                retval = -1;
            }
        }
    }

    BB_FREE (f, BB_SFILE);

    return retval;
}

static int swrite_buffer (BB_SFILE *f)
{
    char *p;
    int nleft;
    int n;

    if (!f) return -1;
    if (f->status != SWRITE && f->status != SRW_WRITE &&
        f->status != SRW_EMPTY) {
        fprintf (stderr, "%s not open for output\n", f->fname);
        return -1;
    }

    p = (char *) f->buffer;
    nleft = f->chars_in_buffer;
    while (nleft) {
        n = (int) write (f->desc, p, nleft);
        if (n == -1) {
            if (errno == EINTR) {
                fprintf (stderr, "swrite_buffer interrupted, retrying\n");
                continue;
            }
            perror ("write");
            fprintf (stderr, "swrite_buffer of %d chars to %s failed\n", nleft,
                     f->fname);
            return -1;
        }
        nleft -= n;
        p += n;
        f->pos += n;
    }
    f->bits_in_last_char = 0;
    f->chars_in_buffer = 0;
    return 0;
}

static int sread_buffer (BB_SFILE *f)
{
    int n;

    if (!f) return -1;
    if (f->status != SREAD && f->status != SRW_READ &&
        f->status != SRW_EMPTY) {
        fprintf (stderr, "%s not open for input\n", f->fname);
        return -1;
    }

    if (f->current_buffer_char + 1 == f->chars_in_buffer) {
        f->chars_in_buffer = 0;
        f->current_buffer_char = -1;
    }
    if (f->chars_in_buffer == BB_SBUFFER_SIZE) {
        fprintf (stderr, "sread_buffer for %s when buffer full\n", f->fname);
        return 0;
    }

  retry:
    n = (int) read (f->desc, (char *) f->buffer + f->chars_in_buffer,
              BB_SBUFFER_SIZE - f->chars_in_buffer);

    if (n == -1) {
        if (errno == EINTR) {
            fprintf (stderr, "sread_buffer interrupted, retrying\n");
            goto retry;
        }
        perror ("read");
        fprintf (stderr, "sread_buffer failed\n");
        return -1;
    }
    if (n == 0) {
        fprintf (stderr, "sread_buffer encountered EOF\n");
        return -1;
    }
    f->pos += n;
    f->chars_in_buffer += n;

    if (f->status == SRW_EMPTY) f->status = SRW_READ;
    
    return 0;
}

static void sinit (BB_SFILE *s)
{
    s->status = 0;
    s->desc = -1;
    s->type = 0;
    s->chars_in_buffer = 0;
    s->current_buffer_char = -1;
    s->bits_in_last_char = 0;
    s->pos = 0;
    s->fname[0] = '\0';
}

int BBsafe_sbits (unsigned int x)
{
    int i;
    int ux = x;
    unsigned int b;

    i = 32;
    b = ((unsigned int) 1) << 31;
    while ((ux & b) == 0 && i > 1) {
        b >>= 1;
        i--;
    }
    return i;
}

static int prepare_write (BB_SFILE *f)
{
    if (!f) return -1;
    if (f->status == SREAD) {
        fprintf (stderr, "%s not open for output\n", f->fname);
        return -1;
    } else if (f->status == SRW_READ) {
        f->chars_in_buffer = 0;
        f->current_buffer_char = -1;
        f->bits_in_last_char = 0;
        f->status = SRW_WRITE;
    } else if (f->status == SRW_EMPTY) {
        f->status = SRW_WRITE;
    } else if (f->status != SWRITE && f->status != SRW_WRITE) {
        fprintf (stderr, "%s has bogus status %d\n", f->fname, f->status);
        return -1;
    }
    
    return 0;
}

static int prepare_read (BB_SFILE *f)
{
    if (!f) return -1;
    if (f->status == SWRITE) {
        fprintf (stderr, "%s not open for input\n", f->fname);
        return -1;
    } else if (f->status == SRW_WRITE) {
        if (BBsafe_sflush (f)) return -1;
        f->chars_in_buffer = 0;
        f->current_buffer_char = -1;
        f->bits_in_last_char = 0;
        f->status = SRW_EMPTY;
    } else if (f->status != SREAD && f->status != SRW_READ &&
               f->status != SRW_EMPTY) {
        fprintf (stderr, "%s has bogus status %d\n", f->fname, f->status);
        return -1;
    }
    
    return 0;
}

BB_SFILE *BBsafe_snet_open (const char *hname, unsigned short p)
{
    struct hostent *h;
    struct sockaddr_in hsock;
    int s;
    BB_SFILE *f = (BB_SFILE *) NULL;

    memset ((void *) &hsock, 0, sizeof (hsock));

    h = gethostbyname (hname);
    if (h == (struct hostent *) NULL) {
        fprintf (stderr, "cannot get host info for %s\n", hname);
        return (BB_SFILE *) NULL;
    }
    memcpy ((void *) &hsock.sin_addr, (void *) h->h_addr, h->h_length);
    hsock.sin_family = AF_INET;
    hsock.sin_port = htons(p);

    s = socket (AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        perror ("socket");
        fprintf (stderr, "Unable to get socket\n");
        return (BB_SFILE *) NULL;
    }
    if (connect (s, (struct sockaddr *) &hsock, sizeof (hsock)) < 0) {
        perror ("connect");
        fprintf (stderr, "Unable to connect to %s\n", hname);
        return (BB_SFILE *) NULL;
    }

    f = sdopen_readwrite (s);
    if (f == (BB_SFILE *) NULL) {
        fprintf (stderr, "sdopen_readwrite failed\n");
        return (BB_SFILE *) NULL;
    }

    return f;
}

BB_SFILE *BBsafe_snet_receive (BB_SPORT *s)
{
    struct sockaddr_in new;
    int l;
    int t;
    BB_SFILE *f = (BB_SFILE *) NULL;

    memset ((void *) &new, 0, sizeof (new));
    new.sin_family = AF_INET;
    new.sin_addr.s_addr = INADDR_ANY;
    new.sin_port = 0;
    l = sizeof (new);

    t = accept (s->t, (struct sockaddr *) &new, (socklen_t *) &l);
    if (t < 0) {
        perror ("accept");
        fprintf (stderr, "accept failed\n");
        return (BB_SFILE *) NULL;
    }

    f = sdopen_readwrite (t);
    if (f == (BB_SFILE *) NULL) {
        fprintf (stderr, "sdopen_readwrite failed\n");
        return (BB_SFILE *) NULL;
    }

    return f;
}

BB_SPORT *BBsafe_snet_listen (unsigned short p)
{
    int s = -1;
    struct sockaddr_in me;
    BB_SPORT *sp = (BB_SPORT *) NULL;
    
    s = socket (AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        perror ("socket");
        fprintf (stderr, "Unable to get socket\n");
        goto FAILURE;
    }

    memset ((void *) &me, 0, sizeof (me));

    me.sin_addr.s_addr = INADDR_ANY;
    me.sin_family = AF_INET;
    me.sin_port = htons (p);

    if (bind (s, (struct sockaddr *) &me, sizeof (me)) < 0) {
        perror ("bind");
        fprintf (stderr, "Cannot bind socket\n");
        goto FAILURE;
    }

    if (listen (s, 100) < 0) {
        perror ("listen");
        fprintf (stderr, "Cannot listen to socket\n");
        goto FAILURE;
    }

    sp = BB_SAFE_MALLOC (1, BB_SPORT);
    if (sp == (BB_SPORT *) NULL) {
        fprintf (stderr, "Out of memory in BBsafe_snet_listen\n");
        goto FAILURE;
    }

    sp->t = s;
    sp->port = p;

    return sp;

  FAILURE:
    if (s >= 0) close (s);
    BB_IFFREE (sp, BB_SPORT);
    return (BB_SPORT *) NULL;
}

void BBsafe_snet_unlisten (BB_SPORT *s)
{
    if (s != (BB_SPORT *) NULL) {
        close (s->t);
        BB_FREE (s, BB_SPORT);
    }
}
