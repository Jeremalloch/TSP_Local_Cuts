# Generated automatically from Makefile.in by configure.
#
#   This file is part of CONCORDE
#
#   (c) Copyright 1995--1999 by David Applegate, Robert Bixby,
#   Vasek Chvatal, and William Cook
#
#   Permission is granted for academic research use.  For other uses,
#   contact the authors for licensing options.
#
#   Use at your own risk.  We make no guarantees about the
#   correctness or usefulness of this code.
#


SHELL = /bin/sh
SRCROOT = ..
BLDROOT = ..
CCINCDIR=$(SRCROOT)/INCLUDE

srcdir = .

CC = gcc
CFLAGS = -O3 -g -ansi -Wall -W -Wstrict-prototypes -Wpointer-arith -Wnested-externs -I$(BLDROOT)/INCLUDE -I$(CCINCDIR)
LDFLAGS = -O3 -g -framework IOKit -framework CoreFoundation
LIBFLAGS = -lm 
RANLIB = ranlib

OBJ_SUFFIX = o
o = $(OBJ_SUFFIX)

THISLIB=bmatch.a
LIBSRCS=bmatch.c christo.c bmatchjump.c bmatchpricing.c bmatchchk.c bmatchpseudo.c bmatchutil.c bmatchexpand.c bmatchancest.c bmatchheaps.c bmatchalloc.c
ALLSRCS=bmatchmain.c bmatchgeom.c $(LIBSRCS)

LIBS=$(BLDROOT)/EDGEGEN/edgegen.a $(BLDROOT)/KDTREE/kdtree.a \
     $(BLDROOT)/FMATCH/fmatch.a $(BLDROOT)/LINKERN/linkern.a \
     $(BLDROOT)/UTIL/util.a

all: bmatch $(THISLIB)

everything: all

bmatch: bmatchmain.$o $(THISLIB) $(LIBS)
	$(CC) $(LDFLAGS) -o $@ bmatchmain.$o $(THISLIB) $(LIBS) $(LIBFLAGS)

bmatchgeom: bmatchgeom.$o $(THISLIB) $(LIBS)
	$(CC) $(LDFLAGS) -o $@ bmatchgeom.$o $(THISLIB) $(LIBS) $(LIBFLAGS)

christo: christomain.$o $(THISLIB) $(LIBS)
	$(CC) $(LDFLAGS) -o $@ christomain.$o $(THISLIB) $(LIBS) $(LIBFLAGS)

clean:
	-rm -f *.$o $(THISLIB) bmatch  bmatchgeom christo

OBJS=$(LIBSRCS:.c=.o)

$(THISLIB): $(OBJS)
	$(AR) $(ARFLAGS) $(THISLIB) $(OBJS)
	$(RANLIB) $(THISLIB)

.PHONY: $(BLDROOT)/concorde.a
$(BLDROOT)/concorde.a: $(OBJS)
	$(AR) $(ARFLAGS) $(BLDROOT)/concorde.a $(OBJS)
	$(RANLIB) $(BLDROOT)/concorde.a

include ../INCLUDE/Makefile.common

# DO NOT DELETE THIS LINE -- make depend depends on it.

I=$(CCINCDIR)
I2=$(BLDROOT)/INCLUDE

bmatchmain.$o: bmatchmain.c      $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/edgegen.h  
bmatchgeom.$o: bmatchgeom.c      $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/edgegen.h  
christomain.$o: christomain.c    $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/edgegen.h  
bmatch.$o:     bmatch.c          $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
christo.$o:    christo.c         $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchjump.$o: bmatchjump.c      $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchpricing.$o: bmatchpricing.c       $(I)/machdefs.h $(I2)/config.h         \
        $(I)/bmatch.h   $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchchk.$o: bmatchchk.c        $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchpseudo.$o: bmatchpseudo.c         $(I)/machdefs.h $(I2)/config.h         \
        $(I)/bmatch.h   $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchutil.$o: bmatchutil.c  $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h     \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchexpand.$o: bmatchexpand.c  $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchancest.$o: bmatchancest.c  $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchheaps.$o: bmatchheaps.c    $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
bmatchalloc.$o: bmatchalloc.c    $(I)/machdefs.h $(I2)/config.h  $(I)/bmatch.h \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
