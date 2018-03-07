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
SRCROOT = .
BLDROOT = .
CCINCDIR=$(SRCROOT)/INCLUDE
BLDINCDIR=$(BLDROOT)/INCLUDE

srcdir = .

CC = gcc
CFLAGS = -O3 -g -ansi -Wall -W -Wstrict-prototypes -Wpointer-arith -Wnested-externs -I$(BLDROOT)/INCLUDE -I$(CCINCDIR)
LDFLAGS = -O3 -g -framework IOKit -framework CoreFoundation
LIBFLAGS = -lm 
RANLIB = ranlib

OBJ_SUFFIX = o
o = $(OBJ_SUFFIX)

DIRS=BIGGUY   BMATCH   LOCALCUT BMATCH   COMBS    CUT      EDGEGEN  \
     FMATCH   HELDKARP KDTREE   LINKERN  LP       PQ       TINY     \
     TOOLS    TSP      UTIL     VERIFY

all: concorde.h concorde.a
	for i in $(DIRS); do ( cd $$i; $(MAKE) $@ ); done

everything: all
	for i in $(DIRS); do ( cd $$i; $(MAKE) $@ ); done

.PHONY: concorde.a
concorde.a:
	for i in $(DIRS); do ( cd $$i; $(MAKE) concorde.a ); done

INC_LIST=$(BLDINCDIR)/config.h  $(CCINCDIR)/machdefs.h $(CCINCDIR)/util.h     \
         $(CCINCDIR)/bigguy.h   $(CCINCDIR)/combs.h    $(CCINCDIR)/cut.h      \
         $(CCINCDIR)/delaunay.h $(CCINCDIR)/edgegen.h  $(CCINCDIR)/fmatch.h   \
         $(CCINCDIR)/heldkarp.h $(CCINCDIR)/kdtree.h   $(CCINCDIR)/linkern.h  \
         $(CCINCDIR)/lp.h       $(CCINCDIR)/tsp.h      $(CCINCDIR)/consec1.h  \
         $(CCINCDIR)/localcut.h $(CCINCDIR)/macrorus.h $(CCINCDIR)/mlinkern.h \
         $(CCINCDIR)/necklace.h $(CCINCDIR)/pqsets.h   $(CCINCDIR)/pq.h       \
         $(CCINCDIR)/cuttree.h  $(CCINCDIR)/bmatch.h                          \
         $(CCINCDIR)/tinytsp.h  $(CCINCDIR)/tsp.h      $(CCINCDIR)/verify.h

concorde.h: $(INC_LIST) Makefile
	cat $(INC_LIST) | grep -v '#include "' > concorde.h

clean:
	-rm -f *.$o concorde.h concorde.a
	for i in $(DIRS); do ( cd $$i; $(MAKE) clean ); done

distclean: clean
	-rm -f Makefile */Makefile INCLUDE/config.h INCLUDE/Makefile.common \
               config.log config.cache config.status

depend.in:
	for i in $(DIRS); do (cd $$i; $(MAKE) depend.in); done
