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
LPSOLVER_LIB = 
DOMINO_LIB = 

OBJ_SUFFIX = o
o = $(OBJ_SUFFIX)

ALLSRCS=fconvert.c tourchk.c  \
        edgunion.c edg2len.c  showres.c  killgrun.c tourlen.c  \
        prob2tsp.c perm2tsp.c

LIBS=$(BLDROOT)/KDTREE/kdtree.a $(BLDROOT)/UTIL/util.a \
     $(BLDROOT)/EDGEGEN/edgegen.a                                 \
     $(BLDROOT)/LINKERN/linkern.a                                 \
     $(BLDROOT)/FMATCH/fmatch.a                                   \
     $(BLDROOT)/KDTREE/kdtree.a                                   \
     $(BLDROOT)/UTIL/util.a

BBLIBS=$(BLDROOT)/TSP/tsp.a                                       \
     $(BLDROOT)/LOCALCUT/localcut.a                               \
     $(BLDROOT)/TINY/tiny.a                                       \
     $(BLDROOT)/EDGEGEN/edgegen.a                                 \
     $(BLDROOT)/LINKERN/linkern.a                                 \
     $(BLDROOT)/CUT/cut.a                                         \
     $(BLDROOT)/COMBS/combs.a                                     \
     $(BLDROOT)/FMATCH/fmatch.a                                   \
     $(BLDROOT)/BIGGUY/bigguy.a                                   \
     $(BLDROOT)/KDTREE/kdtree.a                                   \
     $(BLDROOT)/PQ/pq.a                                           \
     $(BLDROOT)/VERIFY/verify.a                                   \
     $(BLDROOT)/HELDKARP/heldkarp.a                               \
     $(BLDROOT)/LP/lp.a                                           \
     $(BLDROOT)/UTIL/util.a


all: showres tourchk tourlen prob2tsp edg2len edgunion fconvert perm2tsp

everything: all killgrun

fconvert: fconvert.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ fconvert.$o $(LIBS) $(LIBFLAGS)

tourchk: tourchk.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ tourchk.$o $(LIBS) $(LIBFLAGS)

showres: showres.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ showres.$o $(LIBS) $(LIBFLAGS)

killgrun: killgrun.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ killgrun.$o $(LIBS) $(LIBFLAGS)

bbtree: bbtree.$o $(BBLIBS)
	$(CC) $(LDFLAGS) -o $@ bbtree.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

bbproof: bbproof.$o bbcuts.$o bbprice.$o $(BBLIBS)
	$(CC) $(LDFLAGS) -o $@ bbproof.$o bbcuts.$o bbprice.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

bbcheck: bbcheck.$o $(BBLIBS)
	$(CC) $(LDFLAGS) -o $@ bbcheck.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

prob2tsp: prob2tsp.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ prob2tsp.$o $(LIBS) $(LIBFLAGS)

tourlen: tourlen.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ tourlen.$o $(LIBS) $(LIBFLAGS)

tspart: tspart.$o $(BBLIBS)
	$(CC) $(LDFLAGS) -o $@ tspart.$o $(LIBS) $(BLDROOT)/TSP/tsp.a $(BLDROOT)/LINKERN/linkern.a $(BLDROOT)/BMATCH/bmatch.a $(BBLIBS) $(LPSOLVER_LIB) $(LIBFLAGS)

edgunion: edgunion.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ edgunion.$o $(LIBS) $(LIBFLAGS)

edginter: edginter.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ edginter.$o $(LIBS) $(LIBFLAGS)

edg2len: edg2len.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ edg2len.$o $(LIBS) $(LIBFLAGS)

tour2edg: tour2edg.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ tour2edg.$o $(LIBS) $(LIBFLAGS)

perm2tsp: perm2tsp.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ perm2tsp.$o $(LIBS) $(LIBFLAGS)

tdivide: tdivide.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ tdivide.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

ttour: ttour.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ ttour.$o $(LIBS) $(LIBFLAGS)

tcuts: tcuts.$o $(BBLIBS)
	$(CC) $(LDFLAGS) -o $@ tcuts.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

lkhboss: lkhboss.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ lkhboss.$o $(LIBS) $(LIBFLAGS)

lkhgrunt: lkhgrunt.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ lkhgrunt.$o $(LIBS) $(LIBFLAGS)

lpelim: lpelim.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ lpelim.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

cutview: cutview.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ cutview.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)

xcut: xcut.$o $(LIBS)
	$(CC) $(LDFLAGS) -o $@ xcut.$o $(BBLIBS) $(BLDROOT)/TSP/tsp.a $(LPSOLVER_LIB) $(DOMINO_LIB) $(LIBFLAGS)


clean:
	-rm -f *.$o fconvert tourchk \
                    showres \
                    prob2tsp tourlen \
                    edgunion killgrun \
                    edg2len lpelim \
                    cutview

$(BLDROOT)/concorde.a:

include ../INCLUDE/Makefile.common

# DO NOT DELETE THIS LINE -- make depend depends on it.

I=$(CCINCDIR)
I2=$(BLDROOT)/INCLUDE

edg2len.$o:  edg2len.c  $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
edgunion.$o: edgunion.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
fconvert.$o: fconvert.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
killgrun.$o: killgrun.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     \
        $(I)/tsp.h      $(I)/edgegen.h  $(I)/bigguy.h   $(I)/lp.h       \
        $(I)/cut.h      $(I)/kdtree.h   $(I)/combs.h    
prob2tsp.$o: prob2tsp.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
showres.$o:  showres.c  $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     \
        $(I)/tsp.h      $(I)/edgegen.h  $(I)/bigguy.h   $(I)/lp.h       \
        $(I)/cut.h      $(I)/kdtree.h   $(I)/combs.h    
tourchk.$o:  tourchk.c  $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
tourlen.$o:  tourlen.c  $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
tour2edg.$o: tour2edg.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
perm2tsp.$o: perm2tsp.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     
tdivide.$o:  tdivide.c  $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     \ 
        $(I)/tsp.h      $(I)/edgegen.h  $(I)/bigguy.h   $(I)/lp.h       \
        $(I)/cut.h      $(I)/kdtree.h   $(I)/combs.h

