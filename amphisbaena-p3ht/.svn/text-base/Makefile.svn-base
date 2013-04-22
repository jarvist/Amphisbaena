# @file  Makefile
# @brief Makefile
#
# @author Mutsuo Saito (Hiroshima University)
# @author Makoto Matsumoto (Hiroshima University)
#
# Copyright (C) 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
# University. All rights reserved.
#
# The new BSD License is applied to this software.
# see LICENSE.txt
#
# @note
# We could comple test-sse2-Mxxx using gcc 3.4.4 of cygwin.
# We could comple test-sse2-Mxxx using gcc 4.0.1 of Linux.
# We coundn't comple test-sse2-Mxxx using gcc 3.3.2 of Linux.
# We could comple test-alti-Mxxx using gcc 3.3 of osx.
# We could comple test-alti-Mxxx using gcc 4.0 of osx.

WARN = #-Wmissing-prototypes -Wall #-Winline 
#WARN = -Wmissing-prototypes -Wall -W
OPTI = -O9 -finline-functions -fomit-frame-pointer -DNDEBUG \
-fno-strict-aliasing --param max-inline-insns-single=1800 
#OPTI = -O2  #conservative optimisations, JMF 20-8-07
#--param inline-unit-growth=500 --param large-function-growth=900 #for gcc 4
#STD =
#STD = -std=c89 -pedantic
#STD = -std=c99 -pedantic
STD = -std=c99
CC = gcc
CCFLAGS = $(OPTI) $(WARN) $(STD)
ALTIFLAGS = -mabi=altivec -maltivec -DHAVE_ALTIVEC
OSXALTIFLAGS = -faltivec -maltivec -DHAVE_ALTIVEC
SSE2FLAGS = -msse2 -DHAVE_SSE2


# ==========================================================
# comment out or EDIT following lines to get max performance
# ==========================================================
# --------------------
# for gcc 4
# --------------------
#CCFLAGS += --param inline-unit-growth=500 \
#--param large-function-growth=900
# --------------------
# for icl
# --------------------
#CC = icl /Wcheck /O3 /QxB /Qprefetch
# -----------------
# for PowerPC
# -----------------
#CCFLAGS += -arch ppc
# -----------------
# for Pentium M
# -----------------
#CCFLAGS += -march=prescott
# -----------------
# for Athlon 64
# -----------------
CCFLAGS += -march=athlon64

#.PHONY: std-check sse2-check alti-check

all: 
	${CC} ${CCFLAGS} ${SSE2FLAGS} -lm -o amphisbaena amphisbaena.c 

gorgophone:   
	${CC} ${CCFLAGS} ${SSE2FLAGS} -lm -o gorgophone gorgophone.c

gorgophone_steady: gorgophone_steady.c
	${CC} ${CCFLAGS} ${SSE2FLAGS} -lm -o gorgophone_steady gorgophone_steady.c

#profile: SFMT/SFMT.c SFMT/SFMT.h SFMT/SFMT-sse2.h SFMT/SFMT-params2281.h
#	${CC} ${STD} ${SSE2FLAGS} -DMEXP=2281 -lm -o marbles_profile -pg marbles.c


clean:
	rm -f *.o *~
