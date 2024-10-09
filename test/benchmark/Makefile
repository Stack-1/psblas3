INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(INSTALLDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g -fcheck=all -Wall -fbacktrace -O3
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

EXEDIR=./runs


all: runsd dot_vs_conv geaxpby_vs_conv

runsd:
	(if test ! -d runs ; then mkdir runs; fi)


dot_vs_conv: dot_vs_conv.o 
	$(FLINK) dot_vs_conv.o  $(CCOPT) -o dot_vs_conv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv dot_vs_conv $(EXEDIR)

geaxpby_vs_conv: geaxpby_vs_conv.o
	$(FLINK) geaxpby_vs_conv.o $(CCOPT) -o geaxpby_vs_conv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv geaxpby_vs_conv $(EXEDIR)

clean: 
	/bin/rm -f dot_vs_conv.o  *$(.mod) \
	$(EXEDIR)/dot_vs_conv 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



