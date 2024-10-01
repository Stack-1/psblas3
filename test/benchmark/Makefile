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
CCOPT= -g -fcheck=all -Wall -fbacktrace
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).
OBJ = dot_vs_conv.o

EXEDIR=./runs


all: runsd dot_vs_conv

runsd:
	(if test ! -d runs ; then mkdir runs; fi)


dot_vs_conv: $(OBJ)
	$(FLINK) $(OBJ) $(CCOPT) -o dot_vs_conv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv dot_vs_conv $(EXEDIR)


clean: 
	/bin/rm -f dot_vs_conv.o  *$(.mod) \
	$(EXEDIR)/dot_vs_conv 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



