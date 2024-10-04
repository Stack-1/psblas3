INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(INSTALLDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) $(LCUDA) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base -lpsb_ext 
LDLIBS=$(PSBLDLIBS)

LDLIBS=$(PSBGPULDLIBS) 

#
# Compilers and such
#
CCOPT= -g -fcheck=all -Wall -fbacktrace -cpp
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

#
# OBJ files
#
MATRIX_GENERATION 	= ./modules/psb_ds_gen_mod
MIXED_VERSION_1 	= ./modules/psb_ds_cg_1
MIXED_VERSION_2		= ./modules/psb_ds_cg_2
DOUBLE_VERSION		= ./modules/psb_d_cg
SINGLE_VERSION		= ./modules/psb_s_cg
CG_FILE				= psb_dscg

OBJ					= $(MIXED_VERSION_1).o $(MIXED_VERSION_2).o $(DOUBLE_VERSION).o $(SINGLE_VERSION).o 

EXEDIR=./runs


all: --runsd cpu gpu

--runsd:
	(if test ! -d runs ; then mkdir runs; fi)
	(if test ! -d modules ; then mkdir modules; fi)
	(if test ! -d data ; then mkdir data; fi)
	(if test ! -d data/cpu ; then mkdir data/cpu; fi)
	(if test ! -d data/gpu ; then mkdir data/gpu; fi)


cpu: 
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MATRIX_GENERATION).f90 -o $(MATRIX_GENERATION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MIXED_VERSION_1).f90 -o $(MIXED_VERSION_1).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MIXED_VERSION_2).f90 -o $(MIXED_VERSION_2).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(DOUBLE_VERSION).f90 -o $(DOUBLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(SINGLE_VERSION).f90 -o $(SINGLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(CG_FILE).f90 -o $(CG_FILE).o
	$(FLINK) $(LOPT) $(CCOPT) $(OBJ) $(MATRIX_GENERATION).o $(CG_FILE).o -o $(CG_FILE)_cpu $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv $(CG_FILE)_cpu $(EXEDIR)

gpu: 
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c psb_d_pde3d_mod.f90 -o psb_d_pde3d_mod.o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c psb_s_pde3d_mod.f90 -o psb_s_pde3d_mod.o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(MIXED_VERSION_1).f90 -o $(MIXED_VERSION_1).o
	g++ -O3 -c modules/single_to_double.cpp -o modules/single_to_double.o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(MIXED_VERSION_2).f90 -o $(MIXED_VERSION_2).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(DOUBLE_VERSION).f90 -o $(DOUBLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(SINGLE_VERSION).f90 -o $(SINGLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(CG_FILE)_gpu.f90 -o $(CG_FILE)_gpu.o
	$(FLINK) $(LOPT) $(CCOPT) $(OBJ) psb_s_pde3d_mod.o psb_d_pde3d_mod.o modules/single_to_double.o $(CG_FILE)_gpu.o -o psb_dscg_gpu $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv $(CG_FILE)_gpu $(EXEDIR)

gpu-special:
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(SINGLE_VERSION).f90 -o $(SINGLE_VERSION).o
	$(FLINK) $(LOPT) $(CCOPT) $(OBJ) modules/single_to_double.o $(CG_FILE)_gpu.o -o psb_dscg_gpu $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv $(CG_FILE)_gpu $(EXEDIR)

clean: 
	rm -f  modules/*.o *.o  *$(.mod) 
	rm -f $(EXEDIR)/psb_dscg_cpu
	rm -f $(EXEDIR)/psb_dscg_gpu 



