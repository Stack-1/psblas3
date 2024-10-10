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
MATRIX_GENERATION_D 	= modules/matgen/psb_d_pde2d_mod
MATRIX_GENERATION_S 	= modules/matgen/psb_s_pde2d_mod
MIXED_VERSION_1 		= modules/psb_ds_cg_1
MIXED_VERSION_2			= modules/psb_ds_cg_2
DOUBLE_VERSION			= modules/psb_d_cg
SINGLE_VERSION			= modules/psb_s_cg
CG_FILE					= psb_dscg

CPU_OBJ					= $(MATRIX_GENERATION_D).o $(MATRIX_GENERATION_S).o $(MIXED_VERSION_1).o $(DOUBLE_VERSION).o $(SINGLE_VERSION).o $(CG_FILE).o
GPU_OBJ					= $(MATRIX_GENERATION_D).o $(MATRIX_GENERATION_S).o $(MIXED_VERSION_1).o $(MIXED_VERSION_2).o $(DOUBLE_VERSION).o $(SINGLE_VERSION).o 

EXEDIR=./runs


all: --runsd cpu gpu

--runsd:
	(if test ! -d runs ; then mkdir runs; fi)
	(if test ! -d modules ; then mkdir modules; fi)
	(if test ! -d data ; then mkdir data; fi)
	(if test ! -d data/cpu ; then mkdir data/cpu; fi)
	(if test ! -d data/gpu ; then mkdir data/gpu; fi)


cpu: --runsd
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MATRIX_GENERATION_D).f90 -o $(MATRIX_GENERATION_D).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MATRIX_GENERATION_S).f90 -o $(MATRIX_GENERATION_S).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MIXED_VERSION_1).f90 -o $(MIXED_VERSION_1).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MIXED_VERSION_2).f90 -o $(MIXED_VERSION_2).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(DOUBLE_VERSION).f90 -o $(DOUBLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(SINGLE_VERSION).f90 -o $(SINGLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(CG_FILE).f90 -o $(CG_FILE).o
	$(FLINK) $(LOPT) $(CCOPT) $(CPU_OBJ) $(MIXED_VERSION_2).o -o $(CG_FILE)_cpu $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv $(CG_FILE)_cpu $(EXEDIR)

gpu: --runsd
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MATRIX_GENERATION_D).f90 -o $(MATRIX_GENERATION_D).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c $(MATRIX_GENERATION_S).f90 -o $(MATRIX_GENERATION_S).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(MIXED_VERSION_1).f90 -o $(MIXED_VERSION_1).o
	g++ -O3 -c modules/single_to_double.cpp -o modules/single_to_double.o
	g++ -O3 -c modules/double_to_single.cpp -o modules/double_to_single.o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(MIXED_VERSION_2).f90 -o $(MIXED_VERSION_2).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(DOUBLE_VERSION).f90 -o $(DOUBLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(SINGLE_VERSION).f90 -o $(SINGLE_VERSION).o
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -cpp -I../../modules -I. -c -DHAVE_CUDA $(CG_FILE)_gpu.f90 -o $(CG_FILE)_gpu.o
	$(FLINK) $(LOPT) $(CCOPT) $(GPU_OBJ) modules/double_to_single.o modules/single_to_double.o $(CG_FILE)_gpu.o -o psb_dscg_gpu $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv $(CG_FILE)_gpu $(EXEDIR)

gpu-special:
	$(FLINK) $(LOPT) $(CCOPT) $(OBJ) modules/single_to_double.o $(CG_FILE)_gpu.o -o psb_dscg_gpu $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv $(CG_FILE)_gpu $(EXEDIR)

clean: 
	rm -f  modules/*.o *.o  *$(.mod) 
	rm -f $(EXEDIR)/psb_dscg_cpu
	rm -f $(EXEDIR)/psb_dscg_gpu 



