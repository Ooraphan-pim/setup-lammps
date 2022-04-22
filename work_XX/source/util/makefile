#***************************************************************
#
# make file for directory EAM/util/build_mpi on fenrir
#
#***************************************************************
#

INTEL_C_ROOT=/opt/intel/Compiler/11
INTEL_C_BIN=$(INTEL_C_ROOT)/bin/intel64/icpc
INTEL_C_LIB=$(INTEL_C_ROOT)/lib/intel64
INTEL_C_INC=$(INTEL_C_ROOT)/include

INTEL_MKL_ROOT=/opt/intel/Compiler/11/mkl
INTEL_MKL_LIB=$(INTEL_MKL_ROOT)/lib/em64t
INTEL_MKL_INC=$(INTEL_MKL_ROOT)/include

INTEL_F_ROOT=/opt/intel/Compiler/11
INTEL_F_BIN=$(INTEL_F_ROOT)/bin/intel64/ifort
INTEL_F_LIB=$(INTEL_F_ROOT)/lib/intel64
INTEL_F_INC=$(INTEL_F_ROOT)/include

CUDA_ROOT=/opt/nvidia/cuda/3.0
CUDA_LIB=$(CUDA_ROOT)/lib64

#MPI_ROOT=/opt/mpich2/ch3_mrail_gen2-intel11
#MPI_LIB_NAME=-lmpich

MPI_ROOT=/opt/openmpi
MPI_LIB_NAME=-lmpi

MPI_C=$(MPI_ROOT)/bin/mpicxx
MPI_LIB=$(MPI_ROOT)/lib
MPI_INCLUDE=$(MPI_ROOT)/include


GET = bad_juju_here

%.o : %.f90
	$(F90) -c $(FFLAGS) $< -o $@
%.o : %.F90
	$(F90) -c $(FFLAGS) $< -o $@


NR_DIR = /usr2/olmsted/EAM/nrecipes
SVD_OBJS = $(NR_DIR)/svdfit.o $(NR_DIR)/svdcmp.o $(NR_DIR)/svbksb.o

OBJS = d_plane.o e_of_r.o f_to_fdat.o f_to_grad.o

LIB_OBJS = eam_util.o my_util.o repeat_dyn.o ai_gen_data.o read_r_file.o \
           eam_file_util.o crystal.o scatter_randomly.o
# compilation flags for all cases
FC = mpif90
F90C = mpif90
F90=OMPI_FC=/storage/home/singhanat.rer/intel/Compiler/11.0/084/bin/intel64/ifort mpif90

## alpha flags
# -tune generic is on by  default
## ALL_FLAGS = -i4 -r8 -tune generic -warn argument_checking -warn declarations
## DEBUG_FLAGS = -g -O0 -check bounds -check underflow -check overflow -ladebug
## PROD_FLAGS = -O -check bounds
## MATH_LIB = -ldxml

# SGI flags 
## ALL_FLAGS = -i4 -r8 
## DEBUG_FLAGS = -g -check_bounds 
## PROD_FLAGS = -O -check_bounds
## MATH_LIB = -lscs

## ALL_FLAGS = -Wall -Wimplicit-none -i4 -r8 
## DEBUG_FLAGS = -g -fbounds-check
## PROD_FLAGS = -O -fbounds-check

ALL_FLAGS = -assume nounderscore  -warn all  -warn noerrors \
            -check bounds -check format -check uninit       \
            -real-size 64 -DDGELSADDU
DEBUG_FLAGS = -O0 -g ${ALL_FLAGS}
PROD_FLAGS = ${ALL_FLAGS}



FFLAGS = $(DEBUG_FLAGS) 
PREF_PROD_FFLAGS = $(PROD_FLAGS) 
OFLAG =  
KFLAG = 

# flags for link 
LDFLAGS = $(FFLAGS)
PREF_PROD_LDFLAGS = $(PREF_PROD_FFLAGS)

my_util.o : my_util.F90 makefile
	$(FC) -c $(FFLAGS) my_util.F90

d_plane.o : d_plane.f90 makefile my_util.o 
	$(FC) -c $(FFLAGS) d_plane.f90

d_plane : d_plane.o my_util.o makefile
	$(FC) $(LDFLAGS) d_plane.o my_util.o $(MATH_LIB) -o d_plane

bin_energy.o : bin_energy.F90 makefile
	$(FC) -c $(FFLAGS) bin_energy.F90

bin_energy : bin_energy.F90 my_util.o eam_file_util.o makefile
	$(FC) $(LDFLAGS) bin_energy.F90 -DMAIN -o bin_energy

e_of_r.o : e_of_r.F90 makefile
	$(FC) -c $(FFLAGS) e_of_r.F90

e_of_r_full.o : e_of_r.F90 makefile
	$(FC) -c $(FFLAGS) -DFULL e_of_r.F90 -o e_of_r_full.o

e_of_r : e_of_r.o 
	$(FC) $(LDFLAGS) e_of_r.o -o e_of_r

e_of_r_full : e_of_r_full.o 
	$(FC) $(LDFLAGS) e_of_r_full.o -o e_of_r_full

f_to_fdat.o : f_to_fdat.f90 makefile
	$(FC) -c $(FFLAGS) f_to_fdat.f90

f_to_fdat : f_to_fdat.o 
	$(FC) $(LDFLAGS) f_to_fdat.o -o f_to_fdat

f_to_grad.o : f_to_grad.f90 makefile
	$(FC) -c $(FFLAGS) f_to_grad.f90

f_to_grad : f_to_grad.o 
	$(FC) $(LDFLAGS) f_to_grad.o -o f_to_grad

ridge_to_dyn.o : ridge_to_dyn.f90 makefile
	$(FC) -c $(FFLAGS) ridge_to_dyn.f90

ridge_to_dyn : ridge_to_dyn.o 
	$(FC) $(LDFLAGS) ridge_to_dyn.o -o ridge_to_dyn

ridge_to_dyn2.o : ridge_to_dyn2.f90 makefile
	$(FC) -c $(FFLAGS) ridge_to_dyn2.f90

ridge_to_dyn2 : ridge_to_dyn2.o 
	$(FC) $(LDFLAGS) ridge_to_dyn2.o -o ridge_to_dyn2

dyn_to_dyn2.o : dyn_to_dyn2.f90 makefile
	$(FC) -c $(FFLAGS) dyn_to_dyn2.f90

dyn_to_dyn2 : dyn_to_dyn2.o 
	$(FC) $(LDFLAGS) dyn_to_dyn2.o -o dyn_to_dyn2

mix_dyn_to_dyn2.o : mix_dyn_to_dyn2.f90 makefile
	$(FC) -c $(FFLAGS) mix_dyn_to_dyn2.f90

mix_dyn_to_dyn2 : mix_dyn_to_dyn2.o 
	$(FC) $(LDFLAGS) mix_dyn_to_dyn2.o -o mix_dyn_to_dyn2

make_saddle_dyn.o : make_saddle_dyn.f90 makefile
	$(FC) -c $(FFLAGS) make_saddle_dyn.f90

make_saddle_dyn : make_saddle_dyn.o 
	$(FC) $(LDFLAGS) make_saddle_dyn.o -o make_saddle_dyn

f_of_r.o : f_of_r.f90 makefile
	$(FC) -c $(FFLAGS) f_of_r.f90

f_of_r : f_of_r.o 
	$(FC) $(LDFLAGS) f_of_r.o -o f_of_r

extract_eam_alats.o : extract_eam_alats.f90 makefile eam_util.o
	$(FC) -c $(FFLAGS) extract_eam_alats.f90

extract_eam_alats : extract_eam_alats.o eam_util.o 
	$(FC) $(LDFLAGS) extract_eam_alats.o eam_util.o -o extract_eam_alats

eam_util.o : eam_util.F90 makefile
	$(FC) -c $(FFLAGS) eam_util.F90

repeat_dyn.o : repeat_dyn.f90 eam_util.o makefile ../dyn87/libdyn87.a 
	$(FC) -c $(FFLAGS) -I../dyn87 repeat_dyn.f90

fit_util.o : fit_util.f90 makefile
	$(FC) -c $(FFLAGS) fit_util.f90

fit_one_file.o : fit_one_file.f90 makefile fit_util.o
	$(FC) -c $(FFLAGS) fit_one_file.f90

fit_one_file : fit_one_file.o fit_util.o makefile $(SVD_OBJS)
	$(FC) $(LDFLAGS) fit_one_file.o fit_util.o  $(SVD_OBJS) -o fit_one_file

read_r_file.o : read_r_file.f90 makefile
	$(FC) -c $(FFLAGS) read_r_file.f90

ai_gen_data.o : ai_gen_data.f90 makefile
	$(FC) -c $(FFLAGS) ai_gen_data.f90

crystal.o : crystal.f90 makefile
	$(FC) -c $(FFLAGS) crystal.f90

scatter_randomly.o : scatter_randomly.f90 makefile
	$(FC) -c $(PREF_PROD_FFLAGS) scatter_randomly.f90

## read_r_file2.o : read_r_file2.f90 makefile
## 	$(FC) -c $(FFLAGS) read_r_file2.f90

eam_file_util.o : eam_file_util.f90 makefile
	$(FC) -c $(FFLAGS) eam_file_util.f90

read_trans_file.o : read_trans_file.f90 makefile
	$(FC) -c $(FFLAGS) read_trans_file.f90

libmy_utils.a : $(LIB_OBJS)
	rm -f libmy_utils.a
	ar -r libmy_utils.a $(LIB_OBJS) 

#  clean up the object files
clean: 
	-rm -f $(OBJS)






