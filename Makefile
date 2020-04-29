NVCC      = nvcc
CC = gcc
GF = gfortran
GFFLAGS = -lgsl -lgslcblas -O3 -llapack -lblas -g -cpp -Wall -Wtabs
CFLAGS =  -lgsl -lgslcblas -O3  -static -g 
FFTLIB = -I/usr/include -lfftw3 -llapack -lblas -lm  


LDOPTS    = -fopenmp
OPTS      = -O3 -DADD_ -x f95-cpp-input  -Dmagma_devptr_t="integer(kind=8)" 
F77OPTS   = -O3 -DADD_ -fPIC
FOPTS     = -O3 -DADD_ -x f95-cpp-input -fPIC
NVOPTS    = -O3 -DADD_ -Xcompiler -fno-strict-aliasing	
LDOPTS    = -fopenmp

# Depending on how ATLAS and LAPACK were compiled, you may need one or more of:
# -lifcore -ldl -lf2c -lgfortran
LIB       = -llapack -lf77blas -latlas -lcblas -lcublas -lcudart -lstdc++ -lm -lgfortran -lf2c -fPIC
LIB_linux = -llapack -lf77blas -latlas -lcblas -lstdc++ -lm -lgfortran -lf2c -fPIC

# define library directories here or in your environment
MAGMADIR  = $(HOME)/local/magma
LAPACKDIR = /usr/lib64
ATLASDIR  = /usr/lib64/atlas
CUDADIR   = /usr/local/cuda

LIBDIR    = -L$(LAPACKDIR) \
            -L$(ATLASDIR)/lib \
            -L$(CUDADIR)/lib64

INC       = -I$(CUDADIR)/include \
	    -I$(HOME)/local/magma/include	

ifeq ($(FORT), pgfortran)
	FOBJ  = fortran_thunking.o
else
	FOBJ  = fortran.o
endif
#SOURCE_DIR = /home/german/projects/ED-EntanglementSpectrum-Floquet/src/
SOURCE_DIR = /home/german/Dropbox/Programs/ED-EntanglementSpectrum-FloquetV3/src/
#SOURCE_DIR = /home/gs249/Dropbox/Programs/ED-EntanglementSpectrum-FloquetV3/src/

#all: ESFS_GPU ESFS ESFS_THUNKING

all:ESFS_linux ESFS
# $(SOURCE_DIR)fortran_thunking.o: $(SOURCE_DIR)fortran_thunking.c
#	nvcc -O3 -DCUBLAS_USE_THUNKING -DCUBLAS_GFORTRAN -I/usr/local/cuda/include  -c  $(SOURCE_DIR)fortran_thunking.c -o $(SOURCE_DIR)fortran_thunking.o


# $(SOURCE_DIR)fortran.o: $(SOURCE_DIR)fortran.c
#	nvcc -O3 -DCUBLAS_GFORTRAN -I/usr/local/cuda/include  -c  $(SOURCE_DIR)fortran.c -o $(SOURCE_DIR)fortran.o


ESFS_linux: $(SOURCE_DIR)Modules.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)quick-sort-index-table.o  $(SOURCE_DIR)lapack_routines.o  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB_.o $(SOURCE_DIR)density.o  $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 
	$(GF)  $(LDOPTS) -o $@  $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.o $(SOURCE_DIR)lapack_routines.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB_.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o  $(SOURCE_DIR)density.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 $(LIBDIR)	$(LIB_linux)


ESFS: $(SOURCE_DIR)fortran.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)quick-sort-index-table.o  $(SOURCE_DIR)lapack_routines_GPU.o  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB_.o $(SOURCE_DIR)density.o  $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 
	$(GF)  $(LDOPTS) -o $@  $(SOURCE_DIR)fortran.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.o $(SOURCE_DIR)lapack_routines_GPU.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB_.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o  $(SOURCE_DIR)density.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 -L$(MAGMADIR)/lib -lmagma \
	$(LIBDIR) \
	$(LIB)


ESFS_THUNKING: $(SOURCE_DIR)fortran_thunking.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)quick-sort-index-table.o  $(SOURCE_DIR)lapack_routines_GPU.o  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB_TH.o $(SOURCE_DIR)density.o  $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 
	$(GF)  $(LDOPTS)  -o $@  $(SOURCE_DIR)fortran_thunking.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.o $(SOURCE_DIR)lapack_routines_GPU.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB_TH.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o  $(SOURCE_DIR)density.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 -L$(MAGMADIR)/lib -lmagma \
	$(LIBDIR) \
	$(LIB)  -Wall


ESFS_GPU: $(SOURCE_DIR)fortran.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)quick-sort-index-table.o  $(SOURCE_DIR)lapack_routines_GPU.o  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.o $(SOURCE_DIR)density.o  $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 
	$(GF)  $(LDOPTS)    -o $@  $(SOURCE_DIR)fortran.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.o $(SOURCE_DIR)lapack_routines_GPU.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)a_dagger.o $(SOURCE_DIR)ManyBodySpectrum.o $(SOURCE_DIR)TwoBodyEnergySpectrum.o  $(SOURCE_DIR)density.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 -L$(MAGMADIR)/lib -lmagma \
	$(LIBDIR) \
	$(LIB)  -Wall

$(SOURCE_DIR)density.o: $(SOURCE_DIR)density.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)density.f90 $(GFFLAGS)  


$(SOURCE_DIR)Modules.o: $(SOURCE_DIR)Modules.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)Modules.f90 $(GFFLAGS)  

$(SOURCE_DIR)HilbertDimension.o: $(SOURCE_DIR)HilbertDimension.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)HilbertDimension.f90 $(GFFLAGS)  

$(SOURCE_DIR)STATEStoPARTITIONmap.o: $(SOURCE_DIR)STATEStoPARTITIONmap.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)STATEStoPARTITIONmap.f90 $(GFFLAGS)  

$(SOURCE_DIR)quick-sort-index-table.o: $(SOURCE_DIR)quick-sort-index-table.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)quick-sort-index-table.f90 $(GFFLAGS)  

$(SOURCE_DIR)lapack_routines.o: $(SOURCE_DIR)lapack_routines.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)lapack_routines.f90 $(GFFLAGS)  

$(SOURCE_DIR)lapack_routines_GPU.o: $(SOURCE_DIR)lapack_routines_GPU.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)lapack_routines_GPU.f90 $(GFFLAGS)  

$(SOURCE_DIR)INTERFACE.o:$(SOURCE_DIR)INTERFACE.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)INTERFACE.f90 $(GFFLAGS)  

$(SOURCE_DIR)Basis.o:$(SOURCE_DIR)Modules.o $(SOURCE_DIR)Basis.f90 
	$(GF)  -c -o $@ $(SOURCE_DIR)Basis.f90 $(GFFLAGS)  

$(SOURCE_DIR)FLOQUET_OPERATOR_SUB.o:$(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90
	$(GF) -c -o $@ $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90 -DCUBLAS $(GFFLAGS)  

$(SOURCE_DIR)FLOQUET_OPERATOR_SUB_.o:$(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90 $(GFFLAGS)  

$(SOURCE_DIR)FLOQUET_OPERATOR_SUB_TH.o:$(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90 -DCUBLAS -DCUBLAS_USE_THUNKING $(GFFLAGS)  

$(SOURCE_DIR)subset.o:$(SOURCE_DIR)subset.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)subset.f90

$(SOURCE_DIR)ManyBodySpectrum.o:$(SOURCE_DIR)ManyBodySpectrum.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)ManyBodySpectrum.f90 -Wall -Wtabs

$(SOURCE_DIR)TwoBodyEnergySpectrum.o:$(SOURCE_DIR)TwoBodyEnergySpectrum.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)TwoBodyEnergySpectrum.f90 -Wall -Wtabs

$(SOURCE_DIR)a_dagger.o:$(SOURCE_DIR)a_dagger.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)a_dagger.f90 

clean:
	rm $(SOURCE_DIR)*.o *.mod ESFS_GPU *~ ESFS
