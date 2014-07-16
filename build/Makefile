CC = gcc
GF = gfortran
GFFLAGS = -lgsl -lgslcblas -O3 -llapack -lblas -g
CFLAGS =  -lgsl -lgslcblas -O3  -static -g
FFTLIB = -I/usr/include -lfftw3 -llapack -lblas -lm 


SOURCE_DIR = /home/german/projects/ExactDiagonalizationBH/src/
#SOURCE_DIR = /home/gs249/Dropbox/Programs/ExactDiagonalizationBH/src/

all: Boson2DLattice floquetspectrum entanglementBH1D ES2D ESFS

#ESFS: $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90  $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)quick-sort-index-table.f90  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)Modules.o  $(SOURCE_DIR)zufall_sub.o $(SOURCE_DIR)STATEStoPARTITIONmap.f90
#	$(GF)  -o $@   $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.f90 $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)STATEStoPARTITIONmap.f90 $(SOURCE_DIR)zufall_sub.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 $(GFFLAGS) 

ESFS: $(SOURCE_DIR)Modules.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)quick-sort-index-table.o  $(SOURCE_DIR)lapack_routines.o  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)zufall_sub.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 
	$(GF)  -o $@   $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.o $(SOURCE_DIR)lapack_routines.o $(SOURCE_DIR)zufall_sub.o $(SOURCE_DIR)Basis.o $(SOURCE_DIR)INTERFACE.o $(SOURCE_DIR)STATEStoPARTITIONmap.o $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.o $(SOURCE_DIR)subset.o $(SOURCE_DIR)EntanglementSpectrumFloquetStates.f90 $(GFFLAGS) 



ES2D: $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)EntanglementSpectrum2D.f90  $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)quick-sort-index-table.f90  $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)Modules.o
	$(GF)  -o $@   $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)quick-sort-index-table.f90 $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)EntanglementSpectrum2D.f90 $(GFFLAGS) 

entanglementBH1D: $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)fftw_subroutine.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)ENTANGLEMENT_BH_1D.f90
	$(GF)  -o $@ $(SOURCE_DIR)fftw_subroutine.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)ENTANGLEMENT_BH_1D.f90 $(FFTLIB)


floquetspectrum: $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)fftw_subroutine.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)FLOQUET_OPERATOR.f90
	$(GF)  -o $@ $(SOURCE_DIR)fftw_subroutine.o $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)FLOQUET_OPERATOR.f90 $(FFTLIB)

#BosonLattice: $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)BosonLattice1D.f90 $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)quick-sort-index-table.f90 $(SOURCE_DIR)Hamiltonian.f90 $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)Modules.o
#	$(GF)  -o $@   $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)Hamiltonian.f90 $(SOURCE_DIR)quick-sort-index-table.f90 $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)BosonLattice1D.f90 $(GFFLAGS) 

Boson2DLattice: $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)Boson2DLattice.f90  $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)quick-sort-index-table.f90 $(SOURCE_DIR)Hamiltonian.f90 $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)Modules.o
	$(GF)  -o $@   $(SOURCE_DIR)Modules.o $(SOURCE_DIR)HilbertDimension.o $(SOURCE_DIR)Hamiltonian.f90 $(SOURCE_DIR)quick-sort-index-table.f90 $(SOURCE_DIR)Basis.f90 $(SOURCE_DIR)INTERFACE.f90 $(SOURCE_DIR)Boson2DLattice.f90 $(GFFLAGS) 

$(SOURCE_DIR)fftw_subroutine.o: $(SOURCE_DIR)fftw_subroutine.f90
	$(GF) -c -o $@ $(SOURCE_DIR)fftw_subroutine.f90 $(FFTLIB)

$(SOURCE_DIR)Modules.o: $(SOURCE_DIR)Modules.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)Modules.f90 $(GFFLAGS)  

$(SOURCE_DIR)HilbertDimension.o: $(SOURCE_DIR)HilbertDimension.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)HilbertDimension.f90 $(GFFLAGS)  


$(SOURCE_DIR)STATEStoPARTITIONmap.o: $(SOURCE_DIR)STATEStoPARTITIONmap.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)STATEStoPARTITIONmap.f90 $(GFFLAGS)  

$(SOURCE_DIR)zufall_sub.o: $(SOURCE_DIR)zufall_sub.f
	$(GF)  -c -o $@ $(SOURCE_DIR)zufall_sub.f $(GFFLAGS)  


$(SOURCE_DIR)quick-sort-index-table.o: $(SOURCE_DIR)quick-sort-index-table.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)quick-sort-index-table.f90 $(GFFLAGS)  

$(SOURCE_DIR)lapack_routines.o: $(SOURCE_DIR)lapack_routines.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)lapack_routines.f90 $(GFFLAGS)  

$(SOURCE_DIR)INTERFACE.o:$(SOURCE_DIR)INTERFACE.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)INTERFACE.f90 $(GFFLAGS)  

$(SOURCE_DIR)Basis.o:$(SOURCE_DIR)Modules.o $(SOURCE_DIR)Basis.f90 
	$(GF)  -c -o $@ $(SOURCE_DIR)Basis.f90 $(GFFLAGS)  

$(SOURCE_DIR)FLOQUET_OPERATOR_SUB.o:$(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90
	$(GF)  -c -o $@ $(SOURCE_DIR)FLOQUET_OPERATOR_SUB.f90 $(GFFLAGS)  

$(SOURCE_DIR)subset.o:$(SOURCE_DIR)subset.f90
	$(GF) -c -o $@ $(SOURCE_DIR)subset.f90

clean:
	rm $(SOURCE_DIR)*.o *.mod BosonLattice floquetspectrum Boson2DLattice ES2D entanglementBH1D ESFS