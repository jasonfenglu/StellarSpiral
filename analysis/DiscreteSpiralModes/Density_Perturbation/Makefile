F90=ifort
MPIF77=mpif77
FCC=icc

FFLAGS=  -O2 -fp-model strict  -prec-div -prec-sqrt -fp-stack-check -limf  -openmp -module obj -heap-arrays 0 -fpp

LIB   = -L. -lslatecomp -limf -lpthread -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -liomp5  -L$(HDF5HOME)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5
INCL  = -I$(HDF5_HOME)/include/

UTILS  = plot.o rk4.o k3sqrt.o io.o galaxy.o math.o
EXE    = test.exe Density.exe Force.exe Kendall.exe Corotation.exe Rot.exe cDensity.exe RelativeForce.exe StellarGas.exe GasDensity.exe Ecg.exe Follow.exe
MPIEXE = SearchAll.exe
TMPEXE = FindOne.exe  2Density.exe 

OUTILS := $(addprefix obj/,$(UTILS))

vpath obj/ .

all:    $(OUTILS) $(EXE)

mpi:
	$(MPIF77) $(FFLAGS) SearchAll.f90 -c -o obj/SearchAll.o $(LIB) $(INCL)
	$(MPIF77) $(FFLAGS) -o SearchAll.exe $(OUTILS) obj/SearchAll.o  $(LIB)

%.exe:  obj/%.o $(OUTILS)
	$(F90) $(FFLAGS) -o $@ $(OUTILS) $< $(LIB) 

obj/%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $@ $(INCL)
	
obj/%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $@ 

test:	$(OUTILS) test.exe

clean:
	rm -rf *.o *.exe *.mod bc.* obj/*
