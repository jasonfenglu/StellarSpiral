F90=ifort
MPIF77=mpif77
FCC=icc

FFLAGS=  -O3 -fp-model strict  -prec-div -prec-sqrt -fp-stack-check -limf  -openmp -g

LIB   = -L. -lslatecomp -limf -lpthread -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -liomp5 

UTILS  = plot.o rk4.o k3sqrt.o
EXE    = FindOne.exe test.exe 2Density.exe 3Density.exe
MPIEXE = SearchAll.exe

OUTILS := $(addprefix obj/,$(UTILS))

vpath obj/ .

all:    $(OUTILS) mpi $(EXE) $(MPIEXE)

mpi:
	$(MPIF77) $(FFLAGS) SearchAll.f90 -c -o obj/SearchAll.o $(LIB)
	$(MPIF77) $(FFLAGS) -o SearchAll.exe $(OUTILS) obj/SearchAll.o  $(LIB)

%.exe:  obj/%.o $(OUTILS)
	$(F90) $(FFLAGS) -o $@ $(OUTILS) $< $(LIB) 

obj/%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $@  $(LIB)
	
obj/%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $@ 

clean:
	rm -rf *.o *.exe *.mod bc.* obj/*
