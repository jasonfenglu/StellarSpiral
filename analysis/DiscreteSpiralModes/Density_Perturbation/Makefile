F90=ifort
FCC=icc

FFLAGS=  -O2 -g -lm  -L. -lslatec -limf -fp-model strict  -prec-div -prec-sqrt -fp-stack-check
INCL  = -I$(HDF5HOME)/include 

HDF5HOME=$(HDF5_HOME)


OBJECT = plot.o rk4.o k3sqrt.o Density.o
OBJECT2 = plot.o rk4.o k3sqrt.o force.o

all: $(OBJECT) $(OBJECT2)
	$(F90) -o Density.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(F90) -o force.exe $(OBJECT2) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@ 
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe
