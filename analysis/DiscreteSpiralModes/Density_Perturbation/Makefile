F90=ifort
FCC=icc

FFLAGS=  -O3  -g -lm  -L. -lslatec -limf -L$(FFTW_HOME)/lib -lfftw3
INCL  = -I$(HDF5HOME)/include -I$(FFTW_INCL)

HDF5HOME=$(HDF5_HOME)


OBJECT = rk4.o k3sqrt.o Density.o

all: $(OBJECT) $(OBJECT2)
	$(F90) -o Density.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@  $(INCL) 
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe
