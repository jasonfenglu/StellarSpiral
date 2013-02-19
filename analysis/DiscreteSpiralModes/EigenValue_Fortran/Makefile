F90=ifort
FCC=icc

FFLAGS=  -L. -lslatec -openmp -fp-model precise -O3
#INCL  = -I$(HDF5HOME)/include -I$(MPI_HOME)/include

HDF5HOME=$(HDF5_HOME)


OBJECT = plot.o rk4.o k3sqrt.o SearchAll.o
OBJECT2 = rk4.o k3sqrt.o FindOne.o
OBJECT3 = rk4.o k3sqrt.o test_openmp.o

all: $(OBJECT) $(OBJECT2) $(OBJECT3)
	$(F90) -o SearchAll.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o FindOne.exe $(OBJECT2) $(FFLAGS)
	$(F90) -o Test_Openmp.exe $(OBJECT3) $(FFLAGS)
%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@  $(INCL)
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe *.mod
