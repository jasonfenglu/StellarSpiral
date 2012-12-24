F90=ifort
FCC=icc

FFLAGS=  -O2  -L. -lslatec  $(shell gsl-config --cflags)
INCL  = -I$(HDF5HOME)/include

HDF5HOME=$(HDF5_HOME)


#OBJECT = rk4.o k3sqrt.o SearchAll.o
OBJECT2 = rk4.o k3sqrt.o FindOne.o wrapper.o

all: $(OBJECT) $(OBJECT2)
	#$(F90) -o SearchAll.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(F90) -o FindOne.exe $(OBJECT2) $(FFLAGS) $(shell gsl-config --libs)

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@  $(INCL) 
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe
