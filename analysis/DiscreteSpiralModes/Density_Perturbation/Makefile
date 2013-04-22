F90=ifort
MPIF77=mpif77
FCC=icc

FFLAGS=  -O3 -L. -lslatecomp -fp-model strict  -prec-div -prec-sqrt -fp-stack-check -limf -openmp

UTILS  = plot.o rk4.o k3sqrt.o
OBJECT =  $(UTILS) Density.o 
OBJECT3=  $(UTILS) SearchAll.o 
OBJECT4=  $(UTILS) FindOne.o 
OBJECT5=  $(UTILS) test.o
OBJECT7=  $(UTILS) 2Density.o


all: $(UTILS) mpi $(OBJECT) $(OBJECT2) $(OBJECT3) $(OBJECT4) $(OBJECT5) $(OBJECT6) $(OBJECT7)
	$(F90) -o Density.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(MPIF77) -o SearchAll.exe $(OBJECT3) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o FindOne.exe $(OBJECT4) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o test.exe $(OBJECT5) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o 2Density.exe $(OBJECT7) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  

mpi:
	$(MPIF77) $(FFLAGS) SearchAll.f90 -c $(INCL)

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@  $(INCL)
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe *.mod bc.*
