F90=ifort
FCC=icc

FFLAGS=  -O2 -L. -lslatec -fp-model strict  -prec-div -prec-sqrt -fp-stack-check -openmp -g

OBJECT = plot.o rk4.o k3sqrt.o Density.o
OBJECT2= plot.o rk4.o k3sqrt.o force.o
OBJECT3= plot.o rk4.o k3sqrt.o SearchAll.o
OBJECT4= plot.o rk4.o k3sqrt.o FindOne.o

all: $(OBJECT) $(OBJECT2) $(OBJECT3) $(OBJECT4) $(OBJECT5) 
	$(F90) -o Density.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(F90) -o force.exe $(OBJECT2) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(F90) -o SearchAll.exe $(OBJECT3) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o FindOne.exe $(OBJECT4) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@  $(INCL)
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe *.mod
