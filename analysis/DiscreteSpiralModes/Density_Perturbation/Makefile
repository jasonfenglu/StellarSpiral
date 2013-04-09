F90=ifort
FCC=icc

FFLAGS=  -O3 -L. -lslatec -fp-model strict  -prec-div -prec-sqrt -fp-stack-check -openmp  -limf

UTILS  = plot.o rk4.o k3sqrt.o
OBJECT =  $(UTILS) Density.o 
OBJECT2=  $(UTILS) force.o 
OBJECT3=  $(UTILS) SearchAll.o 
OBJECT4=  $(UTILS) FindOne.o 
OBJECT5=  $(UTILS) test.o

all: $(OBJECT) $(OBJECT2) $(OBJECT3) $(OBJECT4) $(OBJECT5) 
	$(F90) -o Density.exe $(OBJECT) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(F90) -o force.exe $(OBJECT2) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf -openmp
	$(F90) -o SearchAll.exe $(OBJECT3) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o FindOne.exe $(OBJECT4) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  
	$(F90) -o test.exe $(OBJECT5) $(FFLAGS) -lm -L$(PGPLOT_HOME) -lpgplot -lX11 -lpng -I$(INCLUDE) -L$(LIB) -liomp5 -lpthread -limf  

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@  $(INCL)
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe *.mod
