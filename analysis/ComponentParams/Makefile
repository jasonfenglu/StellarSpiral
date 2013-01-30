F90=ifort
FCC=icc

FFLAGS =-O2 -L. -lslatec -g

PGPLOT =-L$(PGPLOT_DIR) -lpgplot  -limf -lm -lX11 -lpng


OBJECT2 =  rk4.o k3sqrt.o rc.o
all: $(OBJECT) $(OBJECT2)
	$(F90) -o rc.exe $(OBJECT2) $(FFLAGS) $(PGPLOT)

%.o:%.f90
	$(F90) $(FFLAGS)  -c $< -o $@ 
	
%.o:%.f
	$(F90) $(FFLAGS)  -c $< -o $@ 

%.o:%.c
	$(FCC) $(FFLAGS)  -c $< -o $@ 

clean:
	rm -rf *.o *.exe *.mod
