#################################
###### LAM/MPI
################################# 
F90=ifort
MPIHOME=$(MPI_HOME)
FFTWHOME=$(FFTW_HOME)
HDF5HOME=$(HDF5_HOME)
PGPLOTHOME=$(PGPLOT_HOME)

bin_prefix=

FFLAGS=  -O3 

LIB = -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -fPIC -lpgplot -lX11 -limf -lm -lpng
INCL=-I$(MPIHOME)/include -I$(FFTWHOME)/include -I$(HDF5HOME)/include -L$(PGPLOT_HOME)


OBJECT = \
	main.o

plot.exe: $(addprefix $(bin_prefix),$(OBJECT))
	$(F90) -o plot.exe $(addprefix $(bin_prefix),$(OBJECT)) $(LIB) $(INCL)

$(bin_prefix)%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $@ $(INCL) $(LIB)

clean: 
	rm $(bin_prefix)*.o $(EXEC)*d $(bin_prefix)*.mod  *.exe

rm:
	rm *.h5 
