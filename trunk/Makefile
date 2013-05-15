#################################
#  Antares Makefile 9/4/2011 
#  Hsiang-Hsu Wang
#  hhwang@asiaa.sinica.edu.tw
#################################
NDIM=2
Rie_type=ExactRie
Gravity=
EXEC=antares
AntaresHOME=$(PWD)

#################################
###### LAM/MPI
################################# 
F90=mpif77
MPIHOME=$(MPI_HOME)
FFTWHOME=$(FFTW_HOME)
HDF5HOME=$(HDF5_HOME)
AntaresHOME=$(PWD)
 
FFLAGS= -cpp -O2 -DNDIM=$(NDIM) -DFORCE -DVERBOSE -DISOTHERMAL -D$(Rie_type) -module bin -DRK2  -prec-div -prec-sqrt 

LIB = -L$(MPIHOME)/lib -lmpi -lz -limf -lm -L$(FFTWHOME)/lib -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw -L$(HDF5HOME)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -fPIC -lm -limf -L./src -lslatecomp -L./bin/ -lk3
INCL=-I$(MPIHOME)/include -I$(FFTWHOME)/include -I$(HDF5HOME)/include

vpath *.o bin

VPATH = ./src ./lib
OBJECT= $(k3) common_params.o exceptions.o find_dt.o simcontroll.o main.o init2d.o output1d.o output2d.o bdry2d.o write_collective2d.o dt_loc2d.o dt_all.o sendbdry2d.o sloplimiter.o riemann2d.o bndflux2d.o fluxRieIso2d.o fluxHLLC2d.o force2d.o sgkernel2d.o selfgravity2d.o write_gravity2d.o read_collective2d.o read1d.o read2d.o

k3 = rk4.o k3sqrt.o  galaxy.o 

antares.exe: $(addprefix bin/,$(OBJECT))
	$(F90) -o $(EXEC)$(NDIM)d $(addprefix bin/,$(OBJECT)) $(LIB)

bin/%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $@ $(INCL)

clean: 
	-rm bin/*.o $(EXEC)*d bin/*.mod 

rm:
	rm *.h5 
