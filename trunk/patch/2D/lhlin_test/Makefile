#################################
#  Antares Makefile 9/4/2011 
#  Hsiang-Hsu Wang
#  hhwang@asiaa.sinica.edu.tw
#################################
NDIM=2
Rie_type=ExactRie
Gravity=
PATCH=$(PWD)
EXEC=antares
#################################
###### LAM/MPI
################################# 
F90=mpif77

MPIHOME=$(MPI_HOME)
FFTWHOME=$(FFTW_HOME)
HDF5HOME=$(HDF5_HOME)
AntaresHOME=$(PWD)/../../..

FFLAGS= -cpp -DNDIM=$(NDIM) -DVERBOSE -DISOTHERMAL -DFORCE -D$(Rie_type) -O3

LIB = -L$(MPIHOME)/lib -lmpi -lz -L$(FFTWHOME)/lib -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw -L$(HDF5HOME)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -fPIC
INCL=-I$(MPIHOME)/include -I$(FFTWHOME)/include -I$(HDF5HOME)/include
VPATH = $(PATCH):$(AntaresHOME)/hydro:$(AntaresHOME)/gravity
OBJECT=common_params.o main.o init2d.o output1d.o output2d.o bdry2d.o write_collective2d.o dt_loc2d.o dt_all.o sendbdry2d.o sloplimiter.o riemann2d.o bndflux2d.o fluxRieIso2d.o fluxHLLC2d.o force2d.o sgkernel2d.o selfgravity2d.o write_gravity2d.o read_collective2d.o read1d.o read2d.o

antares.exe: $(OBJECT)
	$(F90) -o $(EXEC)$(NDIM)d $(OBJECT) $(LIB)

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@ $(INCL)

clean: 
	rm *.o $(EXEC)*d *.mod *.h5
