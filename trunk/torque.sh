#!/bin/bash
###### Job name ######
#PBS -N antares2d
###### Output files ######
#PBS -o antares2d.out
#PBS -e antares2d.err
###### Number of nodes and cores ######
#PBS -l nodes=2:ppn=8
###### Queue name ######
#PBS -q small
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque icc/13.0 ifc/13.0  openmpi/1.6.3_ic13.0  fftw/2.1.5_ic13.0_openmpi_1.6.3  HDF/5-1.8.10_ic13.0_openmpi_1.6.3  torque/2.5.3

###### Run your jobs with parameters ######
if [ -n "$PBS_NODEFILE" ]; then
  if [ -f $PBS_NODEFILE ]; then
    NPROCS=`wc -l < $PBS_NODEFILE`
  fi
fi

rm -rf parallel*

$OPENMPI_HOME/bin/mpirun -v -machinefile $PBS_NODEFILE -np $NPROCS ./antares2d > log
