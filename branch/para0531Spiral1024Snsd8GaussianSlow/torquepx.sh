#!/bin/bash
###### Job name ######
#PBS -N antares2d
###### Output files ######
#PBS -o antares2d.out
#PBS -e antares2d.err
###### Number of nodes and cores ######
#PBS -l nodes=2:ppn=16
###### Queue name ######
#PBS -q small
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque fftw/2.1.5_ic13.0_lam_7.1.4 HDF/5-1.8.10_ic13.0_lam_7.1.4 lam/7.1.4_ic13.0 torque

rm antares2d.*

###### Run your jobs with parameters ######
# Start a LAM multicomputer
$LAM_HOME/bin/lamboot $PBS_NODEFILE

$LAM_HOME/bin/mpiexec C ./antares2d  > log

# Shutdown the LAM/MPI run-time environment.
$LAM_HOME/bin/lamhalt
