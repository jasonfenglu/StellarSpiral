#!/bin/bash
###### Job name ######
#PBS -N antares
###### Output files ######
#PBS -e parallel.err
#PBS -o parallel.log
###### Queue name #######
#PBS -q debug
###### Number of nodes and cores ######
#PBS -l nodes=4:ppn=8:dl
###### Sends mail to yourself when the job begins and ends ######
#PBS -M ccfeng@asiaa.sinica.edu.tw
#PBS -m be
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque lam fftw/2.1.5_ic11.0_mpich_1.2.7p1 HDF/5-1.8.7_ic11.0_lam_7.1.4 ifc  icc

rm -rf parallel.err parallel.log log

###### Run parallel jobs ######
$LAM_HOME/bin/lamboot $PBS_NODEFILE
$LAM_HOME/bin/mpiexec C ./antares2d > log
$LAM_HOME/bin/lamhalt
