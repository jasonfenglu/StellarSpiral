#!/bin/bash
###### Job name ######
#PBS -N antares
###### Output files ######
#PBS -e parallel.err
#PBS -o parallel.log
###### Queue name #######
#PBS -q small
###### Number of nodes and cores ######
#PBS -l nodes=2:ppn=16:px
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
module add intel fftw/2.1.5_ic11.0_mpich_1.2.7p1 HDF/5-1.8.7_ic11.0_lam_7.1.4  lam cmake pgplot torque
#module add HDF/5-1.8.10_ic13.0_lam_7.1.4  fftw/2.1.5_ic13.0_lam_7.1.4 pgplot torque lam/7.1.4_ic13.0

rm -f parallel.err
rm -f parallel.log

###### Run parallel jobs ######
$LAM_HOME/bin/lamboot $PBS_NODEFILE
$LAM_HOME/bin/mpiexec C ./antares2d  > log
$LAM_HOME/bin/lamhalt

