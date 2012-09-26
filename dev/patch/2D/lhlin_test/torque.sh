#!/bin/bash
###### Job name ######
#PBS -N spiral_potential
###### Output files ######
#PBS -e parallel_demo3.err
#PBS -o parallel_demo3.log
###### Queue name #######
#PBS -q medium
###### Number of nodes and cores ######
#PBS -l nodes=4:ppn=8:dl
###### Sends mail to yourself when the job begins and ends ######
#PBS -M ccfeng@tiara.sinica.edu.tw
#PBS -m be
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque lam HDF ifc fftw icc

###### Run parallel jobs ######
$LAM_HOME/bin/lamboot $PBS_NODEFILE
$LAM_HOME/bin/mpiexec C ./antares2d > log
$LAM_HOME/bin/lamhalt

