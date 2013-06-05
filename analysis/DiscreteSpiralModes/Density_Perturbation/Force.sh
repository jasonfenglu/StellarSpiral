#!/bin/bash
###### Job name ######
#PBS -N force
###### Output files ######
#PBS -o force/force.out
#PBS -e force/force.err
###### Number of nodes and cores ######
#PBS -l nodes=1:ppn=8
###### Queue name ######
#PBS -q small
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque fftw/2.1.5_ic13.0_lam_7.1.4 HDF/5-1.8.10_ic13.0_lam_7.1.4 lam/7.1.4_ic13.0
module add pgplot

###### Run your jobs with parameters ######
# Set the number of OpenMP threads to share the work.
export OMP_NUM_THREADS=12

./Force.exe 
