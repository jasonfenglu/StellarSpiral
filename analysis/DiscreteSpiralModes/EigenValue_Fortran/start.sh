#!/bin/bash
###### Job name ######
#PBS -N bc
###### Output files ######
#PBS -o bc.out
#PBS -e bc.err
###### Number of nodes and cores ######
#PBS -l nodes=1:ppn=12:dl
###### Queue name ######
#PBS -q serial
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque intel pgplot
OMP_NUM_THREADS=$PBS_NUM_PPN
export OMP_NUM_THREADS

rm -rf bc.log bc.err
###### Run your jobs with parameters ######
./SearchAll.exe 
#./FindOne.exe 44.68 -1.784

