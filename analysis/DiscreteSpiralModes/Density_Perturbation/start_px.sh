#!/bin/bash
###### Job name ######
#PBS -N bc
###### Output files ######
#PBS -o bc.out
#PBS -e bc.err
###### Number of nodes and cores ######
#PBS -l nodes=1:ppn=16:px
###### Queue name ######
#PBS -q small
###### Specific the shell types ######
#PBS -S /bin/bash
#PBS -M ccfeng@tiara.sinica.edu.tw

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh
module purge
module load torque ifc/13.0 icc/13.0 pgplot
OMP_NUM_THREADS=$PBS_NUM_PPN
export OMP_NUM_THREADS

rm -rf bc.log bc.err
rm -rf pgplot.png
###### Run your jobs with parameters ######
#./Density.exe 
./SearchAll.exe > searchall.log
#./FindOne.exe
