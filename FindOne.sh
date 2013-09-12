#!/bin/bash
###### Job name ######
#PBS -N Find_One 
###### Output files ######
#PBS -e FindOne.mpi.err
#PBS -o FindOne.mpi.log
###### Queue name ######
#PBS -q serial
###### Number of nodes and cores(ppn), time(walltime) and memory size(mem) ######
#PBS -l nodes=1:ppn=12:dl
###### Sends mail to yourself when the job begins and ends ######
#PBS -m be
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh

. ~/.module

rm FindOne.mpi.err
rm FineOne.mpi.log
rm FindOne.log

# Set the number of OpenMP threads to share the work.
# This is actually the same as the default value as each node has 8 cores
# and we have asked PBS to use all of them, so we don't need to set it.  It is
# here only for demonstration purposes.
export OMP_NUM_THREADS=12

./FindOne.exe $arg1 $arg2 > FindOne.log
