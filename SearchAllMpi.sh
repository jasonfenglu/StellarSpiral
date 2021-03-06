#!/bin/bash
###### Job name ######
#PBS -N Search_All 
###### Output files ######
#PBS -e searchall.mpi.err
#PBS -o searchall.mpi.log
###### Queue name ######
#PBS -q medium
###### Number of nodes and cores(ppn), time(walltime) and memory size(mem) ######
#PBS -l nodes=6:ppn=12:dl
###### Sends mail to yourself when the job begins and ends ######
#PBS -m be
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh

. ~/.module

rm searchall.mpi.err
rm searchall.mpi.log
rm searchall.log

# Set the number of OpenMP threads to share the work.
# This is actually the same as the default value as each node has 8 cores
# and we have asked PBS to use all of them, so we don't need to set it.  It is
# here only for demonstration purposes.
#export OMP_NUM_THREADS=12

###### Run parallel jobs ######
cat $PBS_NODEFILE | uniq > LAMHOST
$LAM_HOME/bin/lamboot -v LAMHOST
$LAM_HOME/bin/lamclean -v
$LAM_HOME/bin/mpiexec C ./SearchAll.exe $rstart $rend  $istart $iend> searchall.log
$LAM_HOME/bin/lamhalt
rm LAMHOST
