#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile.err
#PBS -l walltime=24:00:00
#PBS -l select=2:ncpus=20
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .
mpif90 inputfile.f
mpiexec.hydra -np 40 -genv OMP_NUM_THREADS=1 -genv I_MPI_PIN=1 -genv I_MPI_FABRICS=shm:ofi -hostfile $PBS_NODEFILE ./a.out >output.txt
rm a.out
mv ../job$tpdir $PBS_O_WORKDIR/.