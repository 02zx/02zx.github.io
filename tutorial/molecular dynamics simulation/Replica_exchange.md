#!/bin/bash
#SBATCH -n 64
#SBATCH -N 2
#SBATCH -J Rep-76-30
#SBATCH -e error.err
#SBATCH -o output.out
module load GROMACS/2020.6-plumed-gpu

mpirun -n 38 gmx_mpi mdrun -s Replica.tpr -deffnm Replica -v -multidir `seq 1 2 76` -replex 100
