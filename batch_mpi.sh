#!/bin/bash

# Request resources:
#SBATCH --time=12:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=16G      # 1 GB RAM
#SBATCH -p shared

#SBATCH -n 8                    # number of MPI ranks
#SBATCH -c 16                   # number of threads per rank (one thread per CPU core)
#SBATCH --ntasks-per-socket=1   # number of MPI ranks per CPU socket
#SBATCH -N 8                    # number of compute nodes. 

module load gcc
module load intelmpi
module load aocl

echo "Running code"
rm output/*

#sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit

#cp ~/quicklisp/local-projects/cl-mpm-worker/mpi-worker ./

#echo $OMP_NUM_THREADS
mpirun ./mpi-worker --dynamic-space-size 16000

