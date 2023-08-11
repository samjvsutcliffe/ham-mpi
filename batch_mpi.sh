#!/bin/bash

# Request resources:
#SBATCH -c 1     # 1 entire node
#SBATCH --time=12:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=64G      # 1 GB RAM
#SBATCH --gres=tmp:1G
#SBATCH -p shared
#SBATCH -n 2
#SBATCH -N 2

module load gcc
#module load openmpi
module load intelmpi
module load aocl

echo "Running code"
rm output/*

sbcl --dynamic-space-size 4000  --disable-debugger --load "build_step.lisp" --quit

cp ~/quicklisp/local-projects/cl-mpm-worker/mpi-worker ./

mpirun ./mpi-worker
#mpirun ~/quicklisp/local-projects/cl-mpm-worker/mpi-worker

