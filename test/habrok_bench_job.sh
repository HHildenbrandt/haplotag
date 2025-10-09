#!/bin/bash

#SBATCH --job-name h4_bench
#SBATCH --output=%j.log
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00

module load GCC/13      # c++ stdlib
export PATH=$PATH:~/haplotag/bin

srun ~/haplotag/test/habrok_bench.sh
