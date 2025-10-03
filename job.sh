#!/bin/bash

#SBATCH --mem=8GB
#SBATCH --cpus-per-task=32
#SBATCH --partition=short
#SBATCH --time=00:30:00


srun ~/haplotag/bin/fastq_h4 ~/haplotag/src/H4.json -f --replace '{"/pool_threads": 32}' --replace '{"/range": "0-1000000"}' --replace '{"/output/root": "/scratch/hb-1000G_str/pilot/raw_fastq/out"}'
