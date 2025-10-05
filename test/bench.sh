#!/bin/bash

export PATH=$PATH:~/haplotag/bin

echo && echo H4_demultfastq_with_clipping_8bp-plateBC scaling
echo && echo 1.000.000 sequences
time H4_demult_fastq_with_clipping_8bp-plateBC Pilot-1/reads/ Pilot-1/reads/out/ Pilot-1/ 1000000
echo && echo 10.000.000 sequences
time H4_demult_fastq_with_clipping_8bp-plateBC Pilot-1/reads/ Pilot-1/reads/out/ Pilot-1/ 10000000

echo fastq_ha scaling, max. pool threads
echo && echo 1.000.000 sequences
time fastq_h4 src/H4.json -f --replace '{"/range": "0-1000000"}'
echo && echo 10.000.000 sequences
time fastq_h4 src/H4.json -f --replace '{"/range": "0-10000000"}'
echo && echo 100.000.000 sequences
time fastq_h4 src/H4.json -f --replace '{"/range": "0-100000000"}'
echo && echo 1.000.000.000 sequences
time fastq_h4 src/H4.json -f --replace '{"/range": "0-1000000000"}'

echo && echo fastq_ha scaling with pool threads, 10.000.000 sequences
echo && echo pool_threads 4
time fastq_h4 src/H4.json -f --replace '{"/range": "0-10000000"}' --replace '{"/pool_threads": 4}'
echo && echo pool_threads 8
time fastq_h4 src/H4.json -f --replace '{"/range": "0-10000000"}' --replace '{"/pool_threads": 8}'
echo && echo pool_threads 16
time fastq_h4 src/H4.json -f --replace '{"/range": "0-10000000"}' --replace '{"/pool_threads": 16}'
echo && echo pool_threads 24
time fastq_h4 src/H4.json -f --replace '{"/range": "0-10000000"}' --replace '{"/pool_threads": 24}'
echo && echo pool_threads 32
time fastq_h4 src/H4.json -f --replace '{"/range": "0-10000000"}' --replace '{"/pool_threads": 32}'
