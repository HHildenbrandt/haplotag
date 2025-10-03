#!/bin/bash

set -e # exit on first error

Red='\033[0;31m'
Green='\033[0;32m'
NOCOLOR='\033[0m'

export PATH=$PATH:~/haplotag/bin

root=~/haplotag/Pilot-1
out=~/haplotag/test/out

# run new implementation
fastq_h4 ~/haplotag/test/test_H4.json -f --replace '{"/range": "0-1000"}'

# run old implementation
mkdir -p ${out}
H4_demult_fastq_with_clipping_8bp-plateBC ${root}/subset/ ${out}/ ${root}/ 1000 > /dev/null

# concat result files
fastq_cat ${out}/R1.fastq.gz ${out}/R2.fastq.gz > ${out}/new.txt
fastq_cat ${out}/R1_001.fastq.gz ${out}/R2_001.fastq.gz > ${out}/old.txt

# compare
if [[ -n $(diff ${out}/new.txt ${out}/old.txt) ]]; then
    echo -e "${Red}test_pilot-1 failed!${NOCOLOR}"
else 
    echo -e "${Green}test_pilot-1 passed!${NOCOLOR}"
fi

#clean up
rm -rf ${out}
