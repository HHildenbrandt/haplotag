#!/bin/bash

Red='\033[0;31m'
Green='\033[0;32m'
NOCOLOR='\033[0m'

export PATH=$PATH:~/haplotag/bin

root=~/haplotag/Pilot-1/H4

# run new implementation
fastq_H4 ~/haplotag/test/test_H4.json -f --replace '{"/range": "0-1000"}'

# run old implementation
mkdir -p ${root}/subset/orig_test_out
H4_demult_fastq_with_clipping_8bp-plateBC ${root}/subset/ ${root}/subset/orig_test_out/ ${root}/ 1000 > /dev/null

# concat result files
fastq_cat ${root}/subset/orig_test_out/R1_001.fastq.gz ${root}/subset/orig_test_out/R2_001.fastq.gz > ${root}/test_pilot_orig.txt
fastq_cat ${root}/subset/test_out/R1_001.fastq.gz ${root}/subset/test_out/R2_001.fastq.gz > ${root}/test_pilot.txt

# compare
if [[ -n $(diff ${root}/test_pilot_orig.txt ${root}/test_pilot.txt) ]]; then
    echo -e "${Red}test_pilot-1.sh failed!${NOCOLOR}"
else 
    echo -e "${Green}test_pilot-1.sh passed!${NOCOLOR}"
    exit 1
fi

# clean up
rm -rf ${root}/subset/orig_test_out
rm -rf ${root}/subset/test_out
rm ${root}/test_pilot*
