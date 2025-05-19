# haplotag

All you need to know about zlib: [Mark Adler's zlib repository](https://github.com/madler/zlib/tree/develop)

```
git clone git@github.com:HHildenbrandt/haplotag.git

cd haplotag
git submodule update --init --recursive --remote
cd vcpkg
./bootstrap-vcpkg.sh -disableMetrics
cd ..

mkdir build && cd build
cmake ..
cmake --build . --target Release

# binaries can be found in ./bin
```

### generate bigger fastq files

blow-up every *.gz file found in ../data:

```
bin$ ./fastq_gen 100000
bin$ time ./fastq_gen 100000
found '../data/R2_001.fastq.gz'
found '../data/test_H4_R2_001.fastq.gz'
found '../data/I1_001.fastq.gz'
found '../data/test_R3_001.fastq.gz'
found '../data/test_R1_001.fastq.gz'
found '../data/test_H4_R4_001.fastq.gz'
found '../data/test_R2_001.fastq.gz'
found '../data/R1_001.fastq.gz'
found '../data/test_H4_I1_001.fastq.gz'
found '../data/I2_001.fastq.gz'
found '../data/test_I2_001.fastq.gz'
found '../data/test_H4_R1_001.fastq.gz'
found '../data/test.fastq.gz'
found '../data/test_I1_001.fastq.gz'
generating 14 '_gen_*.gz' files...
2.10772 MB inflated
177.016 GB deflated
```

## demult_fastq vs. fastq_fuzzy

### License

`demult_fastq` links against `gzStream`, available under LGPL license.<br>
`fastq_fuzzy` links against `zlib-ng`, available under MIT license.

### bench

```
bin& ./fastq_gen 100000

bin$ time ./demult_fastq 
loaded barcodes: 96 A, 96 B, 96 C, 96 D 
real    36m89.157s
user    35m46.333s
sys     0m2.273s

bin$ time ./fastq_fuzzy 
starting
*  50000000 sets processed -
*  45241 MB decompressed
*  41723 MB compressed

real    0m43.564s
user    17m2.988s
sys     0m11.086s
```

### verify

The following diffs should return nothing:

```
bin$ gzip -d -k -f ../data/_R1_001.fastq.gz
bin$ gzip -d -k -f ../data/_R2_001.fastq.gz
bin$ gzip -d -k -f ../data/_fuzzy_R1_001.fastq.gz
bin$ gzip -d -k -f ../data/_fuzzy_R2_001.fastq.gz
bin$ diff ../data/_R1_001.fastq ../data/_fuzzy_R1_001.fastq
bin$ diff ../data/_R2_001.fastq ../data/_fuzzy_R2_001.fastq
bin$ diff ../data/_clearBC.log ../data/_fuzzy_clearBC.log 
bin$ diff ../data/_unclearBC.log ../data/_fuzzy_unclearBC.log 
```

### clean up

```
bin$ rm ../data/_*
```

## Compute nodes on Hábrók

Local storage is on `$TMPDIR`

119 standard nodes with the following components:
* 128 cores @ 2.45 GHz (two AMD 7763 CPUs)
* 512 GB memory
* 3.5 TB internal SSD disk space

24 nodes for multi-node jobs with the following components:
* 128 cores @ 2.45 GHz (two AMD 7763 CPUs)
* 512 GB memory
* 3.5 TB internal SSD disk space
* 100 Gbps Omni-Path link

4 big memory nodes with the following components:
* 80 cores @ 2.3 GHz (two Intel Xeon Platinum 8380 CPUs)
* 4096 GB memory
* 14 TB internal SSD disk space

2 Interactive GPU nodes (Delivered by Fujitsu in an earlier purchase) with the following components:
* 24 cores @ 2.4 GHz (two Intel Xeon Gold 6240R CPUs)
* 768 GB memory
* 1 Nvidia L40s GPU accelerator card with 48GB RAM
* 6 GPU nodes with the following components:
* 64 cores @ 2.6 GHz (two Intel Xeon Platinum 8358 CPUs)
* 512 GB memory

4 Nvidia A100 GPU accelerator cards with 40 GB RAM
* 12 TB internal SSD NVMe disk space
* 100 Gbps Omni-Path link

19 GPU nodes (Delivered by Fujitsu in an earlier purchase) with the following components:
* 36 cores @ 2.7 GHz (two Intel Xeon Gold 6150 CPUs)
* 768 GB memory (621 GB used for temporary disk space)
* 2 Nvidia V100 GPU accelerator cards each with 32 GB RAM
* 621 GB RAM disk

15 nodes with the following components:
* 64 cores @ 2.2 GHz (two AMD EPYC 7601 CPUs)
* 512 GB memory
* 16 TB internal disk space
Only accessible by GELIFES users, see GELIFES Partition

1 node with the following components:
* 64 cores @ 2.1 GHz (two Intel Xeon Gold 6448Y CPUs)
* 1 TB memory
* 440 GB internal disk space
* 4 Nvidia H100 GPU accelerator cards with 80 GB RAM

Only accessible for education purposes in the scope of the Digital Lab project (employee login required)

## Network Hábrók

A 100 Gbps low-latency non-blocking Omni-Path network for 24 compute and 6 GPU nodes<br>
High bandwidth (100 Gigabit per second)<br>
Low latency (few microseconds delay before a client starts receiving the message)<br>
Useful for parallel processing over multiple computers<br>
Two 25 Gbps Ethernet networks<br>
Used for accessing the storage areas and for job communication<br>
Can also be useful to access remote data more quickly<br>


## Storage Hábrók

The cluster has 2.5 PB (2562 TB) of formatted storage available. This scratch storage is set up using the Lustre parallel file system.
50 GB of home directory storage per user

## getCode variants

### min_edit_distance

returns<br>
* 'correct' on first exact match (handled by caller), assumes unique barcodes.
* 'corrected' if the best ed is unique
* 'unclear' otherwise

### demult_fastq

`getCode(I1,codeA,codeC,RX1,QX1,read_type1,13, "A", "C", bc_A, bc_C);`<br>
code_total_length = 13 =? |bc_A| + 1 + |bc_B|, unique barcode lengths |bc_A| == |bc_C| == 6<br>
1st min_ed: (I1.RX[|bc_A| + 1, ], bc_A)    // unbounded <br>
2nd min_ed: (I1.RX[0, |bc_C|], bc_C) <br>
one base skipped.

returns<br>
* read_type = "unclear", codeA = codeB = "00" if |I1.RX| < code_total_length
* read_type = "unclear" if any match is 'unclear', codeA, codeB from `med`
* read_type = "corrected" if 1st or 2nd match is "corrected", codeA, codeB from `med`
* read_type = "correct" if 1st and 2nd match are exact

failed code_total_length - test doesn't fit very well into "unclear".
read_type and codeA, codeB are not always consistent.
The bc_X maps are used for a fancy way to test 'code == barcode'.

Same for getCode(I2, ...)

