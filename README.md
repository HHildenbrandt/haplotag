# haplotag

All you need to know about zlib: [Mark Adler's zlib repository](https://github.com/madler/zlib/tree/develop)

## Build

```bash
git clone git@github.com:HHildenbrandt/haplotag.git

cd haplotag
git submodule update --init --recursive --remote
cd vcpkg
./bootstrap-vcpkg.sh -disableMetrics
cd ..

# we need a c++23 compiler. On Habrok:
module load GCC/13

mkdir build && cd build
cmake ..
cmake --build . --config Release
cd ..

# binaries can be found in ~/haplotag/bin:
bin
├── fastq_cat
├── fastq_h4
├── fastq_paste
├── H4_demult_fastq_with_clipping_7bp-plateBC
├── H4_demult_fastq_with_clipping_8bp-plateBC
└── H4_demult_fastq_with_clipping_noPlateBC
```

For the remainder of this README, we assume that we are im `~/haplotag` and ~/haplotag/bin 
is included in $PATH:

`export PATH=$PATH:~/haplotag/bin`

## fastq_cat

A little tool that behaves similar to `cat`:

```bash
fastq_cat --help
Usage: fastq_cat [OPTIONS] [FILE] ...
Concatenate ranges from fastq[.gz] files.

With no FILE, or when FILE is -, read standard input.

  -f: force overwrite of output file.
  -m <mssk>: only output unmasked lines (max. 64Bit)
    Ex: -m 0010, outputs 2nd line of every 4-line block.
  -o <FILE>: compressed output file
    If not given, writes uncomressed to standard output.
  -r <line range>: only output lines in given range.
    Ex: -r 0-10; -r 10:3
```

Main purpose is to extract subsets and to support UNIX-style piping:

```bash
# show reads from the first 10 sequences
fastq_cat Pilot-1/H4/subset/R1_001.fastq.gz -m 0010 -r 0-40
```

```bash
# show 1st entry in name field
fastq_cat Pilot-1/H4/subset/R1_001.fastq.gz -m 0001 -r 0-40 | awk '{print $1}'
```

```bash
# create subset
fastq_cat Pilot-1/H4/complete/R1_001.fastq.gz -r 0-4000 -o Pilot-1/H4/subset/R1_001.fastq.gz
```

## fastq_paste

A little tool that behaves similar to `paste`:

```bash
fastq_paste --help
Usage: fastq_paste [OPTIONS] [FILE] ...
paste line rangess from fastq[.gz] files.

  -f: force overwrite of output file.
  -m <mssk>: only output unmasked lines (max. 64Bit)
    Ex: -m 0010, outputs 2nd line of every 4-line block.
  -o <FILE>: compressed output file
    If not given, writes to standard output.
  -r <line range>: only output lines in given range.
    Ex: -r 0-10; -r 10:3
  -d: delimiter string
```

## fastq_h4

`H4_demult_fastq_[...]` replacement.

```bash
fastq_h4 --help
Usage: fastq_h4 JSON_FILE [OPTIONS]...
  -h, --help: show this message.
  -f, --force: force overwrite of output directory.
  -v, --verbose: verbose output.
  --replace '{"json_pointer": value}'.
    Ex: --replace '{"/range": "0-1000"}' --replace '{"/barcode/plate/file": "Plate_BC_7.txt"}'
  --dry: dry-run.
```

You can find an example `JSON_FILE` in `~\haplotag\src`.<br>
Note that the comments are *not* part of the json.

```json
{
    "range": "0-1000000", // sequence range, everythin if empty
    "pool_threads": 32,   // number of cores used in thread-pool, -1 for all available cores
    "barcodes": {
        "root": "~/haplotag/Pilot-1",
        "A": {
            "file": "BC_A_H4.txt",
            "unclear_tag": "A00"
        },
        "B": {
            "file": "BC_B.txt",
            "unclear_tag": "B00"
        },
        "C": {
            "file": "BC_C_H4.txt",
            "unclear_tag": "C00"
        },
        "D": {
            "file": "BC_D.txt",
            "unclear_tag": "D00"
        },
        "plate": {
            "file": "Plate_BC_8.txt", // could be empty (no plate)
            "unclear_tag": "P000"
        },
        "stagger": {
            "file": "BC_ME.txt",
            "unclear_tag": "S00",
            "sort_by_tag": true     // this is required for now
        }
    },
    "reads": {
        "root": "~/haplotag/Pilot-1/reads", // root data directory 
        "R1": "R1_001.fastq.gz",
        "R2": "R2_001.fastq.gz",
        "R3": "R3_001.fastq.gz",
        "R4": "R4_001.fastq.gz",
        "I1": "I1_001.fastq.gz"   // ignored if "/barcodes//plate/file" is empty
    },
    "output": {
        "root": "~/haplotag/Pilot-1/reads/out",
        "R1": "R1_001.fastq.gz",  // could be empty (constructable from /reads/R1 and /output/R2)
        "R2": "R2_001.fastq.gz"   // could be empty (no clipping)
    }
}
```

## Pilot-1

```
tree Pilot-1/
Pilot-1/
├── BC_A_H4.txt
├── BC_B.txt
├── BC_C_H4.txt
├── BC_D.txt
├── BC_ME.txt
├── Plate_BC_7.txt
├── Plate_BC_8.txt
├── reads
│   ├── I1_001.fastq.gz -> /scratch/hb-1000G_str/pilot/raw_fastq/Pilot-1_I1_001.fastq.gz
│   ├── R1_001.fastq.gz -> /scratch/hb-1000G_str/pilot/raw_fastq/Pilot-1_R1_001.fastq.gz
│   ├── R2_001.fastq.gz -> /scratch/hb-1000G_str/pilot/raw_fastq/LL-GoNL-Pilot-1_I1_001.fastq.gz
│   ├── R3_001.fastq.gz -> /scratch/hb-1000G_str/pilot/raw_fastq/Pilot-1_R3_001.fastq.gz
│   └── R4_001.fastq.gz -> /scratch/hb-1000G_str/pilot/raw_fastq/Pilot-1_R4_001.fastq.gz
└── subset
    ├── I1_001.fastq.gz
    ├── R1_001.fastq.gz
    ├── R2_001.fastq.gz
    ├── R3_001.fastq.gz
    └── R4_001.fastq.gz
```

Note thet the files under `Pilot-1/reads/' are symbolic links to files that
are *not* included* into the repository. 

### Test

```
./test/test_pilot-1.sh
```

### Mini bench (./test/bench.sh)

Reads from 20GiB USB nvme drive
Writes to local PCIe 4.0 nvme drive

```
H4_demult_fastq_with_clipping_8bp_plateBC scaling
  sequences       time
  1.000.000       1m1.468s
  10.000.000      10m14.738s

fastq_ha scaling, (32 pool threads)
  sequences       time
  1.000.000       0m1.549s
  10.000.000      0m13.176s
  100.000.000     2m12.807s
  1.000.000.000   22m20.853s

fastq_ha scaling with pool threads
  pool_threads    time
  4               1m7.126s
  8               0m33.827s
  16              0m20.724s
  24              0m15.354s
  32              0m13.387s
```

