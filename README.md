# haplotag

All you need to know about zlib: [Mark Adler's zlib repository](https://github.com/madler/zlib/tree/develop)

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
├── fastq_H4
├── fastq_paste
├── H4_demult_fastq_with_clipping_7bp-plateBC
├── H4_demult_fastq_with_clipping_8bp-plateBC
└── H4_demult_fastq_with_clipping_noPlateBC
```

The binaries are statically linked. You can move/copy them to other places.<br>
For the remainder of this README, we assume ~/haplotag/bin in included in $PATH:

`export PATH=$PATH:~/haplotag/bin`

## fastq_cat

A little tool that behaves similar to `cat`:

```bash
fastq_cat --help
Usage: fastq_cp [OPTIONS] [FILE] ...
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

## fastq_H4

`H4_demult_fastq_[...]` replacement.

```bash
fastq_H4 --help
Usage: fastq_H4 JSON_FILE [OPTIONS]...
  -h, --help: show this message.
  -f, --force: force overwrite of output directory.
  -v, --verbose: verbose output.
  --replace '{"json_pointer": value}'.
    Ex: --replace '{"/range": "0-1000"}' --replace '{"/barcode/plate/file": "Plate_BC_7.txt"}'
  --dry: dry-run.
```

You can find an example `JSON_FILE` in `~\haplotag\src.<br>
Note that the comments are *not* part of the json.

```json
{
    "root": "/home/hanno/haplotag/Pilot-1/H4/", // root data directory 
    "range": "",          // if not empty, sequence range e.g. "0-10000", 
    "pool_threads": -1,   // -1: number of cores used in thread-pool
    "barcodes": {
        "root": "",       // directory for the barcode files, same as '/root' if empty
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
            "file": "Plate_BC_8.txt",
            "unclear_tag": "P000"
        },
        "stagger": {
            "file": "BC_ME.txt",
            "unclear_tag": "S00",
            "sort_by_tag": true     // this is required for now
        }
    },
    "reads": {
        "subdir": "complete",       // path below /root
        "R1": "R1_001.fastq.gz",
        "R2": "R2_001.fastq.gz",
        "R3": "R3_001.fastq.gz",
        "R4": "R4_001.fastq.gz",
        "I1": "I1_001.fastq.gz"
    },
    "output": {
        "subdir": "complete/out",   // path below /root
        "clipping": true,
        "R1": "R1_001.fastq.gz",
        "R2": "R2_001.fastq.gz"
    }
}
```

