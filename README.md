# demux10xHTO

A tiny project specific 10x demultiplexer as I did not get any data regarding demultiplexed info. Annoying.


# Test case

```
./target/debug/demux10x -r testData/n10000_HTO_S1_L001_R1_001.fastq.gz -f testData/n10000_HTO_S1_L001_R2_001.fastq.gz -b testData/HTOs.csv -o testData/outpath

# or 

./target/release/demux10x -r testData/n10000_HTO_S1_L001_R1_001.fastq.gz -f testData/n10000_HTO_S1_L001_R2_001.fastq.gz -b testData/HTOs.csv -o testData/outpath

```

# Build

```
source ~/.cargo/env
cargo build --release
```

# Install

```
cp target/release/demux10x /usr/bin/
```

# Usage

```
demux10x -h


demux10x 0.1.0
Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>
Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs

USAGE:
    demux10x --reads <READS> --file <FILE> --bc <BC> --outpath <OUTPATH>

OPTIONS:
    -b, --bc <BC>              the barcodes table name<tab>bc
    -f, --file <FILE>          the input R2 genes file
    -h, --help                 Print help information
    -o, --outpath <OUTPATH>    the outpath
    -r, --reads <READS>        the input R1 reads file
    -V, --version              Print version information


```

The program will create a tab separted table in the outpath containing all sample tags, the cell names and the number of detected fastq reads for the combinations.
