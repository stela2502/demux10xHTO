[![Rust](https://github.com/stela2502/demux10xHTO/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/demux10xHTO/actions/workflows/rust.yml)
# demux10xHTO

This program can be used to quantify some specific transcripts in a 10x data.
It was developed to demultiplex sample reads and antibody reads from 10x data.

Input is not a complex index, but a simple two columns table defining the search name and search sequence.
The tool will 'only' look for these sequences, not the reverse ones.

It has been tested with HTO and A

# Test case

```
./target/debug/demux10x -r testData/n10000_HTO_S1_L001_R1_001.fastq.gz -f testData/n10000_HTO_S1_L001_R2_001.fastq.gz -b testData/HTOs.csv -o testData/outpath

# or 

./target/release/demux10x -r testData/n10000_HTO_S1_L001_R1_001.fastq.gz -f testData/n10000_HTO_S1_L001_R2_001.fastq.gz -b testData/HTOs.csv -o testData/outpath

```

I have not started to implement a test based on this, but you can check the results e.g. with:
```
## first 'gene'
cut -f 2 testData/outpath/Cell2Sample.n10000_HTO_S1_L001_R1_001.fastq.gz.tsv | sort |uniq  -c
## last 'gene'
cut -f 7 testData/outpath/Cell2Sample.n10000_HTO_S1_L001_R1_001.fastq.gz.tsv | sort |uniq  -c

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
Stefan L. <stefan.lang@med.lu.se>
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
