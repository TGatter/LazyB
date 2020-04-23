# Lazy Bastard


## Introduction

This is a preliminary version of the tool.


## Requirements

Following tools must be installed:

- jellyfish 2
- bbduk
- Abyss 2 (using 'abyss-pe')
- minimap2
- modified version of racon


## Usage

call the script lazyB_pipeline.sh as follows:

/path/lazyB\_pipeline.sh \[k-mer-size-filter\] \[k-mer-size-assembly\] \[name\] \[illumina-inputfile-1\] \[illumina-inputfile-2\]] \[nanopore-inputfile\] \[output-folder\]


example:

/path/lazyB\_pipeline.sh 50 90 Scerevisiae illumina\_paired\_end\_1.fa illumina\_paired\_end\_2.fa nanopore\_data.fastq output-folder


Note:
Parts of the pipeline run on 8 cores. To change the number of available cores, set the value manually in lazyB\_pipeline.sh.


\[k-mer-size-filter\] specifies the k-mer size for k-mer counting in raw illumina data. Reads with highly abundant k-mers are removed from the data. Starting at k=50 is recommended.

\[k-mer-size-assembly\] specifies the k-mer size during illumina assembly (here using Abyss). Starting at k=90 is recommended.


# Contact