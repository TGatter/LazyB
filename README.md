# Lazy Bastard


## Introduction

This is a preliminary version of the tool under rapid development.

Economic Genome Assembly from Low Coverage Illumina and Nanopore Data
Thomas Gatter, Sarah von Löhneysen, Polina Drozdova, Tom Hartmann, Peter F. Stadler
bioRxiv 2020.02.07.939454; doi: https://doi.org/10.1101/2020.02.07.939454 

Please contact us with any issues and suggestions.

## Requirements

Following tools must be installed:

- jellyfish 2
- bbduk
- Abyss 2 (using 'abyss-pe')
- minimap2
- modified version of racon (racon_mod in this project, compiled for generic linux)


## Usage

call the script lazyB_pipeline.sh as follows:

/path/lazyB\_pipeline.sh \[k-mer-size-filter\] \[k-mer-size-assembly\] \[name\] \[illumina-inputfile-1\] \[illumina-inputfile-2\]] \[nanopore-inputfile\] \[output-folder\]


example:

/path/lazyB\_pipeline.sh 50 90 Scerevisiae illumina\_paired\_end\_1.fa illumina\_paired\_end\_2.fa nanopore\_data.fastq output-folder


Note:
Parts of the pipeline run on 8 cores. To change the number of available cores, set the value manually in lazyB\_pipeline.sh.

\[k-mer-size-filter\] specifies the k-mer size for k-mer counting in raw illumina data. Reads with highly abundant k-mers are removed from the data. Starting at k=50 is recommended.

\[k-mer-size-assembly\] specifies the k-mer size during illumina assembly (here using Abyss). Starting at k=90 is recommended.

Comment out the following line in lazyB\_pipeline.sh to disable the experimental polishing step.

```
$SCRIPTPATH/racon_prokrast -u -t $CORES $TMP/temp_1.query.fa $TMP/temp_1.align.paf $TMP/temp_1.target.fa > $OUT/03.assembly.fa
```

