#!/bin/bash


SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

K_MER_JELLY=$1
K_MER_ABYSS=$2
NAME=$3
ILLUMINA_RAW_1=$4
ILLUMINA_RAW_2=$5
NANO=$6
OUT=$7

CORES=8                     #number of cores can be set manually here
MINLENGTH=500
ABYSS_MODE=unitigs

mkdir -p $OUT               #create output folder if it doesn't already exist
TMP=$(mktemp -d -p $OUT)    #create a temporary folder - deleted in the end
BASE=$(basename $NANO .fastq)

echo ">>>> K-mer Filtering of Illumina Reads"

jellyfish count -m $K_MER_JELLY -s 100M -t $CORES -C $ILLUMINA_RAW_1 $ILLUMINA_RAW_2 -o ${TMP}/"jelly_count_k${K_MER_JELLY}.jf"  
jellyfish histo ${TMP}/"jelly_count_k${K_MER_JELLY}.jf" > ${TMP}/"jelly_histo_k${K_MER_JELLY}.histo"
TOTAL_NON_UNIQUE_KMERS=$(awk '{if($1 != "1") s += $2} END{print s}' ${TMP}/"jelly_histo_k${K_MER_JELLY}.histo")
ABUNDANCE_THRESHOLD=$($SCRIPTPATH/setAbundanceThresholdFromHisto.py ${TMP}/"jelly_histo_k${K_MER_JELLY}.histo" $TOTAL_NON_UNIQUE_KMERS)
echo "abundance threshold for k-mer filtering: " $ABUNDANCE_THRESHOLD > $OUT/report.txt
jellyfish dump -L $ABUNDANCE_THRESHOLD ${TMP}/"jelly_count_k${K_MER_JELLY}.jf" >  ${TMP}/"filtered_kmers_${K_MER_JELLY}_${ABUNDANCE_THRESHOLD}.fa"
bbduk.sh in1=$ILLUMINA_RAW_1 in2=$ILLUMINA_RAW_2 out1=$TMP/illu_filtered.1.fastq out2=$TMP/illu_filtered.2.fastq ref=${TMP}/"filtered_kmers_${K_MER_JELLY}_${ABUNDANCE_THRESHOLD}.fa" k=$K_MER_JELLY hdist=0

echo ">>>> Illumina Assembly"

mkdir -p $OUT/ABYSS   #create folder "ABYSS" for ABYSS results
abyss-pe -C $OUT/ABYSS np=$CORES name=$NAME k=$K_MER_ABYSS in="../../${TMP}/illu_filtered.1.fastq ../../${TMP}/illu_filtered.2.fastq" ${ABYSS_MODE} 2>&1 | tee $OUT/ABYSS/abyss.log
awk -v min="$MINLENGTH" 'BEGIN {RS = ">" ; ORS = ""} $2 >= min {print ">"$0}' $OUT/ABYSS/"${NAME}-${ABYSS_MODE}.fa"  > $OUT/ABYSS/"${NAME}-${ABYSS_MODE}.l${MINLENGTH}.fa"


echo ">>>> Unitig Filter"

minimap2 -k15 -DP --dual=yes --no-long-join -w5 -m100 -g10000 -r2000 --max-chain-skip 25 --split-prefix foo $NANO $OUT/ABYSS/"${NAME}-${ABYSS_MODE}.l${MINLENGTH}.fa" > $OUT/01_unitigs.to_$BASE.paf
$SCRIPTPATH/unitig_filter.py $OUT/01_unitigs.to_$BASE.paf $OUT/ABYSS/"${NAME}-${ABYSS_MODE}.l${MINLENGTH}.fa" $OUT/report.txt $TMP/unitigs_corrected.fa


echo ">>>> Scrubbing"

minimap2 -k15 -DP --dual=yes --no-long-join -w5 -m100 -g10000 -r2000 --max-chain-skip 25 --split-prefix foo $NANO $NANO > $OUT/01_a_$BASE.to_self.paf
minimap2 -k15 -DP --dual=yes --no-long-join -w5 -m100 -g10000 -r2000 --max-chain-skip 25 --split-prefix foo $NANO $TMP/unitigs_corrected.fa > $OUT/01_b_contigs_corrected.to_$BASE.paf
$SCRIPTPATH/scrubber.py $OUT/01_a_$BASE.to_self.paf $OUT/01_b_contigs_corrected.to_$BASE.paf $NANO $OUT/02_$BASE.scrubbed.fa


echo ">>>> Anchor Mapping"

minimap2  -k15 -DP --dual=yes --no-long-join -w5 -m100 -g10000 -r2000 --max-chain-skip 25 --split-prefix foo $OUT/02_$BASE.scrubbed.fa $TMP/unitigs_corrected.fa > $OUT/02_contigs_corrected.to_$BASE.scrubbed.paf



echo ">>>> Lazy Bastard"

$SCRIPTPATH/prokrastinator.py $OUT/02_contigs_corrected.to_$BASE.scrubbed.paf $TMP/unitigs_corrected.fa $OUT/02_$BASE.scrubbed.fa $TMP
cp $TMP/temp_1.target.fa $OUT/03.assembly.unpolished.fa

echo ">>>> Racon"
$SCRIPTPATH/racon_mod -u -t $CORES $TMP/temp_1.query.fa $TMP/temp_1.align.paf $TMP/temp_1.target.fa > $OUT/03.assembly.fa

rm -rf $TMP























