#!/usr/bin/env bash
set -e
set -u
set -o pipefail

THR=8 # Number of threads to use (check how many are in your computer!)
FASTQDIR="/data/exome"

# List of fastq.gz files found in a specific directory
FASTQS=$( find $FASTQDIR -maxdepth 2 -type f -name "*fastq.gz" -exec ls {} + )

# Sample substrings
SAMPLES=$( ls $FASTQS | sed 's!.*/!!' | awk '{m=split($0,a,"_L"); for(i=1;i<m;i++)printf a[i]}; {print ""}' | sort | uniq )

dirA="/data/preprocess_exo/03_alignment"
dirB="/data/preprocess_exo/04_quant"


for SAMP in $( echo $SAMPLES );
do
    
    echo "Processing of: $SAMP"
    
        if [[ -e "${dirB}/${SAMP}_gene_abundances.tsv" ]]; then
        continue
    fi
    
    # Measure execution time
    START=$(date +%s.%N)
    
    stringtie \
        "${dirA}/${SAMP}_sorted.bam" \
        -G "/data/resources/refgenome/ucsc_hg19/hg19.ensGene.gtf" \
        -e \
        -B \
        -o "${dirB}/${SAMP}_transcripts.gtf" \
        -A "${dirB}/${SAMP}_gene_abundances.tsv" \
        -p 1
    
    # Measure execution time
    END=$(date +%s.%N)
    DIFF=$(echo "$END - $START" | bc)
    echo "$SAMP execution time:"
    echo $DIFF
    
done
