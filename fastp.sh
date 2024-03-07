#!/bin/bash
set -e
set -u
set -o pipefail

THR=8 # Number of threads to use (check how many are in your computer!)
FASTQDIR="/data/exome"

# List of fastq.gz files found in a specific directory
FASTQS=$( find $FASTQDIR -maxdepth 2 -type f -name "*fastq.gz" -exec ls {} + )

# Sample substrings
SAMPLES=$( ls $FASTQS | sed 's!.*/!!' | awk '{m=split($0,a,"_L"); for(i=1;i<m;i++)printf a[i]}; {print ""}' | sort | uniq )

dir01="/data/preprocess_exo/01_cat"
dir02="/data/preprocess_exo/02_fastp_2"

echo $SAMPLES

for SAMP in $( echo $SAMPLES );
do
    
    echo "Processing of: $SAMP"
    
    # Measure execution time
    START=$(date +%s.%N)

    /data/genotools/fastp\
        --in1 ${dir01}/${SAMP}_R1.fastq.gz \
        --in2 ${dir01}/${SAMP}_R2.fastq.gz \
        --out1 ${dir02}/${SAMP}_R1_fastp.fastq.gz \
        --out2 ${dir02}/${SAMP}_R2_fastp.fastq.gz \
        --json ${dir02}/${SAMP}_fastp.json \
        --html ${dir02}/${SAMP}_fastp.html \
        --thread $THR \
        > "${dir02}/${SAMP}_fastp.stdout.txt" \
        2> "${dir02}/${SAMP}_fastp.stderr.txt"
    
    # Measure execution time
    END=$(date +%s.%N)
    DIFF=$(echo "$END - $START" | bc)
    echo "$SAMP execution time:"
    echo $DIFF
    
done
