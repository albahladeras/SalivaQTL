#!/usr/bin/env bash
set -e
set -u
set -o pipefail

THR=16 # Number of threads to use (check how many are in your computer!)
FASTQDIR="/data/exome/"

# List of fastq.gz files found in a specific directory
FASTQS=$( find $FASTQDIR -maxdepth 2 -type f -name "*fastq.gz" -exec ls {} + )

# Sample substrings
SAMPLES=$( ls $FASTQS | sed 's!.*/!!' | awk '{m=split($0,a,"_L"); for(i=1;i<m;i++)printf a[i]}; {print ""}' | sort | uniq )

cat_samples() {
    
    FASTQDIR="/data/exome/"

    # List of fastq.gz files found in a specific directory
    FASTQS=$( find $FASTQDIR -maxdepth 2 -type f -name "*fastq.gz" -exec ls {} + )
    
    dir01="/data/preprocess_exo/01_cat"
    
    R1s=$( ls $FASTQS | grep ${1} | grep "R1")
    R2s=$( ls $FASTQS | grep ${1} | grep "R2")

    # Concatenate sample lanes
    cat $R1s > ${dir01}/${1}_R1.fastq.gz
    cat $R2s > ${dir01}/${1}_R2.fastq.gz

}

# Export the function so it can be used by parallel
export -f cat_samples

# Use parallel to execute the tasks in parallel
parallel -j 16 cat_samples ::: $SAMPLES
