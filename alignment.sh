#!/usr/bin/env bash
set -e
set -u
set -o pipefail

THR=8 # Number of threads to use (check how many are in your computer!)
FASTQDIR="/data/exome"

# Create index (ONLY NEED TO RUN ONCE)

STAR \
    --runThreadN $THR \
    --runMode genomeGenerate \
    --genomeDir 'STAR_index_ensGene' \
    --genomeFastaFiles 'hg19.p13.plusMT.no_alt_analysis_set.fa' \
    --sjdbGTFfile 'ucsc_hg19/hg19.ensGene.gtf' \
    --sjdbOverhang 99
    
# List of fastq.gz files found in a specific directory
FASTQS=$( find $FASTQDIR -maxdepth 2 -type f -name "*fastq.gz" -exec ls {} + )

# Sample substrings
SAMPLES=$( ls $FASTQS | sed 's!.*/!!' | awk '{m=split($0,a,"_L"); for(i=1;i<m;i++)printf a[i]}; {print ""}' | sort | uniq )

dirA="02_fastp"
dirB="03_alignment"

echo $dirA
echo $dirB

for SAMP in $( echo $SAMPLES );
do
    
    echo "Processing of: $SAMP"
    
    # Measure execution time
    START=$(date +%s.%N)
    
    if [[ -e "${dirB}/${SAMP}_sorted.bam.bai" ]]; then
        continue
    fi

    STAR \
	    --runThreadN $THR \
	    --genomeDir 'STAR_index_ensGene' \
	    --readFilesIn "${dirA}/${SAMP}_R1_fastp.fastq.gz" "${dirA}/${SAMP}_R2_fastp.fastq.gz" \
	    --readFilesCommand zcat \
	    --outFileNamePrefix ${dirB}/${SAMP}

	samtools view \
        -bS "${dirB}/${SAMP}Aligned.out.sam" \
        --threads $THR \
        > "${dirB}/${SAMP}.bam"

    samtools sort \
        "${dirB}/${SAMP}.bam" \
        --threads $THR \
        -o "${dirB}/${SAMP}_sorted.bam"
    
    samtools index -@ 4 "${dirB}/${SAMP}_sorted.bam"
    
    # Measure execution time
    END=$(date +%s.%N)
    DIFF=$(echo "$END - $START" | bc)
    echo "$SAMP execution time:"
    echo $DIFF
    
done
