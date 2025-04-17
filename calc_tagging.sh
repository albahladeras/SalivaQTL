#!/bin/bash

#run ldak to obtain weights
ldak --cut-weights weights \
     --bfile data/genotypes/ref_clean \
     --window-cm 0.1 \
     --section-cm 0.1 \
     --no-thin DONE

ldak --calc-weights-all weights \
     --bfile data/genotypes/ref_clean \
     --power -0.25

#rename SNP-list files
mv snps_mqtls.txt  QTL1
mv snps_eqtls.txt  QTL2

# Compute tagging files using LDAK
ldak --calc-tagging results/saliva_tagging \
     --bfile data/genotypes/ref_clean \
     --weights weights/weights.all \
     --power -0.25 \
     --annotation-number 2 \
     --annotation-prefix QTL
