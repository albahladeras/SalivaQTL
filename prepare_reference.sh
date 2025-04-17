#!/bin/bash

# Create PLINK fileset for children only (from imputed genotypes)

awk '$2 !~ /M/ { print $1, $2 }' data/genotypes/all_rsq09.fam > data/genotypes/keep_ids.txt

plink --bfile data/genotypes/all_rsq09 \
      --keep data/genotypes/keep_ids.txt \
      --make-bed \
      --out data/genotypes/filtered

plink --bfile data/genotypes/filtered \
      --maf 0.01 --geno 0.05 \
      --make-bed \
      --out data/genotypes/reference_panel

# Add recombination map
plink --bfile data/genotypes/reference_panel \
      --cm-map data/genetic_map/genetic_map_chr@_combined_b37.txt \
      --make-bed \
      --out data/genotypes/reference_panel_with_cm

# Remove SNPs with missing or invalid cM values
awk '$3 > 0.001' data/genotypes/reference_panel_with_cm.bim > data/genotypes/good_snps.bim
cut -f2 data/genotypes/good_snps.bim > data/genotypes/snps_to_keep.txt

plink --bfile data/genotypes/reference_panel_with_cm \
      --extract data/genotypes/snps_to_keep.txt \
      --make-bed \
      --out data/genotypes/ref_clean

#rename snps as chr:pos (yhey are rsIDs now)
awk '{print $2, $1":"$4}' ref_clean.bim > rsid_to_chrpos.txt
awk '{$2 = $1":"$4; print}' ref_clean.bim > temp && mv temp ref_clean.bim

#clean-up duplicates
#find duplicates
awk '{print $2}' ref_clean.bim | sort | uniq -d > duplicates.txt
plink --bfile ref_clean --exclude duplicates.txt --make-bed --out ref_clean_no_duplicates
