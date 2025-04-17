# Saliva as potential diagnostic medium: DNA methylation biomarkers for disorders beyond the oral cavity
 
In this GitHub repository, you will find the codes developed by the Immunogenetics Research Lab (IRLab) from the University of the Basque Country (UPV/EHU), for the research presented in the article titled: Saliva as diagnostic medium: DNA methylation biomarkers for disorders beyond the oral cavity


## Quality control of genotype data 
### 1. Variant QC
-MAF 0.01  
-HWE 1e-06  
-Call rate 0.05
```
plink
  --bfile saliva_samples \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-06 \
  --make-bed \
  --out saliva_samples_maf001_hwe005_hwe106
```

### 2. Sample QC
-More than 4SD from the mean of heterozigosity  
-Call rate 10%

1. Heterozigosity
1.1 Create the .imiss and .het files with plink
```
plink
  --bfile saliva_samples_maf001_hwe005_hwe106 \
  --missing \
  --het \
  --out saliva_samples_maf001_hwe005_hwe106
```
1.2 Run the [calculate_heterozigosity.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/calculate_heterozigosity.R) script in terminal
```
Rscript calculate_heterozigosity.R saliva_samples_maf001_hwe005_hwe106.imiss \
saliva_samples_maf001_hwe005_hwe106.het
```
1.3 Remove samples with plink
```
plink
  --bfile saliva_samples_maf001_hwe005_hwe106 \
  --remove filter-het.txt \
  --make-bed \
  --out saliva_samples_maf001_hwe005_hwe106_het
```
2. Call rate
```
plink
  --bfile saliva_samples_maf001_hwe005_hwe106_het \
  --mind 0.1 \
  --make-bed \
  --out saliva_samples_maf001_hwe005_hwe106_het_miss010
```

### 3. Prepare files for imputation
1. Create the .frq file with plink
```
plink
  --bfile saliva_samples_maf001_hwe005_hwe106_het_miss010 \
  --freq \
  --out saliva_samples_maf001_hwe005_hwe106_het_miss010
```
2. Execute [Will Rayner preimputation script](https://www.chg.ox.ac.uk/~wrayner/tools/)
```
perl HRC-1000G-check-bim.pl \
  -b saliva_samples_maf001_hwe005_hwe106_het_miss010.bim \
  -f saliva_samples_maf001_hwe005_hwe106_het_miss010.frq \
  -r 100GP_Phase3_combined.legend \
  -g \
  -p EUR \
  -o $dir_preimp \
  -v
```
```
sh Run-plink.sh
```
3. Rename and compress files 
```
for i in {1..23};
do
  vcf-sort saliva_samples_maf001_hwe005_hwe106_het_miss010-updated-chr${i}.vcf | bgzip -c > $i.vcf.gz ;
done
```
4. Impute data using 1000G as reference panel

### 4. Quality control of imputed data
-RsQ > 0.9  
-MAF > 0.01  
-HWE 0.05

1. Concat the imputed files
```
for i in {1..22};
do
echo chr${i}.dose.vcf.gz >>  tmp-concat.txt
done
```
```
echo  chrX.dose.vcf.gz >>  tmp-concat.txt
```
```
bcftools concat -f tmp-concat.txt \
--threads 16 \
-o all_chr_imputed_raw.vcf.gz \
-Oz
```
2. Filter SNPs by rsq09
```
bcftools filter -i 'R2>0.9' \
-o all_chr_imputed_rsq09.vcf \
-Ov \
all_chr_imputed_raw.vcf.gz
```
3. Filter by MAF 0.01 and HWE 0.05
```
plink
  --vcf all_chr_imputed_rsq09.vcf \
  --maf 0.01 \
  --hwe 0.05 \
  --keep-allele-order \
  --make-bed \
  --out all_chr_imputed_rsq09_maf001_hwe005
```
4. Remove sex chromosomes
```
plink \
  --bfile all_chr_imputed_rsq09_maf001_hwe005 \
  --chr 1-22 \
  --keep-allele-order \
  --make-bed \
  --out imputed_rsq09_maf001_hwe005_chr1_22
```

## Methylation data 
### Quality control  
[qc_methylation_data.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/qc_methylation_data.R).

#### 1. CpGs QC:
-Detection P-value < 0.01  
-Bead number 3 in at least 5% of samples per probe  
-All non-CpG probes are removed  
-SNP-related probes are removed: [Zhou et al 2016](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw967) and [McCartney et al 2016](https://www.sciencedirect.com/science/article/pii/S221359601630071X?via%3Dihub)  
-Multi-hit probes are removed: [Nordlund et al 2013](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r105) and [McCartney et al 2016](https://www.sciencedirect.com/science/article/pii/S221359601630071X?via%3Dihub)  
-Probes located in chromosome X and Y are removed  

#### 2. Sample QC: 
-Samples with a fraction of failed probes higher tan 0.05 are removed  

#### 3. Normalization methods: 
  -Functional and Noob Normalization  
  -BMIQ Normalization  
  
#### 4. Batch correction:
  -Combat  
  
#### 5. Outliers treatment:  
  -Winsorization: pct = 0.005  

## mQTL analysis  

### 1. Create [BED File](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
[create_methylation_bedfile.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/create_methylation_bedfile.R)

### 2. Make sure the same samples are present in the methylation and genotype files
Subset from BED and plink files  

### 3. Make sure the samples are in the same order in both files

### 4. Inverse normal transformation of methylation data  
[rnt_methylation.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/rnt_methylation.R)  

### 5. Covariates 
  #### 5 genetic principal components.  
  [pc_air_pc_relate.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/pc_air_pc_relate.R)  
  #### Age  
  #### Sex  
  #### 20 Methylation principal components 
  1. Multiple linear regression with the methylation values as the outcome, and the known confounders as predictors.  
  [pca_residuals_methylation.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/pca_residuals_methylation.R)  
  2. PCA with the residuals.  
  [pca_residuals_methylation.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/pca_residuals_methylation.R)  
  3. Remove methylation PCs associated with the genotype.  
  [gwas_mpcs.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/gwas_mpcs.R) 
  #### Genetic Relatedness Matrix.  
  [pc_air_pc_relate.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/pc_air_pc_relate.R)

### 6. Run Linear mixed model using [Apex](https://github.com/corbinq/apex?tab=readme-ov-file)  
Genotype and bed file must be indexed
```
apex cis \
--vcf imputed_rsq09_maf001_hwe005_chr1_22_common_methyl.vcf.gz \
--bed bed_sorted.bed.gz \
--cov covariates.txt \
--kin kinship.txt \
--window 1000000 \
--threads 16 \
--long \
--prefix MQTL_LMM
```

## Quality control of expression data

### 1. Concat FASTQ files from different lanes:
[concat_lanes.sh](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/concat_lanes.sh)

### 2. Preprocessing with Fastp
[fastp.sh](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/fastp.sh)

### 3. Splice-aware mapping to genome
We used STAR (Spliced Transcripts Alignment to a Reference) software and UCSC hg19 reference genome  
[alignment.sh](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/alignment.sh)  

### 4. Transcript assembly and quantification
We used StringTie software
[stringtie.sh](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/stringtie.sh)

### 5. Sample and gene filtering and TMM transformation
-Samples with more than 5M uniquely mapped reads were selected.  
-Genes with ≥ 6 reads in ≥ 20% of samples were selected.  
-Read counts were normalized between samples using TMM  
[gene_sample_qc_tmm.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/gene_sample_qc_tmm.R)

## eQTL analysis 

### 1. Create BED File (columns: chr, start, end, gene_id, sample_1...sample_n)

### 2. Make sure the same samples are present in the expression and genotype files
Subset from BED and plink files  

### 3. Make sure the samples are in the same order in both files

### 4. Inverse normal transformation of expression data  
[rnt_transformation_exp.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/rnt_transformation_exp.R)

### 5. Covariates 
  #### 5 genetic principal components.  
  Execute the [pc_air_pc_relate.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/pc_air_pc_relate.R) script using the genotype data of samples for which expression data is available.  
  #### Age  
  #### Sex  
  #### Date of RNA-seq
  #### 45 PEER factors 
  [peer_calculation.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/peer_calculation.R)
  #### Genetic Relatedness Matrix.  
  Execute the [pc_air_pc_relate.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/pc_air_pc_relate.R) script using the genotype data of samples for which expression data is available.  
  
### 6. Run Linear mixed model using [Apex](https://github.com/corbinq/apex?tab=readme-ov-file)  
Genotype and bed file must be indexed
```
apex cis \
--vcf imputed_rsq09_maf001_hwe005_chr1_22_common_exp.vcf.gz \
--bed bed_sorted.bed.gz \
--cov covariates.txt \
--kin kinship.txt \
--window 1000000 \
--threads 16 \
--long \
--prefix EQTL_LMM
```
## eQTM analysis  
To perform the [eQTM analyses](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/matrixeqtl.R), we adapted a script provided by Corina Lesseur that utilizes [MatrixeQTL](https://github.com/andreyshabalin/MatrixEQTL) R package. MatrixEQTL was designed to perform linear regressions between genotype and expression data. In our case, we substituted the genotype data with methylation data, and we applied the following model:

```
gene expression ~ methylation + sex + age + time of RNA-seq + 5gPCs + 20mPCs + 45 expression PEER factors + genetic relatedness matrix
transcript expression ~ methylation + sex + age + time of RNA-seq + 5gPCs + 20mPCs + 45 expression PEER factors + genetic relatedness matrix
```
## Multi-SNP Summary-based Mendelian Randomization (SMR) analyses 

### 1. Converting cis-m/eQTL data to BESD format.
We constructed a txt file in [fastqtl format](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/fastqtl_format_for_SMR_and_epi.R). According to the [SMR webpage](https://yanglab.westlake.edu.cn/software/smr/#MakeaBESDfile) this file has no header and contains five columns: gene/CpG site, SNP, distance in bp between the SNP and the gene/CpG, p-value and beta (i.e. SNP effect on gene expression/CpG methylation). We then executed the following command to obtain the BESD format:

```
smr --eqtl-summary fastqtlnomi.txt --fastqtl-nominal-format --make-besd --out mybesd
```
### 2. Updating BESD files
As it is described by the developers, The SNP and probe information in the SMR output files (.esi and .epi) converted from FastQTL output are not complete and need to be updated using the options in Update a BESD file. We first created manually the [esi](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/create_esi.R) and [epi](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/fastqtl_format_for_SMR_and_epi.R) files and we then executed the following smr command:

```
smr --beqtl-summary  mybesd --update-esi mymqtl.esi
smr --beqtl-summary  mybesd --update-esi mymqtl.epi
```
### 3. Run Multi-SNP SMR with the harmonized GWAS summary statistics

```
 smr --bfile imputed_rsq09_maf001_hwe005_chr1_22_common_methyl \
        --gwas-summary gwas.txt \
        --beqtl-summary mybesd \
        --smr-multi \
        --thread-num 8 \
        --out gwas_smr_results
```

## Partitioning Heritability Using LDAK with mQTL and eQTL Annotations 

LDAK (Linkage Disequilibrium Adjusted Kinships) is a widely used toolkit for estimating SNP heritability. Unlike S-LDSC, it accommodates realistic LD patterns and allows fine-tuned annotation models. This makes it ideal for exploring the contribution of biologically defined SNP sets (e.g., mQTLs or eQTLs) to disease risk.
LDAK’s --sum-hers function implements the SumsHer method, which uses summary statistics from GWAS and tagging information from a reference panel to estimate how much heritability is attributable to each SNP category.

Key Concepts  
    • Heritability (h²) is the proportion of phenotypic variance explained by genetic factors.  
    • Partitioned heritability assigns portions of h² to specific SNP sets (e.g., mQTLs, eQTLs).  
    • Annotation: A file labeling SNPs as belonging to a category (e.g., mQTL = 1 or 0).  
    • Tagging: LD-adjusted measure of how well SNPs in each annotation predict genome-wide variation.  

  #### Prepare the reference.  
  [prepare_reference.sh](https://github.com/albahladeras/SalivaQTL/blob/main/prepare_reference.sh)

  ```
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

```

  #### Calculate the tagging file
  [calc_tagging.sh](https://github.com/albahladeras/SalivaQTL/blob/main/calc_tagging.sh)

 ```
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

  ```
  #### Run sum hers 
  [run_sum_hers.sh](https://github.com/albahladeras/SalivaQTL/blob/main/run_sum_hers.sh)

  ```
#!/bin/bash

tagging_file="results/saliva_tagging.tagging"
results_folder="results/outputs"

while IFS=$'\t' read -r file_name ascertainment prevalence; do
    [[ "$file_name" == "file_name" ]] && continue

    if [[ -f "GWAS/${file_name}" ]]; then
        ldak --sum-hers "${results_folder}/output_${file_name}" \
             --tagfile "$tagging_file" \
             --summary "GWAS/${file_name}" \
             --ascertainment "$ascertainment" \
             --prevalence "$prevalence" \
             --check-sums NO
    else
        echo "File GWAS/${file_name} not found. Skipping..."
    fi
done < "GWAS_data.txt"

```
  #### Extract the results
  [extract_results.R](https://github.com/albahladeras/SalivaQTL/blob/main/extract_results.R)
