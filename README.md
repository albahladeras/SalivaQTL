## Quality control of genotype data for the SalivaQTL project
In this GitHub repository, you will find the protocol elaborated by the Immunogenetics Research Lab (IRLab) from the University of the Basque Country (UPV/EHU), to perform the quality control of genotype, methylation and expression data used in the SalivaQTL project.

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
## Quality control of methylation data 
[qc_methylation_data.R](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/qc_methylation_data.R).

### 1. CpGs QC:
-Detection P-value < 0.01  
-Bead number 3 in at least 5% of samples per probe  
-All non-CpG probes are removed  
-SNP-related probes are removed: [Zhou et al 2016](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw967) and [McCartney et al 2016](https://www.sciencedirect.com/science/article/pii/S221359601630071X?via%3Dihub)  
-Multi-hit probes are removed: [Nordlund et al 2013](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r105) and [McCartney et al 2016](https://www.sciencedirect.com/science/article/pii/S221359601630071X?via%3Dihub)  
-Probes located in chromosome X and Y are removed  

### 2. Sample QC: 
-Samples with a fraction of failed probes higher tan 0.05 are removed  

### 3. Normalization methods: 
  -Functional and Noob Normalization  
  -BMIQ Normalization  
  
### 4. Batch correction:
  -Combat  
  
### 5. Outliers treatment:  
  -Winsorization: pct = 0.005

## Quality control of expression data

### 1. Concat FASTQ files from different lanes:
[concat_lanes.sh](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/concat_lanes.sh)

### 2. Preprocessing with Fastp
[fastp.sh](https://github.com/albahladeras/SalivaQTL_QC_genotype/blob/main/fastp.sh)

### 3. Splice-aware mapping to genome
We used STAR (Spliced Transcripts Alignment to a Reference) software and UCSC hg19 reference genome

