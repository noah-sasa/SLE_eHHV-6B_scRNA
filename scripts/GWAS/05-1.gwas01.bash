#!/bin/sh


# before VQSR
# beagle filtering
mkdir -p _log
qsub -o _log -e _log -q long@node0[2-9] -q long@node1[0-9] -q long@node2[0-8] -pe OpenMP 8 -l s_vmem=3G -l m_mem_free=3G Beagle5.4_filtering_AutosomalPAR.sh
qsub -o _log -e _log -q long@node0[2-9] -q long@node1[0-9] -q long@node2[0-8] -pe OpenMP 8 -l s_vmem=3G -l m_mem_free=3G Beagle5.4_filtering_XnonPAR.sh
qsub -o _log -e _log -q long@node0[2-9] -q long@node1[0-9] -q long@node2[0-8] -pe OpenMP 8 -l s_vmem=3G -l m_mem_free=3G Beagle5.4_filtering_YnonPAR.sh

#for i in $(seq 1 300); do
#    echo -e "vcf/JointGenotyping.${i}.annotated.imputed.filtered.vcf.gz"
#done > vcf/AutosomalPAR.list
#source /work23/home/nsasa/conda-pack/gatk4.2.6.1/bin/activate
#REF=/work23/home/nsasa/data/Reference/CHM13/v2.0/chm13v2.0_maskedY.fa
#GATK_DIR=/work23/home/nsasa/data/Script/hg38_for_01/tools/gatk-4.2.6.1
#${GATK_DIR}/gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
#   GatherVcfsCloud\
#    --ignore-safety-checks\
#    --gather-type BLOCK\
#    --create-output-variant-index\
#    --input vcf/AutosomalPAR.list\
#    --output vcf/JointGenotyping.AutosomalPAR.annotated.imputed.filtered.vcf.gz
#Beagleなどでいろいろヘッダー情報など消えているのでもうmergeできない。。。

### vcf to bed,bim,fam  ALT が複数ある行は最初の１個だけ読み込まれるようです→ multiallelicは予め除去済
# --split-x　chm13に合わせてpseudo-autosomal regionの分離（maleのhet.を区別するため）
for i in $(seq 1 299)
do
    vcf=vcf/JointGenotyping.${i}.annotated.imputed.filtered.vcf.gz
    plink --threads 8 --memory 10000\
     --make-bed --vcf ${vcf}\
     --pheno original/pheno.modsex.txt\
     --out original/joint_22_220.AutosomalPAR${i}\
     --double-id\
     --allow-no-sex\
     --keep-allele-order\
     --pheno-name eHHV-6B\
     --update-sex original/pheno.modsex.txt 1
done
vcf=vcf/JointGenotyping.300.annotated.imputed.filtered.vcf.gz
plink --threads 8 --memory 10000\
 --make-bed --vcf ${vcf}\
 --pheno original/pheno.modsex.txt\
 --out original/joint_22_220.AutosomalPAR300\
 --double-id\
 --allow-no-sex\
 --split-x 2394410 153925835\
 --keep-allele-order\
 --pheno-name eHHV-6B\
 --update-sex original/pheno.modsex.txt 1



for i in $(seq 1 50); do
    vcf=vcf/JointGenotyping.XnonPAR_Females_${i}.annotated.imputed.filtered.vcf.gz
    plink --threads 8 --memory 10000\
     --make-bed --vcf ${vcf}\
     --pheno original/pheno.modsex.txt\
     --out original/joint_22_220.XnonPAR_Females${i}\
     --double-id\
     --allow-no-sex\
     --keep-allele-order\
     --pheno-name eHHV-6B\
     --update-sex original/pheno.modsex.txt 1

    vcf=vcf/JointGenotyping.XnonPAR_Males_${i}.annotated.imputed.filtered.vcf.gz
    plink --threads 8 --memory 10000\
     --make-bed --vcf ${vcf}\
     --pheno original/pheno.modsex.txt\
     --out original/joint_22_220.XnonPAR_Males${i}\
     --double-id\
     --allow-no-sex\
     --keep-allele-order\
     --pheno-name eHHV-6B\
     --update-sex original/pheno.modsex.txt 1
done
plink --threads 8 --memory 10000\
 --make-bed --vcf vcf/JointGenotyping.YnonPAR_Males.annotated.imputed.filtered.vcf.gz\
 --pheno original/pheno.modsex.txt\
 --out original/joint_22_220.YnonPAR_Males\
 --double-id\
 --allow-no-sex\
 --keep-allele-order\
 --pheno-name eHHV-6B\
 --update-sex original/pheno.modsex.txt 1


for i in $(seq 2 300); do
    echo -e "original/joint_22_220.AutosomalPAR${i}"
done > original/AutosomalPAR.list
plink --threads 8 --memory 100000\
 --bfile original/joint_22_220.AutosomalPAR1\
 --merge-list original/AutosomalPAR.list\
 --make-bed\
 --out original/joint_22_220.AutosomalPAR.merge



# XnonPAR merge
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220.XnonPAR_Females1\
 --bmerge original/joint_22_220.XnonPAR_Females2\
 --make-bed\
 --out original/joint_22_220.XnonPAR_Females.merge2
rm original/joint_22_220.XnonPAR_Females1.*
rm original/joint_22_220.XnonPAR_Females2.*
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220.XnonPAR_Males1\
 --bmerge original/joint_22_220.XnonPAR_Males2\
 --make-bed\
 --out original/joint_22_220.XnonPAR_Males.merge2
rm original/joint_22_220.XnonPAR_Males1.*
rm original/joint_22_220.XnonPAR_Males2.*
for i in $(seq 2 49); do
    n=$(echo ${i} | awk '{print $1+1}')
    plink --threads 8 --memory 10000\
     --bfile original/joint_22_220.XnonPAR_Females.merge${i}\
     --bmerge original/joint_22_220.XnonPAR_Females${n}\
     --make-bed\
     --out original/joint_22_220.XnonPAR_Females.merge${n}
    rm original/joint_22_220.XnonPAR_Females.merge${i}.*
    rm original/joint_22_220.XnonPAR_Females${n}.*
    plink --threads 8 --memory 10000\
     --bfile original/joint_22_220.XnonPAR_Males.merge${i}\
     --bmerge original/joint_22_220.XnonPAR_Males${n}\
     --make-bed\
     --out original/joint_22_220.XnonPAR_Males.merge${n}
    rm original/joint_22_220.XnonPAR_Males.merge${i}.*
    rm original/joint_22_220.XnonPAR_Males${n}.*
done
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220.XnonPAR_Females.merge50\
 --bmerge original/joint_22_220.XnonPAR_Males.merge50\
 --make-bed\
 --out original/joint_22_220.XnonPAR
#Error: 143 variants with 3+ alleles present.
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220.XnonPAR_Females.merge50\
 --exclude original/joint_22_220.XnonPAR-merge.missnp\
 --make-bed\
 --out original/tmp.FemaleXnonPAR
rm original/joint_22_220.XnonPAR_Females.merge50.*
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220.XnonPAR_Males.merge50\
 --exclude original/joint_22_220.XnonPAR-merge.missnp\
 --make-bed\
 --out original/tmp.MaleXnonPAR
rm original/joint_22_220.XnonPAR_Males.merge50.*
plink --threads 8 --memory 10000\
 --bfile original/tmp.FemaleXnonPAR\
 --bmerge original/tmp.MaleXnonPAR\
 --make-bed\
 --out original/joint_22_220.XnonPAR
rm original/tmp.FemaleXnonPAR.*
rm original/tmp.MaleXnonPAR.*




plink --threads 8 --memory 100000\
 --bfile original/joint_22_220.AutosomalPAR.merge\
 --bmerge original/joint_22_220.XnonPAR\
 --make-bed\
 --out original/joint_22_220.AutosomalPARXnonPAR
rm original/joint_22_220.AutosomalPAR.merge.*
rm original/joint_22_220.XnonPAR.*

plink --threads 8 --memory 100000\
 --bfile original/joint_22_220.AutosomalPARXnonPAR\
 --bmerge original/joint_22_220.YnonPAR_Males\
 --make-bed\
 --out original/joint_22_220
rm original/joint_22_220.AutosomalPARXnonPAR.*
rm original/joint_22_220.YnonPAR_Males.*






## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/
## A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis
##############################################################
############### START ANALISIS ###############################
##############################################################

### Step 1 ### 

# Investigate missingness per individual and per SNP and make histograms.
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220\
 --missing\
 --out original/joint_22_220

# Generate plots to visualize the missingness results.
Rscript script/hist_miss.R

# Delete SNPs and individuals with high levels of missingness, explanation of this and all following steps can be found in box 1 and table 1 of the article mentioned in the comments of this script.
# The following two QC commands will not remove any SNPs or individuals. However, it is good practice to start the QC with these non-stringent thresholds.  
# Note, SNP filtering should be performed before individual filtering.

# Delete SNPs with missingness >0.2.
mkdir -p data
plink --threads 8 --memory 10000\
 --bfile original/joint_22_220\
 --geno 0.2\
 --make-bed --out data/joint_22_220_02
#   432360 variants removed due to missing genotype data (--geno).

# Delete individuals with missingness >0.2.
plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_02\
 --mind 0.2\
 --make-bed --out data/joint_22_220_03
#   0 people removed due to missing genotype data (--mind).

# Delete SNPs with missingness >0.02.
plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_03\
 --geno 0.02\
 --make-bed --out data/joint_22_220_04
#   0 variants removed due to missing genotype data (--geno).

# Delete individuals with missingness >0.02.
plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_04\
 --mind 0.02\
 --make-bed --out data/joint_22_220_05
#   0 people removed due to missing genotype data (--mind).



###################################################################
### Step2 ####

# Check for sex discrepancy.    between sex of the individuals recorded in the dataset and their sex based on X chromosome heterozygosity/homozygosity rates.
# Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_05\
 --check-sex\
 --out data/joint_22_220_05


# Generate plots to visualize the sex-check results.
Rscript --no-save script/gender_check.R
# These checks indicate that there is one woman with a sex discrepancy, F value of 0.99. (When using other datasets often a few discrepancies will be found). 


###################################################
### Step 3 ### 

# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' data/joint_22_220_05.bim > data/joint_22_220_05.bim.snp_1_22.txt
plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_05\
 --extract data/joint_22_220_05.bim.snp_1_22.txt\
 --make-bed --out data/joint_22_220_07
#   --extract: 23316505 variants remaining.

# Generate a plot of the MAF distribution.
plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_07\
 --freq\
 --out data/joint_22_220_07
Rscript --no-save script/MAF_check.R

plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_07\
 --freq case-control\
 --out data/joint_22_220_07

plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_07\
 --freq counts\
 --filter-cases\
 --out data/joint_22_220_07_cases

plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_07\
 --freq counts\
 --filter-controls\
 --out data/joint_22_220_07_ctrls

plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_07\
 --freqx\
 --filter-cases\
 --out data/joint_22_220_07_cases

plink --threads 8 --memory 10000\
 --bfile data/joint_22_220_07\
 --freqx\
 --filter-controls\
 --out data/joint_22_220_07_ctrls


