#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -t 1-50:1
#$ -V


#-pe OpenMP 8 -l s_vmem=6G

Female_vcf=/work23/home/nsasa/data/CHM13v2.0_hhv6/haplotypecaller_jointgenotyping02_sampleploidy_modsexLCAC0489/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.vcf.gz
Male_vcf=/work23/home/nsasa/data/CHM13v2.0_hhv6/haplotypecaller_jointgenotyping02_sampleploidy_modsexLCAC0489/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.vcf.gz

mkdir -p vcf

if [[ ! -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
if [ ! -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz" ]; then
  source /work23/home/nsasa/conda-pack/gatk4.1.9.0/bin/activate
  bcftools annotate\
   --threads 7\
   -a /work23/home/nsasa/data/Reference/CHM13/CHM13v2.0_GATK_Resource/chm13v2.0_dbSNPv155.vcf.gz\
   -c ID\
   --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT'\
   -O z\
   -o vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz\
   ${Female_vcf}
   #-i 'FILTER="PASS"'\
fi
fi
if [ -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz" ]; then
if [ ! -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz.tbi" ]; then
  tabix -p vcf vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz
fi
fi

if [[ ! -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
if [ ! -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz" ]; then
  source /work23/home/nsasa/conda-pack/gatk4.1.9.0/bin/activate
  bcftools annotate\
   --threads 7\
   -a /work23/home/nsasa/data/Reference/CHM13/CHM13v2.0_GATK_Resource/chm13v2.0_dbSNPv155.vcf.gz\
   -c ID\
   --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT'\
   -O z\
   -o vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz\
   ${Male_vcf}
   #-i 'FILTER="PASS"'\
fi
fi
if [ -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz" ]; then
if [ ! -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz.tbi" ]; then
  tabix -p vcf vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz
fi
fi


##############################################
### LD-based genotype refinement by Beagle ###
##############################################

### LD-based genotype refinement for low-confidence genotypes and missing sites in WGS data using BEAGLE v5.4 with default settings
# reference panelを指定した場合はreferenceにあって、VCFにないmarkerはimputeされる。今回のようにSporadic missing genotypes are imputed during phasing.
### https://www.nature.com/articles/s41531-023-00456-6#Sec10
### before excluding multialellic
# https://kgnmg.com/wp/2022/09/04/pop_preparation/

if [[ ! -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
if [[ ! -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.vcf.gz" ]]; then
    source /work23/home/nsasa/conda-pack/beagle5.4/bin/activate
    unset JAVA_TOOL_OPTIONS
    
    java -Xmx6g -Xms6g -jar /work23/home/nsasa/conda-pack/beagle5.4/share/beagle-5.4_22Jul22.46e-0/beagle.jar\
     gt=vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz\
     out=vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed\
     nthreads=6
fi
fi

if [[ ! -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
if [[ ! -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.vcf.gz" ]]; then
    source /work23/home/nsasa/conda-pack/beagle5.4/bin/activate
    unset JAVA_TOOL_OPTIONS
    
    java -Xmx6g -Xms6g -jar /work23/home/nsasa/conda-pack/beagle5.4/share/beagle-5.4_22Jul22.46e-0/beagle.jar\
     gt=vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz\
     out=vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed\
     nthreads=6
fi
fi


### filter multiple alleles and spanning deletion
### https://www.nature.com/articles/s41598-019-53111-7

if [[ ! -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
    source /work23/home/nsasa/conda-pack/gatk4.1.9.0/bin/activate
    bcftools view\
     --threads 7\
     -m2 -M2\
     -e 'ALT="*"'\
     -O z -o vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz\
     vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.vcf.gz
fi

if [[ ! -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
    source /work23/home/nsasa/conda-pack/gatk4.1.9.0/bin/activate
    bcftools view\
     --threads 7\
     -m2 -M2\
     -e 'ALT="*"'\
     -O z -o vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz\
     vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.vcf.gz
fi

if [[ -f "vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
    rm vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.imputed.vcf.gz
    rm vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz
    rm vcf/JointGenotyping.XnonPAR_Females_${SGE_TASK_ID}.annotated.vcf.gz.tbi
fi

if [[ -f "vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.filtered.vcf.gz" ]]; then
    rm vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.imputed.vcf.gz
    rm vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz
    rm vcf/JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.annotated.vcf.gz.tbi
fi
