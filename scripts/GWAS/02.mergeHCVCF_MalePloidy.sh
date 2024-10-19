#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V

SampleID=$1
sex=$2
old_haplotypecaller_dir=$3
#-pe OpenMP 2 -l s_vmem=6G

mkdir -p ${SampleID}
cd ${SampleID}


mkdir -p log

if [ ! -f "${SampleID}.g.vcf.gz" ]; then
  
  ### extract AutosomalPAR
  if [ ! -f "${SampleID}.AutoPAR.g.vcf.gz" ]; then
    old_merged_vcf=${old_haplotypecaller_dir}/${SampleID}/${SampleID}.g.vcf.gz
  
    source /work23/home/nsasa/conda-pack/gatk4.2.6.1/bin/activate
    REF=/work23/home/nsasa/data/Reference/CHM13/v2.0/chm13v2.0_maskedY.fa
    GATK_DIR=/work23/home/nsasa/data/Script/hg38_for_01/tools/gatk-4.2.6.1
    interval=/work23/home/nsasa/data/Script/chm13_hhv6/haplotypecaller/VQSR02/chm13v2.0_Autosomal_PAR.gatk_intervals.list
    /usr/bin/time -f "Memory:%M KB time:%E" -o log/extract_Autosomal_PAR.txt \
    ${GATK_DIR}/gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
     SelectVariants\
      -R $REF\
      -V ${old_merged_vcf}\
      -L $interval\
      -O ${SampleID}.AutoPAR.g.vcf.gz\
      > log/extract_Autosomal_PAR.log 2>&1
  fi
  
  
  ### MergeVCFs
  if [ ! -f "${SampleID}.g.vcf.gz" ]; then
    source /work23/home/nsasa/conda-pack/gatk4.2.6.1/bin/activate
    GATK_DIR=/work23/home/nsasa/data/Script/hg38_for_01/tools/gatk-4.2.6.1
    mkdir -p MergeVcfsTmp
  
    /usr/bin/time -f "Memory:%M KB time:%E" -o log/MergeVcfs.txt \
    ${GATK_DIR}/gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
     MergeVcfs\
     -I ${SampleID}.AutoPAR.g.vcf.gz\
     -I scattered-HaplotypeCaller/chrXnonPAR_HC.g.vcf.gz\
     -I scattered-HaplotypeCaller/chrYnonPAR_HC.g.vcf.gz\
     -O ${SampleID}.g.vcf.gz\
     --TMP_DIR MergeVcfsTmp\
     --MAX_RECORDS_IN_RAM 1000000\
     > log/MergeVcfs.log 2>&1
  
    rm -rf MergeVcfsTmp
    echo "MergeVcfs Done"
  else
    echo "MergeVcfs Already Done"
  fi
  
  if [ -f "${SampleID}.g.vcf.gz" ]; then
    if [ -f "${SampleID}.AutoPAR.g.vcf.gz" ]; then
      rm ${SampleID}.AutoPAR.g.vcf.gz
    fi
  fi
fi



#if [ ! -f "${ID}N.HaplotypeCaller.decomposed.vcf.gz" ]; then
#  source /work23/home/nsasa/conda-pack/pcgr/bin/activate
#  /usr/bin/time -f "Memory:%M KB time:%E" -o log/vt.txt \
#  vt decompose -s ${ID}N.HaplotypeCaller.vcf.gz -o ${ID}N.HaplotypeCaller.decomposed.vcf.gz\
#   > log/vt.log 2>&1
#  tabix -p vcf ${ID}N.HaplotypeCaller.decomposed.vcf.gz
#fi