#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -t 1-50:1
#$ -V

# -pe OpenMP 2 -l s_vmem=24G -l m_mem_free=24G


source /work23/home/nsasa/conda-pack/gatk4.2.6.1/bin/activate
REF=/work23/home/nsasa/data/Reference/CHM13/v2.0/chm13v2.0_maskedY.fa
GATK_DIR=/work23/home/nsasa/data/Script/hg38_for_01/tools/gatk-4.2.6.1


male_haplotypecaller_dir=/work23/home/nsasa/data/CHM13v2.0_hhv6/haplotypecaller_sampleploidy

gvcf_files_males=""
gvcf_files_females=""
for SampleID in $(tail -n +2 /work23/home/nsasa/data/CHM13v2.0_hhv6/haplotypecaller_sampleploidy/all_sex.modsex_LCAC0489.tsv | awk -F "\t" '{print $1}')
do
  sex=$(tail -n +2 /work23/home/nsasa/data/CHM13v2.0_hhv6/haplotypecaller_sampleploidy/all_sex.modsex_LCAC0489.tsv | awk -v SampleID=${SampleID} '$1==SampleID{print $3}')
  haplotypecaller_dir=$(tail -n +2 /work23/home/nsasa/data/CHM13v2.0_hhv6/haplotypecaller_sampleploidy/all_sex.modsex_LCAC0489.tsv | awk -v SampleID=${SampleID} '$1==SampleID{print $5}')

  if [[ "${sex}" == "M" ]]; then
    gvcf_file=${male_haplotypecaller_dir}/${SampleID}/${SampleID}.g.vcf.gz
    gvcf_files_males=${gvcf_files_males}"-V ${gvcf_file} "
  else
    gvcf_file=${haplotypecaller_dir}/${SampleID}/${SampleID}.g.vcf.gz
    gvcf_files_females=${gvcf_files_females}"-V ${gvcf_file} "
  fi
done


interval_N=$(seq -f '%04g' 0 49 | head -n ${SGE_TASK_ID} | tail -n 1)
interval_list=/work23/home/nsasa/data/Script/chm13_hhv6/haplotypecaller/chm13v2.0_X_nonPAR.interval_list_split50/${interval_N}-scattered.interval_list




mkdir -p log

### XnonPAR Males
#XnonPAR_interval_list=/work23/home/nsasa/data/Script/chm13_hhv6/haplotypecaller/chm13v2.0_X_nonPAR.gatk_intervals.list

### GenomicsDBImport
mkdir -p tmp
mkdir -p tmp/tmp_genomicsdbimport_XnonPAR_Males_${SGE_TASK_ID}
mkdir -p gvcfs_db
/usr/bin/time -f "Memory:%M KB time:%E" -o log/GenomicsDBImport_XnonPAR_Males_${SGE_TASK_ID}.txt \
${GATK_DIR}/gatk --java-options "-Xmx24g -Xms24g -XX:+UseSerialGC"\
 GenomicsDBImport\
 --reference ${REF} ${gvcf_files_males}\
 --intervals ${interval_list}\
 --genomicsdb-workspace-path gvcfs_db/gvcfs_db_contig_XnonPAR_Males_${SGE_TASK_ID}\
 --overwrite-existing-genomicsdb-workspace true\
 --tmp-dir tmp/tmp_genomicsdbimport_XnonPAR_Males_${SGE_TASK_ID}\
 > log/GenomicsDBImport_XnonPAR_Males_${SGE_TASK_ID}.log 2>&1

rm -rf tmp/tmp_genomicsdbimport_XnonPAR_Males_${SGE_TASK_ID}

### GenotypeGVCFs
mkdir -p tmp/tmp_genotypegvcfs_XnonPAR_Males_${SGE_TASK_ID}
/usr/bin/time -f "Memory:%M KB time:%E" -o log/GenotypeGVCFs_XnonPAR_Males_${SGE_TASK_ID}.txt \
${GATK_DIR}/gatk --java-options "-Xmx24g -Xms24g -XX:+UseSerialGC"\
 GenotypeGVCFs\
 --reference ${REF}\
 --variant gendb://gvcfs_db/gvcfs_db_contig_XnonPAR_Males_${SGE_TASK_ID}\
 -L ${interval_list}\
 --output JointGenotyping.XnonPAR_Males_${SGE_TASK_ID}.vcf.gz\
 --tmp-dir tmp/tmp_genotypegvcfs_XnonPAR_Males_${SGE_TASK_ID}\
 > log/GenotypeGVCFs_XnonPAR_Males_${SGE_TASK_ID}.log 2>&1

rm -rf tmp/tmp_genotypegvcfs_XnonPAR_Males_${SGE_TASK_ID}

