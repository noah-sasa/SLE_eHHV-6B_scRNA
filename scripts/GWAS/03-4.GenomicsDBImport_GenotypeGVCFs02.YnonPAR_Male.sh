#!/bin/bash

#$ -S /bin/bash
#$ -cwd
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




mkdir -p log




### YnonPAR Males
YnonPAR_interval_list=/work23/home/nsasa/data/Script/chm13_hhv6/haplotypecaller/chm13v2.0_Y_nonPAR.gatk_intervals.list

### GenomicsDBImport
mkdir -p tmp
mkdir -p tmp/tmp_genomicsdbimport_YnonPAR_Males
mkdir -p gvcfs_db
/usr/bin/time -f "Memory:%M KB time:%E" -o log/GenomicsDBImport_YnonPAR_Males.txt \
${GATK_DIR}/gatk --java-options "-Xmx24g -Xms24g -XX:+UseSerialGC"\
 GenomicsDBImport\
 --reference ${REF} ${gvcf_files_males}\
 --intervals ${YnonPAR_interval_list}\
 --genomicsdb-workspace-path gvcfs_db/gvcfs_db_contig_YnonPAR_Males\
 --overwrite-existing-genomicsdb-workspace true\
 --tmp-dir tmp/tmp_genomicsdbimport_YnonPAR_Males\
 > log/GenomicsDBImport_YnonPAR_Males.log 2>&1

rm -rf tmp/tmp_genomicsdbimport_YnonPAR_Males

### GenotypeGVCFs
mkdir -p tmp/tmp_genotypegvcfs_YnonPAR_Males
/usr/bin/time -f "Memory:%M KB time:%E" -o log/GenotypeGVCFs_YnonPAR_Males.txt \
${GATK_DIR}/gatk --java-options "-Xmx24g -Xms24g -XX:+UseSerialGC"\
 GenotypeGVCFs\
 --reference ${REF}\
 --variant gendb://gvcfs_db/gvcfs_db_contig_YnonPAR_Males\
 -L ${YnonPAR_interval_list}\
 --output JointGenotyping.YnonPAR_Males.vcf.gz\
 --tmp-dir tmp/tmp_genotypegvcfs_YnonPAR_Males\
 > log/GenotypeGVCFs_YnonPAR_Males.log 2>&1

rm -rf tmp/tmp_genotypegvcfs_YnonPAR_Males


