#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -t 24-25:1
#$ -V

# -pe OpenMP 1 -l s_vmem=12G -l m_mem_free=12G

SampleID=$1
sex=$2
mapping_dir=$3


mkdir -p ${SampleID}
cd ${SampleID}

NORMAL_BAM=${mapping_dir}/data_processing/ApplyBQSR/${SampleID}/${SampleID}.aligned.rg.marked.cs.fixed.bqsr.cram

#REF=/work23/home/nsasa/data/Reference/hg38/hg38/Homo_sapiens_assembly38.fasta
REF=/work23/home/nsasa/data/Reference/CHM13/v2.0/chm13v2.0_maskedY.fa

CHM13_intervals_dir=/work23/home/nsasa/data/Reference/CHM13/v2.0/haplotypecaller_intervals
interval_list=$(cat ${CHM13_intervals_dir}/intervals.txt | head -n ${SGE_TASK_ID} | tail -n 1)
INTERVAL=${CHM13_intervals_dir}/${interval_list}

chromosome=$(echo $interval_list | awk -F ".interval_list" '{print $1}' | awk -F "chm13v2.0_maskedY_" '{print $2}')

if [[ "${sex}" == "M" ]]; then
	if [[ "${SGE_TASK_ID}" -eq 24 ]]; then
		# chrX nonPAR
		sample_ploidy=1
	elif [[ "${SGE_TASK_ID}" -eq 25 ]]; then
		# chrY nonPAR
		sample_ploidy=1
	else
		sample_ploidy=2
	fi
else
	if [[ "${SGE_TASK_ID}" -eq 25 ]]; then
		# chrY nonPAR
		exit 0
	else
		sample_ploidy=2
	fi
fi

source /work23/home/nsasa/conda-pack/gatk4.2.6.1/bin/activate
GATK_DIR=/work23/home/nsasa/data/Script/hg38_for_01/tools/gatk-4.2.6.1

mkdir -p scattered-HaplotypeCaller
mkdir -p log_HCscatter

### Haplotypecaller subinterval (gvcf -G AS_StandardAnnotation -ERC GVCF)
if [ ! -f "scattered-HaplotypeCaller/${chromosome}_HC.g.vcf.gz" ]; then
#	mkdir -p M2Tmp_${ID}	 --tmp-dir M2Tmp_${ID}\
	/usr/bin/time -f "Memory:%M KB time:%E" -o log_HCscatter/HC_${chromosome}.txt \
	${GATK_DIR}/gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 HaplotypeCaller\
	 -R ${REF}\
	 -I ${NORMAL_BAM}\
	 -contamination 0.0\
	 -G StandardAnnotation -G StandardHCAnnotation\
	 -G AS_StandardAnnotation\
	 -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90\
	 -ERC GVCF\
	 --native-pair-hmm-threads 1\
	 -L ${INTERVAL}\
	 --sample-ploidy ${sample_ploidy}\
	 -O scattered-HaplotypeCaller/${chromosome}_HC.g.vcf.gz\
	 > log_HCscatter/HC_${chromosome}.log 2>&1
#	rm -rf M2Tmp_${ID}
	echo "HaplotypeCaller_${chromosome} Done"
else
	echo "HaplotypeCaller_${chromosome} Already Done"
fi

# ### extract het + merge somatic SNP -> Germline-het-SNP.vcf.gz
# if [ ! -f "scattered-HaplotypeCaller/${ID}_germline-het-SNP.vcf.gz" ]; then
# 	cd scattered-HaplotypeCaller
# 	bash /work23/home/nsasa/data/Script/hg38/CNV/snpfilter_mod.sh ${PRIMARY_BAM} ${ID}_HC.vcf.gz /work23/home/nsasa/data/Reference/hg38/fasta_gz/Homo_sapiens_assembly38.fasta.gz ${ID}
# 	cd ..
# 	echo "Germline-het_${ID} Done"
# else
# 	echo "Germline-het_${ID} Already Done"
# fi
# 
# source /work23/home/nsasa/conda-pack/gatk4.1.9.0/bin/deactivate

