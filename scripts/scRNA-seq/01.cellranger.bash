#!/bin/bash

#$ -S /bin/bash
#$ -t 1-16:1
#$ -cwd
#$ -V

FASTQLIST=/work23/home/nsasa/data/scRNA/SLE/FASTQLIST

SampleInfo=$(head -n ${SGE_TASK_ID} ${FASTQLIST} | tail -n 1)
SampleID=$(echo ${SampleInfo} | awk '{print $1}')
FQ_R1=$(echo ${SampleInfo} | awk '{print $2}')
FQ_R2=$(echo ${SampleInfo} | awk '{print $3}')

export PATH=/work23/home/nsasa/tools/cellranger-7.0.1:$PATH


### count with HHV-6A&B

if [[ ! -d count_${SampleID} ]]; then
    ln -s ${FQ_R1} ${SampleID}_5DE_S1_L001_R1_001.fastq.gz
    ln -s ${FQ_R2} ${SampleID}_5DE_S1_L001_R2_001.fastq.gz

    cellranger count\
     --id count_${SampleID}\
     --fastqs .\
     --sample ${SampleID}_5DE\
     --transcriptome /work23/home/nsasa/data/gencode.v32_HHV6AB\
     --include-introns true\
     --chemistry SC5P-R2\
     --localcores 8
fi
