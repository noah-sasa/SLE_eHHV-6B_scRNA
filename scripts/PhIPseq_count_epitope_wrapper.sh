#!/bin/bash

# Author: Shohei Kojima @ RIKEN (shohei.kojima@riken.jp)
# Version:
#   GNU bash 5.1.16(1)-release (x86_64-pc-linux-gnu)
#   Python 3.7.4
#   fastp 0.20.1
#   pigz 2.6


sample=sample_name
fq1=sample_name_L001_R1_001.fastq.gz
fq2=sample_name_L001_R2_001.fastq.gz
trimmed_fq1=./QCed_fq/${sample}_1.fq
trimmed_fq2=./QCed_fq/${sample}_2.fq

mkdir -p QCed_fq
mkdir -p fastp_html
mkdir -p results


# QC of reads
fastp \
-i ${fq1} \
-o ${trimmed_fq1} \
-I ${fq2} \
-O ${trimmed_fq2} \
-l 6 \
-3 -W 4 -M 20 \
-t 1 -T 1 \
-x \
--compression 1 \
--thread 4 \
--html ./fastp_html/${sample}.html

pigz -p 4 -1 ${trimmed_fq1} ${trimmed_fq2}


# Count detected epitopes
fq1=${trimmed_fq1}.gz
fq2=${trimmed_fq2}.gz
outfile=./results/${sample}.oligo_count.tsv

python PhIPseq_count_epitope.py ${fq1} ${fq2} ${outfile}
