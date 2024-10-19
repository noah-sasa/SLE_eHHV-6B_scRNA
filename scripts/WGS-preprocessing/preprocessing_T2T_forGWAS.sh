#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V


ref="/work23/home/nsasa/data/Reference/CHM13/v2.0/chm13v2.0_maskedY.fa"
interval="/work23/home/nsasa/data/Reference/CHM13/v2.0/interval_list/chm13v2.0_maskedY_Autosomal.interval_list"
intervals_XYM="/work23/home/nsasa/data/Reference/CHM13/v2.0/interval_list/chm13v2.0_maskedY_chrXYM.interval_list"
dbsnp="/work23/home/nsasa/data/Reference/CHM13/CHM13v2.0_GATK_Resource/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.t2t-chm13-v2.0.vcf.gz"
known_indels="/work23/home/nsasa/data/Reference/CHM13/CHM13v2.0_GATK_Resource/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.t2t-chm13-v2.0.vcf.gz"
Mills_1000G_indels="/work23/home/nsasa/data/Reference/CHM13/CHM13v2.0_GATK_Resource/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.t2t-chm13-v2.0.vcf.gz"


### BAM2FASTQ
if [[ ! -f "${SampleID}/${SampleID}_R1.fastq.gz" ]]; then
    if [[ ! -f "${SampleID}/${SampleID}.sorted.bam" ]]; then
        samtools sort\
         -o ${SampleID}/${SampleID}.sorted.bam\
         -n\
         -@ 8\
         ${BAM}
    fi

    samtools fastq\
     -1 ${SampleID}/${SampleID}_R1.fastq.gz\
     -2 ${SampleID}/${SampleID}_R2.fastq.gz\
     -O\
     ${SampleID}/${SampleID}.sorted.bam
     # output quality in the OQ tag if present
fi


if [[ -f "${SampleID}/${SampleID}_R1.fastq.gz" ]]; then
    if [[ -f "${SampleID}/${SampleID}.sorted.bam" ]]; then
        rm ${SampleID}/${SampleID}.sorted.bam
    fi
fi

### FASTQ2SAM
${FASTQ2SAM} \
  -1 ${SampleID}/${SampleID}_R1.fastq.gz \
  -2 ${SampleID}/${SampleID}_R2.fastq.gz \
  -o ${SampleID}/${SampleID}.unmapped.rg.preSort.bam \
  -n ${SampleID} \
  -p ${platform} \
  --hash ${BamHash_fastq} \
  --hash-no-quality ${BamHash_fastq_noQuality}

### SortSAM
mkdir -p ${tmpdir}
${GATK} --java-options "-Xms4000m -Xmx4000m -XX:NewSize=3000m -XX:+UseSerialGC" SortSam \
  -I ${SampleID}/${SampleID}.unmapped.rg.preSort.bam \
  -O ${SampleID}/${SampleID}.unmapped.rg.bam \
  --SORT_ORDER "queryname" \
  --CREATE_INDEX false \
  --CREATE_MD5_FILE false \
  --MAX_RECORDS_IN_RAM 1000000 \
  --TMP_DIR ${tmpdir} 2>&1 | tee ${SortSam}


### BwaMemAndMba
${GATK} --java-options "-Xms330m -Xmx330m -XX:NewSize=300m -XX:+UseSerialGC" SamToFastq \
  --INPUT ${SampleID}/${SampleID}.unmapped.rg.bam \
  --FASTQ /dev/stdout \
  --INTERLEAVE true \
  --NON_PF true | 
  ${BWA} mem -K 100000000 -p -t ${threads} -Y ${ref} /dev/stdin - 2> >(tee ${log.bwa} >&2) |
  ${GATK} --java-options "-Xms250m -Xmx250m -XX:NewSize=200m -XX:+UseSerialGC" MergeBamAlignment \
  --VALIDATION_STRINGENCY SILENT \
  --EXPECTED_ORIENTATIONS FR \
  --ATTRIBUTES_TO_RETAIN X0 \
  --ALIGNED_BAM /dev/stdin \
  --UNMAPPED_BAM ${SampleID}/${SampleID}.unmapped.rg.bam \
  --OUTPUT ${SampleID}/${SampleID}.aligned.rg.bam \
  --REFERENCE_SEQUENCE ${ref} \
  --SORT_ORDER "unsorted" \
  --IS_BISULFITE_SEQUENCE false \
  --ALIGNED_READS_ONLY false \
  --CLIP_ADAPTERS false \
  --MAX_RECORDS_IN_RAM 2000000 \
  --ADD_MATE_CIGAR true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --PROGRAM_RECORD_ID "bwamem" \
  --PROGRAM_GROUP_VERSION ${bwa_version} \
  --PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -t ${threads} -Y ${ref}" \
  --PROGRAM_GROUP_NAME "bwamem" \
  --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
  --ALIGNER_PROPER_PAIR_FLAGS true \
  --UNMAP_CONTAMINANT_READS true 2>&1 | tee ${log.mba}

### MarkDuplicates
${GATK} --java-options "-Xms12000m -Xmx12000m -XX:NewSize=8000m -XX:+UseSerialGC" MarkDuplicates \
  -I ${SampleID}/${SampleID}.aligned.rg.bam \
  -O ${SampleID}/${SampleID}.aligned.rg.marked.bam \
  -M ${SampleID}/${SampleID}.marked_dup_metrics.txt \
  --ASSUME_SORT_ORDER "queryname" \
  --VALIDATION_STRINGENCY SILENT \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_MD5_FILE false \
  --MAX_RECORDS_IN_RAM 1000000 \
  2>&1 | tee ${log}

### Sort_SetTags
mkdir -p ${tmpdir}
${GATK} --java-options "-Xms4000m -Xmx4000m -XX:NewSize=3000m -XX:+UseSerialGC" SortSam \
  -I ${SampleID}/${SampleID}.aligned.rg.marked.bam \
  -O /dev/stdout \
  --SORT_ORDER "coordinate" \
  --CREATE_INDEX false \
  --CREATE_MD5_FILE false \
  --MAX_RECORDS_IN_RAM 1000000 \
  --TMP_DIR ${tmpdir} 2> >(tee ${SortSam} >&2) |
${GATK} --java-options "-Xms4000m -Xmx4000m -XX:NewSize=3000m -XX:+UseSerialGC" SetNmMdAndUqTags \
  -R ${ref} \
  -I /dev/stdin \
  -O ${SampleID}/${SampleID}.aligned.rg.marked.cs.fixed.bam \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE false 2>&1 | tee ${SetNmMdAndUqTags}

### BQSR
${GATK} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms300m -Xmx300m -XX:NewSize=180m -XX:+UseSerialGC" BaseRecalibrator \
  -I ${SampleID}/${SampleID}.aligned.rg.marked.cs.fixed.bam \
  -R ${ref} \
  -L ${interval} \
  --use-original-qualities \
  --known-sites ${dbSNP} \
  --known-sites ${known_indels} \
  --known-sites ${Mills_1000G_indels} \
  -O ${recal_file} 2>&1 | tee ${BaseRecalibrator}
  
${GATK} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms110m -Xmx110m -XX:NewSize=55m -XX:+UseSerialGC" ApplyBQSR \
  -I ${SampleID}/${SampleID}.aligned.rg.marked.cs.fixed.bam \
  -R ${ref} \
  --bqsr-recal-file ${recal_file} \
  -O ${SampleID}/${SampleID}.aligned.rg.marked.cs.fixed.bqsr.bam \
  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
  --add-output-sam-program-record \
  --use-original-qualities \
  --emit-original-quals 2>&1 | tee ${ApplyBQSR}
