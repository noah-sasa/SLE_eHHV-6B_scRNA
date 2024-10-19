#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V


RGTAG="@RG\tID:id\tCN:OSAKA\tLB:Macrogen\tPL:ILLUMINA\tPU:$FCBARCODELANE\tSM:$SAMPLEID"

### BWA MEM
/work1/app/bwa/0.7.8/bwa mem -t ${NP} -M /work1/home/ysuenari/WGSdata/app/bwa/hs37d5.fa -R “$RGTAG”  \
 $TRIMED_PAIRED_READ1 $TRIMED_PAIRED_READ2 2>$BWALOG | \
/work1/app/samtools/1.9/bin/samtools view -@ $NP -b - | \
/work1/app/samtools/1.9/bin/samtools sort -@ $NP - -o $BWABAM > $SAMTOOLSLOG 2>&1
/work1/app/samtools/1.9/bin/samtools index $BWABAM >> $SAMTOOLSLOG 2>&1


### MarkDuplicates
/usr/bin/java -XX:ParallelGCThreads=6 -jar /work1/app/picard/2.17.2/picard.jar \
MarkDuplicates \
INPUT=${BWABAM_L1}  \
INPUT=${BWABAM_L2}  \
INPUT=${BWABAM_L3}  \
OUTPUT=${MDBAM} \
ASSUME_SORTED=false \
REMOVE_DUPLICATES=false \
METRICS_FILE=${METRICSFILE}  \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT > ${LOGFILE} 2>&1
