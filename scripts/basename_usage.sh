#!/usr/bin/bash

#cat fastq files


for file in ./*_L001_R1_001.fastq.gz
do
    name=$(basename ${file} _L001_R1_001.fastq.gz)
    echo ${name}
    cat ${name}*R1* > ${name}_cat.R1.fastq.gz
    cat ${name}*R2* > ${name}_cat.R2.fastq.gz
done
 

