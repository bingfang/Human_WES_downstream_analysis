#!/usr/bin/bash

## move file out of folders
#mv ./Sample*/*_001.fastq.gz ./

#rename
for file in *_passed_AFDP_filtered_vep.txt;
do
    echo $file
    name=$(basename ${file} _passed_AFDP_filtered_vep.txt)
    echo ${name}
    cp ${file} ${name}_FINALmutect2_vep.txt
done

for file in *_passed_AFDP_filtered_vep.txt_summary.html
do
    name=$(basename ${file} _passed_AFDP_filtered_vep.txt_summary.html)
    cp $file ${name}_FINALmutect2_vep.txt_summary.html
done


