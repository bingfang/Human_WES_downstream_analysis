#! /bin/bash

module load VEP

for f in *_passed_AFDP_filtered.txt

do 
    echo $f
    # get basename of files
    bse=$(basename ${f} _passed_AFDP_filtered.txt)
    
    #VEP commend line on Biowulf
    vep -i ./${f} -o ${bse}_passed_AFDP_filtered_vep.txt --offline --cache --dir_cache $VEP_CACHEDIR --species mouse --fasta $VEP_CACHEDIR/mouse.fa --everything --af_exac --vcf
    
    
done