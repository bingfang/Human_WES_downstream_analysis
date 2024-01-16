#! /bin/bash

module load VEP

files="""4_RPZ27422_D001_2_vs_12_RPZ0113.vcf
6_RPZ27424_D001_2_vs_12_RPZ0113.vcf
RPZ26425_D005_2_vs_12_RPZ0113.vcf
3_RPZ27421_D001_2_vs_12_RPZ0113.vcf
5_RPZ27423_D001_2_vs_12_RPZ0113.vcf
RPZ26198_D005_2_vs_12_RPZ0113.vcf
"""

for file in files

do
    vep -i ./file -o file_vep.txt --offline --cache --dir_cache $VEP_CACHEDIR --species mouse --fasta $VEP_CACHEDIR/mouse.fa --everything --af_exac --vcf

echo "DONE"


