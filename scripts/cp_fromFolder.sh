#!/usr/bin/bash



for i in {4..12}
do 
  cp Pair_$i*/*_vs_1_RPZ27911.hard-filtered.vep.vcf.gz ./
done


for i in {4..12}
do 
  gunzip $i*_vs_1_RPZ27911.hard-filtered.vep.vcf.gz 
done



