#!/usr/bin/bash

####  multiple 20bp from blast P2A Flag Kras_exon2
#$Ex2_tail="ACGAATATGATCCAACAATA"
#$Ex2_3_11 = "TAAACTTGTGGTAGTTGGAG"
#${Flag_tail} = "acaaggatgacgatgacaag"
#${Flag_head} = "gactacaaagaccatgacgg"
#${p2A} = "aggctggagacgtggaggag"
#${blast_mid}="attgctgccctctggttatg"
#${blast_start}="aagcctttgtctcaagaaga"
#echo "tggagctgttggcgtaggcaag"| awk '{print toupper($0)}' | echo $seq

#### two 20bp around ATG of kras exon2
# $recom="AAGCGGAGGTATGACTGAAT"
# $WT="CCTGCTGAAAATGACTGAAT"

#lst=("ACGAATATGATCCAACAATA" "TAAACTTGTGGTAGTTGGAG" "acaaggatgacgatgacaag" "gactacaaagaccatgacgg" "aggctggagacgtggaggag" "attgctgccctctggttatg" "aagcct$
lst=("AAGCGGAGGTATGACTGAAT" "CCTGCTGAAAATGACTGAAT")

# To upper case of sequences
for seq in ${lst[@]};
    do
       echo ${seq} | tr '[:lower:]' '[:upper:]'

       SEQ=$(echo ${seq} | tr '[:lower:]' '[:upper:]')

       LST+=($SEQ)
done

# grep 
for SEQ in ${LST[@]};
do
    for file in *RPZ*_R1_001.fastq.gz;
    do
        echo $file
        name=$(basename ${file} _R1_001.fastq.gz)
        echo ${name}
        zgrep -B 1 ${SEQ} ${file} >> ${name}_${SEQ}.txt
        sed 's/@VH00249:1:AAANHJWM5:1:/>/g' ${name}_${SEQ}.txt >> ${name}_${SEQ}_1.txt
        sed 's/:/_/g' ${name}_${SEQ}_1.txt >> ${name}_${SEQ}_2.txt
    done
done

for file in *RPZ*_2.txt;
do
    echo $file >> summary.txt
    wc $file >> summary.txt
done