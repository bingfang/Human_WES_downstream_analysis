#!/usr/bin/bash

#linux use MD5sum
#MAC use MD5



# generate file token
for file in *.fastq.gz;
do
    echo $file
    
    md5sum ${file} > file_MD5.txt
done


