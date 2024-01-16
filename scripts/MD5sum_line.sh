#!/usr/bin/bash



# generate file token
for file in *.fastq.gz;
do
    echo $file
    
    md5sum ${file}
done


