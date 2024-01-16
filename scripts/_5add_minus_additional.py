#!/usr/local/bin/python3.6

import glob
from operator import itemgetter
import re

### check sample name in fieid[-4]

def main():
    with open("IGV_Final.txt", "r") as f:
        data = f.read().split("\n")

    for name in glob.glob("./*annotated.txt"): 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_additional.txt"
        clone=name[2:11]
        print(clone)
        with open(inputfile,'r') as f:
            data_clone = f.read().rstrip().split("\n")
            print(len(data_clone))
        added = add_germline_variants(data1=data_clone, data2=data, clone=clone)
        print(len(added))
        with open(outputfile,'w') as f:
            for line in added:
                f.write(str(line)+'\n')
            
def add_germline_variants(data1, data2, clone): 
    for line in data2:   
        field=line.rstrip().split('\t')
        if clone in field[-4] and field[-2]=="add":  
            data1.append(line.rstrip())
        elif clone in field[-4] and field[-1]=="remove": 
            data1.append(line.rstrip())
        else: 
            data1=data1               
    return data1	

main()    










