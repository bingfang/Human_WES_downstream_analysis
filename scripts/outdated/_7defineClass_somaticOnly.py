#!/usr/local/bin/python3.6

import glob
from operator import itemgetter
import re


def main():

    for name in glob.glob("./*annotated.txt"): 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_Final.txt"
        clone=name[2:9]
        print(clone)
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split("\n")
        print(len(data_in))

        data= defineClass(data_in)
        data=reClass(data)
        data_out=sort_ByAF(data) 
        with open(outputfile,'w') as f:
            f.write(data_in[0].rstrip() + '\tCLASS\n')
            for line in data_out:
                f.write(str(line)+'\n')


### define class based on AF and impact    
def defineClass(data_in):
    data = []
    for line in data_in: 
        field=line.split('\t')
        if "chr" not in field[1]:
            print("not a data line")
        elif float(field[12]) >= 0.9:
            field.insert(16,"I-A")
        elif float(field[12]) < 0.9 and float(field[12]) > 0.2 and field[13] == "HIGH":
            field.insert(16,"I-B")
        elif float(field[12]) < 0.9 and float(field[12]) > 0.2 and field[13] == "MODERATE":
            field.insert(16,"II-A")
        else:
            field.insert(16,"CHECK!")      
        line = "\t".join(field)
        data.append(line)
    return data

### re define class based on amino acid changes
def reClass(data):
    data_out = []
    for line in data: 
        field=line.split('\t')
        if field[16] == "II-A" or field[16] == "I-A" or field[16] == "I-B":
            m=re.search(r'(Ala[0-9]*Gly$)|(Gly[0-9]*Ala$)|(Ser[0-9]*Thr$)|(Thr[0-9]*Ser$)|(Ile[0-9]*Leu$)|(Leu[0-9]*Ile$)',field[15])
            #m=re.search((r'^Ala[0-9]*Gly$'),field[15])
            if m:
                print(m)
                field.insert(16,"II-B")
                line = "\t".join(field)
                data_out.append(line)
            else:    
                data_out.append(line)
        else:
            data_out.append(line)                  
    return data_out
    
### sort list by AF then impact        
def sort_ByAF(data):    
    data_list = []
    for line in data:  
        field=line.split('\t')
        data_list.append(field)
    data_list=sorted(data_list, key =itemgetter(12), reverse = True)
    data_list=sorted(data_list, key =itemgetter(16))    
    data_list=sorted(data_list, key =itemgetter(0))
    data_list=sorted(data_list, key =itemgetter(13))
    data_out=[]
    for item in data_list:
        line="\t".join(item)
        data_out.append(line)
    return data_out
                    

main()    








