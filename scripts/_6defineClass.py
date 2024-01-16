#!/usr/local/bin/python3.6

# check clone name

import glob
from operator import itemgetter
import re

def main():

    with open("IGV_Final.txt",'r') as f:
        sample_info = f.read().rstrip().split("\n")
    samples=find_sample_index(sample_info)


    for name in glob.glob("./*annotated_additional.txt"): 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_Final.txt"
        clone=name[2:13]
        print(clone)
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split("\n")
        print(len(data_in))
        data = change_format(data_in, clone, samples)
        data= defineClass(data)
        data=reClass(data)
        data_out=sort_ByAF(data) 
        with open(outputfile,'w') as f:
            f.write(data_in[0].rstrip() + '\tCLASS\n')
            for line in data_out:
                f.write(str(line)+'\n')

### build a dictionary of sample name with sample index               
def find_sample_index(sample_info):
    samples={}
    for line in sample_info:
        if "GENE" in line and "#CHROM" in line:
            field=line.split("\t")
            sample_index_start=field.index("FORMAT")+1
            sample_index_end=field.index("AF")-1
            for i in range(sample_index_start, (sample_index_end + 1)):  
                samples[field[i]] = i
            break     
    print(samples)
    return samples  

###recalculate AF 
### keep only sample data with variants called           
def change_format(data_in, clone, samples): 
    data = []
    for line in data_in: 
        if "chr" in line:  
            field=line.split('\t')
            if "germline" in field[-5] :
                for item in samples:
                    if clone in item:
                        index=samples[item]
                        print(index)
                sub_field=field[index].split(":")
                reads=sub_field[1].split(",")
                total = 0
                for read in reads:
                    total += int(read)
                AF=""
                for read in reads[1:]:
                    af = round(float(read)/(total+0.00001),3)
                    AF += str(af)
                print(AF)               
                new_field = field[0:10]+["",field[index], AF] + field[-8:]  ### check here
                print(new_field[12:15])
                line = "\t".join(new_field)
                data.append(line)
            else:
                data.append(line)
        else:
            print("invalid line")
    print(len(data))  
    return data

### mark variants to be removed.
### define class based on AF and impact    
def defineClass(data):
    data_out = []
    for line in data:  
        field=line.split('\t')
        if field[-1] == "remove" or "ras" in field[0]:
            field.insert(16,"CHECK!")
        elif len(field[12]) > 5 or "." not in field[12]:  #field[12]is AF
            field.insert(16,"CHECK!")
        elif float(field[12]) >= 0.9:
            field.insert(16,"I-A")
        elif float(field[12]) < 0.9 and float(field[12]) > 0.2 and field[13] == "HIGH":
            field.insert(16,"I-B")
        elif float(field[12]) < 0.9 and float(field[12]) > 0.2 and field[13] == "MODERATE":
            field.insert(16,"II-A")
        else:
            field.insert(16,"CHECK!")      
        line = "\t".join(field)
        data_out.append(line)
    return data_out

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

### sort list by AF, class, gene name, impact    
### define splice variant.    
def sort_ByAF(data):     
    data_list = []
    for line in data:  
        field=line.split('\t')
        if field[-2]=="" and "splice" in field[-3]:
            field[-2] = "Splice variant"
        elif field[15]=="unknown" and "splice" in field[14]:  
            field[15] = "Splice variant"
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








