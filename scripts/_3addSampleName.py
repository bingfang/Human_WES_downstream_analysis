#!/usr/local/bin/python3.6
# check Format of GT

import glob
from operator import itemgetter

def main():
    for name in glob.glob("./combined_germline_somatic_marked.txt"):
        inputfile=str(name)

        outputfile = name[:-4] + "_addSampleName.txt"
        with open (inputfile,"r") as f:
            data=f.read().rstrip()  # remove the empty line at the end
            data=data.split('\n')
        samples,sample_index_start= build_sample_list(data)   
        addSampleName=mark_samples(data, samples,sample_index_start)
        data_3col = add_3column(addSampleName, samples)
        data_out = sort_data(data_3col, samples)
        with open (outputfile,"w") as f:
            for item in data_out:
                f.write(item +"\n")

### build a list of sample name.              
def build_sample_list(data):
    samples=[]
    for line in data:
        if "#CHROM" in line and "clone" not in line: ### The header for germline doesn't have clone.
            field=line.split("\t")
            sample_index_start=field.index("FORMAT")+1
            sample_index_end=field.index("AF")-1
            for i in range(sample_index_start, (sample_index_end + 1)):  
                samples.append(field[i])
            break     
    print(samples)
    return samples,sample_index_start 

### add sample name.
def mark_samples(data, samples,sample_index_start):
    addSampleName = []
    for line in data:
        field=line.rstrip().split('\t')
        if len(field) > 18 and "#CHROM" not in field[1] : ### variants from germline calling have more fields
            name =""
            for j in range(len(samples)):
                sample_name=samples[j]
                sample= str(field[(sample_index_start +j)])

                count=  sample.split(":")
                GT=count[0]

                if GT != '0/0' and GT != './.':
                    name = sample_name + ';'+ name              
            line = line.rstrip() + '\t' +name
            addSampleName.append(line)
        else:      
            addSampleName.append(line)
    print(len(addSampleName))        
    return addSampleName

### add three column, total column =len(samples) + 19
### flag germline_only high_impact variants
### flag variant which is germline_only, moderate_impact and in multiple samples
def add_3column(addSampleName, samples):
    data_3col=[]
    for line in addSampleName:
        line=line.rstrip()
        field=line.split("\t")

        if "Both" in line and len(field) < 20:
            line=line +len(samples)*"\t"+"\tremove"
        elif "Both" in line and len(field) > 20:  
            line = line + "\tnone\tnot\tnot"
        elif "#CHROM" in line and len(field) < 20:
            line=line +"\tremove" 
        elif "#CHROM" in line and len(field) > 20:
            line="\t".join(field[:-1]) +"\tCaller\tsample_name\tNotes\tAdd_to_somatic_list\tRemove_from_somatic_list"    
            print("header line: ",line)
        elif "germline_only" in line and "HIGH" in line:
            line=line+ "\tCHECK\tadd\tnot"
        elif "germline_only" in line and "MODERATE" in line:
            l=(field[-1]).split(';')
            if len(l) < 3 :
                line=line+ "\tnone\tadd\tnot"
            else:
                line=line+ "\tnone\tcheck\tnot"
        elif "somatic_only" in line:
            line = "\t".join(field[0:12])+(len(samples)-1)*"\t" +"\t".join(field[12:16]) + "\t"+ field[-1] + "\t"+ field[-2]+ "\tnone\tnot\tnot"
        else:
            line =line+"\tcheck this line"
        data_3col.append(line)
    return data_3col


### sort by AF,impact, caller         
def sort_data(data_3col, samples):
    data_list =[]
    for line in data_3col:
        field=line.split("\t")
        if len(field) == len(samples) + 19:
            data_list.append(field)
        else:
            pass
    data_list=sorted(data_list, key =itemgetter(-9), reverse = True)
    data_list=sorted(data_list, key =itemgetter(-8))
    data_list=sorted(data_list, key =itemgetter(-5))
    data_out=[]
    for item in data_list:
        line="\t".join(item)
        data_out.append(line)
    return data_out
main()           
            

