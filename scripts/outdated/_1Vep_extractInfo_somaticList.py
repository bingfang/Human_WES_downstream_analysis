#!/usr/local/bin/python3.6

import glob
import re
import operator 
from operator import itemgetter

AFmin = float(input("Enter minimal AF: "))
DPmin = float(input("Enter minimal DP: "))
### check if the vep file name started with the clone name.
def main():
    for name in glob.glob("./*hard-filtered.vep.vcf"): ### vep files generated from dragen
        inputfile = str(name)
        clone=inputfile[2:13]  #### some names have a "_" at the end
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip()  ### remove the empty line at the end.
            data_in =data_in.split('\n') ### read the file and create a list split by newline
        new_header,filtered = filter(data_in,AFmin,DPmin)
        print(len(filtered))
        somatic_high,somatic_moderate=label_clone_gene_mutation(filtered,clone)
        with open("somatic_high_moderate.txt",'a') as f:
            f.write(new_header+'\n') ### write the new header
            for item in somatic_high:
                item="\t".join(item)
                f.write(item+'\n')
            for item in somatic_moderate:
                item="\t".join(item)
                f.write(item+'\n')

def filter(data_in,AFmin,DPmin):
    filtered=[]
    for i in range(len(data_in)):   
        if '#CHROM' in data_in[i][0:6]:
            header_line_number = i
            new_header="Gene\t" + data_in[header_line_number] + '\tAF\tHIGH\tMUTATION\tAMINO ACID CHANGES\tclone' 
            break  
    for line in data_in[(header_line_number + 1):]:
        field=line.split('\t')
        if field[6]=="PASS":
            reads=field[10].split(':')
            AF=reads[3]
            AD=reads[2]
            AD=AD.split(",")
            AD=[float(i) for i in AD]
            DP=sum(AD)
            if float(AF) > AFmin and float(DP) > DPmin:
                line = line + "\t"+ AF 
                filtered.append(line)
    return  new_header, filtered           
            
            
def label_clone_gene_mutation(filtered,clone):
    somatic_high = []
    somatic_moderate = []
    for line in filtered:
        field=line.split('\t')         
        info_index = 7 ####field.index("PASS") + 1
        info = field[info_index]
        sub_field = info.split('|')
        if "HIGH" in str(field[info_index]):
            mutation_index = sub_field.index("HIGH") - 1
            gene_index = mutation_index + 2
            mutation = sub_field[mutation_index]
            gene = sub_field[gene_index]
            for item in sub_field[(gene_index+5):(gene_index+9)]:
                if "ENSMUSP" in item and ":p." in item:  ### mouse protein
                    item = item.split(":p.")
                    AA=item[1]
                elif "ENSP" in item and ":p." in item:   ### human protein
                    item = item.split(":p.")
                    AA=item[1]
                else:
                    AA="unknown"
            line=gene+"\t"+line + "\tHIGH\t"+ mutation + "\t" + AA +"\t"+ clone
            line=line.split('\t')
            somatic_high.append(line)
        elif "MODERATE" in field[info_index] and "deleterious" in str(field[info_index]): 
            sub_field = info.split('|')
            mutation_index = sub_field.index("MODERATE") - 1
            gene_index = mutation_index + 2
            mutation = sub_field[mutation_index]
            gene = sub_field[gene_index]
            for item in sub_field[(gene_index+5):(gene_index+9)]:
                if "ENSMUSP" in item and ":p." in item:  ### mouse protein
                    item = item.split(":p.")
                    AA=item[1]
                elif "ENSP" in item and ":p." in item:  ### Human protein
                    item = item.split(":p.")
                    AA=item[1]
                else:
                    AA="unknown"
            line=gene+"\t"+line + "\tMODERATE\t"+ mutation + "\t" + AA +"\t"+ clone
            line=line.split('\t')
            somatic_moderate.append(line)
    ### sort by gene name               
    somatic_high=sorted(somatic_high, key=itemgetter(0))  
    somatic_moderate=sorted(somatic_moderate, key=itemgetter(0))  
    return somatic_high,somatic_moderate
            
main()













