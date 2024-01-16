#!/usr/local/bin/python3.6

import glob
from operator import itemgetter
import re

# default setting :float(AF) > 0.20 and float(DP) > 30

def main():
    for name in glob.glob("./*_vs_1_RPZ27911.hard-filtered.vep.vcf"): #input is the vep file 
        inputfile = str(name)
        outputfile = str(name)[2:13]+"_somatic_annotated.txt"
        clone=name[2:13]
        print(clone)
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        header_line, passed = filter_pass_AF(data_in)
        high,moderate=label_gene_mutation(passed)
        with open(outputfile,'w') as f:
            f.write("GENE\t"+header_line + '\tAF\tIMPACT\tMUTATION\tAMINO ACID CHANGE\n')
            for line in high:
                line = '\t'.join(line)
                f.write(str(line)+'\n')
            for line in moderate:
                line = '\t'.join(line)
                f.write(str(line)+'\n')

# Select lines with quality "PASS". 
# Extract AF for variants.     
# Append lines with AF>20%, DP>10 to the passed list. 
def filter_pass_AF(data_in):
    passed=[] 
    for i in range(len(data_in)-1):   
        if '#CHROM' in data_in[i][0:6]:
            header_line_number = i
            header_line = data_in[i].split("\t")
            header_line = "\t".join(header_line)
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
            if float(AF) > 0.20 and float(DP) > 30:
                line = line +'\t' + str(AF)  
                passed.append(line)                
    return header_line, passed

# Divide data into "HIGH" and "MODERATE, deleterious" groups.
# Extract inforation about mutation and amino acid change. 
# sort by gene names 
def label_gene_mutation(passed):
    high =[]
    moderate =[]
    for line in passed:
        field=line.split('\t')
        info_index = field.index("PASS") + 1
        info=field[info_index]
        sub_field = info.split('|')
        if "HIGH" in info:
            mutation_index = sub_field.index("HIGH") - 1
            gene_index = mutation_index + 2
            mutation = sub_field[mutation_index]
            gene = sub_field[gene_index]
            for item in sub_field[(gene_index+5):(gene_index+9)]:
                if "ENSMUSP" in item and ":p." in item: ### mouse protein
                    item = item.split(":p.")
                    AA=item[1]
                elif "ENSP" in item and ":p." in item:  ### human protein
                    item = item.split(":p.")
                    AA=item[1]
                else:
                    AA=str()
            line=gene+"\t"+line + "\tHIGH\t"+mutation + "\t" + AA
            line=line.split('\t')
            high.append(line)
        elif "MODERATE" in field[info_index] and "deleterious" in str(field[info_index]): 
            sub_field = info.split('|')
            mutation_index = sub_field.index("MODERATE") - 1
            gene_index = mutation_index + 2
            mutation = sub_field[mutation_index]
            gene = sub_field[gene_index]
            for item in sub_field[(gene_index+5):(gene_index+9)]:
                if "ENSMUSP" in item and ":p." in item:   ### mouse protein
                    item = item.split(":p.")
                    AA=item[1]                
                elif "ENSP" in item and ":p." in item:    ### human protein
                    item = item.split(":p.")
                    AA=item[1]
                else:
                    AA=str()
            line=gene+"\t"+line + "\tMODERATE\t"+mutation + "\t" + AA
            line=line.split('\t')
            moderate.append(line)       
    high=sorted(high, key=itemgetter(0))
    print(len(high))  
    moderate=sorted(moderate, key=itemgetter(0))
    print(len(moderate))  
    return high,moderate    

main()











