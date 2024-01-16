#!/usr/local/bin/python3.6

import glob

def main():
    with open("./500_genes.txt",'r') as f:
        gene500 = f.read().rstrip().split('\n') 
    genelist = make_genelist(gene500) 
    for name in glob.glob("./*_somatic_annotated_additional_Final_fixed.txt"):

        with open(name,'r') as f:
            data_in = f.read().rstrip().split('\n')
   
    
        select500, not_selected =  select_line(data_in, genelist)
        print(name, "skipped variants:",not_selected, "selected variants:", len(select500))     
        with open("tru500_detected.txt",'a') as f:
            f.write(name+'\n')
            for line in select500:
                f.write(str(line)+'\n')
            
# read in gene list form tru500 and make a list
def make_genelist(gene500):
    genelist =[]
    for line in gene500:
        lst = line.split('\t')
        
        print(lst)
        genelist = genelist + lst
    print(genelist)

    return genelist
    

    
def select_line(data_in, genelist):
    select500=[]

    not_selected =0
    for line in data_in:
        field=line.split('\t')
        gene=field[0].upper()
        if gene in genelist:
            select500.append(line)
        else:
            not_selected += 1

    return select500, not_selected

main()
    