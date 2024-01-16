#!/usr/local/bin/python3.6

import glob
from operator import itemgetter
import re

def main():

    with open("IGV_Final.txt",'r') as f:
        data = f.read().rstrip().split("\n")
    removed=build_removeSet(data)


    for name in glob.glob("./*annotated_additional_Final.txt"): 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_fixed.txt"
        clone=name[2:13]
        print(clone)
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split("\n")
        print(len(data_in))
        
        removed=build_removeSet(data_in)
        clean_list, removed_num=remove_falseSomatic(data_in, removed)
 
        with open(outputfile,'w') as f:
            for line in clean_list:
                f.write(str(line)+'\n')
                
                
                

### build a set with variants marked as remove   
def build_removeSet(data_in):
    removed = set()
    for line in data_in:  
        field=line.split('\t')
        if field[-1] == "remove":
            POSID=field[1].strip('\"') + str(field[2])
            removed.add(POSID)
        else:
            pass
    return removed


### remove false somatic variants    
def remove_falseSomatic(data_in, removed):
    clean_list=[]
    removed_num = 0
    for line in data_in:
        field=line.split('\t')        
        POSID=field[1].strip('\"') + str(field[2])
        if POSID in removed:
            removed_num += 1
        else:
            clean_list.append(line)            
    print("removed variants:", removed_num, "final_variants:", len(clean_list))       
    return clean_list, removed_num
    


main()    








