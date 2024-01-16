#!/usr/local/bin/python3.6
from operator import itemgetter
import vep_clade
import glob

def main():
    ## enter one VEP file for analysis
    for name in glob.glob("../dragen_germline/germline_joint_genotyping/Ontarget*nonRPZ113.vep.vcf"):
        inputfile=str(name)
        with open (inputfile,"r") as f:
            data=f.read().rstrip().split('\n')
        print("input lines:", len(data))    
        ## enter first filter condition   
        AF_min=float(input("Enter minimal AF: "))
        DP_min=float(input("Enter minimal DP: "))
        
        Unique,clade1, clade2, clade3,clade4,clades, head_line=filter(data,AF_min, DP_min)
        

        #print("header line number:", head_line)
        print("Total passed filter:",len(Unique), "\nRPZ27911_specific:",len(clade3),"\nclade1_specific:",len(clade1), "\nclade2_specific:", len(clade2), "\nClone_specific and others: ",len(clade4))       
        print("clades_2_7:",len(clades))

        with open ("../clade1_specific.txt","w") as f:
            f.write(data[head_line]+ "\t" + "AF\n")
            for line in clade1:
                f.write(line+'\n')           
        with open ("../clade2_specific.txt","w") as f:
            f.write(data[head_line]+ "\t" + "AF\n")
            for line in clade2:
                f.write(line+'\n')
        with open ("../RPZ27911_specific.txt","w") as f:
            f.write(data[head_line]+ "\t" + "AF\n")
            for line in clade3:
                f.write(line+'\n')
        with open ("../clone_specific_others.txt","w") as f:
            f.write(data[head_line]+ "\t" + "AF\n")
            for line in clade4:
                f.write(line+'\n')
        with open ("../all_clades.txt","w") as f:
            f.write(data[head_line]+ "\t" + "AF\n")
            for line in clades:
                f.write(line+'\n')

## filter by minimal AF and DP, extract gene information, and divide to two groups (high, moderate)
def filter(data,AF_min, DP_min):
    vepA=vep_clade.Vep(data) #build a Vep instance
    Unique, head_line = vepA.filter_AF_DP(AF_min, DP_min) #filter Vep instance
    vepB=vep_clade.Vep(Unique) #build a Vep instance
    clade1, clade2, clade3,clade4=vepB.get_clades() # extract clade1 and 2 gene list

    clades=vepB.detect_clades()
    return Unique,clade1, clade2, clade3,clade4,clades, head_line

## sort by AF    
def sort_sites(vep_list_set):
    list_list=[]
    for line in vep_list_set:
        l = line.split('\t') # itemgetter only work for list of list.
        if l[0] not in ["Skint6","Skint5","Sfi1","Mrgpra3","Mrgpra9","Mrgpra5","Skint11","Mrgprb5","Nlrp1b","Gm10110","Ighg2c", "Gm8909"]:
            list_list.append(l)
    list_list=sorted(list_list, key=itemgetter(0)) # sorted by gene name 
    return list_list      
    
main()    

       

        
         
