#!/usr/local/bin/python3.6


# class VEP handles VEP data 
class Vep:
    def __init__(self,vep_data):  #vep_data is a list from vep file splited by line
        self.vep_data = vep_data

    #select variant with AF>AF_min and DP> DP_min in at least one of the sample. Extract maximum AF number for later sorting
    def filter_AF_DP(self, AF_min, DP_min):
        for i in range(len(self.vep_data)):   
            if '#CHROM' in self.vep_data[i][0:6]:
                head_line = i
                break
        headline=self.vep_data[head_line]
        field=headline.split('\t')
        field_len=len(field)
        sample_index=field.index("FORMAT") + 1        
        filtered = []
        for line in self.vep_data[(head_line + 1):len(self.vep_data)]:
            field=line.split('\t')
            Max_AF_allele = 0
            for sample in field[sample_index:field_len]:
                count=sample.split(":")
                if "GT:AD:AF:DP:" in line:
                    AF=count[2]
                    DP=count[3]
                    if DP != "." and float(DP) > DP_min:
                        AF_allele=AF.split(',')
                        for allele in AF_allele:
                            if allele != "." and float(allele) > Max_AF_allele:
                                Max_AF_allele = float(allele)
                elif "GT:AD:DP:GQ" in line:
                    if len(count)>3:
                        AD=count[1]
                        DP=count[2]
                        if DP != "." and float(DP) > DP_min:
                            AD_allele=AD.split(",")
                            for j in range(1,len(AD_allele)): # first AD_allele is WT
                                AF_allele=float(AD_allele[j])/float(DP)
                                if float(AF_allele) > Max_AF_allele:
                                    Max_AF_allele = float(AF_allele)
                else:
                    print("unusual VEP file",line)
            if Max_AF_allele > float(AF_min):
                line =line + "\t" +  str(Max_AF_allele)
                filtered.append(line)              
        Unique=set(filtered)
        Unique=list(Unique)
        return Unique, head_line


    def get_gene_list(self):
        for i in range(len(self.vep_data)):   
            if '#CHROM' in self.vep_data[i][0:6]:
                header_line_num = i
                headline=self.vep_data[header_line_num]       
                info_index=headline.index("INFO")
                break
            else:
                header_line_num = 0
                info_index= 7
       
        moderate=[]
        high=[]
        for line in self.vep_data[header_line_num:len(self.vep_data)]: ## check if there is a empty line in the end.
            field=line.split('\t')
            info=field[info_index]
            sub_field = info.split("|")
            if "HIGH" in field[info_index]:
                mutation_index = sub_field.index("HIGH") - 1
                gene_index = mutation_index + 2
                mutation = sub_field[mutation_index]
                gene = sub_field[gene_index]
                protein = sub_field[(gene_index+8)]
                if "ENSMUSP" in protein and ":p." in protein:
                    protein = protein.split(":p.")
                    AA=protein[1]
                else:
                    AA="unknown"
                line=gene+"\t"+line + "\tHIGH\t"+ mutation + "\t" + AA
                high.append(line)
            elif "MODERATE" in info and "deleterious" in info:
                mutation_index = sub_field.index("MODERATE") - 1
                gene_index = mutation_index + 2
                mutation = sub_field[mutation_index]
                gene = sub_field[gene_index]
                protein = sub_field[(gene_index+8)]
                if "ENSMUSP" in protein and ":p." in protein:
                    protein = protein.split(":p.")
                    AA=protein[1]
                else:
                    AA="unknown"
                line=gene+"\t"+line + "\tMODERATE\t"+mutation + "\t" + AA
                moderate.append(line)          
        return high, moderate   
    
    
    # divided variants into four groups, PRZ27911_specific, clade1_specific, clade2_specific and others.    
    def get_clades(self):
       
        clade1=[]
        clade2=[]
        clade3=[]
        clade4=[]
        for line in self.vep_data[0:len(self.vep_data)]: 
            field=line.split('\t')
            if "0/0" in field[9] and "0/0" in field[10] and "0/0" in field[11] and "0/0" in field[12] and "0/0" in field[13] and "0/0" in field[14] and "0/0" in field[17] and "0/1" in field[15] and "0/1" in field[16]:
                clade1.append(line)
            elif "0/1" in field[9] and "0/1" in field[10] and "0/1" in field[11] and "0/1" in field[12] and "0/1" in field[13] and "0/1" in field[14] and "0/1" in field[17] and "0/0" in field[15] and "0/0" in field[16]:
                clade2.append(line)
            elif "0/1" in field[9] and "0/1" in field[10] and "0/1" in field[11] and "0/1" in field[12] and "0/1" in field[13] and "0/1" in field[14] and "0/1" in field[17] and "0/1" in field[15] and "0/1" in field[16]:
                clade3.append(line)
            else:
                clade4.append(line)
        return clade1, clade2,clade3,clade4
        
    # detected clades.    
    def detect_clades(self):
       
        clades=[]
        
        for line in self.vep_data[0:len(self.vep_data)]: 
            field=line.split('\t')
            count_het=0


            for info in field[9:]:
                if "0/1" in info:
                    count_het += 1
                elif "./." in info or "0|1" in info or "0/2" in info or "1/" in info:
                    count_het = 0
                    break
                else:
                    pass
            if count_het > 1 and count_het < 8 and "rs"  not in field[2]:
                if len(field[3])==1 and len(field[4])==1:
                    new_line= line + "\t" + str(count_het)
                    clades.append(new_line)
            
        return clades
            
            
  