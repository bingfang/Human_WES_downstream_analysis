#!/usr/local/bin/python3.6

import glob

#save only "passed" and AF>0.2, DP>10 variants


def main():
    for name in glob.glob("./*hard-filtered.pass.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_passed_AFDP_filtered.txt"
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        header_line_number, passed = find_header(data_in)
        passed=filter_PASS(data_in,header_line_number, passed)
        AF_filtered=filter_AF(passed,header_line_number)
        with open(outputfile,'w') as f:
            for line in AF_filtered:
                f.write(str(line)+'\n')

# append header to the passed list.    
def find_header(data_in):
    passed=[] 
    for i in range(len(data_in)):   
        if '#CHROM' not in data_in[i][0:6]:
            passed.append(data_in[i])
        else:
            passed.append(data_in[i])
            header_line_number = i
            break  
    return header_line_number, passed

# append lines with quality "PASS" to the filtered list.   
def filter_PASS(data_in,header_line_number,passed):    
    for line in data_in[(header_line_number + 1):(len(data_in))]:
        field=line.split('\t')
        if field[6]=="PASS":
            passed.append(line)
    print(len(passed))
    return passed
    
# append lines with AF>20%, DP>10 to the filtered list.    
def filter_AF(passed,header_line_number):
    AF_filtered=[]    
    for line in passed[(header_line_number + 1):(len(passed))]:
        field=line.split('\t')
        reads=field[10].split(':')
        AF=reads[2]
        DP=reads[5]
        if float(AF) > 0.20 and float(DP) > 10:
            AF_filtered.append(line)
    print(len(AF_filtered))        
    return AF_filtered

main()
