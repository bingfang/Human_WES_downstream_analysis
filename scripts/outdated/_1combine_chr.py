#!/usr/local/bin/python3.6

import glob

samples="""12_RPZ0113+3_RPZ28460_S7
12_RPZ0113+4_RPZ28465_S8
12_RPZ0113+5_RPZ28441_D001_2_S9
12_RPZ0113+6_RPZ28442_D001_2_S10
12_RPZ0113+7_RPZ28443_D001_2_S11
12_RPZ0113+8_RPZ27082_D002_2_S12"""

samples=samples.split('\n')


dic={}
for sample in samples:
    dic[sample] = ''
    for name in glob.glob("./*RPZ*RPZ*.vcf"): #input is the prefiltered file 
        output_name= sample+"_chrom_all.vcf"
        if sample in name:
            with open(name,'r') as f:
                data_in = f.read()
                dic[sample] += data_in
    with open(output_name,'w') as f: 
        f.write(dic[sample])
