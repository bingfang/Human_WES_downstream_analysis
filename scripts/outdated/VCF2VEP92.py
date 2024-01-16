#!/usr/local/bin/python3.6

import glob
import os

#convert vcf files to vep files


def main():
    for name in glob.glob("./*RP*RPZ0113*passed_AFDP_filtered.txt"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_vep.txt"
        os.system("/opt/nasapps/development/vep/ensembl-vep-release-92/vep --everything --af_exac --species mus_musculus --fork 24 --cache --port 3337 --dir /SeqIdx/vepdb/caches --vcf -i /scratch/cluster_tmp/xubr/" + inputfile + " -o " + " /is2/projects/RAS-intl/static/Genomics/" + outputfile)

main()


"""
### transfer filtered vcf files to moab xubr account, convert vcf to VEP, transfer vep files back to local.
scp _1VCF2VEP.py xubr@moab.ncifcrf.gov:/scratch/cluster_tmp/xubr
scp RP*2*RPZ113*passed_AFDP_filtered.txt xubr@moab.ncifcrf.gov:/scratch/cluster_tmp/xubr

os.system("/opt/nasapps/development/vep/80/bin/variant_effect_predictor.pl --input_file /scratch/cluster_tmp/xubr/" + inputfile + " --species  mus_musculus --format vcf --fork 8  --output_file /scratch/cluster_tmp/xubr/" + outputfile + " --html --sift b --regulatory --protein --symbol  --ccds --canonical --biotype --domains --vcf --check_existing --check_alleles --check_svs  --stats_file $out\.vep.html -cache --dir_cache /SeqIdx/vepdb/caches --offline --force_overwrite")

scp xubr@moab.ncifcrf.gov:/scratch/cluster_tmp/xubr/*RPZ2*RPZ113*passed_AFDP_filtered_vep.txt ./

os.system("vep -i " + inputfile + "-o " + outputfile + "--offline --cache --dir_cache $VEP_CACHEDIR --species mouse --fasta $VEP_CACHEDIR/mouse.fa --everything --af_exac --vcf")
""" 

