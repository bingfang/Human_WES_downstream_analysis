sinteractive --mem=49g --cpus-per-task=32

./run.sh --sourcefq "/data/xubr/fastq_RPZ" --outdir "/data/xubr/exome_RPZ/output" --targets "/data/xubr/S0276129_ID_Regions_mm10.bed" --pairs "/data/xubr/exome_RPZ/pipeline/pairs" --cnv 'true' --dryrun 1


## New folder where we want to run it
NEW_FOLDER="/data/$USER/exome_XXX"

## Create the new folder
if [ ! -d $NEW_FOLDER ]; then mkdir -p $NEW_FOLDER; fi

## Create a subdirectory to store the pipeline in the new folder
cd $NEW_FOLDER
mkdir pipeline

## Copy over the pipeline skeleton
cp -r $PIPE_DIR/* pipeline/

cd pipeline

## Can run this to see all options for the run command
./run.sh --help


## Set up variables for input params for a run
# Folder of fastq.gz files
fq_dir="/data/tandonm/pl_test_data/human/fastq"

# Location of subdirectory to hold pipeline output
output_dir="$NEW_FOLDER/output"

# Location of targets bed file (hg38)
# This is the default file used in the config file, so no need to specify unless it's different
bed_file="/data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed"

# Pairs file; must contain column header "Tumor" and "Normal" (in any order) listing sample IDs
pairs_file="/data/tandonm/pl_test_data/human/pairs"


## Try a dry run to see if it all rules are compiled successfully
./run.sh --sourcefq "$fq_dir" --outdir "$output_dir" --targets "$bed_file" --pairs "$pairs_file" --ffpe "True" --dryrun 1

## Submit to the cluster by excluding the --dryrun flag
./run.sh --sourcefq "$fq_dir" --outdir "$output_dir" --targets "$bed_file" --pairs "$pairs_file" --ffpe "True"

Normal	Tumor
12_RPZ0113	1_RPZ28754_D001pt2_S1
12_RPZ0113	2_RPZ28755_D001pt2_S2
12_RPZ0113	3_RPZ28699_F18_S3
12_RPZ0113	4_RPZ28756_D001pt2_S4
12_RPZ0113	5_RPZ28701_F20_S5
12_RPZ0113	6_RPZ287011_K14_S6
12_RPZ0113	7_RPZ28703_B23_S7
12_RPZ0113	8_RPZ28703_F23_S8
12_RPZ0113	9_RPZ28705_O19_S9
12_RPZ0113	10_RPZ26216_D007pt2_S10
12_RPZ0113	11_RPZ26216_D006pt2_S11
12_RPZ0113	12_RPZ28705_L12_S12
cd 