
library(tidyverse)
#install.packages("officer")
library(officer)  ## for report writing
#install.packages("flextable")
library(flextable)
library(plyr)     ## for join function
#install.packages("readxl")
library("readxl")

############# PCGR filter
# input *.snvs_indels.tiers.tsv file
# select Tier1,2,3

data_files <- list.files("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/PCGR_out/", pattern = "*.snvs_indels.tiers.tsv")  
for (ff in data_files){
  print(ff)
  base <- str_split(ff, ".snvs" )
  print(base[[1]][1])
  maf <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/PCGR_out/",ff), header = TRUE, sep = "\t")
  print(colnames(maf))
  
  
  # PCGR rely on CLINVAR_CLNSIG
  narrow_maf <-  maf %>% select(`SYMBOL`,`CHROM`, `POS`, `REF`, `ALT`, `DBSNPRSID`, `GLOBAL_AF_1KG`, `CLINVAR_CLNSIG`, `CONSEQUENCE`,`PROTEIN_CHANGE`,`TIER`,`TIER_DESCRIPTION`)
  colnames(narrow_maf)
  short_maf <- narrow_maf %>%
    filter(`TIER`== "TIER 1" |`TIER`=="TIER 2"|`TIER`=="TIER 3") 

  write.table(short_maf, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/PCGR_out/",base[[1]][1], "_filter.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  
}


############# cpsr filter
# input *.snvs_indels.tiers.tsv file
# FINAL_CLASSIFICATION is "Pathogenic" or `Likely_Pathogenic`, or "VUS", class5,4,3

data_files <- list.files("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/CPSR_out/", pattern = "*.snvs_indels.tiers.tsv")  
for (ff in data_files){
  print(ff)
  base <- str_split(ff, ".snvs" )
  print(base[[1]][1])
  maf <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/CPSR_out/",ff), header = TRUE, sep = "\t")
  print(colnames(maf))
  
  # CPSR rely on CLINVAR_CLASSIFICATION
  narrow_maf <-  maf %>% select(`SYMBOL`,`VAR_ID`, `DBSNP`, `GENOTYPE`, `CLINVAR_CLASSIFICATION`, `CONSEQUENCE`,`PROTEIN_CHANGE`,`FINAL_CLASSIFICATION`)
  colnames(narrow_maf)
  short_maf <- narrow_maf %>%
    filter(`FINAL_CLASSIFICATION`=="Pathogenic"|`FINAL_CLASSIFICATION`=="VUS" | `FINAL_CLASSIFICATION`=='Likely_Pathogenic') 
  
  ID <- str_split_fixed(short_maf$`VAR_ID`, "_",4)
  colnames(ID) <- c("CHROM", "POS",	"REF","ALT")
  germline_var <- cbind(short_maf,ID)
  write.table(germline_var, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/CPSR_out/",base[[1]][1], "_filter.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
}




############# maf filter
# input *.maf
# relexed_filter selected PASS, High or moderate variants, calculated AF

data_files <- list.files("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/maf/", pattern = "*.maf")  
for (ff in data_files){
  print(ff)
  base <- str_split(ff, ".maf" )
  print(base[[1]][1])
  maf <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/maf/",ff), skip=1, header = TRUE, sep = "\t")
  print(colnames(maf))
  narrow_maf <-  maf %>% select(c(`Hugo_Symbol`,`Chromosome`, `Start_Position`, `Variant_Classification`,`Tumor_Sample_Barcode`,`HGVSp`,`t_depth`,`t_ref_count`,`t_alt_count`,`SIFT`,`IMPACT`,`VARIANT_CLASS`,`FILTER`))%>% 
                        mutate(AF=round((`t_alt_count`/`t_depth`), digits = 2)) 
  narrow_maf$`HGVSp`=substring(narrow_maf$`HGVSp`, first=3) 
  
  short_maf <- narrow_maf %>%
    filter(`IMPACT`=="HIGH"|`IMPACT`=="MODERATE") %>%
    #filter(!(grepl('tolerate',`SIFT`))) %>%
    filter(`FILTER`=="PASS")%>%
    #filter(`AF`>0.2 & t_depth > 30) %>%
    mutate(HGVSp=if_else(HGVSp=='',"Splicing variant", HGVSp)) %>%
    mutate(CLASS=case_when(AF > 0.9 | AF == 0.9 ~ "I-A",
                           AF > 0.2 & AF < 0.9 & `IMPACT`=="HIGH" ~ "I-B",
                           AF > 0.2 & AF < 0.9 & `IMPACT`=="MODERATE" ~ "II-A")) %>%
    mutate(CLASS=if_else(str_detect(HGVSp,"(Ala[0-9]*Gly$)|(Gly[0-9]*Ala$)|(Ser[0-9]*Thr$)|(Thr[0-9]*Ser$)|(Ile[0-9]*Leu$)|(Leu[0-9]*Ile$)"),"II-B", CLASS)) %>%
    arrange(`IMPACT`,`Hugo_Symbol`,`Chromosome`, `Start_Position`)
  
  
  write.table(narrow_maf, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/maf/",base[[1]][1], "_no_filter.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  
  write.table(short_maf, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/maf/",base[[1]][1], "_relexed_filter.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  
}


######## cravat filter
setwd("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/cravat_out")



cravat <- read_excel("DanielSoppet_CS034244_12exome_041723_HTCJ7DRX2_tumor-only_cravat.xlsx", sheet = "Variant")

  cravat_narrow <- cravat %>% select(c(`...2`, `...3`, `...4`, `...5`,`...8`,`CIViC`,`...34` ,`Cancer Gene Census`, `ClinGen Gene`,`ClinVar`, `gnomAD3`,`...14`))
  colnames(cravat_narrow) <-c('Chrom', 'Position', 'Ref Base', 'Alt Base', 'Hugo', 'CIVIC Description', 'COSMIC Variant Count', 'CancerGene Census Driver class', 'ClinGen Gene Disease', 'ClinVar Significance', 'Global AF', 'Samples')


  write.table(cravat_narrow, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/cravat_out/", "cravat_combined_filter.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  


######## merge
setwd("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge")

data_files <- list.files("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/", pattern = "*.cpsr.grch38_filter.txt")  
for (ff in data_files){
  print(ff)
  base <- str_split(ff, ".cpsr.grch38_filter.txt" )
  print(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".cpsr.grch38_filter.txt"))
  
  ## read in data
  cpsr_class345 <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".cpsr.grch38_filter.txt"),header = TRUE, sep = "\t")
  pcgr_Tie123 <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".pcgr_acmg.grch38_filter.txt"),header = TRUE, sep = "\t")
  maf_no_filter <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".hard-filtered.vep_no_filter.txt"), header = TRUE, sep = "\t")
  maf_filtered <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".hard-filtered.vep_relexed_filter.txt"), header = TRUE, sep = "\t")
  cravat_12col <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/", "cravat_combined_filter.txt"), header = TRUE, sep = "\t")
  
  # get the same format of ID
  cpsr_class345 <- cpsr_class345 %>% 
    mutate(Pos1=ifelse(str_count(REF)>1, as.double(`POS`)+1, `POS`))%>%
    mutate(ID=str_c("chr",`CHROM`,"_",`Pos1`))
  pcgr_Tie123 <- pcgr_Tie123 %>% 
    mutate(Pos1=ifelse(str_count(REF)>1, as.double(`POS`)+1, `POS`))%>%
    mutate(ID=str_c("chr",`CHROM`,"_",`Pos1`))
  

  maf_no_filter <- maf_no_filter %>% mutate(ID=str_c(`Chromosome`,"_",`Start_Position`))
  maf_filtered <- maf_filtered %>% mutate(ID=str_c(`Chromosome`,"_",`Start_Position`))
  cravat_narrow_anno <-cravat_12col %>% mutate(ID=str_c(`Chrom`,"_",`Position`)) %>% select(-c(`Chrom`, `Position`, `Ref.Base`, `Alt.Base`, `Hugo`, `Samples`)) 
  cravat_narrow_anno <- cravat_narrow_anno[-c(1),]
  
  
  
  #  join
  
  pcgr_cpsr <- pcgr_Tie123 %>% full_join(cpsr_class345, by="ID")
  
  pcgr_cpsr_anno <- pcgr_cpsr %>% left_join(maf_no_filter, by="ID") %>%
    filter(AF>0.1)
  pcgr_cpsr_maf <- pcgr_cpsr_anno %>% full_join(maf_filtered, by="ID") 
    
  


pcgr_cpsr_maf_cravat <- pcgr_cpsr_maf %>% left_join(cravat_narrow_anno, by="ID") %>% 
    mutate(Gene=ifelse(is.na(Hugo_Symbol.y), `Hugo_Symbol.x`, `Hugo_Symbol.y`)) %>% 
    mutate(Chrome=ifelse(is.na(Hugo_Symbol.y), `Chromosome.x`, `Chromosome.y`)) %>% 
    mutate(Pos=ifelse(is.na(Hugo_Symbol.y), `Start_Position.x`, `Start_Position.y`)) %>% 
    mutate(Hgvsp=ifelse(is.na(Hugo_Symbol.y), `HGVSp.x`, `HGVSp.y` )) %>% 
    mutate(IMPACT=ifelse(is.na(Hugo_Symbol.y), `IMPACT.x`, `IMPACT.y` )) %>%  
    mutate(AF=ifelse(is.na(Hugo_Symbol.y), `AF.x`, `AF.y` )) %>%  
    mutate(Consequence=ifelse(is.na(Hugo_Symbol.y), `Variant_Classification.x`, `Variant_Classification.y`)) %>%
    mutate(SIFT=ifelse(is.na(Hugo_Symbol.y), `SIFT.x`, `SIFT.y`)) %>%
    mutate(Class_cpsr_pcgr_SIFT=paste0(`FINAL_CLASSIFICATION`,',',`TIER_DESCRIPTION`,',', `SIFT`)) %>%
    mutate(Final_class=case_when(str_detect(`Class_cpsr_pcgr_SIFT`,"Pathogenic") ~ "class1",
                                 str_detect(`Class_cpsr_pcgr_SIFT`,'significance') & str_detect(IMPACT,"HIGH")~ "class2",
                                 str_detect(`Class_cpsr_pcgr_SIFT`,'significance') & str_detect(Class_cpsr_pcgr_SIFT,"deleterious") ~ "class3",
                                 str_detect(`Class_cpsr_pcgr_SIFT`,'deleterious')|str_detect(`IMPACT`,'HIGH') ~ "class4",
                                 str_detect(`Class_cpsr_pcgr_SIFT`,'significance')|str_detect(`Class_cpsr_pcgr_SIFT`,'VUS') ~ 'class5',
                                 TRUE~"class6")) %>%
    arrange(`Final_class`,`IMPACT`, `Gene`) 
  
  col_8 <- pcgr_cpsr_maf_cravat %>% select(c(`Gene`,`Chrome`,`Pos`,`Hgvsp`,`AF`,`IMPACT`,`Final_class`,`Consequence`,`Class_cpsr_pcgr_SIFT`))
   col_16<- pcgr_cpsr_maf_cravat %>% select(c(`Gene`,`Chrome`,`Pos`,`Hgvsp`,`AF`,`IMPACT`,`Final_class`,`Global.AF`,`Consequence`,`Class_cpsr_pcgr_SIFT`,`CIVIC.Description`,`COSMIC.Variant.Count`,`CancerGene.Census.Driver.class`,`ClinGen.Gene.Disease`,`ClinVar.Significance`,`SIFT`))
   

  write.table(col_8, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".col_8.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  
  write.table(col_16, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".col_16.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
  } 
  
########## merge with Depmap
  
setwd("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge")

data_files <- list.files("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/", pattern = "*.cpsr.grch38_filter.txt")  
for (ff in data_files){
  print(ff)
  base <- str_split(ff, ".cpsr.grch38_filter.txt" )
  print(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".cpsr.grch38_filter.txt"))
  
  ## read in data
  depmap <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],"_depMap_mutations.csv"),header = TRUE, sep = ",")
  colnames(depmap)
  col16 <- read.delim(paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".col_16.txt"),header = TRUE, sep = "\t")
  colnames(col16)
  
  depmap_tem <- depmap %>% 
    mutate(Pos1=ifelse(str_count(Ref.Allele)>1, as.double(`Position`)+1, `Position`))%>%
    mutate(ID=str_c(`Chromosome`,"_",`Pos1`))
  col16_tem <-col16 %>% mutate(ID=str_c(`Chrome`,"_",`Pos`))
  print("number of variants in our List: ")
  print(nrow(col16_tem))
  
  # add column In_depmap
  col16_depmap <- col16_tem %>% left_join(depmap_tem, by="ID") %>%
    mutate(In_depMap=ifelse(is.na(`Gene.y`),"NO","Yes")) %>%
    select(c(`Gene.x`,`Chrome`,`Pos`,`Hgvsp`,`AF`,`IMPACT`,`Final_class`,`Global.AF`,`Consequence`,`Class_cpsr_pcgr_SIFT`,`CIVIC.Description`,`COSMIC.Variant.Count`,`CancerGene.Census.Driver.class`,`ClinGen.Gene.Disease`,`ClinVar.Significance`,`SIFT`,`Cscape.Score`,`DANN.Score`,`In_depMap`))
  
  #number of depmap
  print("number of variants in depmap: ")
  print(nrow(depmap_tem))
  depmap_add <- depmap_tem %>% left_join(col16_tem, by="ID") %>%
    filter(!(str_detect(`Filter`, "PASS")))
  print("number of variants pass filter: ")
  print(nrow(depmap_add))
  
  # number not in my list
  depmap_na <- depmap_add[is.na(depmap_add$`Gene.y`),]
  print("number of variants not in our list: ")
  print(nrow(depmap_na))
  
  # number has DP>30
  depmap_dp <- depmap_na %>%
    mutate(DP= as.double(`Ref.Count`)+as.double(`Alt.Count`)) %>%
    filter(DP>30)
  print("number of variants having DP>30: ")
  print(nrow(depmap_dp))
  
  #number has AF>0.3
  depmap_af <- depmap_dp %>% 
    filter(`Allele.Fraction` > 0.3) 
  print("number of variants having AF>0.3: ")
  print(nrow(depmap_af))
  
  # number of having coding changes
  depmap_noSilent <- depmap_af %>% 
    filter(!(str_detect(`Variant.Info`, "SILENT")))
  print("number of variants having coding changes: ")
  print(nrow(depmap_noSilent))
  
  depmap_anno<- depmap_noSilent %>% left_join(maf_no_filter, by="ID") 

  write.table(col16_depmap, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".col_16_depmap.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  write.table(depmap_anno, paste0("/Users/xubr/local_projects/20221223_3XFLAG_SW_G12V/Parental_tumor_only/merge/",base[[1]][1],".depmap_no_myList.txt"), sep = "\t",row.names = FALSE, col.names = TRUE)
  
  
  
}  
  
  