library(tidyverse)
#install.packages("officer")
library(officer)
#install.packages("flextable")
library(flextable)
library(dplyr)
setwd("/Users/xubr/local_projects/20230428_9MEF/dragen_combined")
data_files <- list.files("/Users/xubr/local_projects/20230428_9MEF/dragen_combined/", pattern = "*_annotated_additional_Final_fixed.txt")  

for (ff in data_files) {
  base <- str_split(ff, "_somatic" )
  print(base)
  
  maf_temp <- read.delim(paste0("/Users/xubr/local_projects/20230428_9MEF/dragen_combined/",ff), header = FALSE, sep = "\t") %>% select(-c(`V7`,`V8`,`V9`,`V10`,`V11`,`V12`,`V18`,`V19`,`V20`,`V21`,`V22`)) %>% filter(`V1`!="not")
  colnames(maf_temp) <-maf_temp[1,]
  
  maf_temp<-maf_temp[-c(1),] 
  maf_temp<-maf_temp %>% mutate(`VARIANT`=ifelse((nchar(maf_temp$`REF`)==1 & nchar(maf_temp$`ALT`)==1), "SNV", "INDEL"))
  colnames(maf_temp)
  # obtain sample name, file name, date
  sample_name <- base[[1]][1]
  print(sample_name)
  file_name=paste0(base[[1]][1],".docx",collapse ="")
  
  today <- Sys.Date()
  format(today, format="%B %d %Y")
  print(today)
  # open an empty docx file, write sample name
  doc = read_docx()
  prop1=fp_text(font.size = 12, bold = TRUE, font.family = "Times")
  title1=fpar(ftext(sample_name, prop=prop1))
  body_add_fpar(doc, title1)
  
  # write subtitle and date
  body_add_par(doc, "Whole Exome report ")
  body_add_par(doc, today)
  body_add_par(doc, " ")
  
  # write method summary
  body_add_par(doc, "Method Summary: ")
  body_add_par(doc, "DNA was extracted using Qiagen all prep protocol. gDNA libraries from 500ng DNA was prepared and exons were captured using Agilent SureselectXT mouse all exon protocol. Libraries were quantitated via Agilent TapeStation. Sequencing was performed on the Illumina NovaSeq6000. Raw FastQ files were mapped using the Dragen v3.10 protocols for mutation detection. Single Nucleotide variants and insertions or deletions were identified, and a subset of calls filtered by quality, depth of coverage, and allele frequency were compiled into a variant list. A subset of the compiled variants was reviewed manually using the Integrated Genome Viewer software from the Broad Institute. ")
  body_add_par(doc, " ")
  
  body_add_par(doc, "FastQ and BAM files are stored at smb://at-s-is2.ncifcrf.gov/ras-intl/static/Genomics/20230428_9diploidMEF/")
  body_add_par(doc, " ")
  
  # insert a table
  body_add_par(doc, "Summary of findings: ")
  count_impact <- maf_temp %>% dplyr::count(IMPACT,VARIANT) 
  
  
  f.table=qflextable(count_impact)
  f.table=font(f.table,  fontname = "Times", part = "all")
  f.table=fontsize(f.table, size = 10, part = "all")
  doc <- flextable::body_add_flextable(doc, value = f.table, align = "left" )
  body_add_par(doc, " ")
  
  
  
  
  # insert gene list
  body_add_par(doc, "Gene list: ")
  maf1 <- maf_temp %>% select(c(`GENE`,`#CHROM`,`POS`,`AF`,`IMPACT`,`AMINO ACID CHANGE`,`CLASS`))
  f.table=qflextable(maf1)
  set_table_properties(f.table, width = .5, layout = "autofit")
  f.table=font(f.table,  fontname = "Times", part = "all")
  f.table=fontsize(f.table, size = 10, part = "all")
  # also set the table's header font as bold
  f.table=bold(f.table, part = "header")
  doc <- flextable::body_add_flextable(doc, 
                                       value = f.table, 
                                       align = "center" )
  body_add_par(doc, " ")
  body_add_par(doc, " ")
  body_add_par(doc, "I-A    Mutation call is well supported and there is a high probability that the variant will impact RAS dependent proliferation in this cell (e.g. Homozygous mutation of Trp53).")
  body_add_par(doc, " ")
  body_add_par(doc, "Variants called by one of variant calling methods and visualized by IGV. Single variant on a gene, AF>0.9 with high impact. Or multiple variants on a gene, 0.9>AF>0.2 with high impact.")
  body_add_par(doc, " ")
  body_add_par(doc, "Variants called by one of variant calling methods and visualized by IGV. Single variant on a gene, AF>0.9 with moderate impact. Or multiple variants on a gene, 0.9>AF>0.2 with moderate impact.")
  body_add_par(doc, " ")
  body_add_par(doc, "I-B    Mutation call is well supported, the variant is probably damaging, the variant may impact cell growth or genetic stability.")
  body_add_par(doc, " ")
  body_add_par(doc, "Variants called by one of variant calling methods and visualized by IGV. Single variant on a gene, 0.9>AF>0.2 with high impact. ")
  body_add_par(doc, " ")
  body_add_par(doc, " II-A    Call is well supported, variant is non-synonymous and likely impact protein function.") 
  body_add_par(doc, " ")               
  body_add_par(doc, "Variants called by one of variant calling methods and visualized by IGV. Single variant on a gene, 0.9>AF>0.2 with moderate impact.") 
  body_add_par(doc, " ")
  body_add_par(doc, "II-B    Call is well supported, non-synonymous change in amino acids with minor change in amino acid properties (e.g.  Serine to Threonine, Ile to Leu, or Glycine to Alanine) that may alter protein function.")
  body_add_par(doc, " ")
  body_add_par(doc, "Variants called by one of variant calling methods and visualized by IGV. Single variant on a gene, AF>0.2 with moderate impact and with minor change in amino acid properties (e.g.  Serine to Threonine, Ile to Leu, or Glycine to Alanine). ")
  body_add_par(doc, " ")
  body_add_par(doc, "II-C    Mutation call is well supported but the change is unlikely to impact protein function.")
  
  
  # output docx file
  print(doc, target=file_name)
  
}
  