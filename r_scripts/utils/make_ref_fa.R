## generate small nuclear rna reference
if(T){
  rm(list = ls())
  library(data.table)
  library(dplyr)
  library(stringr)
  gtf = fread("./ref/GRCh38/gencode.v42.chr_patch_hapl_scaff.annotation.gtf")
  gtf = as.data.frame(gtf)
  gtf = dplyr::filter(gtf,V3 == "gene")
  ensg = str_split(gtf$V9,";",simplify = T)
  gtf$ensg = str_sub(ensg[,1],10,-4)
  
  urna = fread("./ref/GRCh38/small nuclear RNAs.txt")
  urna_gtf = filter(gtf,ensg %in% urna$`Ensembl gene ID`)
  urna_gtf = dplyr::select(urna_gtf,-c("ensg"))
  
  write.table(urna_gtf,file = "./ref/GRCh38/small_nuclear_RNAs.gtf",sep = "\t",quote = F,row.names = F)
  
  ## 
  # RUN 
  # cat small_nuclear_RNAs.gtf |convert2bed -i gtf > small_nuclear_RNAs.bed
  
  bed = fread("./ref/GRCh38/small_nuclear_RNAs.bed")
  gene_name = str_split(bed$V10,";",simplify = T)
  gene_name = gene_name[,3]
  gene_name = str_sub(gene_name,13,-2)
  rgene_name = paste0("Repeat_regions:",gene_name)
  
  bed$V4 = rgene_name
  
  write.table(bed,"./ref/GRCh38/small_nuclear_RNAs.bed",quote = F,row.names = F,sep = "\t",col.names = F)
  
  # make new reapeat bed
  
  len = bed$V3 - bed$V2
  new_bed = data.frame(bed$V4,
                       0,
                       len,
                       gene_name,
                       0,
                       "+")
  write.table(new_bed,"./ref/GRCh38/small_nuclear_RNAs_append.bed",quote = F,row.names = F,sep = "\t",col.names = F)
  ## 
  # RUN
  # bedtools getfasta -fi GRCh38.fa -bed small_nuclear_RNAs.bed -fo small_nuclear_RNAs.fa -nameOnly
}
