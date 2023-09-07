if(T){
  rm(list = ls())
  library(ballgown)
  library(dplyr)
  library(RSkittleBrewer) 
  library(genefilter)
  library(devtools)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(stringr)
  source("./r_scripts/stringtie_func.R")
  library(snowfall)
  library(parallel)
  library(circlize)
  library(yaml)
  library(rtracklayer)
  # 用连接ensembl数据库作为例子
  library("biomaRt")
  run_ballgown_all = function(data_path,input_file,gene_data){
    ## run ballgown firstly to get dir order -> dir
    bg = ballgown(dataDir=data_path, samplePattern='-')
    dir = bg@dirs
    
    ## load meta data
    input_data = read.csv(input_file)
    # pdata = mutate(input_data,ids = paste0("sample",group,"_",type))
    pdata = dplyr::select(input_data,c("sample","type","group"))
    
    pdata$group = as.character(pdata$group)
    pdata = pdata[match(names(dir),pdata$sample),]
    
    pdata$type = factor(pdata$type,levels = c("input","enrich"))
    
    ## run ballgown secondly to get data
    bg = ballgown(dataDir=data_path, samplePattern='-',pData = pdata)
    # use func get_TPM else
    # bg_fpkm = dplyr::select(bg@expr$trans,-c(contains("cov")))
    # write.csv(bg_fpkm,file.path(data_path,"all_sample_fpkm.csv"))
    bg <<- bg
    bg_filter = subset(bg,"rowVars(texpr(bg)) >10",genomesubset=TRUE)
    bg <<- NA
    results_transcripts = stattest(bg_filter, 
                                   feature="transcript",covariate="type",adjustvars = 
                                     c("group"), getFC=TRUE, meas="FPKM")
    results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filter), 
                                     geneIDs=ballgown::geneIDs(bg_filter),
                                     transcript_id = ballgown::transcriptNames(bg_filter),
                                     results_transcripts)
    res = mutate(results_transcripts ,log2fc = log2(fc)) %>% 
      dplyr::select(c("transcript_id","log2fc","pval","qval"))
    
    colnames(res) = c("transcript_id","all_fc","pval","qval")
    # results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filter), 
    #                                  geneIDs=ballgown::geneIDs(bg_filter), results_transcripts)
    # 
    # results_transcripts = arrange(results_transcripts,fc,pval)
    
    res = add_position(res,gene_data)
    return(res)
    # cut_off_logFC = 2
    # cut_off_pvalue = 0.05
    # results_transcripts = mutate(results_transcripts,gene = if_else(geneNames == ".",geneIDs,geneNames))
    # results_transcripts[,"logfc"] = log2(results_transcripts[,"fc"])
    # results_transcripts$change = ifelse(results_transcripts$pval < cut_off_pvalue & abs(results_transcripts$logfc) >= cut_off_logFC, 
    #                                     ifelse(results_transcripts$fc> cut_off_logFC ,'enrich',' '),
    #                                     ' ')
  }
  
  plot_heatmap = function(){
    
  }
  get_TPM = function(data_path,input_file){
    input_data = read.csv(input_file)
    # sfInit(parallel = TRUE, cpus = min(detectCores() - 1,nrow(input_data)))
    sfInit(parallel = TRUE, cpus = 8)
    sfLibrary(dplyr)
    sfLibrary(rtracklayer)
    sfExport("data_path")
    res = sfLapply(input_data$sample,function(each){
      gtf_file = file.path(data_path,each,"map_genome_rmdup_mapmerged.gtf")
      print(gtf_file)
      gtf_data = rtracklayer::import(gtf_file)
      gtf_data = as.data.frame(gtf_data) %>% 
        filter(type == "transcript") %>%
        select(c("seqnames","start","end","width","transcript_id","TPM"))
      colnames(gtf_data) = c("chr","start","end","width","transcript_id",each)
      return(gtf_data)
    })
    sfStop()
    res = Reduce(left_join,res)
    return(res)
  }
  annote_description = function(gene_data){
    my_ensembl_gene_id = gene_data$known_gene
    my_ensembl_gene_id = my_ensembl_gene_id[grepl("ENS",my_ensembl_gene_id)]
    my_ensembl_gene_id = unique(my_ensembl_gene_id)
    mart <- useMart(biomart="ensembl")
    
    #查看有什么dataset可以选择
    # listDatasets(mart)
    
    #选择一个人的基因注释，加上host连接，不加host需要选择国外mirror
    mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl",host="https://www.ensembl.org")
    # a = listAttributes(mart)
    # getBM 获取注释
    description<-getBM(attributes=c('ensembl_gene_id',"description"), filters= 'ensembl_gene_id', values = my_ensembl_gene_id, mart)
    annote_gene_data = left_join(gene_data,description,by = c("known_gene"="ensembl_gene_id"))
    return(annote_gene_data)
  }
}
# debugonce(run_ballgown_all)

config = read_yaml("./config.yaml")
root_data_path= file.path(config$output$root,"stringtie")
input_file = file.path(config$input$root,"fastq_list.csv")
pipe_list = config$config$stringtie$pipeline

for (each_pipe in pipe_list){
  data_path = file.path(root_data_path,each_pipe)
  if (config$pipeline_config$pipe[[each_pipe]]$stringtie != ""){
    combined_gtf = file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$stringtie) 
  }else{
    combined_gtf = file.path(data_path,"stringtie_merged.gtf")
    
  }
  ## runing ballgown
  cache_gene_data_path = file.path(data_path,"gene_data.rdata")
  if(file.exists(cache_gene_data_path)){
    load(cache_gene_data_path)
    print("loading Rdata for gene_data")
  }else{
    print("construct gene_data")
    gene_data = get_gene_pos(combined_gtf)
    gene_data = annote_description(gene_data)
    save(gene_data,file = cache_gene_data_path)
  }

  
  tpm = get_TPM(data_path,input_file)
  tpm = left_join(gene_data,tpm,by = c("chr","start","end","transcript_id","width"))
  write.csv(tpm,file = file.path(data_path,"all_sample_tpm.csv"),quote = T,row.names = F)
  # debugonce(run_ballgown_all)
  res = run_ballgown_all(data_path,input_file,gene_data)
  # write.csv(res,file.path(data_path,"all_sample_diff_trans.csv"),row.names = F)
  
  #if need to merge TPM and RES files
  # res_tpm = left_join(res,tpm,by = join_by(transcript_id, gene, known_gene,
  #                                          chr, start, end, width))
  
  # enrich_input = lapply(c("enrich","input"),function(each){
  #   c = tpm %>% dplyr::select(contains(each))
  #   c = as.data.frame(lapply(c,as.numeric))
  #   return(rowMeans(c))
  # })
  # enrich_input = as.data.frame(enrich_input)
  # colnames(enrich_input) = c("enrich_tpm","input_tpm")
  # res_tpm = cbind(tpm,enrich_input) %>% dplyr::select(c("transcript_id","enrich_tpm","input_tpm"))
  # res_tpm = right_join(res_tpm,res)
  res_tpm = res
  write.csv(res_tpm,file.path(data_path,"all_sample_diff_trans.csv"),row.names = F)
  
  print("okk")
  ## ploting volcano
  plot_data = prepare_plot_valcano(res)
  p = plot_volcano(plot_data)
  ggsave(p,filename = file.path(data_path,"all_sample_volcano.pdf"),width = 20,height = 16)

  ## ploting heatmap
  data_content = tpm
  data = dplyr::select(data_content,-c("gene","known_gene","chr","start","end","width","strand","description"))
  # colnames(data) = data_content$gene_id
  data_enrich = data %>% dplyr::select(-c("transcript_id"))
  data_enrich=as.data.frame(lapply(data_enrich,as.numeric))
  data_enrich_rowsum = rowSums(data_enrich)
  select_index = data_enrich_rowsum > 10
  filter_data = data[select_index,]

  plot_data = filter_data %>% dplyr::select(-c("transcript_id"))
  plot_data_clust = hclust(dist(plot_data))


  meta_data = read.csv(input_file)
  plot_data_colname = colnames(plot_data)
  sample_name  = meta_data$sample
  col_order = match(sample_name,plot_data_colname)

  col_fun = colorRamp2(c(-5, 0, 5), c( "#0361b0","white","#e60027"))

  pdf(file.path(data_path,"all_sample_heatmap.pdf"),width = 16,height = 20)

  # ComplexHeatmap::Heatmap(plot_data, name = "logFC", col = col_fun,show_row_names = F,na_col = "#aaaaaa",row_order = plot_data_clust$order,column_order = col_order)

  ComplexHeatmap::Heatmap(plot_data, name = "logFC", col = col_fun,show_row_names = F,na_col = "#aaaaaa",row_order = plot_data_clust$order)
  dev.off()

  write.csv(filter_data,file.path(data_path,"all_sample_heatmap.csv"))
}
