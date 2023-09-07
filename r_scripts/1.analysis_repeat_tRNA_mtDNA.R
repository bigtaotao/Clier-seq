if(T){
  rm(list = ls())
  library(dplyr)
  library(yaml)
  library(data.table)
  library(DESeq2)
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  source("./r_scripts/stringtie_func.R")
}



if(T){
  config = read_yaml("./config.yaml")
  root_data_path= file.path(config$output$root,"repeat_tRNA")
  input_file = file.path(config$input$root,"fastq_list.csv")
  
  ## make output dir
  dir.create(root_data_path,recursive = T)
  input_file = read.csv(input_file)
  ## load
  get_data = function(basenames){
    repeat_res = apply(input_file,MARGIN = 1,function(each_sample){
      sample_name = each_sample["name"]
      print(sample_name)
      file = file.path(config$output$root,each_sample["sample"],basenames)
      data = as.data.frame(fread(file))
      data = data[-1,]
      colnames(data) = c("class","val")
      data$sample = sample_name
      data$group = each_sample["group"]
      data$type = each_sample["type"]
      return(data)
    })
    repeat_res = do.call(rbind,repeat_res)
    return(repeat_res)
  }
  
  
  get_mat_col= function(res){
    mat = dplyr::select(res,c("sample","class","val"))
    mat = dcast(mat,class ~ sample)
    mat[is.na(mat)] = 0
    rownames(mat) = mat$class
    mat = dplyr::select(mat,-"class")
    
    
    condition = dplyr::select(res,c("sample","group","type"))
    condition = unique(condition)
    rownames(condition) = condition$sample
    condition = dplyr::select(condition, -c("sample"))
    condition$type = factor(condition$type,levels = c("input","enrich"))
    
    condition = condition[match(colnames(mat),rownames(condition)),]
    return(list(mat,condition))
  }
  
  run_deseq = function(mat,cond){
    # 构建DESeqDataSet对象
    dds <- DESeqDataSetFromMatrix(
      countData = mat,
      colData = cond,
      design = ~ group + type) # 注意这行代码
    
    dds$type<- relevel(dds$type, ref = "input") 
    
    # 差异检验
    dds = DESeq(dds)
    # 提取差异分析结果
    nrDEG_DESeq2 <- as.data.frame(results(dds))
    
  }
  
  
  
}


repeat_data = get_data("map_repeat_rmdup_anno_count.bed")
trna_data = get_data("map_tRNA_rmdup_anno_count.bed")
mt_data = get_data("map_mtDNA_rmdup_anno_count.bed")
comb_data = rbind(repeat_data,trna_data,mt_data)

comb_data = filter(comb_data,grepl("hr1",sample))
# debugonce(get_mat_col)
comb_mat = get_mat_col(comb_data)
comb_deseq = run_deseq(comb_mat[[1]],comb_mat[[2]])

## output 
write.csv(comb_deseq,file = file.path(root_data_path,"deseq2_res.csv"))
write.csv(comb_mat[[1]],file =file.path(root_data_path,"combined_data.csv"))
write.csv(comb_mat[[2]],file = file.path(root_data_path,"condition.csv"))

comb_deseq$genes = rownames(comb_deseq)
plot_data = prepare_plot_valcano(comb_deseq,colname = c("log2FoldChange","pvalue","genes"))
pdf(file.path(root_data_path,"volcano.pdf"),width = 8,height =8)
print(plot_volcano(plot_data))
dev.off()
                                 