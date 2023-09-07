library("rtracklayer")
get_gene_pos = function(combined_gtf){
  
  combined_gtf = import(combined_gtf)
  combined_gtf = as.data.frame(combined_gtf)
  combined_gtf = filter(combined_gtf, type == "transcript") %>% 
    mutate(c = if_else(is.na(gene_name),ref_gene_id,gene_name)) 
  combined_gtf$tg = sub("\\.[^.]*$", "", combined_gtf$transcript_id)
  combined_gtf = combined_gtf %>% mutate(known_gene = if_else(is.na(c),tg,c)) %>%
    dplyr::select(c("gene_id","known_gene","transcript_id","seqnames","start","end","width","strand"))
  colnames(combined_gtf) = c("gene","known_gene","transcript_id","chr","start","end","width","strand")
  
  ## because transcipt_id is unique skip rm dup func
  # combined_gtf = combined_gtf[order(combined_gtf$transcript_id),]
  # unique_gene = unique(combined_gtf$transcript_id)
  # sfInit(parallel = TRUE, cpus = detectCores() - 1)
  # sfLibrary(dplyr)
  # sfExport("combined_gtf")
  # gene_data = sfLapply(unique_gene,function(each_gene){
  #   c = filter(combined_gtf,transcript_id == each_gene)
  #   gene = c$gene[1]
  #   start = min(c$start)
  #   end = max(c$end)
  #   chr = c$chr[1]
  #   known_gene = unique(c$known_gene)[1]
  #   transcript_id = each_gene
  #   length = end - start
  #   return(data.frame(gene = gene,chr = chr,start = start, end = end,known_gene = known_gene,transcript_id =transcript_id,length = length ))
  # })
  # sfStop()
  # gene_data = do.call(rbind,gene_data)
  gene_data = combined_gtf
  return(gene_data)
}

add_position = function(bg_res,gene_data){
  c = left_join(bg_res,gene_data,)
  c = mutate(c,igv = paste0(chr,":",start))
  return(c)
  
}
plot_volcano = function(plot_data,cut_off_logFC = 2,cut_off_pvalue = 0.05){
  p = ggplot(
    # 数据、映射、颜色
    # limma_res, aes(x = log2FoldChange, y = -log10(pvalue))) +
    plot_data, aes(x = fc, y = -log10(pval))) + 
    # geom_point()+
    geom_point(aes(color = change), size=2) +
    #scale_color_manual(values = c("#375284","grey", "#bc311c")) +
    scale_color_manual(values = c("grey", "#bc311c")) +
    
    # 辅助线
    geom_vline(xintercept=c(-cut_off_logFC, cut_off_logFC),lty=4,col="#666666",lwd=0.8) +
    geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="#666666",lwd=0.8) +
    # 注释
    geom_label_repel(
      data = subset(plot_data, pval < cut_off_pvalue & plot_data$fc >= cut_off_logFC),
      aes(label = gene),
      size = 5,fill = "darkred", color = "white",
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (p-value)",) +
    # 图例
    theme_bw()+
    theme(legend.position = "bottom",axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16))
  
}
prepare_plot_valcano = function(res,cut_off_logFC = 2,cut_off_pvalue = 0.05,colname = c("all_fc","pval","known_gene")){
  data = dplyr::select(res,colname)
  colnames(data) = c("fc","pval","gene")
  
  data$change = ifelse(data$pval < cut_off_pvalue & abs(data$fc) >= cut_off_logFC, 
                       ifelse(data$fc > cut_off_logFC ,'enrich',' '),
                       ' ')
  return(data)
}
