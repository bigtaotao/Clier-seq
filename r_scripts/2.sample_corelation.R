if(T){
  rm(list = ls())
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(yaml)
  library(pheatmap)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(purrr)
  library(xlsx)
  source("./r_scripts/stringtie_func.R")
  if(T){
    ### filter data
    filter_by_sample_expr = function(tpm_data){
      # get zero mask
      mask_data = tpm_data
      mask_data[mask_data != 0] = 1
      select_rows = rowSums(mask_data)
      select_tpm_data = tpm_data[select_rows > 1,]
      # select_tpm_data = tpm_data[select_rows > ncol(tpm_data) - 1,]
      return(select_tpm_data)
    }
    filter_by_tpm_sum = function(tpm_data){
      tpm_genes = tpm_data$known_gene
      
      # select by gene names,exclude YRNA rRNA 7SK URNA 5sRNA 5.8sRNA MT-DNA
      arg = "Y_RNA|^RNY\\d+(?!.*P\\d*$).*$|^RN[UV7](?!.*P\\d*$).*$|^U\\d+$|^7SK$|^5.*rRNA$|^RNA5-8S(?!.*P\\d*$).*$|^RNA5S(?!.*P\\d*$).*$|^MT-"
      select_tpm_data = tpm_data[!grepl(arg,tpm_genes,perl = T),]
      del_data = tpm_data[grepl(arg,tpm_genes,perl = T),]
      # split data 
      # tpm_meta = select_tpm_data[,meta_col] #unused
      tpm_mat = select_tpm_data[,-meta_col]
      tpm_mat = as.data.frame(lapply(tpm_mat,as.numeric))
      
      # filter by at least two column not zero
      num_mat = tpm_mat
      num_mat[num_mat > 0 ] = 1
      num_mat[!num_mat == 1] = 0
      # filter by rowSum >10
      mask = rowSums(num_mat) > 1 & rowSums(tpm_mat) > 10
      select_tpm_data = select_tpm_data[mask,]
      
      
      ## debug 
      # a = select_tpm_data[rowSums(num_mat) > 1,]
      # "NT.41208" %in% a$gene
      # a[a$gene == "NT.41208",]
      # b = select_tpm_data[rowSums(tpm_mat) > 10,]
      # b[b$gene == "NT.41208",]
      # c = select_tpm_data[rowSums(num_mat) > 1 & rowSums(tpm_mat) > 10,]
      # c[c$gene == "NT.41208",]
      # select_tpm_data = na.omit(select_tpm_data) #DO NOT USE it at 'select_tpm_data' because of 'description' column
      return(list(select_tpm_data,del_data))
    }
    plot_pheatmap = function(select_tpm_data,file = "corelation_heatmap.pdf" ){
      
      # znn_data = as.data.frame(fread(file.path(data_path,"all_sample_tpm_10_3.csv")))
      # znn_data = znn_data[,-1]
      # sample_order = colnames(znn_data)
      # writeLines(text = sample_order,con = file.path(data_path,"sample_order.txt"))
      
      sample_order = readLines(file.path(data_path,"sample_order.txt"))
      cor_data = select_tpm_data[,-meta_col]
      cor_res = cor(cor_data)
      # r_squared = cor_res^2
      # 
      # hc = hclust(as.dist(1 - r_squared))
      # r_squared = r_squared[hc$order, hc$order]
      # r_squared[lower.tri(r_squared)] = NA
      # 
      ## use znn filter data
      
      cor_res = cor_res[,match(sample_order,colnames(cor_res))]
      cor_res = cor_res[match(colnames(cor_res),rownames(cor_res)),]
      cor_res[lower.tri(cor_res)] = NA
      pdf(file.path(data_path,file),width = 10,height = 10)
      print(pheatmap(cor_res,
               cluster_cols = F,
               cluster_rows = F,
               display_numbers = F,
               border_color = NA,
               na_col = "#ffffff",
               fontsize = 10, 
               angle_col = 90,
               color = colorRampPalette(colors = c("white","#fb9716"))(100)))
      dev.off()
      
    }
    
    # stat rp func
    get_rp_tpm = function(input_data){
      select_tpm_data = input_data
      if(F){
        rp_data = file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"ribosomal_proteins.csv")
        rp_data = as.data.frame(fread(rp_data))
        tpm_rp = filter(select_tpm_data,known_gene %in% rp_data$`Approved symbol`)
        tpm_remain = anti_join(select_tpm_data,tpm_rp)
      }
      # protein coding and pseudogenes
      if(T){
        tpm_rp_gene = select_tpm_data[grepl("(^M?RP[S|L])|(^FAU$)",select_tpm_data$known_gene,perl = T),]
        tpm_rp_ensg = select_tpm_data[grepl("ribosomal protein",select_tpm_data$description),]
        tpm_rp = full_join(tpm_rp_gene,tpm_rp_ensg)
        tpm_remain = anti_join(select_tpm_data,tpm_rp)
      }
      return(list(tpm_rp,tpm_remain))
    }
    
    # stat snoRNA
    get_snoRNA_tpm = function(input_data){
      sno_ensg = readLines(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"genes.snoRNA"))
      tpm_snoRNA_gene = input_data[grepl("(SNOR[AD])|(SCARNA)|(^sno)",input_data$known_gene),]
      tpm_snoRNA_ensg_1 = filter(input_data,known_gene %in% sno_ensg)
      tpm_snoRNA_ensg_2 = input_data[grepl("small Cajal body-specific RNA",input_data$description),]
      tpm_snoRNA = full_join(tpm_snoRNA_gene,tpm_snoRNA_ensg_1) %>% full_join(tpm_snoRNA_ensg_2)
      tpm_remain = anti_join(input_data,tpm_snoRNA)
      return(list(tpm_snoRNA,tpm_remain))
    }
    
    # stat snRNA
    # WARNING snRNA is removed by repeat.fa 
    get_snRNA_tpm = function(input_data){
      sn_ensg = readLines(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"genes.snRNA"))
      tpm_uRNA_gene = input_data[grepl("(RNV?U[0-9]+[A|F|(ATAC)]?)|RN7SK",input_data$known_gene),]
      tpm_uRNA_ensg = input_data %>% dplyr::filter(known_gene %in% sn_ensg)
      tpm_uRNA = full_join(tpm_uRNA_ensg,tpm_uRNA_gene)
      tpm_remain = anti_join(input_data,tpm_uRNA)
      return(list(tpm_uRNA,tpm_remain))
    }
    
    # stat new transcript
    get_newtrans_tpm = function(input_data){
      tpm_new  = input_data[grepl("MSTRG\\.|^NT\\.",input_data$known_gene),]
      tpm_remain = anti_join(input_data,tpm_new)
      return(list(tpm_new,tpm_remain))
    }
    
    # stat protein coding
    get_protein_tpm = function(input_data){
      protein_ensg = read.csv(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"my","genes.protein_coding"))
      tpm_prot_gene_1 = input_data %>% dplyr::filter(known_gene %in% protein_ensg$Gene.name)
      tpm_prot_gene_2 = input_data[grepl("H4-16|CBWD|KIAA0100",input_data$known_gene),]
      tpm_prot_gene = rbind(tpm_prot_gene_1,tpm_prot_gene_2)
      tpm_prot_ensg = input_data %>% dplyr::filter(known_gene %in% protein_ensg$Gene.stable.ID)
      
      znn_protein = readLines(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"znn","znn.protein"))
      tpm_prot_znn = dplyr::filter(input_data,known_gene %in% znn_protein)
      tpm_prot = full_join(tpm_prot_gene,tpm_prot_ensg) %>% full_join(tpm_prot_znn)
      tpm_remain = anti_join(input_data,tpm_prot)
      return(list(tpm_prot,tpm_remain))
    }
    
    # stat novel transcript and pseudo gene
    get_novel_tpm = function(input_data){
      tpm_novel = input_data[grepl("novel transcript|pseudogene",input_data$description),]
      tpm_pseudo = input_data[grepl("P[0-9]+$|^RN[UV7](.*P\\d*$).*$|ENSG00000283390|ENSG00000253683|ENSG00000240036|ENSG00000283041|ENSG00000283252|ZNF815P|FAM185BP",input_data$known_gene),] ## NOT GOOD,example ELP5
      
      znn_newtrans = readLines(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"znn","znn.newtrans"))
      znn_pseudo = readLines(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"znn","znn.pseudo"))
      tpm_novel_znn = dplyr::filter(input_data,known_gene %in% znn_newtrans)
      tpm_pseudo_znn = dplyr::filter(input_data,known_gene %in% znn_pseudo)
      tpm_novel_pseudo = full_join(tpm_novel,tpm_pseudo) %>% full_join(tpm_novel_znn) %>% full_join(tpm_pseudo_znn)
      tpm_remain = anti_join(input_data,tpm_novel_pseudo)
      return(list(tpm_novel_pseudo,tpm_remain))
    }
    
    # stat lncRNA miRNA and vaultRNA
    get_lncRNA_tpm = function(input_data){
      lnc_ensg = read.csv(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"my","genes.lncRNA"))
      mi_ensg = read.csv(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"my","genes.miRNA"))
      vt_ensg = read.csv(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"my","genes.vaultRNA"))
      znn_nc = readLines(file.path(config$pipeline_config$ref_root,config$pipeline_config$pipe$genome$types_file,"znn","znn.ncrna"))
      lnc_ensg = rbind(lnc_ensg,mi_ensg,vt_ensg)
      tpm_lnc_gene = dplyr::filter(input_data,known_gene %in% lnc_ensg$Gene.name)
      tpm_lnc_gene_grepl = input_data[grepl("^VTRNA|RMRP|RPPH1|BCYRN1|ENSG00000272101",input_data$known_gene),]
      tpm_lnc_gene_znn = dplyr::filter(input_data,known_gene %in% znn_nc)
      tpm_lnc_ensg = dplyr::filter(input_data,known_gene %in% lnc_ensg$Gene.stable.ID)
      tpm_lnc = full_join(tpm_lnc_gene,tpm_lnc_ensg) %>% full_join(tpm_lnc_gene_grepl) %>% full_join(tpm_lnc_gene_znn)
      tpm_remain = anti_join(input_data,tpm_lnc)
      return(list(tpm_lnc,tpm_remain))
    }
    
    # get summary of tpm dataframe
    tpm_summary =  function(input_data,type) {
      data = input_data[,-meta_col]
      # get nozero number
      data_num_mat  = data
      data_num_mat[data_num_mat > 0] = 1
      trans_num = colSums(data_num_mat)
      # get tpm count
      tpm_count = colSums(data)
      output = data.frame(sample_name= colnames(data),trans_num = trans_num,tpm_count = tpm_count,type = type)
      return(output)
    }
    
    gene_type_summary = function(input_data,type){
      res = data.frame(gene = input_data$known_gene,type = type)
      return(res)
    }
    stat_tpm = function(select_tpm_data){
      output_path = file.path(data_path,"tpm_summary")
      dir.create(output_path)
      
      #write_delete_data 
      classify_path = file.path(output_path,"classification")
      dir.create(classify_path)
      write.csv(select_tpm_data[[2]],file.path(classify_path,paste0("delete",".csv")))
      ## stat novel transcript
      novel_res = get_novel_tpm(select_tpm_data[[1]])
      tpm_novel = novel_res[[1]]
      
      ## stat RP 
      rp_res = get_rp_tpm(novel_res[[2]])
      tpm_rp = rp_res[[1]]
      
      ## stat snoRNA
      snoRNA_res = get_snoRNA_tpm(rp_res[[2]])
      tpm_snoRNA = snoRNA_res[[1]]
      
      # snRNA_res = get_snRNA_tpm(snoRNA_res[[2]])
      # tpm_snRNA = snRNA_res[[1]] #deprecated
      
      ## stat new transcript
      newt_res = get_newtrans_tpm(snoRNA_res[[2]])
      tpm_newt = newt_res[[1]]
      
      ## stat protein coding
      prot_res = get_protein_tpm(newt_res[[2]])
      tpm_prot = prot_res[[1]]
      
      ## stat lncRNA
      lnc_res = get_lncRNA_tpm(prot_res[[2]])
      tpm_lnc = lnc_res[[1]]
      
      ## trash
      tpm_others = lnc_res[[2]]
      
      tpm_list = list(tpm_rp,tpm_snoRNA,tpm_newt,tpm_prot,tpm_novel,tpm_lnc,tpm_others)
      tpm_type = c('RPS_RPL',"snoRNA","new transcript","protein coding","novel transcript(ENSEMBL)","ncRNA","Others")
      
      
      ## save to csv
      purrr::map2(tpm_list,tpm_type,function(data,path){
        write.csv(data,file.path(classify_path,paste0(path,".csv")))
      })
      
      
      ## point plot
      tpm_sum = purrr::map2(tpm_list,tpm_type,tpm_summary)
      tpm_sum = do.call(rbind,tpm_sum)
      tpm_sum = dplyr::filter(tpm_sum,type != "Others")
      lapply(unique(tpm_sum$sample_name),function(each_sample){
        plot_data = dplyr::filter(tpm_sum,sample_name == each_sample)
        plot_data$percentage = plot_data$tpm_count/sum(plot_data$tpm_count)
        plot_data$log_num = log10(plot_data$trans_num)
        plot_data$type = as.factor(plot_data$type)
        p = ggplot(plot_data,aes(x = type,y = log_num,size= tpm_count)) + 
          geom_point(color = "#ffc64b")+
          scale_y_continuous(breaks = c(1,2,3,4),
                             labels = c(10,100,1000,10000),limits = c(1,5))+
          scale_size_continuous(range=c(2,30),breaks = c(100000,200000,300000),labels = c(1e5,2e5,3e5),limits = c(0,1000000))+
          ylab("transcripts number")+
          xlab("")+
          guides(size=guide_legend(title="Total TPM"))+
          theme_bw()+
          theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
        pdf(file = file.path(output_path,paste0(each_sample,".pdf")),width = 9,height = 6)
        print(p)
        dev.off()
        write.csv(plot_data,file = file.path(output_path,paste0(each_sample,".csv")),row.names = F)
        
      }
      )
      
      ## vocano plot
      type_sum = purrr::map2(tpm_list,tpm_type,gene_type_summary)
      type_sum = do.call(rbind,type_sum)
      diff_gene_data = read.csv(file.path(data_path,"all_sample_diff_trans.csv"))
      lapply(unique(type_sum$type),function(each_type){
        genes = dplyr::filter(type_sum,type == each_type)
        
        res = dplyr::filter(diff_gene_data,known_gene %in% genes$gene )
        if (nrow(res) == 0){
          print("there is no diff genes")
        }
        else{
          plot_data = prepare_plot_valcano(res)
          p = plot_volcano(plot_data)
          ggsave(p,filename = file.path(output_path,paste0(each_type,"_diff.pdf")),width = 20,height = 16)
          write.csv(res,file.path(output_path,paste0(each_type,"_diff.csv")),row.names = F)
        }
        
      })
      
    }
    
    generate_znn_type_data = function(znn_data,znn_root){
      library(xlsx)
      znn_data = read.xlsx(znn_data,1)
      pseudo = znn_data[grepl("pseudo|Pseudogene|novel gene identicle to IGHV1OR15-1",znn_data$By.ZNN.description),]
      protein = znn_data[grepl("protein coding|novel protein|multidrug resistance-related protein|SMIM11B|beta-thionase|novel histone H2B family protein|U2AF1L5|novel immunoglobulin lambda variable gene|Endosomal Transmembrane Epsin Interactor 3",znn_data$By.ZNN.description),]
      new_trans = znn_data[grepl("TEC|misc_RNA|Novel Transcrip",znn_data$By.ZNN.description),]
      nc_RNA = znn_data[grepl("ncRNA|HSALNG|Telomerase-vert|RF01684",znn_data$By.ZNN.description),]
      
      re = anti_join(znn_data,pseudo) %>% anti_join(protein) %>% anti_join(new_trans) %>% anti_join(nc_RNA)
      writeLines(pseudo$known_gene,file.path(znn_root,"znn.pseudo"))
      writeLines(protein$known_gene,file.path(znn_root,"znn.protein"))
      writeLines(new_trans$known_gene,file.path(znn_root,"znn.newtrans"))
      writeLines(nc_RNA$known_gene,file.path(znn_root,"znn.ncrna"))
      
    }
  }
}

if(T){
  config = read_yaml("./config.yaml")
  root_data_path= file.path(config$output$root,"stringtie")
  input_file = file.path(config$input$root,"fastq_list.csv")
  pipe_list = config$config$stringtie$pipeline
  for(each_pipe in pipe_list){
    ## debug use
    # each_pipe = pipe_list[1]
    
    data_path = file.path(root_data_path,each_pipe)
    tpm_data = file.path(data_path,"all_sample_tpm.csv")
    tpm_data = as.data.frame(fread(tpm_data))
    ## debug use 
    # tpm_data = "/mnt/d/tpm_test/all_sample_tpm(3).csv"
    # tpm_data = read.csv(tpm_data)
    
    #used by filter_by_tpm_sum, plot_pheatmap and tpm_summary
    meta_col = c(1:9) 
    select_tpm_data = filter_by_tpm_sum(tpm_data)
    
    #plot heatmap
    # debugonce(plot_pheatmap)
    plot_pheatmap(select_tpm_data[[1]])
    # debugonce(stat_tpm)
    stat_tpm(select_tpm_data)
  }
}

