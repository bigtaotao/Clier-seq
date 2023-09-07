import yaml
import logging
import os
import subprocess
import pandas as pd
from func import *
from itertools import (takewhile, repeat)
def pharse_mapping_reads(path):
    with open(path,'r') as f:
        data = f.readlines()
    
    ## 4 is the number of mapped reads
    ## 9 is the number of properly mapped reads
    mapped_reads = int(data[9].split()[0])
    return mapped_reads
def main():
    ## read configuration
    config = yaml.safe_load(open("./config.yaml", "r"))
    log_file = os.path.join(config["output"]["root"],"history_readstat.log" )
    input_list = config["input"]["list"]
    meta = pd.read_csv(input_list,header=0)
    print(meta)
    # input_samples = meta["sample"].tolist()
     ## logging configuration
         
    check_dir(log_file)  
    # logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
    #                 level=logging.DEBUG,filename=log_file)
    logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    level=logging.DEBUG,filename=log_file)
    pipe_list = config["config"]["pipeline"]
    res_dic = {}
    for i,each_row in meta.iterrows():
        each_sample = each_row["sample"]
        sample_outroot = os.path.join(config["output"]["root"],each_sample)
        
        res_dic[each_sample] = {}
        ## get raw reads
        sample_fq = each_sample + ".R2.fastq.gz"
        sample_fq = os.path.join(sample_outroot, sample_fq)
        print(sample_fq)
        if os.path.exists(sample_fq):
            output_file = sample_fq + ".count"
            cmd = ["zgrep","@",sample_fq,"|wc -l",">",output_file]
            to_str_cmd(cmd,"fastq_stat")
            with open(output_file,'r') as f:
                data = f.readlines()
            res_dic[each_sample]["raw"] = int(data[0]) *2 ## paired end
        else:
            logging.warning("No fastq file for %s"%each_sample)
        ## get mapping reads
        for each in pipe_list:
            bam = os.path.join(sample_outroot,"map_"+each+".bam")
            stat_file = os.path.join(sample_outroot,"map_"+each+".bam.stat")
            if os.path.exists(bam):
                cmd = ["samtools","flagstat",bam,">",stat_file]
                to_str_cmd(cmd,"bam_stat")
                res_dic[each_sample][each] = pharse_mapping_reads(stat_file)
            else:
                logging.warning("No bam file for %s"%each_sample)

            rmdup_bam = os.path.join(sample_outroot,"map_"+each+"_rmdup.bam")
            rmdup_stat_file = os.path.join(sample_outroot,"map_"+each+"_rmdup.bam.stat")
            if os.path.exists(rmdup_bam):
                cmd = ["samtools","flagstat",rmdup_bam,">",rmdup_stat_file]
                to_str_cmd(cmd,"bam_stat")
                res_dic[each_sample][each+"_rmdup"] = pharse_mapping_reads(rmdup_stat_file)
            else:    
                logging.warning("No rmdup bam file for %s"%each_sample)

        res_df = pd.DataFrame({each_sample:res_dic[each_sample]})
        res_df.to_csv(os.path.join(sample_outroot,"map_reads.csv"))
    res = pd.DataFrame(res_dic)
    res.to_csv(os.path.join(config["output"]["root"],"map_reads.csv"))
            
            
if __name__ == "__main__":  
    main()
