import yaml
import logging
import os
import subprocess
import pandas as pd
from func import *
from itertools import (takewhile, repeat)
def filter_fq(bam_file,umap_file,threads):
    def iter_count(file_name):
        buffer = 1024 * 1024
        with open(file_name) as f:
            buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
            return sum(buf.count("\n") for buf in buf_gen)
    
    umap_r1 = umap_file + ".fq.1.gz"
    umap_r2 = umap_file + ".fq.2.gz"

    filter_r1 = umap_file + ".filter.R1.fastq.gz"
    filter_r2 = umap_file + ".filter.R2.fastq.gz"

    cache_file = bam_file + ".fqid"

    ## bamtobed 
    if os.path.exists(cache_file):
        logging.info("skip bamToBed for {}".format(bam_file))
    else:
        cmd = ["bamToBed","-i",bam_file,"| awk '{print $4}'","| awk -F '/' '{print $1}'","| sort -S 80% --parallel={}".format(threads),"| uniq >",cache_file]
        to_str_cmd(cmd,"bamToBed")

    if iter_count(cache_file) == 0:
        logging.info("no fqid found,skip filter fq")
        return umap_r1,umap_r2

    ## filter fq
    if os.path.exists(filter_r1):
        logging.info("skip filter fq for {}".format(umap_r1))
    else:
        cmd = ["pigz -c -k -p {} -d".format(threads) ,umap_r1,"| seqkit -j {} grep -v -f ".format(threads),cache_file,"-o",filter_r1]
        to_str_cmd(cmd,"seqkit")
    if os.path.exists(filter_r2):
        logging.info("skip filter fq for {}".format(umap_r2))
    else:
        cmd = ["pigz -c -k -p {} -d".format(threads) ,umap_r2,"| seqkit -j {} grep -v -f ".format(threads),cache_file,"-o",filter_r2]
        to_str_cmd(cmd,"seqkit")
    return filter_r1,filter_r2 

def run_mapping(threads,ref,r1,r2,umap,bam,picard_mem,is_utp,is_umi,gencore_fasta):
    

    rmdup_bam = bam[:-4] + "_rmdup.bam"
    metrc = bam[:-4] + "_rmdup_metrc.txt"
    bed = rmdup_bam[:-4] + ".bed"
    
    ## run mapping
    if os.path.exists(bam):
        logging.info("skip bowtie for %s"%bam)
    else:
        #hisat2
        #dUTP   
        if(str(is_utp) == "1"):
            logging.info("is dUTP library,run hisat2 with dUTP")
            cmd1 = ["hisat2","--dta","--rna-strandness RF","-p",str(threads),"-x",ref,"-1",r1,"-2",r2,"--un-conc-gz",umap+".fq.gz"]
        else:
            logging.info("is normal library")
            cmd1 = ["hisat2","--dta","-p",str(threads),"-x",ref,"-1",r1,"-2",r2,"--un-conc-gz",umap+".fq.gz"]
        cmd2 = ["samtools","sort","-@",str(threads),"-o",bam]
        
        cmd = cmd1 + ["|"] + cmd2
        to_str_cmd(cmd,"hisat2")

    ## filter fq    
    next_fq = filter_fq(bam,umap,str(threads))

    ## run rmdup
    if os.path.exists(rmdup_bam):
        logging.info("skip rmdup for %s"%rmdup_bam)
    else:
        html_log = bam[:-4] + "_rmdup.html"
        if(is_umi):
            # update gencore to 0.17.0, bug fixed
            # cmd = ["gencore" ,"-i ",bam, "-r", gencore_fasta, "-u UMI", "-h ",html_log,"-o",rmdup_bam]
            try: ## because of UMI reads may not combine to another reads, so we need to rm these reads
                cmd = ["gencore" ,"-i ",bam, "-r", gencore_fasta, "-u UMI", "-h ",html_log]
                cmd = cmd + ["|"] + ["samtools","sort","-@",str(threads),"-o",rmdup_bam]
                to_str_cmd(cmd,"rmdup")
            except:
                cache_bam = bam[:-4] + "_rmunpair.bam"
                cmd = ["samtools", "view" ,"-@",str(threads),"-bh" ,"-f 2" ,bam ,">" ,cache_bam]
                to_str_cmd(cmd,"samtools")

                cmd = ["gencore" ,"-i ",cache_bam, "-r", gencore_fasta, "-u UMI", "-h ",html_log]
                cmd = cmd + ["|"] + ["samtools","sort","-@",str(threads),"-o",rmdup_bam]
                to_str_cmd(cmd,"rmdup")
        else:
            mem = "-Xmx"+str(picard_mem)+"g"
            cmd = ["picard","MarkDuplicates",mem,"VALIDATION_STRINGENCY=LENIENT","REMOVE_DUPLICATES=true","I=",bam,"O=",rmdup_bam,"M=",metrc]
            # to_str_cmd(cmd,"picard")

            # update gencore to 0.17.0, bug fixed
            # cmd = ["gencore" ,"-i ",bam, "-h ",html_log,"|","samtools","sort","-@",str(threads),"-o",rmdup_bam]
            # cmd = ["gencore" ,"-i ",bam, "-r", gencore_fasta,"-h ",html_log,"-o",rmdup_bam]

            ## use gencore to rm dup of un umi library
            # cmd = ["gencore" ,"-i ",bam, "-r", gencore_fasta,"-h ",html_log]
            # cmd = cmd + ["|"] + ["samtools","sort","-@",str(threads),"-o",rmdup_bam]
            to_str_cmd(cmd,"rmdup")
    
    ## run bamToBed
    if os.path.exists(bed):
        logging.info("skip bed for %s"%bed)
    else:
        cmd = ["bamToBed","-i",rmdup_bam,"|","sort","-k","1,1","-k","2,2n",">",bed]
        to_str_cmd(cmd,"bamToBed")
        
    return bed,rmdup_bam,next_fq


def call_intersect(bed,ref):
    
    res = bed[:-4] + "_anno.bed"
    # return res
    try:
        cmd = ["bedtools intersect","-a",bed,"-b",ref,"-wa","-wb","-s","-sorted",">",res]
        to_str_cmd(cmd,"bedtools")
        # status = subprocess.run(cmd).returncode
        
        other_res = bed[:-4] + "_other.bed"
        cmd = ["bedtools intersect","-a",bed,"-b",ref,"-wa","-v","-s","-sorted",">",other_res]
        to_str_cmd(cmd,"bedtools")
    except:
        logging.error("intersect failed for %s"%bed)
        
    return res

def bed_count(bed,is_utp):
    try:
        data=pd.DataFrame(pd.read_table(bed, index_col=None,header=None,sep='\t',low_memory=False))
        data = data.iloc[:,:12]
        data.columns=['Chr','Start','End','ReadName','Q','Strand','winChr','winStart','winEmd','geneName','winP','winStrand']
        ## need to analysis antisense
        # split_readname = data['ReadName'].str.split('/',expand=True)
        # split_readname.columns = ["rn","R"]
        # data = pd.concat([data, split_readname], axis=1)
        # read2 = data.query("R == '2'")
        # uniq_data = read2.drop_duplicates(subset=['ReadName'],keep='first')
        uniq_data = data.drop_duplicates(subset=['ReadName'],keep='first') #rm dup count
        geneCounts=uniq_data.groupby('geneName').size()
        geneCounts.sort_values(ascending=False)
    except:
        logging.error("bed_count failed for %s"%bed)
        geneCounts = pd.DataFrame(columns=['Chr','Start','End','ReadName','Q','Strand','winChr','winStart','winEmd','geneName','winP','winStrand'])
    outfile = bed[:-4] + "_count.bed" 
    geneCounts.to_csv(outfile)
    
def bed_count_nointersect(bed):
    bf=pd.DataFrame(pd.read_table(bed,header=None))
    bf.columns=['Chr','Start','Stop','name','Q','Strand']
    geneCounts=bf.groupby('Chr').size()
    geneCounts.sort_values(ascending=False)
    outfile = bed[:-4] + "_count.bed"
    geneCounts.to_csv(outfile)


def run_batch(pipe_list,init_r1,init_r2,is_utp,sample_outroot,config,is_umi):
    mapping_threads = config["config"]["mapping"]["threads"]
    picard_mem = config["config"]["picard"]["mem"]
    for i in range(len(pipe_list)):
        if i == 0:
            r1_path = init_r1
            r2_path = init_r2
            next_path = os.path.join(sample_outroot,"umap_"+pipe_list[i])
        else:
            r1_path = next_r1
            r2_path = next_r2
            next_path = next_path_cache + "_" + pipe_list[i]
            # next_path = os.path.join(sample_outroot,"umap_"+pipe_list[i])
        next_path_cache = next_path
        ref = os.path.join(config["pipeline_config"]["ref_root"],config["pipeline_config"]["pipe"][pipe_list[i]]["ref"])
        bam = os.path.join(sample_outroot,"map_"+pipe_list[i]+".bam")
        gencore_fasta = os.path.join(config["pipeline_config"]["ref_root"],config["pipeline_config"]["pipe"][pipe_list[i]]["fasta"])

        bed,rmdup_bam,next_fq = run_mapping(mapping_threads,ref,r1_path,r2_path,next_path,bam,picard_mem,is_utp,is_umi,gencore_fasta)

        intersect = call_intersect(bed,os.path.join(config["pipeline_config"]["ref_root"],config["pipeline_config"]["pipe"][pipe_list[i]]["bed"]))
        bed_count(intersect,is_utp)
        next_r1 = next_fq[0]
        next_r2 = next_fq[1]

def main():
    ## read configuration
    config = yaml.safe_load(open("./config.yaml", "r"))
    log_file = os.path.join(config["output"]["root"],"history.log" )
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

    for i,each_row in meta.iterrows():
        each_sample = each_row["sample"]
        is_utp = each_row["dutp"]
        is_umi = each_row["umi"]
        sample_outroot = os.path.join(config["output"]["root"],each_sample)
        check_dir(sample_outroot,is_dir = True)

       
        if config["config"]["fastp"]["run"]:
            sample_r1 = each_sample + ".R1.fastq.gz"
            sample_r2 = each_sample + ".R2.fastq.gz"
            sample_r1_path = os.path.join(config["input"]["root"],sample_r1)
            sample_r2_path = os.path.join(config["input"]["root"],sample_r2)
            sample_r1_out_path = os.path.join(sample_outroot,sample_r1)
            sample_r2_out_path = os.path.join(sample_outroot,sample_r2)
            sample_html_log = os.path.join(sample_outroot,each_sample + ".html")
            if(os.path.exists(sample_r1_out_path) and os.path.exists(sample_r2_out_path)):
                logging.info("skipping running fastp: %s"%each_sample)
            else:
                if(int(is_umi) == 0):
                # cmd = ["fastp","-w",str(config["config"]["fastp"]["threads"]),"-f","14","-F","14","-i",sample_r1_path,"-I",sample_r2_path,"-o",sample_r1_out_path,"-O",sample_r2_out_path,"-h",sample_html_log]
                    cmd = ["fastp","-w",str(config["config"]["fastp"]["threads"]),"-i",sample_r1_path,"-I",sample_r2_path,"-o",sample_r1_out_path,"-O",sample_r2_out_path,"-h",sample_html_log]
                    to_str_cmd(cmd,"fastp")
                else:
                    logging.info("run fastp as umi mode")
                    ## UMI SKIP is constentially: 4
                    cmd = ["fastp","-w",str(config["config"]["fastp"]["threads"]),"-i",sample_r1_path,"-I",sample_r2_path,"-o",sample_r1_out_path,"-O",sample_r2_out_path,"-h",sample_html_log, "-U --umi_loc per_read --umi_len ",str(is_umi),"--umi_prefix UMI --umi_skip 4"]
                    to_str_cmd(cmd,"fastp")
            mapping_r1_path = sample_r1_out_path
            mapping_r2_path = sample_r2_out_path
        else:
            sample_r1 = each_sample + ".R1.fastq.gz"
            sample_r2 = each_sample + ".R2.fastq.gz"
            sample_r1_path = os.path.join(config["input"]["root"],sample_r1)
            sample_r2_path = os.path.join(config["input"]["root"],sample_r2)
            mapping_r1_path = sample_r1_path
            mapping_r2_path = sample_r2_path
        
        # pipe_list = ["ebv","repeat","rRNA","tRNA","genome"]
        pipe_list = config["config"]["pipeline"]
        run_batch(pipe_list,mapping_r1_path,mapping_r2_path,is_utp,sample_outroot,config,is_umi)
if __name__ == "__main__":  
    main()
# %%
