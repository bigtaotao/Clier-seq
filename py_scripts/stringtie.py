#In[]
import yaml
import logging
import os
import subprocess
import pandas as pd
from func import *
import multiprocessing
## should change the thread number to constant
def run_stringtie(rmdup_bam,ref,threads,is_utp,sample_outroot,pipe_name = ""):
    output_file = os.path.join(sample_outroot,"map_{}_rmdup.gtf".format(pipe_name))
    logging.info("run stringtie firstly")
    ## https://www.biostars.org/p/390873/ DO NOT USE --rf
    cmd = ["stringtie","-p",str(threads),"-G ",ref,"-l","NT","-o",output_file, rmdup_bam]
    to_str_cmd(cmd,"stringtie") 
    return(output_file)
def make_combined(merge_list,ref,threads,outputpath):
    cmd = ["stringtie --merge -p",str(threads),"-G",ref,"-l","NT","-o",outputpath,merge_list]
    to_str_cmd(cmd,"stringtie merge")
def compare_gtf(ref,combined_gtf,output_root):
    cmd = ["gffcompare -r",ref," -G -o",output_root,combined_gtf]
    to_str_cmd(cmd,"gffcompare")
def run_stringtie_eB(rmdup_bam,ref,threads,is_utp,sample_outroot,pipe_name = ""):
    output_file = os.path.join(sample_outroot,"map_{}_rmdup_mapmerged.gtf".format(pipe_name))
    print(sample_outroot)
    logging.info("run stringtie-eB secondly")
    ## https://www.biostars.org/p/390873/ DO NOT USE --rf
    cmd = ["stringtie","-e -B -p",str(threads),"-G ",ref,"-l","NT", "-o",output_file, rmdup_bam]
    to_str_cmd(cmd,"stringtie-eB")
def run_stringtie_each(each_row,config,pipe_outroot,ref,each_pipe,parallel_threads):
                is_utp = each_row["dutp"]
                each_sample = each_row["sample"]   
                sample_outroot = os.path.join(pipe_outroot,each_sample)
                check_dir(sample_outroot,is_dir = True)
                bam_root = os.path.join(config["output"]["root"],each_sample)
                rmdup_bam = os.path.join(bam_root,"map_{}_rmdup.bam".format(each_pipe))
                print("processing sample: {}".format(each_sample))
                res = run_stringtie(rmdup_bam,ref,parallel_threads,is_utp,sample_outroot,each_pipe)
                return(res)
def run_stringtie_eB_each(each_row,config,pipe_outroot,combined_gtf,each_pipe,parallel_threads):
            is_utp = each_row["dutp"]
            each_sample = each_row["sample"]   
            sample_outroot = os.path.join(pipe_outroot,each_sample)

            bam_root = os.path.join(config["output"]["root"],each_sample)
            rmdup_bam = os.path.join(bam_root,"map_{}_rmdup.bam".format(each_pipe))        

            run_stringtie_eB(rmdup_bam,combined_gtf,parallel_threads,is_utp,sample_outroot,each_pipe)
def main():
    ## read configuration
    config = yaml.safe_load(open("./config.yaml", "r"))
    log_file = os.path.join(config["output"]["root"],"stringtie_history.log" )
    check_dir(log_file)
    logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                        level=logging.DEBUG,filename=log_file)

    input_list = config["input"]["list"]
    meta = pd.read_csv(input_list,header=0)
    print(meta)
    stringtie_outroot = os.path.join(config["output"]["root"],"stringtie")
    parallel_threads = config["config"]["stringtie"]["parallel"]
    stringtie_threads = config["config"]["stringtie"]["threads"]
    pipe_list = config["config"]["stringtie"]["pipeline"]
    for each_pipe in pipe_list:
        ref = os.path.join(config["pipeline_config"]["ref_root"],config["pipeline_config"]["pipe"][each_pipe]["gtf"])
        pipe_outroot = os.path.join(stringtie_outroot,each_pipe)
        check_dir(pipe_outroot,is_dir=True)
        
        if config["pipeline_config"]["pipe"][each_pipe]["stringtie"] != "":
            combined_gtf = os.path.join(config["pipeline_config"]["ref_root"],config["pipeline_config"]["pipe"][each_pipe]["stringtie"])
            logging.info("stringtie use gtf file from config.yaml")
        else:
            ## first stringtie
            stringtie_input = []
            for i,each_row in meta.iterrows():
                stringtie_input.append(each_row)
            
            res = []
            pool = multiprocessing.Pool(processes=parallel_threads)
            for each in stringtie_input:
                res.append(pool.apply_async(run_stringtie_each,(each,config,pipe_outroot,ref,each_pipe,stringtie_threads,)))
            pool.close()
            pool.join()
            stringtie_file = [each.get() for each in res]
            ## make merge list
            merge_list_file = os.path.join(pipe_outroot,"mergelist.txt")
            with open(merge_list_file,"w") as f:
                for each in stringtie_file:
                    f.writelines(each + "\n")
            ## make combined gtf
            combined_gtf = os.path.join(pipe_outroot,"stringtie_merged.gtf")
            make_combined(merge_list_file,ref,parallel_threads,combined_gtf)
            ## compare with ref
            compare_outroot = os.path.join(pipe_outroot,"merged")
            compare_gtf(ref,combined_gtf,compare_outroot)
        ## second stringtie
        stringtie_eB_input = []
        for i,each_row in meta.iterrows():
            stringtie_eB_input.append(each_row)
        
        res = []
        pool = multiprocessing.Pool(processes=parallel_threads)
        for each in stringtie_eB_input:
            res.append(pool.apply_async(run_stringtie_eB_each,(each,config,pipe_outroot,combined_gtf,each_pipe,stringtie_threads,)))
        pool.close()
        pool.join()
        logging.info("stringtie pipeline {} finished".format(each_pipe))
        
if __name__ == "__main__":  
    main()
# %%
