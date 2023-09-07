import logging
import os
def  to_str_cmd(cmd,logginginfo):
    str_cmd = ""
    for each in cmd:
        str_cmd += each + " "
    logging.info("running %s: %s"%(logginginfo,str_cmd))
    status = os.system(str_cmd)
    if status != 0:
        logging.error("Failed to run %s: %s"%(logginginfo,str_cmd))
        exit(1)
    return (0)
def check_dir(path,is_dir = False):
    if(is_dir):
        # print(path)
        if not os.path.exists(path):
            # print("make dir")
            os.makedirs(path,exist_ok=True)
    else:
        dir_name = os.path.dirname(path)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
    return 0 