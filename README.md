Fully Automated glycoRNA sequencing analysis pipeline.

This project referenced the following projects:
https://github.com/ChangLab/FAST-iCLIP.git
 

# Usage
1. download ref data from google drive:
  https://drive.google.com/file/d/1GaTAQnNClPnUi9mHN2MftQczP9DSkrDO/view?usp=drive_link
2. edit the metadata file **fastq_list.csv** and the configuration file **config.yaml**
   
   1.  fastq_list.csv:
        - name: sample name
        - sample: fastq file name, if the sample name is `TEST`, then the name of the fastq file in the input folder must be `TEST.R1.fastq.gz` and `TEST.R2.fastq.gz`
        - type: input or enrich type of glycoRNA 
        - group: paired sample numbers
        - dutp: If it is a dUTP library construction, then fill in as 1, otherwise as 0.
        - umi: f it is UMI library construction, then fill in the UMI length, otherwise, fill in as 0.
    2. config.yaml:

        - input:

          - root: fastq file folders
          - list: the path of fastq_list.csv
        - output:
          - root: output path
3. run `python py_scripts/main.py`
4. run `python py_scripts/stringtie.py`

# Dependencies
The version numbers listed have been tested successfully. There can be difficulties if you choose to run updated versions of some of these dependencies.

- python 3.6.11
- bedtools 2.29.2
- fastp 0.23.2
- hisat2 2.2.1
- pandas 1.15
- picard 2.23.4
- samtools 1.10
- stringtie 2.2.1: https://github.com/gpertea/stringtie
- ballgown 2.28.0: https://github.com/alyssafrazee/ballgown
- gencore 0.17.2: https://github.com/OpenGene/gencore