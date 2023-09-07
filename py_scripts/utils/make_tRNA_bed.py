#In[]
fa_file = "/mnt/d/znn/pipeline/glycoRNA/ref/tRNA/hg38_tRNAs.fa"

with open(fa_file,"r") as f:
    fa = f.read().splitlines()

# %%
import re
fa_name = []
for each_line in fa:
    if each_line.startswith(">"):
        fa_name.append(each_line)

#In[]
fa_contents = []
for each in fa_name:
    start_span = 200
    end_span = 200
    contents = each.split(' ')
    name = contents[0].replace('>','')
    position = contents[10]
    position_chr = position.split(':')[0]
    position_start = int(position.split(':')[1].split('-')[0])-start_span
    position_end = int(position.split(':')[1].split('-')[1])+end_span
    strand = contents[11][1]
    fa_contents.append([position_chr,position_start,position_end,name,strand])
# %%
import pandas as pd
df = pd.DataFrame(fa_contents,columns=['chr','start','end','name','strand'])
df.to_csv("/mnt/d/znn/pipeline/glycoRNA/ref/tRNA/hg38_tRNAs.bed",sep='\t',index=False,header=False)
# %%
## 
# then use bedtools getfasta to get the sequence
# bedtools getfasta -nameOnly -fi ../GRCh38/GRCh38.fa -bed hg38_tRNAs.bed -fo tRNAs.fa
# hisat2-build -p 8 tRNAs.fa tRNA
# faidx --transform bed tRNAs.fa > tRNAs.bed
# awk -v OFS="\t" '{print $0,$1,0,"+"}' tRNAs.bed > trnas_add.bed
# mv trnas_add.bed tRNAs.bed