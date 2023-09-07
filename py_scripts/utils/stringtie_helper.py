from func import *
import argparse


def transform(old_gtf,new_gtf,transform_gtf):
    cmd = ['gffcompare', '-r', old_gtf, '-G', '-o', transform_gtf, new_gtf]
    to_str_cmd(cmd, "gffcompare")
parser = argparse.ArgumentParser(description='stringtie transter')
parser.add_argument('-o', help='old gtf', required=True)
parser.add_argument('-n', help='new gtf', required=True)
parser.add_argument('-p', help='output prffix', required=True,default="transform")


if __name__ == '__main__':
    args = parser.parse_args()
    old_gtf = args.o
    new_gtf = args.n
    transform_gtf = args.p
    transform(old_gtf,new_gtf,transform_gtf)