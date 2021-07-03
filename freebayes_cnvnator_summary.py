#!/usr/bin/env python
# @coding: utf-8 
# @Author : masen

import os 
from collections import Counter
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d','--dir', required=True, help='work dir')
parser.add_argument('-l','--list', required=True, help='all samples id  in a file,each one line')
parser.add_argument('-o','--out', required=True, help='out_put dir')
args = parser.parse_args()

dir=args.dir
id_file=args.list
out_put_dir=args.out
def sum_str(str):
    num=0
    if "/" in str:
        alit=str.split("/")[:]
        for i in alit:
            if i==".":
                num +=0
            elif i=="0":
                num +=0
            else :
                num +=1
    else:
        if str==".":
            num +=0
        elif str=="0":
            num +=0
        else :
            num +=1
    return num

def snp_summary():
    all_samples={}
    with open(id_file) as f:
        for i in f.readlines():
            id=i.rstrip()
            snp_file=dir+"/"+id+".freebayes.vcf"
            all_samples[id]=snp_file
    region_list=[]
    dict_id={}
    for id,snp in all_samples.items():
        dict_cnv_limit={}
        if os.path.exists("{dir}/{id}.cnv.txt".format(dir=dir,id=id)) and os.path.getsize("{dir}/{id}.cnv.txt".format(id=id,dir=dir)) !=0:
            with open("{dir}/{id}.cnv.txt".format(dir=dir,id=id))as f:
                for i in f.readlines():
                    cnv_chr,cnv_start,cnv_end,cnv_id,cnv_number=i.rstrip().split("\t")[:]
                    dict_cnv_limit[(cnv_chr,cnv_start,cnv_end)]=cnv_number 
        with open(snp) as f:
            dict_a={}
            for i in f.readlines():
                if not i.startswith("#") and dict_cnv_limit:
                    CHROM ,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,id_genotype=i.rstrip().split("\t")[:]
                    region_list.append((CHROM ,POS))
                    SVTYPE=INFO.split(";")[1].split("=")[1]
                    for strpos,up in dict_cnv_limit.items():
                        if CHROM == strpos[0] and int(POS)<=int(strpos[2]) and  int(POS)>=int(strpos[1]):
                            if "/" in id_genotype.split(":")[0]:
                                if sum_str(id_genotype.split(":")[0]) <= int(up):
                                    dict_a[(CHROM ,POS)]=str(sum_str(id_genotype.split(":")[0]))
                                else:
                                    dict_a[(CHROM ,POS)]=str(up)
                            else:
                                if sum_str(id_genotype.split(":")[0]) <= int(up):
                                    dict_a[(CHROM ,POS)]=str(sum_str(id_genotype.split(":")[0]))
                                else:
                                    dict_a[(CHROM ,POS)]=str(up)
                        else:
                            if "/" in id_genotype.split(":")[0]:
                                dict_a[(CHROM ,POS)]=str(sum_str(id_genotype.split(":")[0]))
                            else:
                                dict_a[(CHROM ,POS)]=str(sum_str(id_genotype.split(":")[0]))
                elif not i.startswith("#") :
                    CHROM ,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,id_genotype=i.rstrip().split("\t")[:]
                    region_list.append((CHROM ,POS))
                    SVTYPE=INFO.split(";")[1].split("=")[1]
                    if "/" in id_genotype.split(":")[0]:
                        dict_a[(CHROM ,POS)]=str(sum_str(id_genotype.split(":")[0]))
                    else:
                        dict_a[(CHROM ,POS)]=str(sum_str(id_genotype.split(":")[0]))
                else:
                    pass
                    #dict_id[i]=str(sum_str(id_genotype.split(":")[0]))
        dict_id[id]=dict_a
    region_list=set(region_list)
    out=open(out_put_dir+"/snp_summary.txt","wt")
    print("chr","pos","\t".join(sorted(all_samples.keys())),sep="\t",file=out)
    for i in region_list:
        chr=i[0];pos=i[1]
        final_line=[dict_id[j].get(i,"0") for j in sorted(all_samples.keys())]
        print(chr,pos,"\t".join(final_line),sep="\t",file=out)
    out.close
def sv_summary():
    all_samples={}
    #{dir}/{id}.cnvnator.txt
    with open(id_file) as f:
        for i in f.readlines():
            id=i.rstrip()
            sv_file=dir+"/"+id+".cnv.txt"
            all_samples[id]=sv_file
    region_list=[]
    dict_id={}
    for id,sv in all_samples.items():
        with open(sv) as f:
            dict_a={}
            for i in f.readlines():
                if i:#3       81596001        81600000        Druze-2 1
                    CHROM,POS,END,ID,CNV=i.rstrip().split("\t")[:]
                    region_list.append((CHROM ,POS,END))
                    dict_a[(CHROM,POS,END)]=CNV
                else:
                    pass
                    #dict_id[i]=str(sum_str(id_genotype.split(":")[0]))
        dict_id[id]=dict_a
    region_list=set(region_list)
    out=open(out_put_dir+"/sv_summary.txt","wt")
    print("chr","pos","end","\t".join(sorted(all_samples.keys())),sep="\t",file=out)
    for i in region_list:
        chr=i[0];pos=i[1];end=i[2]
        final_line=[dict_id[j].get(i,"2") for j in sorted(all_samples.keys())]
        print(chr,pos,end,"\t".join(final_line),sep="\t",file=out)
    out.close
if __name__ == '__main__':
    snp_summary();sv_summary()
