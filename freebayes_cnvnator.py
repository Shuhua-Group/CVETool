#!/usr/bin/env python
# @coding: utf-8 
# @Author : masen

import os 
from collections import Counter
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d','--dir', required=True, help='work dir')
parser.add_argument('-id','--id', required=True, help='sampel id ')
parser.add_argument('-bam','--bam', required=True, help='bam file path')
parser.add_argument('-ref','--ref', required=True, help='ref file path')
parser.add_argument('-sv_region','--sv_region', required=True,help='sv region such as  chr:start-end')
parser.add_argument('-snp_region','--snp_region', required=True, help='snp region such as  chr:position')
parser.add_argument('-cnvnator_loc','--cnvnator_loc', required=True, help='Soft loc of CNVnator')
parser.add_argument('-freebayes_loc','--freebayes_loc', required=True, help='Soft loc of freebayes')
parser.add_argument('-cr','--cnv_root', help='cnvnator root file path')

args = parser.parse_args()

dir=args.dir
id=args.id
bam=args.bam
ref=args.ref
sv_region=args.sv_region
#snp_region=args.snp_region
snp_region=args.snp_region.split(":")[0]+":"+str(int(args.snp_region.split(":")[1])-1)+"-"+args.snp_region.split(":")[1]
cnvnator_loc=args.cnvnator_loc
freebayes_loc=args.freebayes_loc

cnv_root_file=args.cnv_root

chrom=sv_region.split(":")[0]
region_start=sv_region.split(":")[1].split("-")[0]
region_end=sv_region.split(":")[1].split("-")[1]

#def manta(id,bam,ref,dir,sv_region):
#    os.system("/picb/humpopg-bigdata5/masen/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py --bam {bam} --referenceFasta {ref} --runDir {dir}/{id}_manta --region={region} ;python2 {dir}/{id}_manta/runWorkflow.py".format(id=id,bam=bam,dir=dir,ref=ref,region=sv_region))
#    os.system("bgzip -d {dir}/{id}_manta/results/variants/diploidSV.vcf.gz".format(dir=dir,id=id))

def freebayes(id,bam,ref,dir,snp_region,freebayes_loc):
    os.system("{freebayes_loc} -f {ref} --region {region} --cnv-map {dir}/{id}.cnv.txt {bam} --vcf {dir}/{id}.freebayes.vcf".format(id=id,bam=bam,dir=dir,ref=ref,region=snp_region,freebayes_loc=freebayes_loc))

def freebayes2(id,bam,ref,dir,snp_region,freebayes_loc):
    os.system("{freebayes_loc} -f {ref} --region {region} {bam} --vcf {dir}/{id}.freebayes.vcf".format(id=id,bam=bam,dir=dir,ref=ref,region=snp_region,freebayes_loc=freebayes_loc))

def cnvnator(id,bam,dir,sv_region,cnvnator_loc):#ref -d {dir}
    chrom=sv_region.split(":")[0]
    with open("{dir}/{samples_id}.cnv_input.region.txt".format(dir=dir,samples_id=id),"wt")as f:
        f.write(sv_region+"\nexit\nEOF")
    os.system("{cnvnator_loc} -root {dir}/{samples_id}.root -tree {bam} -chrom {chrom};{cnvnator_loc} -root {dir}/{samples_id}.root -his 100 -chrom {chrom} -fasta {ref};{cnvnator_loc} -root {dir}/{samples_id}.root -chrom {chrom} -stat 100;{cnvnator_loc} -root {dir}/{samples_id}.root -chrom {chrom}  -partition 100;{cnvnator_loc} -root {samples_id}.root -genotype 100 < {dir}/{samples_id}.cnv_input.region.txt >{dir}/{samples_id}.cnvnator.txt".format(bam=bam,dir=dir,samples_id=id,chrom=chrom,sv_region=sv_region,cnvnator_loc=cnvnator_loc))
def cnvnator2(id,bam,dir,sv_region,cnvnator_loc,cnv_root_file):#ref -d {dir}
    chrom=sv_region.split(":")[0]
    with open("{dir}/{samples_id}.cnv_input.region.txt".format(dir=dir,samples_id=id),"wt")as f:
        f.write(sv_region+"\nexit\nEOF")
    os.system("{cnvnator_loc} -root {cnv_root_file} -genotype 100 < {dir}/{samples_id}.cnv_input.region.txt >{dir}/{samples_id}.cnvnator.txt".format(bam=bam,dir=dir,samples_id=id,chrom=chrom,sv_region=sv_region,cnv_root_file=cnv_root_file,cnvnator_loc=cnvnator_loc))
    
def file_make(id,dir):
    man_results="{dir}/{id}_manta/results/variants/diploidSV.vcf".format(id=id,dir=dir)
    out=open("{dir}/{id}.txt".format(id=id,dir=dir),"wt")
    with open(man_results) as f:
        for i in f.readlines():
            if not i.startswith("#"):
                CHROM ,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,id_genotype=i.split("\t")[:]
                end=INFO.split(";")[0].split("=")[1]
                SVTYPE=INFO.split(";")[1].split("=")[1]
                type=id_genotype.split(":")[0]
                if SVTYPE=="DEL":
                    print(CHROM,POS,end,id,"1",sep="\t",file=out)
                elif SVTYPE=="DUP":
                    if type=="0/1":
                        print(CHROM,POS,end,id,"3",sep="\t",file=out)
                    elif type=="1/1":
                        print(CHROM,POS,end,id,"4",sep="\t",file=out)
                    else:
                        print(CHROM,POS,end,id,"#",sep="\t",file=out)    
            else:
                pass
    out.close()
def cnv_read(id,dir):
    cnv_results="{dir}/{id}.cnv.vcf".format(id=id,dir=dir)
    out=open("{dir}/{id}.cnv.txt".format(id=id,dir=dir),"wt")
    with open(cnv_results) as f:
        for i in f.readlines():
            if not i.startswith("#"):
                CHROM ,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,id_genotype=i.split("\t")[:]
                end=INFO.split(";")[0].split("=")[1]
                if int(POS) >=int(region_start) and int(POS) <=int(region_end):
                    SVTYPE=INFO.split(";")[1].split("=")[1]
                    count_1=Counter(id_genotype.split(":")[0].split("/")[:])["1"]
                    count_0=Counter(id_genotype.split(":")[0].split("/")[:])["0"]
                    cn=id_genotype.split(":")[1]
                    num_out=int(count_1)*int(cn)+int(count_0)*1
                    if SVTYPE=="DEL" and num_out !=0:
                        print(CHROM,POS,end,id,"1",sep="\t",file=out)
                    elif SVTYPE=="DUP":
                        print(CHROM,POS,end,id,num_out,sep="\t",file=out)
            else:
                pass
    out.close()
def cnv_change(id,dir):
    cnv_results="{dir}/{id}.cnvnator.txt".format(id=id,dir=dir)
    out=open("{dir}/{id}.cnv.txt".format(id=id,dir=dir),"wt")
    with open(cnv_results) as f:
        for i in f.readlines():
            cnv_num=float(i.rstrip().split(" ")[-2])
            print(sv_region.split(":")[0],sv_region.split(":")[1].split("-")[0],sv_region.split(":")[1].split("-")[1],id,round(cnv_num),sep="\t",file=out)
    out.close()
if __name__ == '__main__':
    if cnv_root_file:
        cnvnator2(id,bam,dir,sv_region,cnvnator_loc,cnv_root_file)
    else:
        cnvnator(id,bam,dir,sv_region,cnvnator_loc)
    cnv_change(id,dir)
    if os.path.exists("{dir}/{id}.cnv.txt".format(id=id,dir=dir)) and os.path.getsize("{dir}/{id}.cnv.txt".format(id=id,dir=dir)) !=0:
        freebayes(id,bam,ref,dir,snp_region,freebayes_loc)
    else:
        freebayes2(id,bam,ref,dir,snp_region,freebayes_loc)
        #print("region maybe with no DUP/DEL,please check again!")


