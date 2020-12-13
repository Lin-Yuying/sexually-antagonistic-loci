'''
This script detects male specific SNP from genotype .vcf file.
Author: Y.Lin
Date: Jul 13, 2020
Usage: python sex_specific_SNP.py -h
'''
import gzip
import os
import argparse
import collections
import argparse
import pandas as pd

############ 00. get command line arguments #############
parser = argparse.ArgumentParser()
parser.add_argument('-g','--genotype',
                    help='VCF file',required=True)
parser.add_argument('-pop','--population',
                    help='population .txt file, with tab seperator',required=True)
parser.add_argument('-w','--windows',
                   help = "slide window size")
parser.add_argument('-f','--frequency',
                   help='mini freq of male-specific SNP,default: 0')
parser.add_argument('-o','--output',
                   help='Output .csv name')
parser.add_argument('-c','--chr_len',
                   help='chromosome length .csv, format: chr1\t34115677',required=True)
args=parser.parse_args()

#args
vcf = args.genotype
pop = args.population
freq = float(args.frequency) if args.frequency else 0
chr_len_csv = args.chr_len
out = args.output
windows = int(args.windows)

########## 01. get column info for individuals ##########
def getCol(vcf,pop):
    with open(vcf,'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                pass
            elif line.startswith("#"):
                line = line.lstrip("#")
                chrom,pos,ID,ref,alt,qual,ftr,info,fmt,*indvs=line.split("\t")
            else:
                pass
    col=[]
    with open(pop,'r') as p:
        for i in p:
            i=i.strip("\n")
            col.append(indvs.index(i))
    return col
 
########## calculate prop of male heterozygous#####
def getProp(male_geno,female_geno):
    # males are heterozygous, females are homozygous 
    male = [i for i in male_geno if i not in set(female_geno)]
    return len(male)/len(male_geno)
  
col = getCol(vcf,pop)
######## 02. get genotype data base on col information########
######## 0/0:0,66,255:22:22,0:74 
with open(vcf,'r') as inf,open(out+".prop.csv",'w') as outf:
    for line in inf:
        line = line.rstrip("\n")
        if line.startswith('#'):
            pass
        else:
            chrom,pos,ID,ref,alt,qual,ftr,info,fmt,*indvs=line.split("\t")
            male_geno=[]
            female_geno=[]
            for i in range(len(indvs)):
                if i in col:
                    male_geno.append(indvs[i].split(":")[0])
                else:
                    female_geno.append(indvs[i].split(":")[0])
            if './.' in male_geno:
                male_geno = [i for i in male_geno if i != './.']
            else:
                pass
            if './.' in female_geno:
                female_geno = [i for i in female_geno if i != './.']
            else:
                pass
            if len(set(female_geno)) == 1 and set(female_geno) != {'0/1'} and set(female_geno) != set(male_geno):
                #print(chrom,pos,male_geno,female_geno)
                prop = getProp(male_geno,female_geno)
                if float(prop) > freq:
                    outf.write("{}\t{}\t{}\n".format(chrom,pos,prop))
            else:
                pass

##### 03. slide window based male-specific SNP calculation ########
# import data with pd
df = pd.read_csv(out+".prop.csv", sep = "\t", header = None, names =['chr','pos','prop'])
#chr_list = sorted(list(set(df["chr"])))
chr_list= []
max_pos = {}
with open(chr_len_csv,"r") as infile:
    for line in infile:
        line = line.rstrip("\n")
        chrom, m_pos = line.split("\t")
        max_pos[chrom] = m_pos
        chr_list.append(chrom)
with open(out+"_specific_snp.csv","w") as outfile:
  for chrom in chr_list:
    start = 0
    end = start + windows
    while end <= int(max_pos[chrom]):
        count = len(df[(df["chr"]== chrom) & (df["pos"] <= end) & (df["pos"] > start)])
        outfile.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, count))
        start = start + windows/2
        end += windows/2
    count = len(df[(df["chr"]== chrom) & (df["pos"] <= int(max_pos[chrom])) & (df["pos"] > end)])
    outfile.write("{}\t{}\t{}\t{}\n".format(chrom, start, max_pos[chrom], count))
