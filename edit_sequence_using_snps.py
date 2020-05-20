
import pandas as pd
import sys

inputFile = sys.argv[1]
outFile = sys.argv[2]
refgenome = sys.argv[3]
afcuoff = sys.argv[4]
gapfile = sys.argv[5]


def readFasta(fastafile):
    """
    read fasta seqence
    """
    seq={}
    for line in fastafile:
        if line.startswith('>'):
            name=line.replace('>','').split()[0]
            seq[name]=''
        else:
            seq[name]+=line.replace('\n','').strip()
    return seq

def snpReplace(start,alt,seq):
    """
    replace reference with snp sites
    """
    seq = list(seq)
    seq[start] = alt
    seq = ''.join(seq)
    return seq

# 读取reference
f_o=open(refgenome)
reference = readFasta(f_o)
f_o.close()

seq = reference['BetaCoV/Wuhan-Hu-1/2019|EPI_ISL_402125']

# 读取vcf文件
vcfinfo = pd.read_csv(inputFile,sep='\t',header=None)
vcfinfo.columns = ['libs','Start','End','Ref','Alt','qual','DP','altCount','frequency']
vcfinfo = vcfinfo[vcfinfo.frequency >= float(afcuoff)]
for index,row in vcfinfo.iterrows():
    seq = snpReplace(row.Start,row.Alt,seq)

# 插入coverage=0的位点 并进行改写
gapdf = pd.read_csv(gapfile,sep='\t',header=None)
gapdf.columns = ['libs','position','coverage']
gapdf.head()

for index,row in gapdf.iterrows():
    seq = snpReplace(row.position-1,"N",seq)


# 写出
_str = []
l = open(outFile,'w')
_str.append(">%s\n%s\n"%(outFile,seq))
l.writelines(_str)
l.close()