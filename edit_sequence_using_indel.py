### analysis indel data

import pandas as pd
import sys

inputFile = sys.argv[1]
outFile = sys.argv[2]
refgenome = sys.argv[3]
afcuoff = sys.argv[4]


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

# 读取reference
f_o=open(refgenome)
reference = readFasta(f_o)
f_o.close()

seq = reference['BetaCoV/Wuhan-Hu-1/2019|EPI_ISL_402125']

indeldf = pd.read_csv(inputFile,sep='\t',header=None)
indeldf.columns = ['libs','Start','End','Ref','Alt','qual','DP','altCount','frequency']
indeldf = indeldf[indeldf.frequency >= float(afcuoff)]

seq = list(seq)
for index,row in indeldf.iterrows():
    start = row.Start
    ref = row.Ref
    alt = row.Alt
    if len(alt) > len(ref):
        type = 'insertion'
        seq[start] = alt
        # 将涉及到的位置进行清空
        i = start + 1
        while i < start+len(ref):
            seq[i] = ''
            i += 1
        print(seq[start])
    else:
        type = 'deletion'
        # 添加“-”至等长
        gapnumber = len(ref) - len(alt)
        i = 0
        while i < gapnumber:
            alt = alt + '-'
            i += 1
        
        seq[start] = alt
        i = start + 1
        while i < start+len(ref):
            seq[i] = ''
            i += 1
        print(seq[start])

seq = ''.join(seq)

# 写出
_str = []
l = open(outFile,'w')
_str.append(">%s\n%s\n"%(outFile,seq))
l.writelines(_str)
l.close()