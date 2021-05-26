import numpy as np
from multiprocessing import Pool
import pandas as pd
import sys
import scipy.sparse as sp
from scipy.sparse import save_npz
from util import readFileid


def genkmerSampleAbunDict(kmeridfiledir,filelist,fileidDict,lowfeat,highfeat):
    print('genkmerSampleAbunDict')
    kmerSampleAbunDict = {}
    for infile in filelist:
        f = open(kmeridfiledir + infile + 'id')
        for line in f:
            kmerid, abun = line.strip().split('\t')
            kmerid = int(kmerid)
            if kmerid >= lowfeat and kmerid <= highfeat:
                if kmerid not in kmerSampleAbunDict.keys():
                    kmerSampleAbunDict[kmerid] = {}
                kmerSampleAbunDict[kmerid][fileidDict[infile]] = abun 
    return kmerSampleAbunDict


def kmerAbunbybatch_para(args):
    subsubkmers=args[0]
    filelist=args[1]
    outdir=args[2]
    kmeridfiledir=args[3]
    fileidfile=args[4]
    fileidDict, fileNamedict, filelist = readFileid(fileidfile)
    lowfeat = subsubkmers[0]
    highfeat = subsubkmers[-1]
    outfile=outdir+'kmer_'+str(lowfeat)+'_'+str(highfeat)+'.abun.npz'
    kmerSampleAbunDict=genkmerSampleAbunDict(kmeridfiledir,filelist,fileidDict,lowfeat,highfeat)
    print('finish genkmerSampleAbunDict')

    ROWNUM=len(subsubkmers)
    COLNUM=len(filelist)

    mat = sp.dok_matrix((ROWNUM, COLNUM), dtype=np.float32)
    for kmerid, fileidDict in kmerSampleAbunDict.items():
        for fileid, abun in fileidDict.items():
            mat[kmerid-lowfeat, fileid] = int(abun)
    mat = mat.tocsr()
    save_npz(outfile, mat)


def proKmerBatch(filelist,kmeridfiledir,kmermatrixdir,subkmers,corenum,fileidfile):
    data=[]
    subsubkmers = np.array_split(subkmers, corenum)
    data = []
    for i in range(corenum):
        onedata = [subsubkmers[i],filelist,kmermatrixdir,kmeridfiledir,fileidfile]
        data.append(onedata)
    pool = Pool(corenum)
    pool.map(kmerAbunbybatch_para, data)
    pool.close()


def genKmersparsematrix(KMERNUM,filelist,kmeridfiledir,kmermatrixdir,nodenum,docNum,coreNum,fileidfile):
    kmeridlist=[i for i in range(KMERNUM)]
    subkmeridlist=np.array_split(kmeridlist,nodenum)
    print('number of kmers: '+str(len(kmeridlist)))
    print('number of kmer subsets: '+str(len(subkmeridlist)))
    subkmers=subkmeridlist[docNum]
    print('number of kmers in subsets: '+str(len(subkmers)))
    proKmerBatch(filelist,kmeridfiledir,kmermatrixdir,subkmers,coreNum,fileidfile)

