import sys
import multiprocessing
from multiprocessing import Pool
import numpy as np
import glob
import os
import subprocess
import math
from numpy import linalg as LA


def removeAllfiles(dir):
    filelist=[]
    for filename in glob.iglob(dir + '*', recursive=True):  # attention
        filelist.append(filename)
    for filename in filelist:
        os.remove(filename)

##################################################################################################
def detect_outlier(kmerfreq, threshold):
    outliers = []
    mean = np.mean(kmerfreq)
    std = np.std(kmerfreq)
    for y in kmerfreq:
        if std!=0:
            z_score = (y - mean) / std
            if np.abs(z_score) > threshold:
                outliers.append(y)
        elif std==0:
            outliers=[]
    return outliers

def calMean(kmerfreq,threshold):
    normal = []
    outliers = detect_outlier(kmerfreq, threshold)
    if len(outliers) != 0:
        for i in kmerfreq:
            if i not in outliers:
                normal.append(i)
        normalmean = np.array(normal).mean()
    else:
        normalmean=np.array(kmerfreq).mean()
    return normalmean

def readACodekmers(filename):
    f=open(filename)
    codefreqsDict={}
    kmerfreqs=[]
    codes=[]
    for line in f:
        line=line.strip()
        code,freq=line.split('\t')
        if code not in codefreqsDict.keys():
            codefreqsDict[code]=[]
        codefreqsDict[code].append(int(freq))
    return codefreqsDict

def calfilesCodemean_pro(args):
    subidfiles=args[0]
    threshold=args[1]
    indir=args[2]
    outdir=args[3]

    for filename in subidfiles:
        filename=filename.replace('id','')
        codemeans=[]
        codes=[]
        infile =indir + filename+'code'
        codefreqsDict = readACodekmers(infile)
        outfile=outdir+filename+'codemean'
        output=open(outfile,'w')
        for code,kmerfreq in codefreqsDict.items():
            codemean=calMean(kmerfreq,threshold)
            codemeans.append(codemean)
            codes.append(code)

        for code,codemean in zip(codes,codemeans):
            output.write(code)
            output.write(',')
            output.write(str(round(codemean,2)))
            output.write('\n')

def calfilesCodemean(idfilelist, corenum,threshold, indir, outdir):
    removeAllfiles(outdir)
    subidfilelist = np.array_split(idfilelist, corenum)
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [subidfilelist[i],threshold,indir,outdir]
        data.append(tuple(onedata))

    pool.map(calfilesCodemean_pro, data)
    pool.close()
    pool.join()
###########################################################################################

def readACodekmerscodelist(filename):
    f=open(filename)
    codelist=[]
    for line in f:
        line=line.strip()
        code,freq=line.split('\t')
        codelist.append(code)
    return codelist

def readACodemeancodelist(filename):
    f = open(filename)
    codeset = set()
    for line in f:
        line = line.strip()
        code, mean = line.split(',')
        codeset.add(code)
    return codeset

def checkCodemeancomplete(infilecodefreq,infilecodemean):
    codelistfreq = readACodekmerscodelist(infilecodefreq)
    codesetmean = readACodemeancodelist(infilecodemean)
    for code in codelistfreq:
        if code not in codesetmean:
            print('code '+str(code)+' does not have mean value')
            return False
    return True

def genCodelist(codefile):
    f=open(codefile)
    codelist=[]
    for line in f:
        code=line.strip()
        codelist.append(int(code))
    print('Finish reading code list: ')
    print('# of code list: '+str(len(codelist)))
    print('# of code set: '+str(len(set(codelist))))
    return sorted(codelist)
####################################################

def calCodevalue(codeVallist):
    procodeVallist=[]
    codeValarr=np.array(codeVallist)
    for val in codeVallist:
        procodeVallist.append(str(int(val)))
    return procodeVallist

def genSamplematrix_para(args):
    subfiles=args[0]
    codelist=args[1]
    indir=args[2]
    abundir=args[3]
    outdir=args[4]
    fileSampledict=args[5]

    featstr = 'feat\t' + '\t'.join([str(i) for i in codelist])

    sampleIDs=[]
    codeVallists=[]
    for filename in subfiles:
        filename = filename.replace('id', '')
        print(filename)
        infilecodemean=indir+filename+'codemean'
        infilecodefreq=abundir+filename+'code'
        if checkCodemeancomplete(infilecodefreq, infilecodemean)==False:
            print(filename+' does not have complete mean value')
            return 0

        outfile=outdir+filename+'.sample'
        f = open(infilecodemean)
        print(outfile)
        output=open(outfile,'w')
        print(outfile)


        sampleID = fileSampledict[filename]

        codeValdict = {}
        codeVallist = []
        for line in f:
            code, val = line.strip().split(',')
            if code not in codeValdict.keys():
                codeValdict[code] = float(val)
            else:
                print(code + ' exists more than one in this file')
                return 0
        # same order of feature list
        for code in codelist:
            if str(code) in codeValdict.keys():
                codeVallist.append(codeValdict[str(code)])
            else:
                codeVallist.append(0)
        codeVallist = calCodevalue(codeVallist)

        output.write(featstr)
        output.write('\n')
        output.write(sampleID)
        output.write('\t')
        output.write('\t'.join(codeVallist))
        output.write('\n')

def genSamplematrix(idfilelist,codelist, corenum,abundir, indir, outdir,fileSampledict):
    removeAllfiles(outdir)
    subidfilelist = np.array_split(idfilelist, corenum)
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [subidfilelist[i],codelist,indir,abundir,outdir,fileSampledict]
        data.append(tuple(onedata))
    pool.map(genSamplematrix_para, data)
    pool.close()
    pool.join()



