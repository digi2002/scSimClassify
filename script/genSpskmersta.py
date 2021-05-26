"""
	Function: build a CKG look-up table
	Input:
        task=sys.argv[1](BREAST,BREAST16)
        datatype=sys.argv[2] (all,mapped)
        random=int(sys.argv[3]) (0,1,2,3,4,5,1000)
        fold = int(sys.argv[4]) (0,1,2,3,4)
        nbit=int(sys.argv[5]) (32,16)
        zeroper=float(sys.argv[6])(0.1)
        topper=float(sys.argv[7])(0.05)
        flag=sys.argv[8](selstd,selhash)
    Input Description:
        task: k-mer length
        datatype: read type
        random,fold: determine k-mer matrix from 'random' run, 'fold' fold
        nbit: fingerprint bit size
        zeroper,topper: representative k-mer selection criterion
        flag: representative k-mer selection (selkmer) or build a k-mer group look-up table (selhash)

	Author: Qi Sun
"""


import sys
import multiprocessing
from multiprocessing import Pool
import numpy as np
from simhash import Simhash,SimhashIndex
from util import readFilelist,readReadCount,genReadcountVec,linecount,linecountpnz,readFileid
import util
from scipy.sparse import save_npz,load_npz
import glob
import datetime
import os
import subprocess
from genKmercode import selCode,genKmercodedict





def readFeatName(filename):
    f = open(filename)
    line = f.readline()
    tokens = line.strip().split('\t')
    return tokens[1:]

def cntnonZero(vals):
    cnt=0
    for val in vals:
        if val!=0:
            cnt=cnt+1
    return cnt

def calkmerSta(vals):
    cnt = cntnonZero(vals)
    std = np.std(vals)
    mean = np.mean(vals)
    return round(std,3),round(mean,3),cnt

def checksampleidlist(sampleids,selsampleids):
    if len(sampleids)!=len(selsampleids):
        print('# of train samples are not equal')
        return False
    for m in sampleids:
        if m not in set(selsampleids):
            print('id not the same')
            return False
    return True
#####################################################################################

def checksta(abunfile,stafile):
    flag=False
    if os.path.exists(stafile):
        if linecountpnz(abunfile) == linecount(stafile):
            flag = True
            print(stafile + ' exists')
    return flag

#kmer_28623798_28722500.abun.npz
def getLowfeat(filename):
    tokens=filename.split('_')
    lowfeat=int(tokens[1])
    return lowfeat

def getAbundocs(indir):
    abunfilelist = []
    for filename in glob.iglob(indir + 'kmer_*.abun.npz', recursive=True):  # attention
        abunfile=filename.replace(indir,'')
        abunfilelist.append(abunfile)
    return abunfilelist

def genTrainidselector(filelist,trainrunids):
    trainselecter=[]
    for i in range(len(filelist)):
        if filelist[i] in trainrunids:
            trainselecter.append(i)
    return np.array(trainselecter)

def kmerSta_para(args):
    indir=args[0]
    hashdir=args[1]
    subabunfilelist=args[2]
    trainselecter=args[3]
    readCountarr=args[4]
    trainsampleids=args[5]

    for abunfile in subabunfilelist:
        stafile = abunfile.replace('.abun.npz', '.sta')
        stafile = hashdir + stafile
        flag=checksta(indir + abunfile,stafile)
        if flag==False:
            output = open(stafile, 'w')
            mat = load_npz(indir + abunfile)
            ROW, COL = mat.shape
            print(abunfile)
            print(str(ROW)+','+str(COL))
            newmat = mat.tocsr()[:,trainselecter]
            newreadCountarr=readCountarr[trainselecter]
            for rownum in range(ROW):
                nonzero=newmat.getrow(rownum).count_nonzero()
                if nonzero < 0.005*len(trainsampleids):
                    rownum = rownum + getLowfeat(abunfile)
                    output.write(str(rownum) + ',0,0,0\n')
                else:
                    row = newmat.getrow(rownum)
                    denserow=np.array(row.todense()).flatten()
                    normDenserow=[]
                    for i in range(len(newreadCountarr)):
                        normDenserow.append(round(int(denserow[i]) / newreadCountarr[i]))
                    normDenserow=np.array(normDenserow)
                    std, mean, cnt = calkmerSta(normDenserow)  # 计算是否满足条件
                    rownum = rownum + getLowfeat(abunfile)
                    output.write(str(rownum) + ',' + str(std) + ',' + str(mean) + ',' + str(cnt) + '\n')
###########################################################################################
def mergeDict(dictList):
    mergedDict = {}
    for d in dictList:  # you can list as many input dicts as you want here
        for key, value in d.items():
            if key not in mergedDict.keys():
                mergedDict[key] = 0
            mergedDict[key]=mergedDict[key]+value
    return mergedDict


def selStd(allSelstdcount,selNum):
    count=0
    selstds=[]
    stdList=list(allSelstdcount.keys())
    stdList.sort(reverse=True)
    for std in stdList:
        if count<selNum:
            count=allSelstdcount[std]+count
            selstds.append(std)
        else:
            break
    return selstds,count


def readstafile(filename,samplenum,threshold):
    f=open(filename)
    selStdcount={}
    linecnt=0
    kmerlinecnt=0
    for line in f:
        kmerlinecnt=kmerlinecnt+1
        kmerid,std,mean,cnt=line.strip().split(',')
        std=float(std)
        if int(cnt) >= (threshold*samplenum):
            linecnt = linecnt + 1
            if std not in selStdcount.keys():
                selStdcount[std]=0
            selStdcount[std]=selStdcount[std]+1
    return selStdcount,linecnt,kmerlinecnt


def kmerSel_para(args):
    indir=args[0]
    samplenum=args[1]
    subabunfilelist=args[2]
    threshold=args[3]
    selStdcountlist=[]
    alllinecnt=0
    allkmerlinecnt=0
    for abunfile in subabunfilelist:
        stafile=indir+abunfile.replace('.abun.npz','.sta')
        selStdcount, linecnt,kmerlinecnt=readstafile(stafile,samplenum,threshold)
        selStdcountlist.append(selStdcount)
        alllinecnt=alllinecnt+linecnt
        allkmerlinecnt=allkmerlinecnt+kmerlinecnt
    return selStdcountlist,alllinecnt,allkmerlinecnt

def globSelkmer(hashdir,samplenum,abunfilelist,threshold,top,corenum,outfile,outsta):
    data = []
    subabunfilelist = np.array_split(abunfilelist, corenum)
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [hashdir,samplenum, subabunfilelist[i], threshold]
        data.append(onedata)
    results,alllinecnts,allkmerlinecnts=zip(*pool.map(kmerSel_para, data))
    pool.close()
    pool.join()
    if len(results)!=corenum:
        return 0

    dictList = []
    for dicts in results:
        for dict in dicts:
            dictList.append(dict)
    allSelstdcount = mergeDict(dictList)

    selkmernum=0
    for cnt in alllinecnts:
        selkmernum=cnt+selkmernum
    selNum = round(top * selkmernum)
    selstds, count = selStd(allSelstdcount, selNum)

    allkmer=0
    for cnt in allkmerlinecnts:
        allkmer=cnt+allkmer

    output=open(outfile,'w')
    for selstd in selstds:
        output.write(str(selstd)+'\n')

    outsta=open(outsta,'w')
    outsta.write('all kmers: '+str(allkmer)+'\n')
    outsta.write('train samples: '+str(samplenum)+'\n')
    outsta.write('zero threshold: '+str(threshold)+'\n')
    outsta.write('kmer occurance in more than '+str((threshold*samplenum))+' samples'+'\n')
    outsta.write('all seleted kmers by zero: '+str(selkmernum)+'\n')
    outsta.write('top : '+str(top)+' kmers in selected kmers: '+str(selNum)+'\n')
    outsta.write('acturally seleted: '+str(count) + ' kmers\n')
    outsta.write('smallest std : '+str(selstd)+'\n')

############################################################################################
def getStdsel(filename):
    stdset=set()
    f=open(filename)
    for line in f:
        stdset.add(line.strip())
    return stdset

def genFeat(newmat,rownum,newreadCountarr,newfilelist):
    row = newmat.getrow(rownum)
    denserow = np.array(row.todense()).flatten()
    normDenserow = []
    feat={}
    for i in range(len(newreadCountarr)):
        normDenserow.append(round(int(denserow[i]) / newreadCountarr[i]))
    for val, featname in zip(normDenserow, newfilelist):
        feat[featname] = val
    return feat


def genFeatnonzero(newmat,rownum,newreadCountarr,newfilelist):
    row = newmat.getrow(rownum)
    denserow = np.array(row.todense()).flatten()
    normDenserow = []
    selfilelist=[]
    feat={}
    for i in range(len(newreadCountarr)):
        if denserow[i]!=0:
            normDenserow.append(round(int(denserow[i]) / newreadCountarr[i]))
            selfilelist.append(newfilelist[i])
    for val, featname in zip(normDenserow, selfilelist):
        feat[featname] = val
    return feat

def kmerHashcode_para(args):
    indir=args[0]
    hashdir=args[1]
    subabunfilelist=args[2]
    trainselecter=args[3]
    fileReadcount=args[4]
    trainsampleids=args[5]
    nbit=args[6]
    selstdfile=args[7]
    samplenum = len(trainsampleids)
    zeroper=args[8]
    topper=args[9]
    readCountarr=args[10]
    fileidfile=args[11]
    stdset=getStdsel(selstdfile)
    fileidDict, fileNamedict, filelist = readFileid(fileidfile)

    for abunfile in subabunfilelist:
        abunf = open(indir+abunfile)
        outfile=abunfile.replace('.abun.npz','.hashcode'+str(nbit)+'_' + str(zeroper) + '_' + str(topper))
        if os.path.exists(hashdir+outfile)==False or linecount(hashdir+outfile)==0:
            print(hashdir+outfile)
            stafile = abunfile.replace('.abun.npz', '.sta')
            staf=open(hashdir+stafile)
            output=open(hashdir+outfile,'w')

            mat = load_npz(indir+abunfile)
            ROW, COL = mat.shape
            newmat = mat.tocsr()[:, trainselecter]
            newreadCountarr = readCountarr[trainselecter]
            newfilelist=np.array(filelist)[trainselecter]

            for rownum in range(ROW):
                kmerid=rownum + getLowfeat(abunfile)
                kmeridsta, std, mean, cnt=staf.readline().strip().split(',')
                if kmerid==int(kmeridsta):
                    if std in stdset and int(cnt) >= (zeroper*samplenum):
                        trainfeat = genFeatnonzero(newmat, rownum, newreadCountarr, newfilelist)
                        hashcode = Simhash(trainfeat, f=nbit)
                        output.write(str(kmerid) + '\t' + str(hashcode.value) +'\n')
                else:
                    print(str(kmerid) + ' '+str(kmeridsta)+' not equal')
        else:
            print(outfile+' exists')

def kmerSelection_para(args):
    hashdir = args[0]
    subabunfilelist = args[1]
    selstdfile = args[2]
    samplenum = args[3]
    threshold = args[4]
    kmerlist=[]

    stdset=getStdsel(selstdfile)

    for abunfile in subabunfilelist:
        stafile = abunfile.replace('.abun', '.sta')
        staf=open(hashdir+stafile)
        for staline in staf:
            kmeridsta,std,mean,cnt=staline.strip().split(',')
            if std in stdset and int(cnt) >= (threshold*samplenum):
                kmerlist.append(kmeridsta)
    return kmerlist


#############################################################################################
def getLastline(filename):
    f=open(filename)
    lines=list(f.readlines())
    lastline=lines[-1]
    return lastline

def checkall(allkmerfile,stafile):
    f=open(stafile)
    lines=list(f.readlines())
    lastsecondline=lines[-2]
    tokens=lastsecondline.strip().split(':')
    num,kmer=tokens[1].strip().split(' ')
    allkmer=linecount(allkmerfile)
    if allkmer==int(num):
        print('Passed final exam')
        print(str(allkmer)+' kmers are selected')
        return True,allkmer,num
    else:
        return False,allkmer,num

#############################################################################################
def genoutdir(indir,random,fold):
    if not os.path.exists(indir + 'random' + str(random)):
        os.mkdir(indir + 'random' + str(random))
    if not os.path.exists(indir + 'random' + str(random)+'/fold' + str(fold) + '/'):
        os.mkdir(indir + 'random' + str(random)+'/fold' + str(fold) + '/')
    if not os.path.exists(indir + 'random' + str(random)+'/fold' + str(fold) + '/tmp/'):
        os.mkdir(indir + 'random' + str(random)+'/fold' + str(fold) + '/tmp/')


def main(task,datatype,random,fold,corenum,nbit,zeroper,topper,flag):
    dirinfo='AAAA'
    datadir = dirinfo + task + '/data/'
    #if random==0, inter-dataset experiment
    if task == 'PBMC' and random!=1000:
        datadir = dirinfo+'/PBMC_old' + '/data/'
    elif task=='PBMC' and random==1000:
        datadir=dirinfo+'/PBMCtest/data/'
    elif task == 'PBMC16':
        datadir = dirinfo+'/PBMC_old' + '/data/'

    barcodefile =barcodefiles[task]
    fileidfile = barcodefile.replace('tsv','id')

    if task=='PBMC' and random==1000:
        barcodefile='Cells_labelName_Seurat.csv.inter.barcode'
        fileidfile='Cells_labelName_Seurat.csvinter.barcode.id'
    print(barcodefile)
    print(fileidfile)

    kmeridfiledir = dirinfo + task + '/' + datatype + 'kmerIDfiles_10/'
    readcountfile = 'filereadCount.txt'
    indir = dirinfo +task+'/'+datatype+'kmerMatrix_10/'
    hashdir = indir + 'random' + str(random) + '/fold' + str(fold) +'/'
    #generate output directories
    genoutdir(indir, random, fold)
    if task in PBMCtask:
        if random!=1000:
            trainsampleids,infile=util.getTrainfilelistrandom(random,fold+1)#take a subset of samples as train samples
        else:
            trainsampleids = readFilelist(datadir + barcodefile) # take all samples as train samples
        print('# of trainsampleids: '+str(len(trainsampleids)))

        fileReadcount = readReadCount(kmeridfiledir + readcountfile)
        readCountarr,sampleids = genReadcountVec(fileReadcount, datadir+fileidfile)

    #go through all kmermatrix files, return kmermatrix file list and sampleid list
    abunfilelist=getAbundocs(indir+'fold999/tmp/')

    #generate a trainfile selection array
    print('# of sampleids: '+str(len(sampleids)))
    print('# of trainsampleids: ' + str(len(trainsampleids)))
    trainselecter=genTrainidselector(sampleids,trainsampleids)

    #split all kmer matrixs files into corenum subsets
    subabunfilelist=np.array_split(abunfilelist,corenum)

    # for each k-mer, build a vector of normalized abundences from train samples
    # claculate std and occurrence for each k-mer
    # generate .sta files with multiple threads, check the existence of .sta file, if exists, do not generate
    data=[]
    pool=Pool(corenum)
    for i in range(corenum):
        onedata = [indir+'fold999/tmp/',hashdir+'tmp/',subabunfilelist[i], trainselecter, readCountarr,trainsampleids]
        data.append(onedata)
    pool.map(kmerSta_para,data)
    pool.close()
    pool.join()

#################################################################################################################
    #check whether selstd_* file exists
    selstdfile = 'selstd' + '_' + str(zeroper) + '_' + str(topper)
    selstafile = 'selection.sta' +'_'+ str(zeroper) + '_' + str(topper)
    existflag=False
    if os.path.exists(hashdir+selstdfile)==True and os.path.exists(hashdir+selstafile)==True:
        laststd=getLastline(hashdir+selstdfile).strip()
        line=getLastline(hashdir+selstafile)
        tokens=line.strip().split(':')
        laststaid=tokens[1].strip()
        if laststaid==laststd:
            existflag=True
            print(selstdfile+' and '+selstafile +' exist')
            print('samllest std: '+str(laststd))
    else:
        print(selstafile+' not exist')
        print(selstdfile + ' not exist')

    if (flag=='selstd' or flag=='selhash') and existflag==False:
        globSelkmer(hashdir=hashdir+'tmp/', samplenum=len(trainsampleids),
                    abunfilelist=abunfilelist, threshold=zeroper,  top=topper,
                    corenum=corenum,outfile=hashdir+selstdfile,outsta=hashdir+selstafile)
################################################################################################################
    if flag=='selhash':
        data = []
        pool = Pool(corenum)

        for i in range(corenum):
            onedata = [indir + 'fold999/tmp/', hashdir + 'tmp/', subabunfilelist[i], trainselecter, fileReadcount,
                       trainsampleids, nbit,hashdir+selstdfile,zeroper,topper,readCountarr,datadir+fileidfile]
            data.append(onedata)

        pool.map(kmerHashcode_para, data)
        pool.close()
        pool.join()

        outcodefile = hashdir + str(nbit) + 'bit.codegroup'+ '_' + str(zeroper) + '_' + str(topper)
        outcodefile10 = hashdir + str(nbit) + 'bit.codegroup10'+ '_' + str(zeroper) + '_' + str(topper)
        selCodeset = selCode(nbit,zeroper,topper,hashdir + 'tmp/', abunfilelist, outcodefile, outcodefile10)

        outfile10 = hashdir + str(nbit) + 'bit.codegroup10.kmercode'+ '_' + str(zeroper) + '_' + str(topper)
        outfile = hashdir + str(nbit) + 'bit.codegroup.kmercode'+ '_' + str(zeroper) + '_' + str(topper)
        genKmercodedict(nbit,zeroper,topper,hashdir + 'tmp/', selCodeset, outfile10, outfile)
        outfilesorted = hashdir + str(nbit) + 'bit.codegroup10.kmercode.sorted'+ '_' + str(zeroper) + '_' + str(topper)
        cmd = 'sort -k1 -n ' + outfile10 + ' > ' + outfilesorted
        p = subprocess.Popen(args=cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
        p.wait()
        print(p.stdout.read())
        print(p.stderr.read())

        checkresult, allkmer, num = checkall(outfile, hashdir + selstafile)
        if checkresult==False:
            current_time = datetime.datetime.now()
            outlog = open(dirinfo+'logfile.log', 'a')
            outlog.write("Time now is : " + str(current_time) + '\n')
            outlog.write('genSpskmersta '+task+' '+datatype+' '+str(random)+' '+str(fold)+' '+
                         str(nbit)+' '+str(zeroper)+' '+str(topper)+' '+flag+' : wrong\n')
            outlog.write('all kmer in kmercode: ' + str(allkmer) + '\n')
            outlog.write('kmer num in selection: ' + str(num) + '\n')

if __name__ == "__main__":
    # execute only if run as a script
    #sbatch genKmersta.sh PBMC mapped 5 4 32 0.05 0.3 selhash

    PBMCtask = ['PBMC', 'PBMC16']
    #barcodefiles contains all the barcodes
    barcodefiles = {'PBMC': 'barcodes.tsv', 'PBMC16': 'barcodes.tsv'}

    task=sys.argv[1]#T2D21 T2D
    datatype=sys.argv[2]
    random=int(sys.argv[3])
    fold = int(sys.argv[4])
    nbit=int(sys.argv[5])
    zeroper=float(sys.argv[6])
    topper=float(sys.argv[7])
    flag=sys.argv[8]

    if datatype=='all':
        datatype=''
    corenum=32
    main(task,datatype,random,fold,corenum,nbit,zeroper,topper,flag)


