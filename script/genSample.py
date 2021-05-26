"""
	Function: generate CKG features for samples
	Input:
        task=sys.argv[1](BREAST,BREAST16)
        datatype=sys.argv[2] (all,mapped)
        random=int(sys.argv[3]) (0,1,2,3,4,5,1000(for BREASTtest,fold=0))
        fold = int(sys.argv[4]) (0,1,2,3,4)
        nbit=int(sys.argv[5]) (32,16)
        zeroper=float(sys.argv[6])(0.1)
        topper=float(sys.argv[7])(0.05)
        filter=sys.argv[8](nofilter)

    Input Description:
        task,datatype,random,fold: samples
        task,datatype,random,fold,nbit,zeroper,topper: a k-mer group look-up table
        filter: no use here

	Author: Qi Sun
"""

import sys
import multiprocessing
from multiprocessing import Pool
import numpy as np
import glob
import os
import subprocess
from genSpskmersta import checkall
from calSample import calfilesCodemean,genCodelist,genSamplematrix
from util import readFileDict
import util


def readDict(infile):
    f=open(infile)
    kmerDict={}
    for line in f:
        tokens=line.strip().split('\t')
        if len(tokens)==2:
            kmer=tokens[0]
            index=tokens[1]
            if kmer not in kmerDict.keys():
                kmerDict[kmer] = int(index)
            else:
                print(kmer + ' exists more than once')
                return 0
    print('# of kmers in dict: '+str(len(kmerDict)))
    return kmerDict

def removeAllfiles(dir):
    filelist=[]
    for filename in glob.iglob(dir + '*', recursive=True):
        filelist.append(filename)
    for filename in filelist:
        os.remove(filename)

def removesplitfiles(dir):
    filelist=[]
    for filename in glob.iglob(dir + 'split*', recursive=True):
        filelist.append(filename)
    for filename in filelist:
        os.remove(filename)

def linecount(filename):
    count = 0
    for line in open(filename).readlines(): count += 1
    return count

def splitfile(codedir,codefile):
    subdictlist=[]
    os.chdir(codedir)
    prefix='split'
    removesplitfiles(codedir)
    p = subprocess.Popen(args='split -l 2000000 ' + codefile +' \''+prefix+'\'', stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
    p.wait()
    print(p.stdout.read())
    print(p.stderr.read())
    for filename in glob.iglob(codedir + 'split' + '*', recursive=True):  # attention
        subdictlist.append(filename)
    return subdictlist

def checksubdictlist(codefile,subdictlist):
    subdictcnt=0
    for subdictfile in subdictlist:
        subdictcnt=subdictcnt+linecount(subdictfile)
    dictcnt=linecount(codefile)
    if subdictcnt==dictcnt:
        return True
    else:
        return False


def getSampleidfiles(indir,task):
    idfilelist = []
    if task in PBMCtask:
        filestr='*id'
    else:
        filestr='SRR*id'
    for filename in glob.iglob(indir + filestr, recursive=True):  # attention
        idfile=filename.replace(indir,'')
        idfilelist.append(idfile)
    return idfilelist

def convfile_para(args):
    subdict=args[0]
    kmerCodedict=readDict(subdict)
    indir=args[1]
    subfiles=args[2]
    dictnum=args[3]
    outdir=args[4]

    for idfile in subfiles:
        codefile=idfile.replace('id','code')
        outfile=outdir+codefile+str(dictnum)
        f=open(indir+idfile)
        output=open(outfile,'w')
        for line in f:
            kmer,abun=line.strip().split('\t')
            if kmer in kmerCodedict.keys():
                codeid=kmerCodedict[kmer]
                line=str(codeid)+'\t'+abun+'\n'
                output.write(line)


def convKmer2code(subdict,dictnum,idfilelist,indir,outdir,corenum):
    subidfilelsit=np.array_split(idfilelist,corenum)
    print('len of subidfilelsit: '+str(len(subidfilelsit)))
    pool = Pool(corenum)
    data=[]
    outdir = outdir + str(dictnum) + '/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        removeAllfiles(outdir)
    for i in range(corenum):
        onedata = [subdict,indir,subidfilelsit[i],dictnum,outdir]
        data.append(tuple(onedata))
    pool.map(convfile_para, data)
    pool.close()
    pool.join()


def checkComplete(idfilelist,dictnum,outtmpdir):
    for idfile in idfilelist:
        for i in range(dictnum):
            codefile = idfile.replace('id', 'code')
            subcodefile = outtmpdir+str(i)+'/' + codefile + str(i)
            if not os.path.exists(subcodefile):
                print(subcodefile+' does not exist')
                return False
    return True

def combineCodefile(outtmpdir,outmergedir,idfilelist,dictnum,corenum):
    removeAllfiles(outmergedir)
    commandlist=[]
    for filename in idfilelist:
        command='cat '
        for i in range(dictnum):
            subcodefile=outtmpdir+str(i)+'/'+filename.replace('id','code')+str(i)
            command=command+subcodefile+' '
        command=command+'> '+ outmergedir+filename.replace('id','code')
        commandlist.append(command)

    if len(commandlist)<corenum:
        corenum=len(commandlist)
    subcommandlist = np.array_split(commandlist, corenum)
    data = []
    pool = Pool(corenum)
    pool.map(process, subcommandlist)
    pool.close()
    pool.join()


def process(commandlist):
    for command in commandlist:
        os.system(command)


def genoutdir(outdir,random,fold,suffix):
    if not os.path.exists(outdir + 'random' + str(random)):
        os.mkdir(outdir + 'random' + str(random))
    if not os.path.exists(outdir + 'random' + str(random)+'/fold' + str(fold) + '/'):
        os.mkdir(outdir + 'random' + str(random)+'/fold' + str(fold) + '/')
    if not os.path.exists(outdir + 'random' + str(random)+'/fold' + str(fold) + '/tmp'+suffix):
        os.mkdir(outdir + 'random' + str(random)+'/fold' + str(fold) + '/tmp'+suffix)
    if not os.path.exists(outdir + 'random' + str(random)+'/fold' + str(fold) + '/mergeCodeabunfile'+suffix):
        os.mkdir(outdir + 'random' + str(random)+'/fold' + str(fold) + '/mergeCodeabunfile'+suffix)
    if not os.path.exists(outdir + 'random' + str(random)+'/fold' + str(fold) + '/codeMeanfile'+suffix):
        os.mkdir(outdir + 'random' + str(random)+'/fold' + str(fold) + '/codeMeanfile'+suffix)
    if not os.path.exists(outdir + 'random' + str(random)+'/fold' + str(fold) + '/sampleMatrix'+suffix):
        os.mkdir(outdir + 'random' + str(random)+'/fold' + str(fold) + '/sampleMatrix'+suffix)
    if not os.path.exists(outdir + 'random' + str(random)+'/fold' + str(fold) + '/sampleMatrix'+suffix+'tmp/'):
        os.mkdir(outdir + 'random' + str(random)+'/fold' + str(fold) + '/sampleMatrix'+suffix+'tmp/')


def getKmersel(filename):
    kmerset=set()
    f=open(filename)
    for line in f:
        kmerset.add(line.strip())
    return kmerset

def checkkmercode(dirinfo,nbit,zeroper,topper,datatype,random,fold):
    import datetime
    current_time = datetime.datetime.now()
    indir = dirinfo + task + '/' + datatype + 'kmerMatrix_10/'
    hashdir = indir + 'random' + str(random) + '/fold' + str(fold) + '/'
    outfile = hashdir + str(nbit) + 'bit.codegroup.kmercode' + '_' + str(zeroper) + '_' + str(topper)
    selstafile = 'selection.sta' + '_' + str(zeroper) + '_' + str(topper)
    checkresult, allkmer, num=checkall(outfile, hashdir + selstafile)
    if  checkresult== False:
        outlog = open(dirinfo+'logfile.log', 'a')
        outlog.write("Time now is : " + str(current_time)+'\n')
        outlog.write('in genCodefile genKmersta ' + task + ' ' + datatype + ' ' + str(random) + ' ' + str(fold) + ' ' +
                     str(nbit) + ' ' + str(zeroper) + ' ' + str(topper) + ' : wrong\n')
        outlog.write('all kmer in kmercode: '+str(allkmer)+'\n')
        outlog.write('kmer num in selection: '+str(num)+'\n')
        return False
    else:
        return True

def combinefiles(sampledir,filelist,outfile):
    output=open(outfile,'w')
    for infile in filelist:
        filename=sampledir+infile+'.sample'
        f=open(filename)
        f.readline()
        line=f.readline().strip()
        output.write(line+'\n')


def genidfilelist(filename):
    idfilelist=[]
    f=open(filename)
    for line in f:
        idfilelist.append(line.strip()+'id')
    return idfilelist

def main(task,nbit,zeroper,topper,datatype,random,fold,corenum,filter):
    dirinfo='AAAA'
    if task!='PBMCtest':
        suffix = '_' + str(zeroper) + '_' + str(topper) + '_' + str(nbit)+'/'
        status=checkkmercode(dirinfo, nbit, zeroper, topper, datatype, random, fold)
        if status==False:
            return 0
    else:
        suffix='_'+str(zeroper)+'_'+str(topper)+'_'+str(nbit)+str(filter)+'/'

    if task=='PBMCtest':
        indir = dirinfo + 'COVID'+patient + '/' + datatype + 'kmerIDfiles_10/'
    elif task in PBMCtask:
        indir = dirinfo + task + '/' + datatype + 'kmerIDfiles_10/'
    else:
        indir = dirinfo + task + '/' + datatype + 'kmerIDfiles_10/mergeKmeridfile/'

    if task=='PBMCtest':
        codetask='PBMC'
    else:
        codetask=task

    codedir = dirinfo+ codetask + '/' + datatype + 'kmerMatrix_10/' + 'random' + str(random) + '/fold' + str(fold) + '/'
    kmercodefile=codedir+str(nbit)+'bit.codegroup10.kmercode'+'_' + str(zeroper) + '_' + str(topper)
    print('kmercode file : '+kmercodefile)

    outdir = dirinfo + task + '/'+patient + datatype + 'sampleMatrix_10/'
    print(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    genoutdir(outdir, random, fold,suffix)

    if task=='PBMC' and random==1000:
        datadir = dirinfo + 'PBMCtest' + '/data/'
        barcodefile='Cells_labelName_Seurat.csv.inter.barcode'

    elif task == 'PBMC' or task=='PBMC16':
        datadir = dirinfo+'/PBMC_old' + '/data/'
        barcodefile = barcodefiles[task]
    elif task=='PBMCtest':
        datadir = dirinfo + task + '/data/'
        barcodefile = 'count_w_metadata.' + patient + '.barcode.label.csv.inter.barcode'  # N2
    else:
        datadir = dirinfo + task + '/data/'
        barcodefile = barcodefiles[task]

    idfilelist=genidfilelist(datadir+barcodefile)

    print('len of idfilelist: '+str(len(idfilelist)))

    subdictlist=splitfile(codedir, kmercodefile)
    if checksubdictlist(kmercodefile, subdictlist)==False:
        print('dict num wrong')
        return 0

    outtmpdir = dirinfo + task + '/' + datatype + 'sampleMatrix_10/' \
                + 'random' + str(random) + '/fold' + str(fold) + '/tmp'+suffix
    dictnum=len(subdictlist)

    for i in range(dictnum):
        print('dict num in main: '+str(i))
        subdict=subdictlist[i]
        convKmer2code(subdict, i, idfilelist, indir, outtmpdir, corenum)
    if checkComplete(idfilelist,dictnum,outtmpdir)==False:
        print('uncomplete subcodefile')
        return 0

    outmergedir = outdir + 'random' + str(random) + '/fold' + str(fold) +'/mergeCodeabunfile'+suffix
    combineCodefile(outtmpdir, outmergedir, idfilelist, dictnum, corenum)
    print('Finish conv kmer to code')
    print(outmergedir)
    ################################################################################################################
    threshold = 2
    # 计算meanfile
    abundir = outdir + 'random' + str(random) + '/fold' + str(fold) + '/mergeCodeabunfile' + suffix
    outmeandir = outdir + 'random' + str(random) + '/fold' + str(fold) + '/codeMeanfile' + suffix
    calfilesCodemean(idfilelist, corenum, threshold, abundir, outmeandir)

    print('Finish cal codemean')
    print(outmeandir)

    if task == 'PBMC' and random == 1000:
        datadir = dirinfo + 'PBMCtest' + '/data/'
        barcodefile = 'Cells_labelName_Seurat.csv.inter.barcode'

    elif task == 'PBMC' or task=='PBMC16':
        datadir = dirinfo+'/PBMC_old' + '/data/'
        barcodefile = barcodefiles[task]
    elif task=='PBMCtest':
        datadir = dirinfo + task + '/data/'
        barcodefile = 'count_w_metadata.'+patient+'.barcode.label.csv.inter.barcode'#N2
    else:
        datadir = dirinfo + task + '/data/'
        barcodefile = barcodefiles[task]
    fileSampledict = readFileDict(datadir + barcodefile)

    print('all samples: '+datadir+barcodefile)
    print('len of fileSampledict: '+str(len(fileSampledict)))

    codefile = codedir + str(nbit) + 'bit.codegroup10' + '_' + str(zeroper) + '_' + str(topper)
    print('codefile : '+codefile)
    codelist = genCodelist(codefile)
    outsampledir = outdir + 'random' + str(random) + '/fold' + str(fold) + '/sampleMatrix' + suffix + 'tmp/'
    genSamplematrix(idfilelist, codelist, corenum, abundir, outmeandir, outsampledir, fileSampledict)

    print('Finish gen samples')
    print(outsampledir)
    filelist = util.readFilelist(datadir + barcodefile)
    sampledir = outdir + 'random' + str(random) + '/fold' + str(fold) + '/sampleMatrix' + suffix + 'tmp/'
    outfile = outdir + 'random' + str(random) + '/fold' + str(fold) + '/sampleMatrix' + suffix + task + suffix.replace(
        '/', '') + '.txt'
    combinefiles(sampledir, filelist, outfile)


if __name__ == "__main__":
    # execute only if run as a script
    PBMCtask = ['PBMC', 'PBMC16']
    barcodefiles = {'PBMC': 'barcodes.tsv', 'PBMC16': 'barcodes.tsv'}
    task=sys.argv[1]
    datatype=sys.argv[2]
    random=int(sys.argv[3])
    fold = int(sys.argv[4])
    nbit = int(sys.argv[5])
    zeroper = float(sys.argv[6])
    topper = float(sys.argv[7])
    filter=sys.argv[8]
    print(filter)
    corenum=32
    if datatype=='all':
        datatype=''
    if filter!='nofilter':
        filter=''
    if task=='PBMC':
        patient=''
    elif task=='PBMCtest':
        patient = 'C2'

    main(task, nbit, zeroper, topper, datatype, random, fold, corenum,filter)



