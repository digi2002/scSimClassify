import glob
import os
from multiprocessing import Pool
import numpy as np
from scipy.sparse import load_npz


#jellyfish count -F 2 -m 21 -s 100M -t 10 -C -o ../testresult **-L 10** SRR341634_1.fastq SRR341634_2.fastq
def genCommandlist(outKmerdir,outHistodir,outSRRdir,filelist):
    commandlist=[]
    for runid in filelist:
        filename1=outSRRdir+runid+'_1.fastq'
        filename2=outSRRdir+runid+'_2.fastq'
        outfile=outKmerdir+runid+'.kmer'
        outhishto=outHistodir+runid+'.histo'
        command='jellyfish count -F 2 -m 12 -s 100M -t 10 -C -o '+outfile+' '+filename1+' '+filename2+'\n'
        command=command+'jellyfish histo '+outfile+' -o '+outhishto
        commandlist.append(command)
    return commandlist


def convertFastq_old(filelist):
    outlog = open(taskdir + '/log.txt', 'a+')
    outBarcodedir = taskdir + '/barcodefiles/'
    outfastadir = taskdir + '/fastafiles/'
    for filename in filelist:
        readmapped={}
        duplist=[]
        reads=[]
        readids=[]
        flagset = set()
        outfile = filename + '.fasta'
        f=open(outBarcodedir+filename+'.sam')
        outfasta=open(outfastadir+outfile,'w')
        for line in f:
            tokens=line.strip().split('\t')
            readid=tokens[0]
            flag=tokens[1]
            flagset.add(flag)
            read=tokens[9]
            reads.append(read)
            readids.append(readid)
            dupline = readid
            duplist.append(dupline)

            if dupline not in readmapped.keys():
                readmapped[dupline]=[]
            readmapped[dupline].append(read)

            barcode=tokens[18].replace('CB:Z:','')
            if barcode==filename and len(read)==98:
                outfasta.write('>'+flag+'\n')
                outfasta.write(read+'\n')
            else:
                outlog.write('wrong_'+filename+'_'+barcode+'_'+read+'\n')

        for dupline,flags in readmapped.items():
            if len(set(flags))!=1:
                flags=set(flags)
                flags=list(flags)
                outlog.write(filename+'_'+dupline+'_'+','.join(flags)+'\n')
        outlog.write('\n')


def convertFastq_old(filelist):
    outlog = open(taskdir + '/log.txt', 'a+')
    outBarcodedir = taskdir + '/barcodefiles/'
    outfastadir = taskdir + '/fastafiles/'
    for filename in filelist:
        readidseq={}
        flagset=set()
        outfile = filename + '.fasta'
        f=open(outBarcodedir+filename+'.sam')
        outfasta=open(outfastadir+outfile,'w')
        for line in f:
            tokens=line.strip().split('\t')
            readid=tokens[0]
            flag=tokens[1]
            flagset.add(flag)
            read=tokens[9]
            barcode = tokens[18].replace('CB:Z:', '')
            if readid not in readidseq.keys():
                readidseq[readid]=list()
            readidseq[readid].append((read,barcode,flag))

        readidseqnew={}
        for readid,readidinfos in readidseq.items():
            readdict={}
            for read,barcode,flag in readidinfos:
                if (read,barcode) not in readdict.keys():
                    readdict[(read,barcode)]=[]
                readdict[(read,barcode)].append(flag)
            readidseqnew[readid]=readdict
        for readid,readdict in readidseqnew.items():
            for (read,barcode),flags in readdict.items():
                flagstr='_'.join(flags)
                if barcode==filename and len(read)==98:
                    outfasta.write('>'+flagstr+'\n')
                    outfasta.write(read+'\n')
                else:
                    outlog.write('wrong_'+filename+'_'+barcode+'_'+read+'\n')
        outlog.write('\n')

def linecount(filename):
    count = 0
    for line in open(filename).readlines(): count += 1
    return count


def linecountpnz(filename):
    mat=load_npz(filename)
    ROW, COL = mat.shape
    return ROW


#############################################
def readFilelist(filename):
    filelist=[]
    f=open(filename)
    for line in f:
        filelist.append(line.strip())
    return filelist

def readFileDict(filename):
    fileDict={}
    f=open(filename)
    for line in f:
        fileDict[line.strip()]=line.strip()
    return fileDict


def removeAllfiles(dir):
    filelist=[]
    for filename in glob.iglob(dir + '*', recursive=True):  # attention
        filelist.append(filename)
    for filename in filelist:
        os.remove(filename)

def exeCommand_multi(commandlist,corenum):
    print('exeCommand_multi: '+str(len(commandlist)))
    subcommandlist = np.array_split(commandlist, corenum)
    print('number of subcommandlist in exeCommand: '+str(len(subcommandlist)))
    pool = Pool(corenum)
    pool.map(exeCommand, subcommandlist)
    pool.close()

def exeCommand(commandlist):
    if  isinstance(commandlist, str):
        os.system(commandlist)
    else:
        for command in commandlist:
            os.system(command)

def readReadCount(filename,task=''):
    print('readCount file: '+str(filename))
    fileReadcount = {}
    f = open(filename)
    for line in f:
        infile, sum = line.strip().split('\t')
        sum=int(sum)
        scaler = round(sum / 1000000, 3)
        if scaler < 1 or scaler > 10:
            print(infile+' out of range: '+str(scaler))
        fileReadcount[infile] = scaler
    print('Finish reading file readcount: '+str(len(fileReadcount)))
    return fileReadcount


def genReadcountVec(fileReadcount,fileidfile):
    readCountarr=[]
    fileidDict, fileNamedict,filelist = readFileid(fileidfile)
    for i in range(len(filelist)):
        infile=fileNamedict[i]
        readCountarr.append(fileReadcount[infile])
    return np.array(readCountarr),filelist

def readFileid(filename):
    f=open(filename)
    fileIDdict={}
    fileNamedict={}
    filelist=[]
    for line in f:
        infile, id = line.strip().split('\t')
        fileIDdict[infile]=int(id)
        fileNamedict[int(id)]=infile
        filelist.append(infile)
    return fileIDdict,fileNamedict,filelist


def combineDict(dictlist):
    newdict={}
    for onedict in dictlist:
        for onekey,onevalue in onedict.items():
            if onekey in newdict.keys():
                print(onekey+' has multiple values')
                return False
            newdict[onekey]=onevalue
    return newdict

METADATA='CCCC/PBMC/sharefiles/Cells_labelName_Seurat.txt'
def getTrainfilelistrandom(random,fold,infile=METADATA):
    if fold==999:
        infile=METADATA
    elif fold!=0:
        infile=infile.replace('.txt','_random'+str(random)+'_train'+str(fold)+'.txt')
    f=open (infile)
    filelist=[]
    for line in f:
        infile,label=line.strip().split(',')
        filelist.append(infile)
    print('Finish getting filelist: '+str (len(filelist)))
    return filelist,infile

