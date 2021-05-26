import sys
import multiprocessing
import numpy as np
from multiprocessing import Pool
from util import readFilelist,combineDict,exeCommand_multi,exeCommand
import os
from collections import defaultdict
from genKmersparsematrix import genKmersparsematrix
from infofile import getKMERNUM
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO

def gensplitBarcodeCommandlist(datadir,bamfile,outBarcodedir,filelist):
    commandlist=[]
    for barcode in filelist:
        command=command+'samtools view '+datadir+bamfile+' | grep '+barcode+' > '+ outBarcodedir+barcode+'.sam; '
        command=command+'cat '+datadir + bamfile + '.header '+outBarcodedir+barcode+'.sam > '+outBarcodedir+barcode+'.sam'
        commandlist.append(command)
    return commandlist

def genJellyfishcountCommandlist(filelist,fastadir,kmerdir,k):
    commandlist=[]
    for infile in filelist:
        filename=fastadir+infile+'.fasta'
        outfile=kmerdir+infile+'.kmer'
        command='jellyfish count -m '+str(k)+' -s 100M -t 10 -C -o '+outfile+' '+filename+'\n'
        commandlist.append(command)
    return commandlist

def genJellyfishdumpCommandlist(filelist,kmerdir,kmerdumpdir):
    commandlist=[]
    for infile in filelist:
        filename=kmerdir+infile+'.kmer'
        outfile=kmerdumpdir+infile+'.kmer'
        if os.path.exists(filename) == True:
            command = 'jellyfish dump -c -t -L 5 ' + filename + ' > ' + outfile
            commandlist.append(command)
        else:
            print(filename+' does not exist')
    return commandlist


def combineKmerfile_multi(filelist,corenum,kmerdumpdir,kmerlistdir):
    print('combineKmerfile_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in convertFastq: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [kmerdumpdir, subfilelist[i]]
        data.append(onedata)
    results=pool.map(combineKmerfile, data)
    pool.close()
    allset=set()
    for oneset in results:
        allset=allset.union(oneset)

    outfile=filelist[0]+'_'+filelist[-1]+'.cmb'
    output=open(kmerlistdir+outfile,'w')

    for kmerstr in allset:
        output.write(kmerstr+'\n')


def combineKmerfile(args):
    kmerdumpdir=args[0]
    filelist=args[1]
    kmerset=set()
    for infile in filelist:
        filename=kmerdumpdir+infile+'.kmer'
        f=open(filename)
        for line in f:
            tokens=line.strip().split('\t')
            kmerstr=tokens[0]
            kmerset.add(kmerstr)
    return kmerset


def genKmerIDfile_multi(filelist,corenum,kmerdumpdir,kmeridfiledir,kmerlistdir):
    print('genKmerIDfile_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in convertFastq: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [kmerdumpdir, kmeridfiledir,kmerlistdir,subfilelist[i]]
        data.append(onedata)

    pool.map(genKmerIDfile, data)
    pool.close()


def genKmerIDfile(args):
    kmerdumpdir=args[0]
    kmeridfiledir=args[1]
    kmerlistdir=args[2]
    filelist=args[3]
    kmeridDict=readIDdict(kmerlistdir)
    if kmeridDict==False:
        return False
    for infile in filelist:
        filename=kmerdumpdir+infile+'.kmer'
        outfile=kmeridfiledir+infile+'id'
        if os.path.exists(filename) == True:
            f=open(filename)
            output=open(outfile,'w')
            for line in f:
                kmerstr,cnt=line.strip().split('\t')
                if kmerstr in kmeridDict.keys():
                    kmerid=kmeridDict[kmerstr]
                    output.write(kmerid+'\t'+cnt+'\n')


def calFilereadcount_multi(filelist,corenum,kmeridfiledir,app=''):
    outfile=kmeridfiledir+app+'filereadCount.txt'
    output=open(outfile,'w')
    print('calFilereadcount_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in convertFastq: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)

    for i in range(corenum):
        onedata = [kmeridfiledir,subfilelist[i]]
        data.append(onedata)

    results=pool.map(calFilereadcount, data)
    pool.close()
    for result in results:
        for infile,sum in result.items():
            output.write(infile+'\t'+str(sum)+'\n')

def calFilereadcount(args):
    kmeridfiledir= args[0]
    filelist=args[1]
    fileReadcount={}
    for infile in filelist:
        sum=0
        filename = kmeridfiledir + infile + 'id'
        f=open(filename)
        for line in f:
            kmerid,abun=line.strip().split('\t')
            abun=int(abun)
            sum=sum+abun
        fileReadcount[infile]=sum
    return fileReadcount


def readIDdict(kmerlistdir):
    kmeridDict={}
    filename = kmerlistdir+'allkmerstr.kmerid'
    f=open(filename)
    for line in f:
        kmerstr,kmerid=line.strip().split('\t')
        if kmerstr not in kmeridDict.keys():
            kmeridDict[kmerstr]=kmerid
        else:
            return False
    return kmeridDict



def genKmerid(kmerlistdir):
    filename = kmerlistdir+'allkmerstr'
    outfile=filename+'.kmerid'
    output=open(outfile,'w')
    f=open(filename)
    id=0
    for line in f:
        kmerstr=line.strip()
        output.write(kmerstr+'\t'+str(id)+'\n')
        id=id+1


def genFileid(filelist,outfile):
    output=open(outfile,'w')
    id=0
    for infile in filelist:
        output.write(infile + '\t' + str(id) + '\n')
        id = id + 1



def bamtofasta(datadir,bamfile):
    command = 'samtools fasta ' + datadir + bamfile + ' --threads 10 > ' + datadir + bamfile + '.fasta'
    exeCommand(command)

def checkfasta(args):
    fastadir=args[0]
    fileset = set(args[1])
    for barcode in fileset:
        fastainput = open(fastadir+barcode+'.fasta')

        fastas = SimpleFastaParser(fastainput)
        nameDict = defaultdict(int)
        for (name, seq) in fastas:
            nameDict[name]=nameDict[name]+1

        namelist=[]
        for name,cnt in nameDict.items():
            if cnt>1:
                namelist.append((name,cnt))
        if len(namelist)!=0:
            output = open(fastadir + barcode + '.fasta.sta', 'w')
            for (name,cnt) in namelist:
                output.write(name+','+str(cnt)+'\n')

def checkfasta_multi(filelist,corenum,fastadir):
    print('checkfasta_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in checkfasta: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [fastadir, subfilelist[i]]
        data.append(onedata)
    pool.map(checkfasta, data)
    pool.close()

def genBCQnamefile_multi(filelist,corenum,outBarcodedir,bcqnamefile):
    print('genBCQnamefile_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in genBCQnamefile: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [outBarcodedir, subfilelist[i], bcqnamefile]
        data.append(onedata)
    pool.map(genBCQnamefile, data)
    pool.close()

def genBCQnamefile(args):
    outBarcodedir=args[0]
    fileset=set(args[1])
    bcqnamefile=args[2]
    barcodeDict=defaultdict(list)
    f=open(bcqnamefile)
    for line in f:
        if line.find('CB:Z:')!=-1:
            qname, barcode = line.strip().split('\t')
            barcode = barcode.replace('CB:Z:', '')
            if barcode in fileset:
                barcodeDict[barcode].append(qname)
    for barcode,qnames in barcodeDict.items():
        output=open(outBarcodedir+barcode,'w')
        for qname in qnames:
            output.write(qname+'\n')


def genBCFastafile_multi(filelist,corenum,outfastadir,barcodedir,fastafile):
    print('genBCFastafile_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in genBCFastafile: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [outfastadir, barcodedir,subfilelist[i], fastafile]
        data.append(onedata)
    pool.map(genBCFastafile, data)
    pool.close()

def genSRRBCFastafile_multi(filelist,corenum,outfastadir,fastafile):
    print('genSRRBCFastafile_multi: ' + str(len(filelist)))
    subfilelist = np.array_split(filelist, corenum)
    print('number of subcommandlist in genSRRBCFastafile: ' + str(len(subfilelist)))
    data = []
    pool = Pool(corenum)
    for i in range(corenum):
        onedata = [outfastadir, subfilelist[i], fastafile]
        data.append(onedata)
    pool.map(genSRRfastafile, data)
    pool.close()



#read qname into dictionary with barcodes as values
#Note: duplicated qname will be eliminated
def readBCQname(dir,barcode):
    f=open(dir+barcode)
    qnameDict = defaultdict(str)
    for line in f:
        qnameDict[line.strip()]=barcode
    return qnameDict


def barcodelistTrim(barcodelist):
    newbarcodelist=[]
    for barcodenum in barcodelist:
        barcode,num=barcodenum.split('-')
        newbarcodelist.append(barcode)
    return newbarcodelist,str(num)

def genSRRfastafile(args):
    outfastadir = args[0]
    barcodelist = args[1]
    barcodelist,num=barcodelistTrim(barcodelist)
    fastafile = args[2]
    BCfastaDict = defaultdict(list)
    fastainput1 = open(fastafile+'_1.fasta')
    fastainput2 = open(fastafile+'_2.fasta')

    fastas1 = SimpleFastaParser(fastainput1)
    fastas2 = SimpleFastaParser(fastainput2)

    for (name1, seq1), (name2, seq2) in zip(fastas1, fastas2):
        name1=name1.split(' ')[0]
        name2 = name2.split(' ')[0]
        if name1==name2:
            barcode=seq1[:16]
            if barcode in barcodelist:
                BCfastaDict[barcode].append((name2, seq2))
        else:
            print(name1+'__'+name2)

    for barcode, fastapair in BCfastaDict.items():
        output = open(outfastadir + barcode +'-'+str(num)+ '.fasta', 'w')
        for (name, seq) in fastapair:
            output.write('>' + name + '\n')
            output.write(seq + '\n')
import os
def renamefasta(subfiles, fastadir):
    newbarcodelist=barcodelistTrim(subfiles)
    for barcode in newbarcodelist:
        os.rename(fastadir+barcode+'.fasta', fastadir+barcode+'-1.fasta')


def genBCFastafile(args):
    outfastadir=args[0]
    barcodedir=args[1]
    fileset=set(args[2])
    fastafile=args[3]
    BCfastaDict=defaultdict(list)
    qnameDicts=[]
    for barcode in fileset:
        qnameDicts.append(readBCQname(barcodedir,barcode))
    qnameDict=combineDict(qnameDicts)
    if qnameDict==False:
        return 0
    fastainput=open(fastafile)
    fastas=SimpleFastaParser(fastainput)
    for (name, seq) in fastas:
        if name in qnameDict.keys():
            barcode = qnameDict[name]
            BCfastaDict[barcode].append((name, seq))
    for barcode,fastapair in BCfastaDict.items():
        output=open(outfastadir+barcode+'.fasta','w')
        for (name,seq) in fastapair:
            output.write('>'+name+'\n')
            output.write(seq+'\n')


def prorandomsel(filename):
    outfile=filename.replace('label','tsv')
    outid=filename.replace('label','id')
    output=open(outfile,'w')
    outputid=open(outid,'w')
    filelist=[]
    f=open(filename)
    cnt=0
    for line in f:
        tokens=line.strip().split(',')
        filelist.append(tokens[0])
        output.write(tokens[0]+'\n')
        outputid.write(tokens[0]+'\t'+str(cnt)+'\n')
        cnt=cnt+1
    return filelist


originfiles={'COVID1C':'SRR11680207','PBMC16':'','PBMCtest':'','COVIDC4':'SRR11680216','COVIDC2':'SRR11680208',
             'COVIDN2': 'SRR11680219','PBMC':'','COVIDF4':'SRR11680213','COVIDF2':'SRR11680210','COVIDN4':'SRR11680225'}
barcodefiles={'PBMC':'barcodes.tsv','PBMC16':'barcodes.tsv',
             'COVID1C':'count_w_metadata.C1.barcode.tsv','PBMCtest':'',
              'COVIDN2': 'count_w_metadata.N2.barcode.tsv',
              'COVIDF4':'count_w_metadata.F4.barcode.tsv',
              'COVIDN4':'count_w_metadata.N4.barcode.tsv',
              'COVIDF2':'count_w_metadata.F2.barcode.tsv',
              'COVIDC2':'count_w_metadata.C2.barcode.tsv'}

def main(step,dir,taskdir,datatype,nodenum,docNum,coreNum):
    datadir=taskdir+'/data/'
    barcodefile=barcodefiles[TASKDIR]
    bamfile=originfiles[TASKDIR]
    barcodedir=taskdir+'/barcodefiles/'
    fastqdir=taskdir+'/fastqall/'
    if datatype=='all': datatype=''
    fastadir = taskdir + '/' + datatype + 'fasta/'

    kmerdir = taskdir + '/' + datatype + 'kmer/'
    kmerdumpdir = taskdir + '/' + datatype + 'kmer_dump/'
    kmerlistdir = taskdir + '/' + datatype + 'kmerList/'
    kmeridfiledir = taskdir + '/' + datatype + 'kmerIDfiles_10/'
    kmermatrixdir = taskdir + '/' + datatype + 'kmerMatrix_10/'

    if TASKDIR == 'PBMC' or TASKDIR=='PBMC16':
        barcodedir = dir+'PBMC_old' + '/barcodefiles/'
        datadir = dir+'PBMC_old' + '/data/'
        fastadir = dir+'PBMC_old' + '/' + datatype + 'fasta/'
        if TASKDIR=='PBMC16':
            kmerdir = dir+'PBMC16_old' + '/' + datatype + 'kmer/'
        elif TASKDIR=='PBMC':
            kmerdir = dir+'/PBMC_old' + '/' + datatype + 'kmer/'

    elif TASKDIR=='PBMCtest':
        datadir=dir+'PBMCtest/data/'
        barcodefile='count_w_metadata.'+patient+'.barcode.label.csv.inter.barcode'
        kmeridfiledir = dir + '/' + 'COVID'+patient+'/'+datatype + 'kmerIDfiles_10/'
        kmerdumpdir = dir + '/' + 'COVID'+patient+'/' + datatype+'kmer_dump/'
        kmerlistdir = dir + '/' + 'PBMC/'+ datatype + 'kmerList/'

    filelist=readFilelist(datadir+barcodefile)
    subfilelist=np.array_split(filelist,nodenum)
    print('number of barcodes: '+str(len(filelist)))
    print('number of subsets: '+str(len(subfilelist)))
    subfiles=subfilelist[docNum]
    print('number of barcodes in subsets: '+str(len(subfiles)))

    #taskdir='/scratch/qsu226/PBMC16'#generate 16-mer

    if step=='genbcqname':
        bcqnamefile=datadir+bamfile+'.namebarcode'
        genBCQnamefile_multi(subfiles, coreNum, barcodedir, bcqnamefile)

    elif step=='genbcfasta':
        fastafile=datadir+bamfile+datatype+'.fasta'
        genBCFastafile_multi(subfiles, coreNum, fastadir, barcodedir, fastafile)

    elif step=='genSRRbcfasta':
        fastafile = fastqdir + bamfile
        genSRRBCFastafile_multi(subfiles, coreNum, fastadir, fastafile)

    elif step=='checkfasta':
        checkfasta_multi(subfiles, coreNum, fastadir)

    elif step=='genkmer':
        commandlist = genJellyfishcountCommandlist(subfiles,fastadir,kmerdir,21)
        exeCommand_multi(commandlist, coreNum)

    elif step=='dumpkmer':
        commandlist=genJellyfishdumpCommandlist(subfiles, kmerdir, kmerdumpdir)
        exeCommand_multi(commandlist, coreNum)

    elif step=='combinekmer':
        combineKmerfile_multi(subfiles, coreNum, kmerdumpdir, kmerlistdir)

    elif step=='genkmerid':
        genKmerid(kmerlistdir)

    elif step=='genInterkmeridfile' and TASKDIR=='PBMCtest':
        genKmerIDfile_multi(subfiles, coreNum, kmerdumpdir, kmeridfiledir, kmerlistdir)

    elif step =='genkmeridfile':
        genKmerIDfile_multi(subfiles, coreNum, kmerdumpdir, kmeridfiledir, kmerlistdir)

    elif step=='genreadcount':
        calFilereadcount_multi(filelist, coreNum, kmeridfiledir)

    elif step=='genfileid':
         genFileid(filelist, datadir + barcodefile.replace('tsv','id'))

    elif step=='genkmermatrix':
        kmermatrixdir=kmermatrixdir+'fold999/tmp/'
        fileidfile = datadir + barcodefile.replace('tsv','id')
        KMERNUM=getKMERNUM(datatype,taskdir)
        genKmersparsematrix(KMERNUM, filelist, kmeridfiledir, kmermatrixdir, nodenum, docNum, coreNum,fileidfile)


if __name__ == "__main__":
    # execute only if run as a script
    #TASKDIR = 'PBMCtest'
    TASKDIR='PBMCtest'
    patient='C2'
    dir= 'AAAA'
    taskdir = dir + TASKDIR
    step=sys.argv[1]
    datatype=sys.argv[2]
    nodenum=int(sys.argv[3])
    docNum=int(sys.argv[4])
    coreNum=int(sys.argv[5])

    main(step,dir,taskdir,datatype,nodenum,docNum,coreNum)
