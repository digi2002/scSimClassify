import numpy as np
import keras
import util


from BREASTfileinfo import readFilereadcount_BREAST,readFilereadcount_BREASTunmapped,\
    readFilereadcount_BREASTmapped,readFilereadcount_BREASTtest,readFilereadcount_BREASTtestunmapped,\
    readFilereadcount_BREASTtestmapped
import math

PBMCtask = ['PBMC', 'PBMC16']


def readData(filename,Samplelabelid):
    print('in readData')
    print(len(Samplelabelid))
    X=[]
    labels=[]
    sampleids=[]
    #arr = np.loadtxt(filename, dtype='str', delimiter='\t')
    arr = np.loadtxt(filename, dtype='str', delimiter=',')
    print('arr shape: '+str(arr.shape[0])+' '+str(arr.shape[1]))
    arrtrans = np.transpose(arr)
    geneName=arrtrans[0][1:]
    print(filename)
    featArr=arrtrans[1:]
    for feat in featArr:
        sampleid=feat[0].strip('\"')
        sampleid=sampleid.replace('.','-')
        if sampleid in Samplelabelid.keys():
            label=Samplelabelid[sampleid]
            feat = [float(i) for i in feat[1:]]
            X.append(feat)
            labels.append(label)
            sampleids.append(sampleid)
    y = keras.utils.to_categorical(labels, len(set(Samplelabelid.values())))
    return np.array(X),y,labels,sampleids




def readKmerData(filename,Samplelabelid):
    print('in readData')
    print(len(Samplelabelid))
    X=[]
    labels=[]
    sampleids=[]
    print(filename)
    #arr = np.loadtxt(filename, dtype='str', delimiter='\t')
    arr = np.loadtxt(filename, dtype='str', delimiter=',')
    #arrtrans = np.transpose(arr)
    geneName=arr[0][1:]
    print(filename)
    featArr=arr[1:]
    for feat in featArr:
        sampleid=feat[0].strip('\"')
        if sampleid in Samplelabelid.keys():
            label=Samplelabelid[sampleid]
            feat = [float(i) for i in feat[1:]]
            X.append(feat)
            labels.append(label)
            sampleids.append(sampleid)
    y = keras.utils.to_categorical(labels, len(set(Samplelabelid.values())))
    return np.array(X),y,labels,sampleids

def getKmerdata(indir, filelist,labelDict):
    X = []
    sampleids = []
    labels=[]
    for filename in filelist:
        infile = indir + filename + '.sample'
        f = open(infile)
        line1 = f.readline()
        line2 = f.readline()
        tokens = line2.strip().split('\t')
        sampleid = tokens[0]
        feats = [float(i) for i in tokens[1:]]
        if sampleid in labelDict.keys():
            X.append(feats)
            sampleids.append(sampleid)
            labels.append(labelDict[sampleid])
    X = np.array(X)
    sampleids = np.array(sampleids)
    y = keras.utils.to_categorical(labels, len(set(labelDict.values())))
    return X, y,labels,sampleids


def getReadcount(datatype):
    if datatype=='':
        fileReadcount = readFilereadcount_BREAST()
    elif datatype=='unmapped':
        fileReadcount = readFilereadcount_BREASTunmapped()
    elif datatype=='mapped':
        fileReadcount = readFilereadcount_BREASTmapped()
    return fileReadcount


def getTestreadcount(datatype):
    if datatype == '':
        fileReadcount = readFilereadcount_BREASTtest()
    elif datatype == 'unmapped':
        fileReadcount = readFilereadcount_BREASTtestunmapped()
    elif datatype == 'mapped':
        fileReadcount = readFilereadcount_BREASTtestmapped()
    return fileReadcount


def convXreads(dirinfo,task,feats,sampleids,datatype):
    feats_new=[]
    if task in PBMCtask:
        kmeridfiledir =dirinfo+str(task)+'/'+datatype+'kmerIDfiles_10/'
        readcountfile = 'filereadCount.txt'
        fileReadcount=util.readReadCount(kmeridfiledir+readcountfile,task)
    else:
        fileReadcount=getReadcount(datatype)
    for feat,sampleid in zip(feats,sampleids):
        feat_new=[]
        for x in feat:
            x=round(x/fileReadcount[sampleid],6)
            feat_new.append(x)
        feats_new.append(np.array(feat_new))
    return np.array(feats_new)

def convXtestreads(dirinfo,task,patient,feats,sampleids,datatype):
    feats_new = []
    if task=='PBMC':
        kmeridfiledir = dirinfo + 'COVID'+patient + '/' + datatype + 'kmerIDfiles_10/'
        readcountfile = 'filereadCount.txt'
        fileReadcount = util.readReadCount(kmeridfiledir + readcountfile, task)
    else:
        fileReadcount = getTestreadcount(datatype)
    for feat, sampleid in zip(feats, sampleids):
        feat_new = []
        for x in feat:
            x = round(x / fileReadcount[sampleid], 6)
            feat_new.append(x)
        feats_new.append(np.array(feat_new))
    return np.array(feats_new)

def convnorm2X(feats):
    feats_new=[]
    for feat in feats:
        norm = np.linalg.norm(feat)
        if norm!=0:
            feats_new.append(feat / norm)
        else:
            feats_new.append(feat)
    return np.array(feats_new)
