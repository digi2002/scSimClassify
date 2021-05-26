"""
	Function: cell type classification with four classification methods
	Input:
	    logfile = sys.argv[1] results' filename
        settype = sys.argv[2] ('intra','interdev')
        task=sys.argv[3] ('BREAST')
        datatype = sys.argv[4] ('all','mapped')
        random=int(sys.argv[5]) (0,1,2,3,4) 5 for intra best parameter seeking
        fold = int(sys.argv[6]) (0,1,2,3,4)
        nbit=int(sys.argv[7]) (16,32)
        feattype=sys.argv[8] ('kmer','gene')
        method = sys.argv[9]  ('mlp','rf')
        maxDepth = sys.argv[10] (mlp: layer, rf: maxDepth)
        maxFeat = sys.argv[11] ('sqrt','log','None')
        nEst = int(sys.argv[12])(mlp: dropout, rf: nEst)
    Input Description:
        logfile: to record classification results
        task,datatype,random,fold,nbit: to locate training/test file directory for intra-dataset
        for inter-dataset: file locations are fixed /random1000/fold0, training/test are selected by random, fold
        feattype: compressed features or gene
        method: classifier
        maxDepth,maxFeat,nEst: parameters
	Author: Qi Sun
"""
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from readCellfeat import readData,getKmerdata,convXreads,convXtestreads,convnorm2X
from readCelllabel import getInterDevLabel,getInterDevLabelPBMC,getLabelinfo,getIntraLabelRandom
from outputResult import parse,genlogfile,genInfostr,outlogfile,outrecord,outSingleresult,outmodel,openlogfile

import sys,os,math,pickle
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from mlpclassification import mlpmodel
from sklearn.metrics import f1_score,precision_score,recall_score
from sklearn.metrics import accuracy_score
from BREASTfileinfo import getFilelist_BREAST,getFilelist_BREASTtest
from sklearn.svm import NuSVC
from util import readFilelist
import keras
from mlpclassification import convertPrediction

def classify(modeldir,outdir, modelname,outfilename,maxDepth, maxFeat, nEst, X_train, X_test, y_train, y_test, label_train, label_test, method,
             trainSamplelabelid, testSamplelabelid, label2id):
    output = open(outdir + outfilename,'w')
    print('modelname: '+modeldir + modelname)
    if method == 'rf' or method=='gbm':
        if os.path.exists(modeldir + modelname)==True:
            clf=pickle.load(open(modeldir + modelname, 'rb'))
            print('loaded model')
        else:
            if maxFeat=='None':
                maxFeat=len(X_train[0])
            if method=='rf':
                clf = RandomForestClassifier(max_depth=maxDepth, max_features=maxFeat, n_estimators=nEst,random_state=0)
            elif method=='gbm':
                clf = GradientBoostingClassifier(max_depth=maxDepth, max_features=maxFeat, n_estimators=nEst,
                                                 random_state=0)
            clf.fit(X_train, label_train)
            pickle.dump(clf, open(modeldir + modelname, 'wb'))
        acc = clf.score(X_test, label_test)
        label_pred = clf.predict(X_test)

    elif method=='svm':
        print('Scaling the featuer with StandardScaler')
        scaler = StandardScaler().fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)

        if os.path.exists(modeldir + modelname):
            clf = pickle.load(open(modeldir + modelname, 'rb'))
            print('loaded model')
        else:
            KERNEL=maxDepth
            NU=nEst
            clf = NuSVC(kernel=KERNEL,nu=NU,random_state=0)
            clf.fit(X_train, label_train)
            print(X_train[0])
            print(label_train)
            pickle.dump(clf, open(modeldir + modelname, 'wb'))
        acc = clf.score(X_test, label_test)
        print(label_test)
        label_pred = clf.predict(X_test)
        print(label_pred)

    elif method=='mlp':
        print('Scaling the featuer with StandardScaler')
        scaler = StandardScaler().fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)

        if os.path.exists(modeldir + modelname)==True:
            model = keras.models.load_model(modeldir + modelname)
            print('loaded model')
            predictions = model.predict(X_test)
            label_pred = convertPrediction(set(trainSamplelabelid.values()), predictions)
        else:
            featsize = len(X_train[0])
            lastlayer = len(set(trainSamplelabelid.values()))
            layers = maxDepth
            dropout=nEst
            print('featsize ' + str(featsize))

            if layers == 1:
                layersizes=[]
                first = round((featsize + lastlayer) / 2)
                second = 0
                third = 0
                layersizes.append(first)
            elif layers == 2:
                layersizes=[]
                first = round((featsize + lastlayer) / 2)
                second = round((first + lastlayer) / 2)
                third = 0
                layersizes.append(first)
                layersizes.append(second)
            elif layers == 3:
                layersizes = []
                first = round((featsize + lastlayer) / 2)
                second = round((first + lastlayer) / 2)
                third = round((second + lastlayer) / 2)
                layersizes.append(first)
                layersizes.append(second)
                layersizes.append(third)
            label_pred = mlpmodel(modeldir,modelname,trainSamplelabelid, layersizes,dropout, X_train, X_test, y_train, y_test, label_test)

    acc = accuracy_score(label_test, label_pred)
    f1 = f1_score(label_test, label_pred, average='weighted')
    recall = recall_score(label_test, label_pred, average='weighted')
    precision = precision_score(label_test, label_pred, average='weighted')

    outSingleresult(output, label2id, trainSamplelabelid, testSamplelabelid, label_test, label_pred)

    print(acc)
    print(f1)
    print(recall)
    print(precision)

    id2label = {}
    labelname = []
    for label, id in label2id.items():
        output.write(label + ' : ' + str(id) + '\n')
        id2label[id] = label
    for id in sorted(id2label.keys()):
        labelname.append(id2label[id])

    report = classification_report(label_test, label_pred, target_names=labelname, digits=3)
    print(report)
    output.write(report)
    reportdict = parse(report)
    return acc, f1,precision,recall,reportdict, labelname, str(len(X_train[0]))

def proNAN(X_train,X_test):
    X_test = pd.DataFrame(X_test)
    X_test.fillna(X_train.mean(), inplace=True)
    X_test = np.array(X_test)
    return X_test


def main(zeroper,topper,LOGFILE, settype, task, datatype,random,fold,nbit,feattype,method,maxDepth, maxFeat, nEst, num,protype):
    suffix = '_' + str(zeroper) + '_' + str(topper) + '_' + str(nbit) + '/'
    suffixtest = '_' + str(zeroper) + '_' + str(topper) + '_' + str(nbit) + 'nofilter/'
    if task in taskList:
        dirinfo = 'AAAA'
    else:
        dirinfo = 'BBBB'

    # set up output information
    logfile = settype + method
    if task=='BREAST16':
        modeldir = dirinfo + 'BREAST' + '/RESULT/models/' + settype + '/' + method + '/'
        outdir, logfile = openlogfile(dirinfo, protype, logfile, settype, 'BREAST')
    elif task == 'PBMC16':
        modeldir=dirinfo+'PBMC'+'/RESULT/models/'+settype+'/'+method+'/'
        outdir, logfile = openlogfile(dirinfo, protype, logfile, settype,'PBMC')
    else:
        modeldir = dirinfo + task + '/RESULT/models/' + settype + '/' + method + '/'
        outdir, logfile = openlogfile(dirinfo, protype, logfile, settype, task)
    print('outdir: '+str(outdir))
    print ('logfile: '+str(logfile))

    modelname=outmodel(settype, task, datatype,random,fold,nbit,feattype,method,maxDepth, maxFeat, nEst, num)
    outfilename=outrecord(settype, task, datatype,random,fold,nbit,feattype,method,maxDepth, maxFeat, nEst, num)
    trainSamplelabelid, testSamplelabelid, label2id = getIntraLabelRandom(task, random, fold + 1)

    if settype == 'intra':
        if task in taskList:
            datadir = dirinfo + task + '/data/'
            barcodefile = barcodefiles[task]
            if task == 'PBMC' or task == 'PBMC16':
                datadir = dirinfo+'/PBMC_old' + '/data/'
            filelist = readFilelist(datadir + barcodefile)
        else:
            filelist = getFilelist_BREAST()

        if feattype == 'kmer':
            if datatype=='all':
                datatype=''
            dir = dirinfo + task + '/'  + datatype + 'sampleMatrix_10/'
            sampledir = dir + 'random' + str(random) + '/fold' + str(fold) + '/sampleMatrix' + suffix + 'tmp/'
            print('sample dir: '+sampledir)
            X_train, y_train, labels_train, sampleids_train = getKmerdata(sampledir, filelist, trainSamplelabelid)
            X_test, y_test, labels_test, sampleids_test = getKmerdata(sampledir, filelist, testSamplelabelid)
            X_train= convXreads(dirinfo,task,X_train, sampleids_train, datatype)
            X_test = convXreads(dirinfo,task,X_test, sampleids_test, datatype)
            X_train = np.log2(X_train + 1)
            X_test = np.log2(X_test + 1)

            print(len(X_train[0]))
            print(len(X_test[0]))

        elif feattype == 'gene':
            datatype='none'
            if task=='PBMC':
                filterfeatlog = dirinfo + 'PBMC/data/filtered_gene_bc_matrices/hg19/random' + str(random) + '/fold' + str(
                    fold) + '/matrix_genefiltered.csv'
            else:
                filterfeatlog = dirinfo + 'BREASTtest/geneExpre/random'+str(random)+'/fold'+str(fold)+'/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix_filterlog.txt'
            X_train, y_train, labels_train, sampleids_train = readData(filterfeatlog, trainSamplelabelid)
            X_test, y_test, labels_test, sampleids_test = readData(filterfeatlog, testSamplelabelid)


    if settype=='interdev':
        if task=='BREAST':
            filelist_train = getFilelist_BREAST()
            filelist_test = getFilelist_BREASTtest()
            trainSamplelabelid, testSamplelabelid, label2id = getInterDevLabel(task, random, fold + 1)  # 对于取label用
        elif task=='PBMC':
            filelist_train =readFilelist(dirinfo+'/PBMCtest/data/Cells_labelName_Seurat.csv.inter.barcode')
            filelist_test = readFilelist(dirinfo+'/PBMCtest/data/count_w_metadata.'+patient+'.barcode.label.csv.inter.barcode')
            trainSamplelabelid, testSamplelabelid, label2id = getInterDevLabelPBMC(task, patient,random, fold + 1)

        if feattype == 'kmer':
            if datatype == 'all':
                datatype = ''
            traindir = dirinfo + task + '/' + datatype + 'sampleMatrix_10/'
            trainsampledir = traindir + 'random1000/fold0/sampleMatrix' + suffix + 'tmp/'
            testdir = dirinfo + task + 'test/' + patient + datatype + 'sampleMatrix_10/'

            testsampledir = testdir + 'random1000/fold0/sampleMatrix' + suffixtest + 'tmp/'
            X_train, y_train, labels_train, sampleids_train = getKmerdata(trainsampledir, filelist_train,trainSamplelabelid)
            X_test, y_test, labels_test, sampleids_test = getKmerdata(testsampledir, filelist_test, testSamplelabelid)
            X_train = convXreads(dirinfo,task,X_train, sampleids_train, datatype)
            X_test = convXtestreads(dirinfo,task,patient,X_test, sampleids_test, datatype)
            X_train = np.log2(X_train + 1)
            X_test = np.log2(X_test + 1)
            X_test=proNAN(X_train, X_test)#deal with missing data (C2patient)

        elif feattype == 'gene':
            readtype = 'none'
            if task=='BREAST':
                trainfile = dirinfo+'BREASTtest/geneExpre/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix_filterlog.txt'
                testfile = dirinfo+'BREASTtest/geneExpre/GSE118389_tpm_rsem_added_filterlog.txt'
            elif task=='PBMC':
                trainfile = dirinfo + 'PBMCtest/geneExpre/matrix_genesamplefiltered_inter.csv.genesorted'
                testfile = dirinfo + 'PBMCtest/geneExpre/count_w_metadata.'+patient+'.ge.filteredinter.txt.genesorted'

            X_train, y_train, labels_train, sampleids_train = readData(trainfile, trainSamplelabelid)
            X_test, y_test, labels_test, sampleids_test = readData(testfile, testSamplelabelid)

    print('Train shape')
    print(X_train.shape)
    print('Test shape')
    print(X_test.shape)
    acc, f1,precision,recall,reportdict, labelname, featnum = \
        classify(modeldir,outdir, modelname,outfilename,maxDepth, maxFeat, nEst, X_train, X_test, y_train, y_test, labels_train, labels_test, method,
             trainSamplelabelid, testSamplelabelid, label2id)
    infostr = genInfostr(settype, task, random, fold, nbit,method, feattype, datatype,maxDepth, maxFeat, nEst, num, featnum)
    outlog = open(outdir + logfile, 'a')
    outlogfile(outlog, infostr, acc, f1,precision,recall,reportdict, labelname)

if __name__ == "__main__":
    barcodefiles = {'PBMC': 'barcodes.tsv', 'PBMC16': 'barcodes.tsv'}
    taskList=['PBMC','PBMC16']

    # execute only if run as a script
    logfile = sys.argv[1]
    settype = sys.argv[2]  # 'intra'
    task=sys.argv[3]#BREAST, BREAST16,PBMC，
    datatype = sys.argv[4]  # 'all' 'unmapped','mapped'
    random=int(sys.argv[5])
    fold = int(sys.argv[6])  # 0,1,2,3
    nbit=int(sys.argv[7])
    feattype=sys.argv[8]#kmer ,gene
    method = sys.argv[9]  # gbm,rf
    if method=='svm':
        maxDepth=sys.argv[10] #kernel
        maxFeat=sys.argv[11] #no use
        nEst=float(sys.argv[12]) #C
    elif method=='rf' or method=='gbm':
        if sys.argv[10] != 'None':
            maxDepth = int(sys.argv[10])
        else:
            maxDepth = sys.argv[10]
        maxFeat=sys.argv[11]
        nEst=int(sys.argv[12])
    elif method=='mlp':
        maxDepth=int(sys.argv[10])
        maxFeat=sys.argv[11]
        nEst=float(sys.argv[12])

    num = sys.argv[13]
    protype=sys.argv[14]

    if task in taskList:
        zeroper=0.005
        patient = 'N2' #if PBMC , patient should have value
    else:
        zeroper=0.1
        patient =''
    topper=0.05
    main(zeroper,topper,logfile, settype, task, datatype, random, fold, nbit,feattype, method, maxDepth, maxFeat, nEst, num,protype)

