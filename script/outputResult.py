import os
from readBreastlabel import getLabelinfo
from sklearn.metrics import confusion_matrix

def genlogfile(outdir,logfile,settype):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    logfile = logfile + '.csv'
    if not os.path.exists(outdir + logfile):
        outlog = open(outdir + logfile, 'w')
        if settype=='intra':
            outlog.write(genIntraTitle() + '\n')
        elif settype=='interdev':
            outlog.write(genInterTitle() + '\n')
    return logfile

def openlogfile(dirinfo,protype,logfile,settype,task):
    if protype == 'grid':
        outdir = dirinfo + task+'/RESULT/gridSearch/' + logfile + '/'
    elif protype == 'ave':
        outdir = dirinfo +task+'/RESULT/aveBatch/' + logfile + '/'
    try:
        if not os.path.exists(os.path.dirname(outdir)):
            os.makedirs(os.path.dirname(outdir))
    except OSError as err:
        print(err)

    logfile = logfile + '.csv'
    try:
        if os.path.exists(outdir+logfile)==False:
            outlog = open(outdir + logfile, 'w')
            if settype == 'intra':
                outlog.write(genIntraTitle() + '\n')
            elif settype == 'interdev':
                outlog.write(genInterTitle() + '\n')
    except OSError as err:
        print(err)
    return outdir,logfile


def outrecord(settype, task, datatype,random,fold,nbit,feattype,method,maxDepth, maxFeat, nEst, num):
    outfilename = settype + '_' + task + '_' + datatype + '_' + str(random) + '_' + str(fold) + '_' + \
                str(nbit) + '_' + feattype + '_' + method + '_' + str(maxDepth) + '_' + str(maxFeat) + '_' + str(nEst) + '_' + str(num)
    return outfilename

def outmodel(settype, task, datatype,random,fold,nbit,feattype,method,maxDepth, maxFeat, nEst, num):
    if datatype=='':
        datatype='all'
    modelname='model_'+settype+'_'+task+'_'+datatype+'_'+str(random)+'_'+str(fold)+'_'+\
              str(nbit)+'_'+feattype+'_'+method+'_'+str(maxDepth)+'_'+str(maxFeat)+'_'+str(nEst)+'_'+str(num)
    return modelname


def outSingleresult(output,label2id,trainSamplelabelid,testSamplelabelid,label_test,label_pred):
    trainlabelinfo = getLabelinfo(trainSamplelabelid)
    testlabelinfo = getLabelinfo(testSamplelabelid)

    output.write('Train:\n')
    for label, cnt in trainlabelinfo.items():
        output.write(str(label) + ' : ' + str(cnt) + '\n')
    output.write('Test:\n')
    for label, cnt in testlabelinfo.items():
        output.write(str(label) + ' : ' + str(cnt) + '\n')

    for label, id in label2id.items():
        output.write(label + ' : ' + str(id) + '\n')

    m = confusion_matrix(label_test, label_pred)
    output.write('finish computing confusion matrix\n')
    for i in range(len(m)):
        for j in range(len(m[0])):
            output.write(str(m[i][j]))
            output.write(' ')
        output.write('\n')



def genIntraTitle():
    str = 'settype,task,random,fold,nbit,method,feattype,datatype,maxDepth,maxFeat,sEst,experTime,' \
          'featNum,acc,f1,precision,recall,Bcell_pre,Bcell_recall,Bcell_f1,Her2_pre,Her2_recall,' \
          'Her2_f1,LA_pre,LA_recall,LA_f1,LB_pre,LB_recall,LB_f1,Myeloid_pre,Myeloid_recall,' \
          'Myeloid_f1,Stromal_pre,Stromal_recall,Stromal_f1,TNBC_pre,TNBC_recall,TNBC_f1,Tcell_pre,Tcell_recall,Tcell_f1'
    return str


def genInterTitle():
    str = 'settype,task,random,fold,nbit,method,feattype,datatype,maxDepth,maxFeat,sEst,experTime,' \
          'featNum,acc,f1,precision,recall,Bcell_pre,Bcell_recall,Bcell_f1,Tcell_pre,Tcell_recall,' \
          'Tcell_f1,epithelial_pre,epithelial_recall,epithelial_f1,macrophage_pre,macrophage_recall,' \
          'macrophage_f1,Stromal_pre,Stromal_recall,Stromal_f1'
    return str


def genInfostr(settype, task, random,fold, nbit,method, feattype, datatype,maxDepth, maxFeat, nEst, num, featnum):
    infostr = settype + ',' + task + ','+str(random)+','+str(fold)+ ',' +str(nbit)+','+ method + ',' + feattype + ',' +datatype+','+ str(maxDepth) + ',' + str(
        maxFeat) + ',' + str(nEst) + ',' + num + ',' + str(featnum)
    return infostr

def outlogfile(outlog, infostr, acc, f1,precision,recall,reportdict, labelname):
    outlog.write('\n')
    outlog.write(infostr)
    outlog.write(',')
    outlog.write(str(round(acc, 3)))
    outlog.write(',')
    outlog.write(str(round(f1, 3)))
    outlog.write(',')
    outlog.write(str(round(precision, 3)))
    outlog.write(',')
    outlog.write(str(round(recall, 3)))
    outlog.write(',')


    for label in labelname:
        outlog.write(reportdict[label]['precision'] + ',')
        outlog.write(reportdict[label]['recall'] + ',')
        outlog.write(reportdict[label]['f1'] + ',')


def parse(report):
    reportdict = {}
    lines = report.split('\n')
    for line in lines[2:]:
        if line != '':
            line = line.strip()
            line = line.replace('avg / total', 'avg/total')
            line = '\t'.join(line.split())
            print(line)
            classlabel, precision, recall, f1, support = line.strip().split('\t')
            if classlabel not in reportdict.keys():
                reportdict[classlabel] = {}
            reportdict[classlabel]['precision'] = precision
            reportdict[classlabel]['recall'] = recall
            reportdict[classlabel]['f1'] = f1
        else:
            break
    return reportdict