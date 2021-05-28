from collections import Counter
dirinfo='CCCC'

labelCon = {'TNBC_SC_Tumor_Tumor_Tumor': 'TNBC',
            'nonTumor_Stromal_Stromal': 'Stromal',
            'nonTumor_Immune_Myeloid': 'Myeloid',
            'HER2_SC_Tumor_Tumor_Tumor': 'Her2',
            'nonTumor_Immune_Tcell': 'Tcell',
            'LA_SC_Tumor_Tumor_Tumor': 'LA',
            'LB_SC_Tumor_Tumor_Tumor': 'LB',
            'nonTumor_Immune_Bcell': 'Bcell'}

labelConInter = {'TNBC_SC_Tumor_Tumor_Tumor': 'epithelial',
                 'nonTumor_Stromal_Stromal': 'stroma',
                 'nonTumor_Immune_Myeloid': 'macrophage',
                 'HER2_SC_Tumor_Tumor_Tumor': 'epithelial',
                 'nonTumor_Immune_Tcell': 'Tcell',
                 'LA_SC_Tumor_Tumor_Tumor': 'epithelial',
                 'LB_SC_Tumor_Tumor_Tumor': 'epithelial',
                 'nonTumor_Immune_Bcell': 'Bcell'
                 }

def labelIndex(trainSamplelabel):
    labelset = set()
    label2id = {}
    trainSamplelabelid = {}
    for sampleid, label in trainSamplelabel.items():
        labelset.add(label)
    sortedlabelset = sorted(labelset)
    cnt = 0
    for label in sortedlabelset:
        label2id[label] = cnt
        cnt = cnt + 1
    for sampleid, label in trainSamplelabel.items():
        trainSamplelabelid[sampleid] = label2id[label]
    return trainSamplelabelid, label2id


def getTrainsamplelabel(task,infile, flag):
    f = open(infile)
    sampleLabel = {}
    for line in f:
        if task=='PBMC' or task=='COVIDN2':
            infile,label=line.strip().split(',')
            label=label.replace(' ', '_')
            sampleLabel[infile] = label
        else:
            runid, spot, sampleid, label, newlabel = line.strip().split(',')
            if flag == 'intra':
                sampleLabel[sampleid] = labelCon[newlabel]
            elif flag == 'inter':
                sampleLabel[sampleid] = labelConInter[newlabel]
    return sampleLabel


def getIntratestsamplelabel(task,label2id, infile):
    f = open(infile)
    sampleLabel = {}
    for line in f:
        if task=='PBMC' or task=='COVIDN2':
            infile,label=line.strip().split(',')
            label = label.replace(' ', '_')
            sampleLabel[infile] =label2id[label]
        else:
            runid, spot, sampleid, label, newlabel = line.strip().split(',')
            if labelCon[newlabel] in label2id.keys():
                sampleLabel[sampleid] = label2id[labelCon[newlabel]]
    return sampleLabel


def getIntrawithIntertestsamplelabel(label2id, infile):
    f = open(infile)
    sampleLabel = {}
    for line in f:
        runid, spot, sampleid, label, newlabel = line.strip().split(',')
        if labelConInter[newlabel] in label2id.keys():
            sampleLabel[sampleid] = label2id[labelConInter[newlabel]]
    return sampleLabel


def getIntertestsamplelabel(task,label2id, infile):
    f = open(infile)
    sampleLabel = {}
    for line in f:
        if task == 'PBMC' or task == 'COVIDN2':
            infile, label = line.strip().split(',')
            label = label.replace(' ', '_')
            sampleLabel[infile] = label2id[label]
        else:
            #runid, spot, sampleid, label, newlabel = line.strip().split(',')
            runid, spot, sampleid, label = line.strip().split(',')
            if label in label2id.keys():
                sampleLabel[sampleid] = label2id[label]
    return sampleLabel


def getLabelinfo(sampleLabel):
    labelCount = {}
    for sampleid, label in sampleLabel.items():
        if label not in labelCount.keys():
            labelCount[label] = 0
        labelCount[label] = labelCount[label] + 1
    return labelCount


def getIntraLabelRandom(task,random, fold):
    if task=='PBMC' or task=='PBMC16':
        task='PBMC'
        METADATA_TRAIN = dirinfo + task + '/sharefiles/Cells_labelName_Seurat_random' + str(random) + '_train' + str(fold) + '.txt'
        METADATA_TEST = dirinfo + task + '/sharefiles/Cells_labelName_Seurat_random' + str(random) + '_test' + str(fold) + '.txt'
    elif task=='COVIDN2':
        METADATA_TRAIN = dirinfo + task + '/sharefiles/count_w_metadata.N2.barcode.label_random' + str(random) + '_train' + str(
            fold) + '.txt'
        METADATA_TEST = dirinfo + task + '/sharefiles/count_w_metadata.N2.barcode.label_random' + str(random) + '_test' + str(
            fold) + '.txt'
    elif task=='BREAST' or task=='BREAST16':
        task='BREAST'
        METADATA_TRAIN = dirinfo+task+'/sharefiles/breast_meta_random' + str(random) + '_train' + str(fold) + '.txt'
        METADATA_TEST = dirinfo+task+'/sharefiles/breast_meta_random' + str(random) + '_test' + str(fold) + '.txt'

    trainSamplelabel = getTrainsamplelabel(task,METADATA_TRAIN, 'intra')
    trainSamplelabelid, label2id = labelIndex(trainSamplelabel)
    testSamplelabelid = getIntratestsamplelabel(task,label2id, METADATA_TEST)

    return trainSamplelabelid, testSamplelabelid, label2id


def getInterDevLabel(task,random, fold):
    METADATA_TRAIN = dirinfo+'BREAST/sharefiles/breast_meta.txt'
    METADATA_TEST = dirinfo+'BREASTtest/sharefiles/breasttest_meta' + '_random' + str(random) + '_test' + str(fold) + '.txt'
    trainSamplelabel = getTrainsamplelabel(task,METADATA_TRAIN, 'inter')
    print(Counter(trainSamplelabel.values()))
    trainSamplelabelid, label2id = labelIndex(trainSamplelabel)
    print(Counter(trainSamplelabelid.values()))
    testSamplelabelid = getIntertestsamplelabel(task,label2id, METADATA_TEST)
    print(Counter(testSamplelabelid.values()))
    return trainSamplelabelid, testSamplelabelid, label2id

def getInterDevLabelPBMC(task,patient,random, fold):
    METADATA_TRAIN = dirinfo+'PBMC/sharefiles/Cells_labelName_Seurat.csv.inter'
    METADATA_TEST = dirinfo+'PBMCtest/sharefiles/count_w_metadata.'+patient+'.barcode.label.inter' + '_random' + str(random) + '_test' + str(fold) + '.txt'
    trainSamplelabel = getTrainsamplelabel(task,METADATA_TRAIN, 'inter')
    print('train set')
    print(Counter(trainSamplelabel.values()))
    trainSamplelabelid, label2id = labelIndex(trainSamplelabel)
    print('train set')
    print(Counter(trainSamplelabelid.values()))
    testSamplelabelid = getIntertestsamplelabel(task,label2id, METADATA_TEST)
    print('train test')
    print(Counter(testSamplelabelid.values()))
    return trainSamplelabelid, testSamplelabelid, label2id

def getIntraLabelPBMC():
    METADATA_TRAIN = dirinfo  + 'PBMC/sharefiles/Cells_labelName_Seurat.csv'
    trainSamplelabel = getTrainsamplelabel('PBMC',METADATA_TRAIN, 'intra')
    trainSamplelabelid, label2id = labelIndex(trainSamplelabel)
    return trainSamplelabelid