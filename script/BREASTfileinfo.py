
codedir='CCCC'
METADATA = codedir+'/BREAST/sharefiles/breast_meta.txt'
METADATAtest = codedir+'/BREASTtest/sharefiles/breasttest_meta.txt'


def getFilelist_BREAST(infile=METADATA):
    f=open (infile)
    filelist=[]
    for line in f:
        runid,spot,sampleid,label,newlabel=line.strip().split(',')
        filelist.append(runid)
    print('Finish getting filelist: '+str (len(filelist)))
    return filelist

def getFilelist_BREASTtest(infile=METADATAtest):
    f=open (infile)
    filelist=[]
    for line in f:
        runid,spot,sampleid,label=line.strip().split(',')
        filelist.append(runid)
    print('Finish getting filelist: '+str (len(filelist)))
    return filelist