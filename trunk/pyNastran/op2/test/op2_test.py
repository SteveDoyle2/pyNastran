import os
import sys
from pyNastran.op2.test.test_op2 import get_failed_files,runLotsOfFiles
from pyNastran.general.general import get_files_of_type

def parseSkippedCards(fname):
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    
    results = {}
    for line in lines:
        if 'OES' in line[0:3]:
            (fore,aft) = line.strip().split('->')
            (oes,form,elementTypeNum) = fore.strip().split(' ')
            (elementType,eType) = elementTypeNum.strip().split('=')
            (msg,fpath) = aft.strip().split('-')
            #print "fpath=|%s|" %(fpath)
            fpath = fpath.lstrip()[6:]
            eName = msg.split(' ')[0]
            #print "eName=%s eType=%s form=%s fpath=|%s|" %(eName,eType,form,fpath)
            key = (eName,eType,form)
            if key not in results:
                results[key] = fpath
            ###
        ###
    
    filesToAnalyze = []
    for (key,value) in sorted(results.iteritems()):
        filesToAnalyze.append(value)
    
    f = open('newElements.in','wb')
    for fname in filesToAnalyze:
        f.write(fname+'\n')
    f.close()
    return filesToAnalyze

def get_all_files(foldersFile,fileType):
    f = open(foldersFile,'r')
    lines = f.readlines()
    files2 = []
    for line in lines:
        moveDir = os.path.join('r"'+line.strip()+'"')
        moveDir = line.strip()
        if moveDir and moveDir[0] != '#':
            print("moveDir = %s" %(moveDir))
            assert os.path.exists(moveDir),'%s doesnt exist' %(moveDir)
            files2 += get_files_of_type(moveDir,fileType,maxSize=100.)
        ###
    ###
    return files2

if __name__=='__main__':
    # works
    files = get_files_of_type('tests','.op2')
    
    foldersFile = 'tests/foldersRead.txt'

    iSubcases = []
    debug     = False
    makeGeom  = False
    writeBDF  = False
    writeF06  = True
    writeMatlab = False
    printResults = False

    deleteF06 = True
    saveCases = True
    regenerate = False
    stopOnFailure = False
    getSkipCards = False

    if getSkipCards:
        files2 = parseSkippedCards('skippedCards.out')
    elif regenerate:
        files2 = get_all_files(foldersFile,'.op2')
        files2 += files
    else:
        files2 = get_failed_files('failedCases.in')
    files = files2
    
    skipFiles = ['nltrot99.op2','rot12901.op2','plan20s.op2'] # giant

    nStart = 0
    nStop = 10000
    try:
        os.remove('skippedCards.out')
    except:
        pass

    print("nFiles = %s" %(len(files)))
    runLotsOfFiles(files,makeGeom,writeBDF,writeF06,writeMatlab,deleteF06,
                   printResults,debug,saveCases,skipFiles,stopOnFailure,
                   nStart,nStop)
    sys.exit('final stop...')
    

