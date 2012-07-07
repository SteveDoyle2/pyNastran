import os
import sys
import time
from traceback import print_exc
from pyNastran.op2.op2    import OP2,EndOfFileError
from pyNastran.bdf.errors import *
from pyNastran.op2.op2Errors import *
from pyNastran.op2.test.test_op2 import getFailedFiles,runOP2,runLotsOfFiles
from pyNastran.general.general import getFilesOfType

def parseSkippedCards(fname):
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    
    results = {}
    for line in lines:
        if 'OES' in line[0:3]:
            #print line
            (fore,aft) = line.strip().split('->')
            (oes,form,elementTypeNum) = fore.strip().split(' ')
            (elementType,eType) = elementTypeNum.strip().split('=')
            (msg,fpath) = aft.strip().split('-')
            #print "fpath=|%s|" %(fpath)
            fpath = fpath.lstrip()[6:]
            #print fpath
            #sys.exit()
            eName = msg.split(' ')[0]
            #print "eName=%s eType=%s form=%s fpath=|%s|" %(eName,eType,form,fpath)
            key = (eName,eType,form)
            if key not in results:
                results[key] = fpath
            ###
        ###
    
    filesToAnalyze = []
    for (key,value) in sorted(results.iteritems()):
        #print key,value
        filesToAnalyze.append(value)
    
    f = open('newElements.in','wb')
    for fname in filesToAnalyze:
        f.write(fname+'\n')
    f.close()
    return filesToAnalyze

def getAllFiles(foldersFile,fileType):
    f = open(foldersFile,'r')
    lines = f.readlines()
    files2 = []
    for line in lines:
        moveDir = os.path.join('r"'+line.strip()+'"')
        moveDir = line.strip()
        if moveDir and moveDir[0] != '#':
            print("moveDir = %s" %(moveDir))
            assert os.path.exists(moveDir),'%s doesnt exist' %(moveDir)
            files2 += getFilesOfType(moveDir,fileType,maxSize=100.)
        ###
    ###
    return files2

if __name__=='__main__':
    # works
    files = getFilesOfType('tests','.op2')
    
    #moveDir = r'D:\work\move\hard_demo'
    #moveDir = r'D:\work\move\move_tpl'
    foldersFile = 'tests/foldersRead.txt'
    #moveDir = r'D:\work\move\solid_shell_bar'
    #files2 = ['ann6611.op2']

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
        files2 = getAllFiles(foldersFile,'.op2')
        files2 += files
    else:
        files2 = getFailedFiles('failedCases.in')
    
    #files2 = [r'D:\work\move\move_tpl\ar29sadl.op2']
    #files = files+files2
    files = files2
    #files = [r'D:\work\move\move_tpl\see101hs.op2']
    #print len(files)
    #files = []
    
    #            HIS, R1B        EQEXIN
    #skipFiles = ['accopt3.op2','acms111m.op2','adjoint.op2','aerobeam.op2',] # tpl
    skipFiles = ['nltrot99.op2','rot12901.op2','plan20s.op2'] # giant
    #skipFiles = []
    #print files

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
    #runLotsOfFiles(files,makeGeom,writeBDF,debug,saveCases,stopOnFailure,nStart,nStop)
    sys.exit('final stop...')
    

