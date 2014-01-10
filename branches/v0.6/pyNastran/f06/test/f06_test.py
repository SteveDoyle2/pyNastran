import os
import sys

from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.op2.test.op2_test import get_all_files
from pyNastran.f06.test.test_f06 import run_lots_of_files
from pyNastran.utils import get_files_of_type


def parse_skipped_cards(fname):
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()

    results = {}
    for line in lines:
        if 'OES' in line[0:3]:
            (fore, aft) = line.strip().split('->')
            (oes, form, elementTypeNum) = fore.strip().split(' ')
            (element_type, eType) = elementTypeNum.strip().split('=')
            (msg, fpath) = aft.strip().split('-')
            #print "fpath=|%s|" %(fpath)
            fpath = fpath.lstrip()[6:]
            eName = msg.split(' ')[0]
            #print "eName=%s eType=%s form=%s fpath=|%s|" %(eName,eType,form,fpath)
            key = (eName, eType, form)
            if key not in results:
                results[key] = fpath

    filesToAnalyze = []
    for (key, value) in sorted(results.iteritems()):
        filesToAnalyze.append(value)

    f = open('newElements.in', 'wb')
    for fname in filesToAnalyze:
        f.write(fname + '\n')
    f.close()
    return filesToAnalyze


def main():
    # works
    files = get_files_of_type('tests', '.f06')

    foldersFile = 'tests/foldersRead.txt'
    #files2 = ['ann6611.f06']

    iSubcases = []
    debug = False
    saveCases = True
    regenerate = True
    stopOnFailure = False
    getSkipCards = False

    if getSkipCards:
        files2 = parse_skipped_cards('skippedCards.out')
    elif regenerate:
        files2 = get_all_files(foldersFile, '.f06')
        for fname in files2:
            if 'test_f06' in fname:
                os.remove(fname)
        
        files3 = []
        for fname in files2:
            if 'test_f06' not in fname:
                files3.append(fname)
        #files2 = [fname if 'test_f06' not in fname for fname in files2]
        files2 = files3
        #print files2
        #files2 = []
        files2 += files
        files2.sort()
    else:
        files2 = get_failed_files('failedCases.in')

    #files2 = [r'D:\work\move\move_tpl\ar29sadl.f06']
    #files = files+files2
    files = files2
    #files = [r'D:\work\move\move_tpl\see101hs.f06']
    #print len(files)
    #files = []

    #            HIS, R1B        EQEXIN
    #skipFiles = ['accopt3.f06','acms111m.f06','adjoint.f06','aerobeam.f06',] # tpl
    skipFiles = ['nltrot99.f06', 'rot12901.f06']  # giant
    #print files

    nStart = 0
    nStop = 10000
    try:
        os.remove('skippedCards.out')
    except:
        pass

    print("nFiles = %s" % len(files))
    #print files
    import time
    t0 = time.time()
    run_lots_of_files(files, debug, saveCases, skipFiles,
                      stopOnFailure, nStart, nStop)
    print("dt = %f" %(time.time() - t0))
    sys.exit('final stop...')

if __name__ == '__main__':
    main()
