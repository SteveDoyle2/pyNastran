#pylint: disable=C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
#import time

from pyNastran.bdf.test.test_bdf import run_lots_of_files
from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.op2.test.op2_test import getAllFiles
from pyNastran.general.general import getFilesOfType

if __name__ == '__main__':
    # works
    files  = getFilesOfType('tests', '.bdf')
    files += getFilesOfType('tests', '.dat')
    
    foldersFile = 'tests/foldersRead.txt'

    iSubcases = []
    debug     = False

    saveCases = True
    regenerate = True
    stopOnFailure = False

    if regenerate:
        files2  = getAllFiles(foldersFile, '.bdf')
        files2 += getAllFiles(foldersFile, '.nas')
        files2 += getAllFiles(foldersFile, '.dat')
        files2 += files
    else:
        files2 = get_failed_files('failedCases.in')
    files = files2
    
    skipFiles = [] # giant

    nStart = 0
    nStop = 10000
    try:
        os.remove('skippedCards.out')
    except:
        pass

    print("nFiles = %s" %(len(files)))
    cid = None
    check = True
    xref = False
    debug = False
    failed_files = run_lots_of_files(files, debug=debug, xref=xref,
                                     check=check, cid=cid)
    f = open('failedCases.in','wb')
    for fname in failed_files:
        f.write('%s\n' %(fname))
    f.close()
    sys.exit('final stop...')
    
