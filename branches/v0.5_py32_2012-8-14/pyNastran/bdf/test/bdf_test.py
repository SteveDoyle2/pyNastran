## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
#pylint: disable=C0103


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
    f = open('failedCases.in','w',encoding='utf-8')
    for fname in failed_files:
        f.write('%s\n' %(fname))
    f.close()
    sys.exit('final stop...')
    
