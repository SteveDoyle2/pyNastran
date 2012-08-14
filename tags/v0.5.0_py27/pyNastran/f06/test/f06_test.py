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
import os
import sys
#import time
from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.f06.test.test_f06 import runLotsOfFiles
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

def main():
    # works
    files = getFilesOfType('tests','.f06')
    
    #moveDir = r'D:\work\move\hard_demo'
    moveDir = r'D:\work\move\move_tpl'
    #moveDir = r'D:\work\move\solid_shell_bar'
    #files2 = ['ann6611.f06']

    iSubcases = []
    debug     = False
    saveCases = True
    regenerate = True
    stopOnFailure = False
    getSkipCards = False

    if getSkipCards:
        files2 = parseSkippedCards('skippedCards.out')
    elif regenerate:
        files2 = getFilesOfType(moveDir,'.f06')
        files2 = []
        files2 += files
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
    skipFiles = ['nltrot99.f06','rot12901.f06'] # giant
    #print files

    nStart = 0
    nStop = 10000
    try:
        os.remove('skippedCards.out')
    except:
        pass

    print("nFiles = %s" %(len(files)))
    runLotsOfFiles(files,debug,saveCases,skipFiles,stopOnFailure,nStart,nStop)
    #runLotsOfFiles(files,debug,saveCases,stopOnFailure,nStart,nStop)
    sys.exit('final stop...')

if __name__ == '__main__':
    main()
