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
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
#import resource
#import time

from pyNastran.bdf.test.test_bdf import run_lots_of_files
from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.op2.test.op2_test import get_all_files
from pyNastran.utils import get_files_of_type

def remove_marc_files(files):
    files2 = []
    for f in files:
        if 'marc' not in f:
            files2.append(f)
    return files2

def get_open_fds():
    fds = []
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    for fd in range(0, soft):
        try:
            flags = fcntl.fcntl(fd, fcntl.F_GETFD)
        except IOError:
            continue
        fds.append(fd)
    return fds

def get_file_names_from_file_number(fds):
    names = []
    for fd in fds:
        names.append(os.readlink('/proc/self/fd/%d' % fd))
    return names

if __name__ == '__main__':
    files = get_files_of_type('tests', '.bdf')
    files += get_files_of_type('tests', '.dat')
    foldersFile = 'tests/foldersRead.txt'

    iSubcases = []
    debug = False

    saveCases = True
    regenerate = False
    stopOnFailure = False

    if regenerate:
        files2 = get_all_files(foldersFile, '.bdf')
        files2 += get_all_files(foldersFile, '.nas')
        files2 += get_all_files(foldersFile, '.dat')
        files2 += files
    else:
        files2 = get_failed_files('failedCases.in')

    files = remove_marc_files(files2)

    skipFiles = []  # giant

    nStart = 0
    nStop = 10000
    try:
        os.remove('skippedCards.out')
    except:
        pass

    print("nFiles = %s" % len(files))
    cid = None
    check = True
    xref = True
    debug = False
    failed_files = run_lots_of_files(files, debug=debug, xref=xref,
                                     check=check, cid=cid)
    try:
        f = open('failedCases.in', 'wb')
    except IOError:
        #fds = get_open_fds()
        #print(get_file_names_from_file_number(fds))
        raise
    for fname in failed_files:
        f.write('%s\n' % fname)
    f.close()
    sys.exit('finished...')
