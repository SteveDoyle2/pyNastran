#pylint: disable=C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import PY2
import os
import sys
#import time

from pyNastran.bdf.test.test_bdf import run_lots_of_files
from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.op2.test.op2_test import get_all_files
from pyNastran.utils.dev import get_files_of_type

def remove_marc_files(files):
    files2 = []
    for f in files:
        if 'marc' not in f:
            files2.append(f)
    return files2

def get_open_fds():
    import resource
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

def main():
    # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test
    files = get_files_of_type('tests', '.bdf')
    files += get_files_of_type('tests', '.dat')
    foldersFile = 'tests/foldersRead.txt'

    iSubcases = []
    debug = False

    saveCases = True
    regenerate = False
    stopOnFailure = False
    nastran = r'C:\MSC.Software\MSC.Nastran\bin\nastran.exe scr=yes bat=no old=no '
    nastran = ''

    if regenerate:
        files2 = get_all_files(foldersFile, '.bdf')
        files2 += get_all_files(foldersFile, '.nas')
        files2 += get_all_files(foldersFile, '.dat')
        files2 += files
        files2.sort()
    else:
        files2 = get_failed_files('failedCases.in')

    files = remove_marc_files(files2)
    files = [fname for fname in files
             if not os.path.basename(fname).startswith('out_')]  # removing test output files

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
    size = [8]
    is_double = [False]
    post = -1
    failed_files = run_lots_of_files(files, debug=debug, xref=xref,
                                     check=check, cid=cid,
                                     nastran=nastran,
                                     size=size, is_double=is_double, post=post)
    ntotal = len(files)
    nfailed = len(failed_files)
    npassed = ntotal - nfailed
    sys.stderr.write('%i/%i passed\n' % (npassed, ntotal))
    try:
        if PY2:
            f = open('failedCases.in', 'wb')
        else:
            f = open('failedCases.in', 'w')
    except IOError:
        #fds = get_open_fds()
        #print(get_file_names_from_file_number(fds))
        raise
    for fname in failed_files:
        f.write('%s\n' % fname)
    f.close()
    sys.exit('finished...')

if __name__ == '__main__':  # pragma: no cover
    main()
