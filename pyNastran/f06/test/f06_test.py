from __future__ import print_function
from six import iteritems
import os
import sys

import pyNastran
from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.op2.test.op2_test import get_all_files
from pyNastran.f06.test.test_f06 import run_lots_of_files
from pyNastran.utils.dev import get_files_of_type
pkg_path = pyNastran.__path__[0]


def parse_skipped_cards(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()

    results = {}
    for line in lines:
        if 'OES' in line[0:3]:
            (fore, aft) = line.strip().split('->')
            (oes, form, element_type_num) = fore.strip().split(' ')
            (element_type, etype) = element_type_num.strip().split('=')
            (msg, fpath) = aft.strip().split('-')
            #print("fpath=%r" % fpath)
            fpath = fpath.lstrip()[6:]
            ename = msg.split(' ')[0]
            #print("eName=%s etype=%s form=%s fpath=%r" % (ename, etype, form, fpath))
            key = (ename, etype, form)
            if key not in results:
                results[key] = fpath

    files_to_analyze = []
    for (key, value) in sorted(iteritems(results)):
        files_to_analyze.append(value)

    f = open('new_elements.in', 'wb')
    for fname in files_to_analyze:
        f.write(fname + '\n')
    f.close()
    return files_to_analyze


def main():
    # works
    files = get_files_of_type('tests', '.f06')

    folders_file = os.path.join(pkg_path, 'bdf', 'test', 'tests', 'foldersRead.txt')
    #files2 = ['ann6611.f06']

    isubcases = []
    debug = False
    save_cases = True
    regenerate = True
    stop_on_failure = False
    get_skip_cards = False

    if get_skip_cards:
        files2 = parse_skipped_cards('skipped_cards.out')
    elif regenerate:
        files2 = get_all_files(folders_file, '.f06')
        for fname in files2:
            if 'test_f06' in fname:
                os.remove(fname)

        files3 = []
        for fname in files2:
            if 'test_f06' not in fname:
                files3.append(fname)
        #files2 = [fname if 'test_f06' not in fname for fname in files2]
        files2 = files3
        #print(files2)
        #files2 = []
        files2 += files
        files2.sort()
    else:
        files2 = get_failed_files('failed_cases.in')

    files2 = [fname for fname in files2
              if '.test_op2.f06' not in fname
              and '.test_f06.f06' not in fname]
    #files2 = [r'D:\work\move\move_tpl\ar29sadl.f06']
    #files = files+files2
    files = files2
    #files = [r'D:\work\move\move_tpl\see101hs.f06']
    #print(len(files))
    #files = []

    #            HIS, R1B        EQEXIN
    #skipFiles = ['accopt3.f06','acms111m.f06','adjoint.f06','aerobeam.f06',] # tpl
    skip_files = ['nltrot99.f06', 'rot12901.f06']  # giant
    #print(files)

    nstart = 0
    nstop = 10000
    try:
        os.remove('skipped_cards.out')
    except:
        pass

    print("nfiles = %s" % len(files))
    #print(files)
    import time
    t0 = time.time()
    run_lots_of_files(files, debug, save_cases, skip_files,
                      stop_on_failure, nstart, nstop)
    print("dt = %f" % (time.time() - t0))
    sys.exit('final stop...')

if __name__ == '__main__':  # pragma: no cover
    main()
