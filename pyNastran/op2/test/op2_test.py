from __future__ import print_function
from six import iteritems
import os
import sys

import pyNastran
from pyNastran.op2.test.test_op2 import get_failed_files, run_lots_of_files
from pyNastran.utils.dev import get_files_of_type
pkg_path = pyNastran.__path__[0]

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
            #print("fpath=%r" % fpath)
            fpath = fpath.lstrip()[6:]
            eName = msg.split(' ')[0]
            #print "eName=%s eType=%s form=%s fpath=|%s|" %(eName,eType,form,fpath)
            key = (eName, eType, form)
            if key not in results:
                results[key] = fpath

    filesToAnalyze = []
    for key, value in sorted(iteritems(results)):
        filesToAnalyze.append(value)

    f = open('newElements.in', 'wb')
    for fname in filesToAnalyze:
        f.write(fname + '\n')
    f.close()
    return filesToAnalyze


def get_all_files(folders_file, file_type):
    with open(folders_file, 'r') as f:
        lines = f.readlines()

    files2 = []
    for line in lines:
        move_dir = os.path.join('r"'+line.strip()+'"')
        move_dir = line.strip()
        if move_dir and move_dir[0] != '#':
            if not os.path.exists(move_dir):
                print("***move_dir doesn't exist = %s" % move_dir)
                continue
            print("move_dir = %s" % move_dir)
            #assert os.path.exists(moveDir), '%s doesnt exist' % move_dir
            files2 += get_files_of_type(move_dir, file_type, maxSize=4.2)
    return files2


def main(regenerate=True):
    # works
    files = get_files_of_type('tests', '.op2')

    folders_file = os.path.join(pkg_path, 'bdf', 'test', 'tests', 'foldersRead.txt')

    iSubcases = []
    debug = False
    make_geom = False
    write_bdf = False
    write_f06 = True
    write_op2 = False
    is_vector = [True] # is this vectorized
    vector_stop = [True]  # corresponds to is_vector; stop if case fails=True
    binary_debug = False  # catch any errors

    delete_f06 = True
    saveCases = True
    stop_on_failure = False
    get_skip_cards = False

    if get_skip_cards:
        files2 = parse_skipped_cards('skippedCards.out')
    elif regenerate:
        files2 = get_all_files(folders_file, '.op2')
        files2 += files
        files2.sort()
    else:
        files2 = get_failed_files('failedCases.in')
    files = files2

    skip_files = []
    #skip_files = ['nltrot99.op2', 'rot12901.op2', 'plan20s.op2'] # giant

    nstart = 0
    nstop = 10000
    try:
        os.remove('skippedCards.out')
    except:
        pass

    print("nFiles = %s" % len(files))
    import time
    t0 = time.time()
    run_lots_of_files(files, make_geom=make_geom, write_bdf=write_bdf,
                   write_f06=write_f06, delete_f06=delete_f06,
                   write_op2=write_op2, debug=debug, saveCases=saveCases, skip_files=skip_files,
                   stop_on_failure=stop_on_failure,
                   is_vector=is_vector, vector_stop=vector_stop,
                   nstart=nstart, nstop=nstop, binary_debug=binary_debug,
                   compare=not data['--disablecompare'])
    print("dt = %f" %(time.time() - t0))
    sys.exit('final stop...')

if __name__ == '__main__':
    """the interface for op2_test"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg = "Usage:\n"
    is_release = False
    msg += "op2_test [-r] [-c] [-u]\n"
    msg += "  op2_test -h | --help\n"
    msg += "  op2_test -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP2 will work with pyNastran %s.\n" % ver
    msg += "\n"
    #msg += "Positional Arguments:\n"
    #msg += "  OP2_FILENAME         Path to OP2 file\n"
    #msg += "\n"
    msg += "Options:\n"
    msg += "  -r, --regenerate      Dumps the OP2 as a readable text file\n"
    msg += "  -c, --disablecompare  Doesn't do a validation of the vectorized result\n"
    #msg += "  -z, --is_mag_phase    F06 Writer writes Magnitude/Phase instead of\n"
    #msg += "                        Real/Imaginary (still stores Real/Imag); [default: False]\n"
    #msg += "  -s <sub>, --subcase   Specify one or more subcases to parse; (e.g. 2_5)\n"
    if len(sys.argv) == 0:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    regenerate = data['--regenerate']
    main(regenerate=regenerate)
