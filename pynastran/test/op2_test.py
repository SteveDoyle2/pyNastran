from __future__ import print_function
from six import iteritems
import os
import sys

import pyNastran
from pyNastran.op2.test.test_op2 import get_failed_files, run_lots_of_files
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
            #print "eName=%s eType=%s form=%s fpath=|%s|" %(ename, etype, form, fpath)
            key = (ename, etype, form)
            if key not in results:
                results[key] = fpath

    files_to_analyze = []
    for key, value in sorted(iteritems(results)):
        files_to_analyze.append(value)

    f = open('new_elements.in', 'wb')
    for fname in files_to_analyze:
        f.write(fname + '\n')
    f.close()
    return files_to_analyze


def get_all_files(folders_file, file_type):
    with open(folders_file, 'r') as f:
        lines = f.readlines()

    files2 = []
    for line in lines:
        move_dir = os.path.join('r"' + line.strip() + '"')
        move_dir = line.strip()
        if move_dir and move_dir[0] != '#':
            if not os.path.exists(move_dir):
                print("***move_dir doesn't exist = %s" % move_dir)
                continue
            print("move_dir = %s" % move_dir)
            #assert os.path.exists(move_dir), '%s doesnt exist' % move_dir
            files2 += get_files_of_type(move_dir, file_type, max_size=4.2)
    return files2


def main(regenerate=True, make_geom=False, write_bdf=False, save_cases=True,
         debug=False, write_f06=True):
    # works
    files = get_files_of_type('tests', '.op2')

    folders_file = os.path.join(pkg_path, 'bdf', 'test', 'tests', 'foldersRead.txt')

    isubcases = []
    write_op2 = False
    is_vector = [True] # is this vectorized
    vector_stop = [True]  # corresponds to is_vector; stop if case fails=True
    binary_debug = False  # catch any errors
    quiet = True

    delete_f06 = True
    stop_on_failure = False
    get_skip_cards = False

    if get_skip_cards:
        files2 = parse_skipped_cards('skipped_cards.out')
    elif regenerate:
        files2 = get_all_files(folders_file, '.op2')
        files2 += files
        files2.sort()
    else:
        files2 = get_failed_files('failed_cases.in')
    files = files2

    skip_files = []
    #skip_files = ['nltrot99.op2', 'rot12901.op2', 'plan20s.op2'] # giant

    nstart = 0
    nstop = 10000
    try:
        os.remove('skipped_cards.out')
    except:
        pass

    print("nfiles = %s" % len(files))
    import time
    t0 = time.time()
    short_stats = data['--short_stats']
    run_lots_of_files(files, make_geom=make_geom, write_bdf=write_bdf,
                      write_f06=write_f06, delete_f06=delete_f06,
                      write_op2=write_op2, debug=debug, save_cases=save_cases, skip_files=skip_files,
                      stop_on_failure=stop_on_failure,
                      is_vector=is_vector, vector_stop=vector_stop,
                      nstart=nstart, nstop=nstop, binary_debug=binary_debug,
                      compare=not data['--disablecompare'], short_stats=short_stats,
                      quiet=quiet)
    print("dt = %f" %(time.time() - t0))
    sys.exit('final stop...')

if __name__ == '__main__':
    """the interface for op2_test"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg = "Usage:\n"
    is_release = False
    msg += "op2_test [-r] [-s] [-c] [-u] [-t] [-g] [-n] [-d] [-f]\n"
    msg += "  op2_test -h | --help\n"
    msg += "  op2_test -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP2 will work with pyNastran %s.\n" % ver
    msg += "\n"
    #msg += "Positional Arguments:\n"
    #msg += "  OP2_FILENAME         Path to OP2 file\n"
    #msg += "\n"
    msg += "Options:\n"
    msg += "  -d, --debug           debug logging\n"
    msg += "  -r, --regenerate      Dumps the OP2 as a readable text file\n"
    msg += "  -c, --disablecompare  Doesn't do a validation of the vectorized result\n"
    msg += "  -t, --short_stats     Short get_op2_stats printout\n"
    msg += "  -g, --geometry        Reads the OP2 for geometry, which can be written out\n"
    msg += "  -n, --write_bdf       Writes the bdf to fem.test_op2.bdf (default=False)\n" # n is for NAS
    msg += "  -f, --write_f06       Writes the f06 to fem.test_op2.f06\n"
    msg += "  -s, --save_cases      Disables saving of the cases (default=False)\n"
    #msg += "  -z, --is_mag_phase    F06 Writer writes Magnitude/Phase instead of\n"
    #msg += "                        Real/Imaginary (still stores Real/Imag); [default: False]\n"
    #msg += "  -s <sub>, --subcase   Specify one or more subcases to parse; (e.g. 2_5)\n"
    if len(sys.argv) == 0:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    debug = data['--debug']
    regenerate = data['--regenerate']
    make_geom = data['--geometry']
    write_bdf = data['--write_bdf']
    write_f06 = data['--write_f06']
    save_cases = not data['--save_cases']
    main(regenerate=regenerate, make_geom=make_geom, write_bdf=write_bdf,
         save_cases=save_cases, write_f06=write_f06, debug=debug)
