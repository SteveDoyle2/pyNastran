# coding: utf-8
"""the interface for bdf_test"""
import os
import sys

import pyNastran
from pyNastran.gui.test.test_gui import run_lots_of_files
from pyNastran.op2.test.test_op2 import get_failed_files
from pyNastran.op2.test.op2_test import get_all_files
from pyNastran.utils.dev import get_files_of_type


def remove_marc_files(filenames):
    """Marc files are not supported"""
    filenames2 = []
    for filename in filenames:
        if 'marc' not in filename:
            filenames2.append(filename)
    return filenames2

# def get_open_fds():
    # import resource
    # fds = []
    # soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    # for fd in range(0, soft):
        # try:
            # flags = fcntl.fcntl(fd, fcntl.F_GETFD)
        # except IOError:
            # continue
        # fds.append(fd)
    # return fds

# def get_file_names_from_file_number(fds):
    # names = []
    # for fd in fds:
        # names.append(os.readlink('/proc/self/fd/%d' % fd))
    # return names


def run(regenerate=True):
    """Runs the full BDF test suite"""
    if crash_cards is None:
        crash_cards = []
    # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test
    files = get_files_of_type('tests', '.bdf')
    files += get_files_of_type('tests', '.dat')
    folders_file = 'tests/foldersRead.txt'

    # isubcases = []
    # save_cases = True
    # stop_on_failure = False

    failed_cases_filename = 'failed_cases%s%s.in' % (sys.version_info[:2])
    if regenerate:
        files2 = get_all_files(folders_file, '.bdf')
        files2 += get_all_files(folders_file, '.nas')
        files2 += get_all_files(folders_file, '.dat')
        files2 += files
        files2.sort()
    else:
        print('failed_cases_filename = %r' % failed_cases_filename)
        files2 = get_failed_files(failed_cases_filename)

    skip_files = [
        'mp10a.dat',
        'mp20e.dat',
        'mp30.dat',
        'mp30b.dat',
        'mp60bd.dat',
        'mp60br.dat',
        'mp60cd.dat',
        'mp60cr.dat',
        'mp70a.dat',
        #'heli112em8.dat',  # horrible CORD1x model
    ]

    files = remove_marc_files(files2)
    files = [fname for fname in files
             if not os.path.basename(fname).startswith('out_')
             and '.test_op2.' not in fname  # removing test output files
             and '.test_bdf.' not in fname
             and 'tecplot' not in fname
             and os.path.basename(fname) not in skip_files]

    print("nfiles = %s" % len(files))
    check = True
    debug = False
    failed_files = run_lots_of_files(files, debug=debug,
                                     encoding='latin1',
                                     dev=True)
    ntotal = len(files)
    nfailed = len(failed_files)
    npassed = ntotal - nfailed
    sys.stderr.write('%i/%i passed\n' % (npassed, ntotal))

    write = 'w'
    with open(failed_cases_filename, write) as failed_cases_file:
        for fname in failed_files:
            failed_cases_file.write('%s\n' % fname)
    sys.exit('finished...')


def run_lots_of_files(files, debug=False, encoding='latin1', dev=True):
    """used by gui_test.py to run thousands of files"""
    nfailed = 0
    ntotal = 0
    npassed = 0
    #time0 = time.time()
    failed_files = []
    for filename in files:
        try:
            is_passed = run_()
        except:
            failed_files.append(filename)
        if not is_passed:
            sys.stderr.write('**file=%s vector_failed=%s\n' % (
                op2file, is_vector_failed))
            #failed_cases.append(op2file)
            nfailed += 1
        else:
            npassed += 1
    return failed_files


def main():
    """the interface for gui_test"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg = "Usage:\n"
    is_release = False
    msg += "gui_test [-r]\n"
    msg += "  gui_test -h | --help\n"
    msg += "  gui_test -v | --version\n"
    msg += "\n"
    msg += "Tests to see if many BDFs will work with pyNastran %s.\n" % ver
    msg += "\n"
    #msg += "Positional Arguments:\n"
    #msg += "  OP2_FILENAME         Path to OP2 file\n"
    #msg += "\n"
    msg += "Options:\n"
    msg += "  -r, --regenerate     Resets the tests\n"
    if len(sys.argv) == 0:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    #print(data)
    regenerate = data['--regenerate']

    run(regenerate=regenerate)


if __name__ == '__main__':  # pragma: no cover
    main()
