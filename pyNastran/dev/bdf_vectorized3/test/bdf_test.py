# coding: utf-8
"""the interface for bdf_test"""
import os
import sys

import pyNastran
from pyNastran.dev.bdf_vectorized3.test.test_bdf import run_lots_of_files
from pyNastran.bdf.test.bdf_test import remove_marc_files
from pyNastran.op2.test.op2_test import get_failed_files, get_all_files
from pyNastran.utils.dev import get_files_of_type


def run(regenerate=True, run_nastran=False, debug=False,
        sum_load: bool=True,
        run_nominal: bool=True,
        run_equivalence: bool=True,
        xref: bool=True,
        crash_cards=None):
    """Runs the full BDF test suite"""
    if crash_cards is None:
        crash_cards = []
    # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test
    files = get_files_of_type('tests', '.bdf')
    files += get_files_of_type('tests', '.dat')

    test_dirname = os.path.dirname(pyNastran.bdf.test.__file__)
    folders_file = os.path.join(test_dirname, 'tests', 'foldersRead.txt')

    # isubcases = []
    # save_cases = True
    # stop_on_failure = False
    if run_nastran:
        if os.path.exists(r'C:\MSC.Software\MSC.Nastran\bin\nastran.exe'):
            nastran = r'C:\MSC.Software\MSC.Nastran\bin\nastran.exe scr=yes bat=no old=no '
        else:
            raise RuntimeError('cant find Nastran')
        # elif os.path.exists():
    else:
        nastran = ''

    failed_cases_filename = 'failed_cases%s%s.in' % (sys.version_info[:2])
    if regenerate:
        files2 = get_all_files(folders_file, '.bdf')
        files2 += get_all_files(folders_file, '.nas')
        files2 += get_all_files(folders_file, '.dat')
        files2 = list(set(files2))
        files2.sort()
    else:
        print('failed_cases_filename = %r' % failed_cases_filename)
        files2 = get_failed_files(failed_cases_filename)

    #for filename in files2:
        #print(filename)
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
             and '.test_op2.' not in fname # removing test output files
             and '.test_bdf.' not in fname
             and '.test_bdfv.' not in fname
             and 'tecplot' not in fname
             and os.path.basename(fname) not in skip_files]

    if os.path.exists('skipped_cards.out'):
        os.remove('skipped_cards.out')

    print("nfiles = %s" % len(files))
    check = True
    debug = False
    size = [8]
    is_double = [False]
    post = -1
    failed_files = run_lots_of_files(files, debug=debug, xref=xref,
                                     check=check,
                                     nastran=nastran,
                                     size=size, is_double=is_double, post=post,
                                     encoding='latin1', crash_cards=crash_cards,
                                     run_nominal=run_nominal,
                                     run_equivalence=run_equivalence,
                                     dev=True, pickle_obj=True)
    ntotal = len(files)
    nfailed = len(failed_files)
    npassed = ntotal - nfailed
    sys.stderr.write('%i/%i passed\n' % (npassed, ntotal))

    with open(failed_cases_filename, 'w') as failed_cases_file:
        for fname in failed_files:
            failed_cases_file.write('%s\n' % fname)
    sys.exit('finished...')


def main():
    """the interface for bdf_test"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    #is_release = False
    msg = (
        'Usage:  bdf_test [-r] [-n] [-s S...] [-e E] [-L] [-x] [-c C] [--safe] [--skip_nominal] [--skip_equivalence]\n'
        '        bdf_test -h | --help\n'
        '        bdf_test -v | --version\n'
        '\n'
        "Tests to see if many BDFs will work with pyNastran %s.\n"
        '\n'
        "Options:\n"
        '  -r, --regenerate     Resets the tests\n'
        '  -c C, --crash_cards  Crash on specific cards (e.g. CGEN,EGRID)\n'
        '  -n, --run_nastran    Runs Nastran\n'
        '  -L, --sum_loads      Disables static/dynamic loads sum\n'
        '  -s S, --size S       Sets the field size\n'
        '  -e E, --nerrors E    Allow for cross-reference errors (default=100)\n'
        '  -x, --xref           disables cross-referencing and checks of the BDF.\n'
        '                       (default=False -> on)\n'
        '  --skip_nominal       skips loading & comparison with the unvectorized model\n'
        '  --skip_equivalence   skips nodal equivalencing\n'
        '  --safe               Use safe cross-reference (default=False)\n' % ver
    )
    if len(sys.argv) == 0:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    #print(data)
    regenerate = data['--regenerate']
    run_nastran = data['--run_nastran']
    sum_load = not data['--sum_loads']
    xref = not data['--xref']
    run_nominal = not data['--skip_nominal']
    run_equivalence = not data['--skip_equivalence']

    crash_cards = []
    if data['--crash_cards']:
        crash_cards = data['--crash_cards'].split(',')
    run(regenerate=regenerate, run_nastran=run_nastran,
        sum_load=sum_load,
        run_nominal=run_nominal,
        run_equivalence=run_equivalence,
        xref=xref, crash_cards=crash_cards)

if __name__ == '__main__':  # pragma: no cover
    main()
