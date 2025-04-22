# coding: utf-8
"""the interface for bdf_test"""
import os
import sys

import pyNastran
from pyNastran.bdf.test.test_bdf import run_lots_of_files
from pyNastran.op2.test.op2_test import get_failed_files, get_all_files
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

def run(regenerate: bool=True, run_nastran: bool=False, debug: bool=False,
        sum_load: bool=True, sum_mass: bool=True, run_mcid: bool=True,
        run_skin_solids: bool=True,
        run_export_caero: bool=True,
        xref: bool=True, is_lax_parser: bool=False,
        crash_cards=None):
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
        print(f'failed_cases_filename = {failed_cases_filename!r}')
        files2 = get_failed_files(failed_cases_filename)

    #for filename in files2:
        #print(filename)
    skip_files = [
        'pyNastran_dump.bdf',
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

    print(f'nfiles = {len(files):d}')
    check = True
    debug = False
    size = [8]
    is_double = [False]
    is_lax_parsers = [is_lax_parser]
    post = -1
    failed_files = run_lots_of_files(
        files, debug=debug, xref=xref,
        check=check,
        nastran=nastran,
        size=size, is_double=is_double, post=post,
        is_lax_parser=is_lax_parsers,
        sum_load=sum_load, run_mass=sum_mass,
        run_mcid=run_mcid,
        run_export_caero=run_export_caero,
        run_skin_solids=run_skin_solids,
        encoding='latin1', crash_cards=crash_cards,
        dev=True, run_pickle=True)
    ntotal = len(files)
    nfailed = len(failed_files)
    npassed = ntotal - nfailed
    sys.stderr.write(f'{npassed:d}/{ntotal:d} passed\n')

    with open(failed_cases_filename, 'w') as failed_cases_file:
        for fname in failed_files:
            failed_cases_file.write('%s\n' % fname)
    sys.exit('finished...')


def main():
    """the interface for bdf_test"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    #is_release = False
    skips = '[--skip_loads] [--skip_mass] [--skip_mcid] [--skip_aero] [--skip_skin]'
    msg = (
        f'Usage:  bdf_test [-r] [-n] [-s S...] [-e E] [-x] [-c C] [--safe] [--lax] {skips}\n'
        '        bdf_test -h | --help\n'
        '        bdf_test -v | --version\n'
        '\n'
        f'Tests to see if many BDFs will work with pyNastran {ver}.\n'
        '\n'
        "Options:\n"
        '  -r, --regenerate     Resets the tests\n'
        '  -c C, --crash_cards  Crash on specific cards (e.g. CGEN,EGRID)\n'
        '  -n, --run_nastran    Runs Nastran\n'
        '  -s S, --size S       Sets the field size\n'
        '  -e E, --nerrors E    Allow for cross-reference errors (default=100)\n'
        '  -x, --xref           disables cross-referencing and checks of the BDF.\n'
        '                       (default=False -> on)\n'
        '  --safe               Use safe cross-reference (default=False)\n'
        '  --skip_loads         Disables static/dynamic loads sum\n'
        '  --skip_mass          Disables mass sum\n'
        '  --skip_mcid          Disables MCID checks\n'
        '  --skip_aero          Disables aero exporting\n'
        '  --skip_skin          Disables solid skinning\n'
        '  --lax                Use the lax card parser (default=False)\n'
    )
    if len(sys.argv) == 0:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    #print(data)
    regenerate = data['--regenerate']
    run_nastran = data['--run_nastran']
    sum_load = not data['--skip_loads']
    sum_mass = not data['--skip_mass']
    run_mcid = not data['--skip_mcid']
    run_skin_solids = not data['--skip_skin']
    run_export_caero = not data['--skip_aero']
    xref = not data['--xref']
    is_lax_parser = data['--lax']

    crash_cards = []
    if data['--crash_cards']:
        crash_cards = data['--crash_cards'].split(',')
    run(regenerate=regenerate, run_nastran=run_nastran,
        sum_load=sum_load, sum_mass=sum_mass,
        run_mcid=run_mcid,
        run_export_caero=run_export_caero,
        run_skin_solids=run_skin_solids,
        is_lax_parser=is_lax_parser,
        xref=xref, crash_cards=crash_cards)

if __name__ == '__main__':  # pragma: no cover
    main()
