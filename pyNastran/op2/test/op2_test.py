import os
import sys

import pyNastran
from pyNastran.utils.dev import get_files_of_type
PKG_PATH = pyNastran.__path__[0]

def get_failed_files(filename):
    """Gets the list of failed files"""
    with open(filename, 'r') as infile:
        lines = infile.readlines()

    files = []
    for line in lines:
        line = line.strip()
        if line not in files:
            files.append(line)
    return files


def parse_skipped_cards(fname):
    """
    An outdated method to read the log file created by this script to
    determine what was skipped."""
    with open(fname, 'r') as skip_file:
        lines = skip_file.readlines()

    results = {}
    for line in lines:
        if 'OES' in line[0:3]:
            (fore, aft) = line.strip().split('->')
            (unused_oes, form, element_type_num) = fore.strip().split(' ')
            (unused_element_type, etype) = element_type_num.strip().split('=')
            (msg, fpath) = aft.strip().split('-')
            #print("fpath=%r" % fpath)
            fpath = fpath.lstrip()[6:]
            ename = msg.split(' ')[0]
            #print "eName=%s eType=%s form=%s fpath=%r" % (ename, etype, form, fpath)
            key = (ename, etype, form)
            if key not in results:
                results[key] = fpath

    files_to_analyze = []
    for key, value in sorted(results.items()):
        files_to_analyze.append(value)

    with open('new_elements.in', 'w') as new_elements_file:
        for file_to_analyze in files_to_analyze:
            new_elements_file.write(file_to_analyze + '\n')
    return files_to_analyze


def get_all_files(folders_file, file_type, max_size=4.2):
    """
    Gets all the files in the folder and subfolders.  Ignores missing folders.

    Parameters
    ----------
    folders_file : str
        path to the file with a list of folders
    file_type : str
        a file extension
    max_size : float; default=4.2
        size in MB for max file size

    Returns
    -------
    filenames : List[str]
        a series of filenames that were found
    """
    with open(folders_file, 'r') as file_obj:
        lines = file_obj.readlines()

    files2 = []
    for line in lines:
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        if '"' in line:
            # "C:\Program Files\Siemens\NX 12.0\NXNASTRAN\nxn12\nast"
            line = line.strip('"')
            pthi = line.split('\\')
            unused_pth = os.path.join(*pthi)
            move_dir = os.path.join(line)
        else:
            # C:\MSC.Software\MSC.Nastran\msc20051\nast\doc
            move_dir = os.path.join(line)
        #move_dir = line.strip()
        if move_dir:
            if not os.path.exists(move_dir):
                #print("***move_dir doesn't exist = %r" % move_dir)
                continue
            print("move_dir = %s" % move_dir)
            #assert os.path.exists(move_dir), '%s doesnt exist' % move_dir
            files_in_dir = get_files_of_type(move_dir, file_type, max_size=max_size)
            files2 += files_in_dir
            #print('nfiles = %s/%s' % (len(files_in_dir), len(files2)))
    #print('nfiles = %s' % len(files2))
    return files2

def run(regenerate=True, make_geom=False, combine=True,
        write_bdf=False, build_pandas=True,
        xref_safe=False,
        save_cases=True, debug=False, write_f06=True, write_op2=False,
        compare=True, short_stats=False, write_hdf5=True):
    # works
    files = get_files_of_type('tests', '.op2')

    folders_file1 = os.path.join(PKG_PATH, 'bdf', 'test', 'tests', 'foldersRead.txt')
    folders_file2 = os.path.join(PKG_PATH, 'op2', 'test', 'folders_read.txt')

    unused_isubcases = []
    binary_debug = [True, False]  # catch any errors
    quiet = True

    stop_on_failure = False
    get_skip_cards = False

    max_size = 4000. # MB
    failed_cases_filename = 'failed_cases%s%s.in' % (sys.version_info[:2])
    if get_skip_cards:
        files2 = parse_skipped_cards('skipped_cards.out')
    elif regenerate or not os.path.exists(failed_cases_filename):
        files2 = get_all_files(folders_file1, '.op2', max_size=max_size)
        files2 = get_all_files(folders_file2, '.op2', max_size=max_size)
        files2 += files
        assert len(files2) > 0, files2
    else:
        print('failed_cases_filename = %r' % failed_cases_filename)
        files2 = get_failed_files(failed_cases_filename)
    assert len(files2) > 0, files2
    files = list(set(files2))
    files.sort()
    files = [filename for filename in files if '.test_op2.' not in filename]

    skip_files = []
    #skip_files = ['nltrot99.op2', 'rot12901.op2', 'plan20s.op2'] # giant

    nstart = 0
    nstop = 20000
    if os.path.exists('skipped_cards.out'):
        os.remove('skipped_cards.out')

    #try:
    #os.remove('skipped_cards.out')
    #except F:
        #pass

    print("nfiles = %s" % len(files))
    import time
    time0 = time.time()

    from pyNastran.op2.test.test_op2 import run_lots_of_files
    failed_files = run_lots_of_files(files, make_geom=make_geom, combine=combine,
                                     write_bdf=write_bdf,
                                     xref_safe=xref_safe,
                                     write_f06=write_f06, delete_f06=True,
                                     write_op2=write_op2, delete_op2=True,
                                     write_hdf5=write_hdf5, delete_hdf5=True,
                                     build_pandas=build_pandas,
                                     debug=debug,
                                     skip_files=skip_files, stop_on_failure=stop_on_failure,
                                     nstart=nstart, nstop=nstop, binary_debug=binary_debug,
                                     compare=compare, short_stats=short_stats,
                                     quiet=quiet, dev=True)
    if save_cases:
        with open(failed_cases_filename, 'w') as failed_cases_file:
            for op2file in failed_files:
                failed_cases_file.write('%s\n' % op2file)

    seconds = time.time() - time0
    minutes = seconds / 60.
    print("dt = %.2f seconds = %.2f minutes" % (seconds, minutes))
    ntotal = len(files)
    nfailed = len(failed_files)
    npassed = ntotal - nfailed

    msg = '-----done with all models %s/%s=%.2f%%  nfailed=%s-----' % (
        npassed, ntotal,
        100. * npassed / float(ntotal),
        ntotal - npassed)
    print(msg)
    sys.exit("%s\ndt = %.2f seconds = %.2f minutes" % (msg, seconds, minutes))


def main():
    """the interface for op2_test"""
    from docopt import docopt
    ver = str(pyNastran.__version__)

    msg = "Usage:  "
    #is_release = False
    is_dev = 'dev' in ver
    msg += "op2_test [-r] [-s] [-c] [-u] [-t] [-g] [-n] [-f] [-o] [-h] [-d] [-b] [--safe] [--skip_dataframe] [--nocombine]\n"
    msg += "        op2_test -h | --help\n"
    msg += "        op2_test -v | --version\n"
    msg += "\n"
    msg += "Tests to see if an OP2 will work with pyNastran %s.\n" % ver
    msg += "\n"
    #msg += "Positional Arguments:\n"
    #msg += "  OP2_FILENAME         Path to OP2 file\n"
    #msg += "\n"
    msg += "Options:\n"
    msg += "  -r, --regenerate       Resets the tests\n"
    msg += "  -b, --binary_debug     Dumps the OP2 as a readable text file\n"
    msg += "  -c, --disablecompare   Doesn't do a validation of the vectorized result\n"
    msg += "  -t, --short_stats      Short get_op2_stats printout\n"
    msg += "  -g, --geometry         Reads the OP2 for geometry, which can be written out\n"
    # n is for NAS
    msg += "  -n, --write_bdf        Writes the bdf to fem.test_op2.bdf (default=False)\n"
    msg += "  -f, --write_f06        Writes the f06 to fem.test_op2.f06\n"
    msg += "  -o, --write_op2        Writes the op2 to fem.test_op2.op2\n"
    msg += "  -h, --write_hdf5       Writes the hdf5 to fem.test_op2.h5\n"
    msg += "  --skip_dataframe       Disables pandas dataframe building; [default: False]\n"
    msg += "  --nocombine            Disables case combination\n"
    msg += "  -s, --save_cases       Disables saving of the cases (default=False)\n"
    msg += "  --safe                 Safe cross-references BDF (default=False)\n"
    #msg += "  -z, --is_mag_phase    F06 Writer writes Magnitude/Phase instead of\n"
    #msg += "                        Real/Imaginary (still stores Real/Imag); [default: False]\n"
    #msg += "  -s <sub>, --subcase   Specify one or more subcases to parse; (e.g. 2_5)\n"
    msg += "  -d, --debug            debug logging\n"
    if len(sys.argv) == 0:
        sys.exit(msg)

    data = docopt(msg, version=ver)
    debug = data['--debug']
    #binary_debug = data['--binary_debug']
    regenerate = data['--regenerate']
    make_geom = data['--geometry']
    write_bdf = data['--write_bdf']
    write_f06 = data['--write_f06']
    write_op2 = data['--write_op2']
    write_hdf5 = data['--write_hdf5']
    save_cases = not data['--save_cases']
    short_stats = data['--short_stats']
    compare = not data['--disablecompare']
    build_pandas = not data['--skip_dataframe']
    xref_safe = data['--safe']
    combine = not data['--nocombine']
    run(regenerate=regenerate, make_geom=make_geom,
        combine=combine,
        write_bdf=write_bdf,
        xref_safe=xref_safe,
        save_cases=save_cases, write_f06=write_f06, write_op2=write_op2,
        write_hdf5=write_hdf5, short_stats=short_stats,
        build_pandas=build_pandas, compare=compare, debug=debug)

if __name__ == '__main__':   # pragma: no cover
    main()
