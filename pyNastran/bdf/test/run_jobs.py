import os
from pathlib import Path
import sys

import argparse
import pyNastran

def run_jobs_cmd_line(argv=None, quiet: bool=False):
    """
    run_nastran_job dirname
    run_nastran_job filename -x C:\bin\nastran.exe
    """
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser()

    parser.add_argument("bdf_dirname_filename", help='path to Nastran filename')
    parser.add_argument('-x', '--exe', help='path to Nastran execuable')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    bdf_filename_dirname = Path(args.bdf_dirname_filename)
    nastran_exe = args.exe
    if nastran_exe is None:
        nastran_exe = 'nastran'
    elif '.bat' in nastran_exe or '.exe' in nastran_exe:
        nastran_exe = Path(nastran_exe)
        assert nastran_exe.exists(), nastran_exe
        assert nastran_exe.is_file(), nastran_exe

    assert bdf_filename_dirname.exists(), bdf_filename_dirname

    from pyNastran.utils.nastran_utils import run_nastran
    if bdf_filename_dirname.is_dir():
        dirname = bdf_filename_dirname
        extensions = ['.dat', '.bdf']
        bdf_filenames = [dirname / fname.name for fname in dirname.iterdir()
                         if fname.suffix in extensions and '.test_bdf.' not in fname.name]
        print('bdf_filenames =', bdf_filenames)
        assert len(bdf_filenames) > 0, dirname
        for bdf_filename in bdf_filenames:
            run_nastran(bdf_filename, nastran_cmd=nastran_exe, cleanup=True)


    elif bdf_filename_dirname.is_file():
        run_nastran(bdf_filename_dirname, nastran_cmd=nastran_exe, cleanup=True)
    else:
        raise NotImplementedError(bdf_filename_dirname)

if __name__ == '__main__':
    run_jobs_cmd_line()
