import os
import sys
from pathlib import Path

import argparse
from cpylog import SimpleLogger
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.utils.nastran_utils import run_nastran


def cmd_line_run_jobs(argv=None, quiet: bool=False):
    """
    run_nastran_job dirname
    run_nastran_job filename -x C:\bin\nastran.exe
    """
    FILE = os.path.abspath(__file__)
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [FILE] + sys.argv[2:] # ['run_jobs'] + sys.argv[2:]
    print(f'argv = {argv}')
    parser = argparse.ArgumentParser(prog='run_jobs')

    parser.add_argument("bdf_dirname_filename", help='path to Nastran filename')
    #parser.add_argument('-o', '--overwrite', default=False, help='overwrite files')
    parser.add_argument('-x', '--exe', default='nastran', help='path to Nastran execuable')
    parser.add_argument('-c', '--cleanup', default=True, help='cleanup the junk output files (log, f04, plt)')
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

    cleanup = args.cleanup
    extensions = ['.dat', '.bdf']
    run_jobs(bdf_filename_dirname, nastran_exe,
             extensions=extensions, cleanup=cleanup)


def run_jobs(bdf_filename_dirname: Path, nastran_exe: str | Path,
             extensions: list[str], cleanup: bool=True,
             run: bool=True) -> None:
    """runs a series of jobs in a specific folder"""
    assert bdf_filename_dirname.exists(), bdf_filename_dirname
    if isinstance(extensions, str):
        extensions = [extensions]
    assert isinstance(extensions, list), extensions
    for ext in extensions:
        assert ext.startswith('.'), extensions

    log = SimpleLogger()
    if bdf_filename_dirname.is_dir():
        dirname = bdf_filename_dirname
        #suffixs = [fname.suffix for fname in dirname.iterdir() if '.test_bdf.' not in fname.name]
        bdf_filenames = [dirname / fname.name for fname in dirname.iterdir()
                         if fname.suffix in extensions and '.test_bdf.' not in fname.name]
        assert len(bdf_filenames) > 0, dirname
    elif bdf_filename_dirname.is_file():
        bdf_filenames = [bdf_filename_dirname]
    else:  # pragma: no cover
        raise NotImplementedError(bdf_filename_dirname)

    bdf_filenames_str = [str(bdf_filename) for bdf_filename in bdf_filenames]
    log.info(f'bdf_filenames = {bdf_filenames_str}')
    for bdf_filename in bdf_filenames:
        assert bdf_filename.exists(), print_bad_path(bdf_filename)

    nfiles = len(bdf_filenames)
    for ifile, bdf_filename in enumerate(bdf_filenames):
        base = os.path.splitext(str(bdf_filename))[0]
        op2_filename = base + '.op2'
        if os.path.exists(op2_filename):
            log.warning(f'skipping {str(bdf_filename)} because {op2_filename} already exists')
            continue

        percent = (ifile + 1) / nfiles * 100
        log.info(f'running {ifile+1}/{nfiles}={percent:.0f}%: {str(bdf_filename)}')
        run_nastran(bdf_filename, nastran_cmd=nastran_exe,
                    cleanup=cleanup, run=run)
        log.debug(f'finished {ifile+1}/{nfiles}={percent:.0f}%: {str(bdf_filename)}')
    log.info('done')


if __name__ == '__main__':
    cmd_line_run_jobs()
