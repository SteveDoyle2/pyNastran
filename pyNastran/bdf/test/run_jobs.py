import os
import sys
from pathlib import Path
import time
import argparse

from cpylog import SimpleLogger
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.utils.nastran_utils import run_nastran
from pyNastran.utils.dev import get_files_of_type


def cmd_line_run_jobs(argv=None, quiet: bool=False):
    """
    run_nastran_job dirname
    run_nastran_job filename.bdf -x C:\bin\nastran.exe
    run_nastran_job .            -x C:\bin\nastran.exe --cleanup -r
    """
    FILE = os.path.abspath(__file__)
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [FILE] + sys.argv[2:]  # ['run_jobs'] + sys.argv[2:]
    print(f'argv = {argv}')
    parser = argparse.ArgumentParser(prog='run_jobs')

    parser.add_argument("bdf_dirname_filename", help='path to Nastran filename')
    #parser.add_argument('-o', '--overwrite', default=False, help='overwrite files')
    parser.add_argument('-x', '--exe', default='nastran', help='path to Nastran execuable')
    parser.add_argument('-c', '--cleanup', action='store_true', help='cleanup the junk output files (log, f04, plt)')
    parser.add_argument('-r', '--recursive', action='store_false', help='recursively search for directories')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                        version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    recursive = args.recursive
    bdf_filename_dirname = Path(args.bdf_dirname_filename)
    nastran_exe = args.exe
    # if nastran_exe is None:
    #     nastran_exe = 'nastran'
    if '.bat' in nastran_exe or '.exe' in nastran_exe:
        nastran_exe = Path(nastran_exe)
        assert nastran_exe.exists(), print_bad_path(nastran_exe)
        assert nastran_exe.is_file(), nastran_exe

    cleanup = args.cleanup
    extensions = ['.dat', '.bdf']
    run_jobs(bdf_filename_dirname, nastran_exe,
             extensions=extensions, cleanup=cleanup,
             recursive=True)


def get_bdf_filenames_to_run(bdf_filename_dirname: Path | list[Path],
                             extensions: list[str],
                             recursive: bool=False) -> list[Path]:
    if isinstance(extensions, str):
        extensions = [extensions]
    assert isinstance(extensions, list), extensions
    for ext in extensions:
        assert ext.startswith('.'), extensions

    bdf_filename_dirname_list = bdf_filename_dirname
    if not isinstance(bdf_filename_dirname, list):
        assert isinstance(bdf_filename_dirname, Path), bdf_filename_dirname
        assert bdf_filename_dirname.exists(), bdf_filename_dirname
        bdf_filename_dirname_list = [bdf_filename_dirname]
    del bdf_filename_dirname

    #----------------------------------------
    bdf_filenames: list[Path] = []
    for bdf_filename_dirnamei in bdf_filename_dirname_list:
        assert bdf_filename_dirnamei.exists(), bdf_filename_dirnamei
        if bdf_filename_dirnamei.is_dir():
            dirname = os.path.abspath(bdf_filename_dirnamei)
            if recursive:
                bdf_filenamesi = []
                for ext in extensions:
                    files = get_files_of_type(
                        dirname, extension=ext, max_size=0.)  # no size limit (in MB)
                    bdf_filenamesi += [Path(fname) for fname in files
                                       if '.test_bdf.' not in os.path.basename(fname)]
            else:
                #suffixs = [fname.suffix for fname in dirname.iterdir() if '.test_bdf.' not in fname.name]
                bdf_filenamesi = [dirname / fname.name for fname in dirname.iterdir()
                                  if fname.suffix in extensions and '.test_bdf.' not in fname.name]
            assert len(bdf_filenamesi) > 0, dirname

        elif bdf_filename_dirnamei.is_file():
            bdf_filenamesi = [bdf_filename_dirnamei]
        else:  # pragma: no cover
            raise NotImplementedError(bdf_filename_dirnamei)
        bdf_filenames.extend(bdf_filenamesi)
    for bdf_filename in bdf_filenames:
        assert bdf_filename.exists(), print_bad_path(bdf_filename)
    assert len(bdf_filenames) > 0, bdf_filenames
    return bdf_filenames


def run_jobs(bdf_filename_dirname: Path, nastran_exe: str | Path,
             extensions: list[str], cleanup: bool=True,
             recursive: bool=False,
             run: bool=True) -> int:
    """runs a series of jobs in a specific folder

    Parameters
    ----------
    recursive: bool; default=False
        finds all bdf/dat files in all sub-directories
        NOTE: doesn't apply to files
    """
    bdf_filenames = get_bdf_filenames_to_run(
        bdf_filename_dirname, extensions, recursive=recursive)

    log = SimpleLogger(level='debug')
    bdf_filenames_str = [str(bdf_filename) for bdf_filename in bdf_filenames]
    log.info(f'bdf_filenames = {bdf_filenames_str}')
    for bdf_filename in bdf_filenames:
        assert bdf_filename.exists(), print_bad_path(bdf_filename)

    nfiles = len(bdf_filenames)
    t_run_min = 0.
    t_est_min = 0.
    t_est_hr = 0.
    t0 = time.time()
    for ifile, bdf_filename in enumerate(bdf_filenames):
        base = os.path.splitext(str(bdf_filename))[0]
        op2_filename = base + '.op2'
        if os.path.exists(op2_filename):
            log.warning(f'skipping {str(bdf_filename)} because {op2_filename} already exists')
            continue

        #nfiles_remaining0 = nfiles - ifile
        nfiles_remaining1 = nfiles - (ifile + 1)
        percent0 = ifile / nfiles * 100
        percent1 = (ifile + 1) / nfiles * 100
        log.debug(f'estimated time remaining: {t_est_min:.0f} min = {t_est_hr:.1f} hr; time/run={t_run_min:.1f} min')
        log.info(f'running  {ifile+1}/{nfiles}={percent0:.0f}%: {str(bdf_filename)}')
        run_nastran(bdf_filename, nastran_cmd=nastran_exe,
                    cleanup=cleanup, run=run)
        log.debug(f'finished {ifile+1}/{nfiles}={percent1:.0f}%: {str(bdf_filename)}')

        # if 0:
        dt = time.time() - t0
        # else:
        #     dt = 20. * 60. * (ifile + 1)

        t_run_min = dt / (ifile + 1) / 60
        t_est_sec = dt * nfiles_remaining1 / (ifile + 1)
        t_est_min = t_est_sec / 60.
        t_est_hr = t_est_min / 60.
    log.info('done')
    return nfiles


if __name__ == '__main__':
    cmd_line_run_jobs()
