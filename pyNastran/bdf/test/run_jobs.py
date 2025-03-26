import os
import sys
import datetime
from itertools import count
from pathlib import Path
import time
import argparse
from typing import Optional

from cpylog import SimpleLogger
import pyNastran
from pyNastran.utils import print_bad_path, PathLike
from pyNastran.utils.nastran_utils import run_nastran
from pyNastran.utils.dev import get_files_of_type

def cmd_line_run_jobs(argv=None, quiet: bool=False) -> int:
    """
    run_nastran_job dirname
    run_nastran_job filename.bdf -x C:\bin\nastran.exe
    run_nastran_job .            -x C:\bin\nastran.exe --cleanup -r --test
    run_nastran_job filename.bdf filename2.bdf
    """
    FILE = os.path.abspath(__file__)
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [FILE] + argv[2:]  # ['run_jobs'] + sys.argv[2:]
    if not quiet:
        print(f'argv = {argv}')
    #print(f'argv = {argv}')

    parser = argparse.ArgumentParser(prog='run_jobs')
    parser.add_argument('bdf_dirname_filename', nargs='+', help='path to Nastran filename')
    #parser.add_argument('-o', '--overwrite', default=False, help='overwrite files')
    parser.add_argument('-x', '--exe', default='nastran', help='path to Nastran execuable')
    parser.add_argument('-c', '--cleanup', action='store_true', help='cleanup the junk output files (log, f04, plt)')
    parser.add_argument('-r', '--recursive', action='store_true', help='recursively search for directories')
    parser.add_argument('--nofolder', action='store_true', help='dont show the directory path')
    parser.add_argument('--args', help='additional arguments')

    file_group = parser.add_mutually_exclusive_group(required=False)
    file_group.add_argument('--infile', help='run only files listed in the file; overwrites bdf_dirname_filename')
    file_group.add_argument('--outfile', help='skip run the jobs')

    parser.add_argument('--test', action='store_false', help='skip run the jobs')
    parser.add_argument('--debug', action='store_true', help='more debugging')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                        version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    show_folder = not args.nofolder
    run = args.test
    recursive = args.recursive
    debug = args.debug
    #print(args)
    nastran_args = args.args
    if nastran_args is None or len(nastran_args) == 0:
        keywords = []
    else:
        keywords = nastran_args.split()

    if args.infile is not None:
        assert os.path.exists(args.infile), print_bad_path(args.infile)
        with open(args.infile, 'r') as txt_file:
            lines = [line.strip() for line in txt_file.readlines()]
        bdf_filename_dirname = [Path(filenamei) for filenamei in lines]
    else:
        bdf_filename_dirname = [Path(filenamei) for filenamei in args.bdf_dirname_filename]
    nastran_exe = args.exe
    # if nastran_exe is None:
    #     nastran_exe = 'nastran'
    if '.bat' in nastran_exe or '.exe' in nastran_exe:
        nastran_exe = Path(nastran_exe)
        assert nastran_exe.exists(), print_bad_path(nastran_exe)
        assert nastran_exe.is_file(), nastran_exe

    cleanup = args.cleanup
    extensions = ['.dat', '.bdf']

    level = 'warning' if quiet else 'debug'
    out_filename = '' if args.outfile is None else args.outfile
    nfiles = run_jobs(
        bdf_filename_dirname, nastran_exe,
        extensions=extensions, cleanup=cleanup,
        keywords=keywords, show_folder=show_folder,
        recursive=recursive, run=run,
        out_filename=out_filename,
        debug=debug, log=level)
    return nfiles


def get_bdf_filenames_to_run(bdf_filename_dirname: PathLike | list[PathLike],
                             extensions: str | list[str],
                             recursive: bool=False) -> list[Path]:
    if isinstance(extensions, str):
        extensions = [extensions]
    assert isinstance(extensions, list), extensions
    for ext in extensions:
        assert ext.startswith('.'), extensions

    bdf_filename_dirname_list = bdf_filename_dirname
    if not isinstance(bdf_filename_dirname, list):
        if isinstance(bdf_filename_dirname, str):
            bdf_filename_dirname = Path(bdf_filename_dirname)
        assert isinstance(bdf_filename_dirname, Path), bdf_filename_dirname
        assert bdf_filename_dirname.exists(), bdf_filename_dirname
        bdf_filename_dirname_list = [bdf_filename_dirname]
    elif isinstance(bdf_filename_dirname, str):
        bdf_filename_dirname_list = [Path(bdf_filename_dirname)]
    elif isinstance(bdf_filename_dirname, list):
        bdf_filename_dirname_list = [Path(pathi) for pathi in bdf_filename_dirname]
    #else:
        #print(type(bdf_filename_dirname_list), bdf_filename_dirname_list)
        #adsf
    del bdf_filename_dirname

    #----------------------------------------
    bdf_filenames: list[Path] = []
    for bdf_filename_dirnamei in bdf_filename_dirname_list:
        assert bdf_filename_dirnamei.exists(), bdf_filename_dirnamei
        if bdf_filename_dirnamei.is_dir():
            dirname = Path(os.path.abspath(bdf_filename_dirnamei))
            if recursive:
                bdf_filenamesi = []
                for ext in extensions:
                    files = get_files_of_type(
                        dirname, extension=ext, max_size=0.)  # no size limit (in MB)
                    bdf_filenamesi += [Path(fname) for fname in files
                                       if ('.test_bdf.' not in os.path.basename(fname) and
                                           '.test_op2.' not in os.path.basename(fname))]
            else:
                #suffixs = [fname.suffix for fname in dirname.iterdir() if '.test_bdf.' not in fname.name]
                bdf_filenamesi = [dirname / fname.name for fname in dirname.iterdir()
                                  if (fname.suffix in extensions and
                                      '.test_bdf.' not in fname.name and
                                      '.test_op2.' not in fname.name)]
            assert len(bdf_filenamesi) > 0, dirname

        elif bdf_filename_dirnamei.is_file():
            bdf_filenamesi = [bdf_filename_dirnamei]
        else:  # pragma: no cover
            raise NotImplementedError(bdf_filename_dirnamei)
        bdf_filenames.extend(bdf_filenamesi)

    bdf_filenames_run = []
    for bdf_filename in bdf_filenames:
        assert bdf_filename.exists(), print_bad_path(bdf_filename)
        base, ext = os.path.splitext(str(bdf_filename))
        op2_filename = Path(base + '.op2')
        if op2_filename.exists():
            continue
        bdf_filenames_run.append(bdf_filename)
    assert len(bdf_filenames_run) > 0, bdf_filenames_run
    return bdf_filenames_run



def run_jobs(bdf_filename_dirname: PathLike | list[PathLike],
             nastran_exe: PathLike,
             extensions: str | list[str],
             cleanup: bool=True,
             recursive: bool=False,
             keywords: Optional[str | list[str] | dict[str, str]]=None,
             show_folder: bool=True,
             run: bool=True,
             out_filename: str='',
             debug: bool=False,
             log: SimpleLogger | str='debug') -> int:
    """
    Runs a series of jobs in a specific folder

    Parameters
    ----------
    recursive: bool; default=False
        finds all bdf/dat files in all sub-directories
        NOTE: doesn't apply to files

    Returns
    -------
    nfiles: int
        the number of files
    out_filename: str
        path to file to write list of jobs

    TODO: remove failed jobs from the time estimator
    """
    # log = SimpleLogger(level='debug')
    log = SimpleLogger(log) if isinstance(log, str) else log
    bdf_filenames = get_bdf_filenames_to_run(
        bdf_filename_dirname, extensions, recursive=recursive)
    #print(bdf_filenames, len(bdf_filenames))

    bdf_filenames_str = [str(bdf_filename) for bdf_filename in bdf_filenames]
    bdf_filenames_str_short = [os.path.basename(bdf_filename) for bdf_filename in bdf_filenames]
    if show_folder:
        msg = '\n - '.join(bdf_filenames_str)
    else:
        msg = '\n - '.join(bdf_filenames_str_short)
    log.info(f'bdf_filenames:\n - {msg}')

    widthcases = len(str(len(bdf_filenames))) + 1
    msg2 = ''
    for ifile, bdf_filename, bdf_filenames_str_short in zip(count(), bdf_filenames, bdf_filenames_str_short):
        assert bdf_filename.exists(), print_bad_path(bdf_filename)
        if debug:
            ifile_str = f'{ifile+1}:'
            if show_folder:
                msg2 += f'{ifile_str:{widthcases}s} {bdf_filename}\n'
            else:
                msg2 += f'{ifile_str:{widthcases}s} {bdf_filenames_str_short}\n'
    if msg2:
        log.debug(f'bdf_filenames:\n{msg2}')

    if out_filename:
        with open(out_filename, 'w') as out_file:
            for bdf_filename in bdf_filenames:
                out_file.write(f'{str(bdf_filename)}\n')

    nfiles = len(bdf_filenames)
    eta = 'N/A'
    eta_next = 'N/A'
    t_run_min = 0.
    t_est_min = 0.
    t_est_hr = 0.
    t0 = time.time()
    for ifile, bdf_filename in enumerate(bdf_filenames):
        if not os.path.exists(bdf_filename):
            log.warning(f'skipping {str(bdf_filename)} because {bdf_filename} doesnt exist')
            continue

        base = os.path.splitext(str(bdf_filename))[0]
        op2_filename = base + '.op2'
        if os.path.exists(op2_filename):
            log.warning(f'skipping {str(bdf_filename)} because {op2_filename} already exists')
            continue

        #nfiles_remaining0 = nfiles - ifile
        nfiles_remaining1 = nfiles - (ifile + 1)
        percent0 = ifile / nfiles * 100
        percent1 = (ifile + 1) / nfiles * 100

        log.debug(f'ETA:{eta}; time remaining: {t_est_min:.0f} min = {t_est_hr:.1f} hr; time/run={t_run_min:.1f} min; ETA next:{eta_next}')
        log.info(f'running  {ifile+1}/{nfiles}={percent0:.0f}%: {str(bdf_filename)}')
        return_code, call_args = run_nastran(bdf_filename, nastran_cmd=nastran_exe,
                                             keywords=keywords, cleanup=cleanup, run=run,
                                             debug=debug, log=log)
        log.debug(f'finished {ifile+1}/{nfiles}={percent1:.0f}%: {str(bdf_filename)}; return_code={return_code}')

        # if 0:
        dt = time.time() - t0
        # else:
        #     dt = 20. * 60. * (ifile + 1)

        t_est_next = dt / 60.
        t_run_min = dt / (ifile + 1) / 60.
        t_est_sec = dt * nfiles_remaining1 / (ifile + 1)
        t_est_min = t_est_sec / 60.
        t_est_hr = t_est_min / 60.
        new = datetime.datetime.now() + datetime.timedelta(minutes=t_est_min)
        nexti = datetime.datetime.now() + datetime.timedelta(minutes=t_est_next)
        eta = new.strftime("%Y-%m-%d %I:%M %p")  # '2025-01-29 05:30 PM'
        eta_next = nexti.strftime("%Y-%m-%d %I:%M %p")  # '2025-01-29 05:30 PM'
    log.info('done')
    return nfiles


if __name__ == '__main__':
    cmd_line_run_jobs()
