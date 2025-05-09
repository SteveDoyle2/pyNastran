"""
TODO: support old=no
"""
# pep8: disable=E252
import os
import sys
# import copy
import shlex
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
    r"""
    run_nastran_job dirname
    run_nastran_job fem.bdf -x C:\bin\nastran.exe
    run_nastran_job .       -x C:\bin\nastran.exe --cleanup -r --test
    run_nastran_job .    --outfile files_to_run.out --test -args "old=no mem=10gb"
    run_nastran_job junk --infile  files_to_run.out
    run_nastran_job fem1.bdf fem2.bdf
    run_nastran_job . --skip fem1.bdf fem2.bdf --all
    """
    local_file = os.path.abspath(__file__)
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [local_file] + argv[2:]  # ['run_jobs'] + sys.argv[2:]
    if not quiet:
        print(f'argv = {argv}')
    #print(f'argv = {argv}')

    parser = argparse.ArgumentParser(prog='run_jobs')
    parser.add_argument('bdf_dirname_filename', nargs='+', help='path to Nastran filename/directory')
    #parser.add_argument('-o', '--overwrite', default=False, help='overwrite files')
    parser.add_argument('-x', '--exe', default='nastran', help='path to Nastran executable')
    parser.add_argument('-c', '--cleanup', action='store_true', help='cleanup the junk output files (log, f04, plt)')
    parser.add_argument('-r', '--recursive', action='store_true', help='recursively search for directories')
    parser.add_argument('--skip', nargs='+', help='dont process specific files')
    parser.add_argument('-a', '--all', help='dont skip files that have an op2')
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

    run = args.test
    recursive = args.recursive
    debug = args.debug
    #print(args)
    skip_files = args.skip
    nastran_args = args.args
    process_all = args.all
    if nastran_args is None or len(nastran_args) == 0:
        keywords = []
    else:
        keywords = nastran_args.split()
    extensions = ['.dat', '.bdf']
    if args.infile is not None:
        bdf_filename_dirname, all_keywords_list = load_infile(args.infile, extensions)
        print(f'bdf_filename_dirname = {bdf_filename_dirname}')
    else:
        bdf_filename_dirname = [Path(filenamei) for filenamei in args.bdf_dirname_filename]
        all_keywords_list: list[str] = keywords

    nastran_exe = args.exe
    # if nastran_exe is None:
    #     nastran_exe = 'nastran'
    if '.bat' in nastran_exe or '.exe' in nastran_exe:
        nastran_exe = Path(nastran_exe)
        assert nastran_exe.exists(), print_bad_path(nastran_exe)
        assert nastran_exe.is_file(), nastran_exe

    cleanup = args.cleanup
    level = 'warning' if quiet else 'debug'
    out_filename = '' if args.outfile is None else args.outfile
    nfiles = run_jobs(
        bdf_filename_dirname, nastran_exe,
        extensions=extensions, cleanup=cleanup,
        keywords=all_keywords_list,
        recursive=recursive, run=run,
        out_filename=out_filename,
        skip_files=skip_files,
        process_all=process_all,
        debug=debug, log=level)
    return nfiles


def load_infile(infilename: str,
                extensions: list[str]) -> tuple[list[Path], list[list[str]]]:
    assert os.path.exists(infilename), print_bad_path(infilename)
    with open(infilename, 'r') as txt_file:
        lines = [line.strip() for line in txt_file.readlines()]
    #bdf_filename_dirname = [Path(filenamei) for filenamei in lines]

    all_call_args = []
    all_keywords_list = []
    bdf_filename_dirname = []
    for line in lines:
        line = line.strip()
        # print(line)
        if len(line) == 0:
            continue
        call_args = shlex.split(line, posix=False)
        all_call_args.append(call_args)
        name = ''
        keywords_list = []
        for arg in call_args[1:]:
            base, ext = os.path.splitext(arg)
            ext = ext.lower()
            if ext in extensions:
                assert name == '', (name, arg)
                name = os.path.abspath(arg)
            else:
                keywords_list.append(arg)
        bdf_filename_dirname.append(Path(name))
        all_keywords_list.append(keywords_list)

    return bdf_filename_dirname, all_keywords_list


def get_bdf_filenames_to_run(bdf_filename_dirname: PathLike | list[PathLike],
                             extensions: str | list[str],
                             recursive: bool=False,
                             process_all: bool=False) -> list[Path]:
    if isinstance(extensions, str):
        extensions = [extensions]
    assert isinstance(extensions, list), extensions
    for ext in extensions:
        assert ext.startswith('.'), extensions

    bdf_filename_dirname_list_in = bdf_filename_dirname
    if not isinstance(bdf_filename_dirname, list):
        if isinstance(bdf_filename_dirname, str):
            bdf_filename_dirname = Path(bdf_filename_dirname)
        assert isinstance(bdf_filename_dirname, Path), bdf_filename_dirname
        assert bdf_filename_dirname.exists(), bdf_filename_dirname
        bdf_filename_dirname_list_in = [bdf_filename_dirname]
    elif isinstance(bdf_filename_dirname, str):
        bdf_filename_dirname_list_in = [Path(bdf_filename_dirname)]
    elif isinstance(bdf_filename_dirname, list):
        bdf_filename_dirname_list_in = [Path(pathi) for pathi in bdf_filename_dirname]
    #else:
        #print(type(bdf_filename_dirname_list), bdf_filename_dirname_list)
        #adsf
    del bdf_filename_dirname

    #----------------------------------------
    bdf_filename_dirname_list_out = _deglob(bdf_filename_dirname_list_in, recursive=recursive)
    del bdf_filename_dirname_list_in

    bdf_filenames = _directory_to_files(
        bdf_filename_dirname_list_out, extensions, recursive)
    del bdf_filename_dirname_list_out

    bdf_filenames_run = []
    allow_op2_skip = not process_all
    for bdf_filename in bdf_filenames:
        assert bdf_filename.exists(), print_bad_path(bdf_filename)
        base, ext = os.path.splitext(str(bdf_filename))
        op2_filename = Path(base + '.op2')
        if op2_filename.exists() and allow_op2_skip:
            continue
        bdf_filenames_run.append(bdf_filename)
    assert len(bdf_filenames_run) > 0, bdf_filenames_run
    return bdf_filenames_run


def _deglob(bdf_filename_dirname_list_in: list[Path],
            recursive: bool=False) -> list[Path]:
    """handle globbing"""
    pwd = Path(os.getcwd())
    #print(f"pwd = {str(pwd)}")
    bdf_filename_dirname_list_out: list[Path] = []
    for bdf_filename_dirnamei in bdf_filename_dirname_list_in:
        if '*' in str(bdf_filename_dirnamei):
            print(f'bdf_filename_dirnamei={bdf_filename_dirnamei!r}; recursive={recursive}')
            if recursive:
                outi = pwd.rglob(str(bdf_filename_dirnamei))
            else:
                outi = pwd.glob(str(bdf_filename_dirnamei))
            outi = list(outi)
            for outii in outi:
                print(f'  outi = {outii}')
            bdf_filename_dirname_list_out.extend(outi)
        else:
            bdf_filename_dirname_list_out.append(bdf_filename_dirnamei)
    return bdf_filename_dirname_list_out


def _directory_to_files(bdf_filename_dirname_list: list[PathLike],
                        extensions: list[str], recursive: bool) -> list[Path]:
    bdf_filenames: list[Path] = []
    for bdf_filename_dirnamei in bdf_filename_dirname_list:
        assert bdf_filename_dirnamei.exists(), bdf_filename_dirnamei
        if bdf_filename_dirnamei.is_dir():
            dirname = Path(os.path.abspath(bdf_filename_dirnamei))
            if recursive:
                bdf_filenamesi = []
                for ext in extensions:
                    files = get_files_of_type(
                        dirname, extension=ext, max_size_mb=0.0)  # no size limit (in MB)
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
            #print(f'file {bdf_filename_dirnamei}')
            bdf_filenamesi = [bdf_filename_dirnamei]
        else:  # pragma: no cover
            raise NotImplementedError(bdf_filename_dirnamei)
        bdf_filenames.extend(bdf_filenamesi)
    return bdf_filenames


def run_jobs(bdf_filename_dirname: PathLike | list[PathLike],
             nastran_exe: PathLike,
             extensions: str | list[str],
             cleanup: bool=True,
             recursive: bool=False,
             keywords: Optional[str | list[str] | dict[str, str] | list[list[str]]]=None,
             run: bool=True,
             out_filename: PathLike='',
             skip_files: PathLike | list[PathLike]='',
             process_all: bool=False,
             debug: bool=False,
             log: SimpleLogger | str='debug') -> int:
    """
    Runs a series of jobs in a specific folder with
    specific file extensions

    Parameters
    ----------
    bdf_filename_dirname: PathLike | list[PathLike]
        a series of filenames or directories to process
    nastran_exe: str, list[str]
        the path to Nastran
    extensions: str, list[str]
        the file extensions that will be analyzed (*.dat, *.bdf)
    recursive: bool; default=False
        finds all bdf/dat files in all sub-directories
        NOTE: doesn't apply to files
    keywords : str, list[str], dict[str, Any]
        keywords = 'old=no'
        keywords = ['old=no', 'parallel=8']
        keywords = {'old': 'no', 'parallel': 8}
    out_filename: PathLike; default=''
        file to load previously saved jobs
    skip_files: PathLike | list[skip_files]; default=''
        files that should be skipped
    process_all: bool; default=False
        True:  process all bdf/dat files in all sub-directories
               other than those specified by skip_files
        False: skip files that have op2s
    run: bool; default=True
        lets you disable actually running Nastran to test out code/get the call arguments
    cleanup: bool; default=False
        removes the *.asg, *.asm, *.log, *.f04, *.mon1, *.mon2 and *.plt files
    debug: bool; default=False
        print the call args

    Returns
    -------
    nfiles: int
        the number of files

    TODO: remove failed jobs from the time estimator
    TODO: add out_filenames to output
    """
    if skip_files == '':
        skip_files = []
    if isinstance(skip_files, PathLike):
        skip_files = [skip_files]
    skip_files_path: list[Path] = [Path(pth).absolute() for pth in skip_files]

    #print(f'skip_files = {skip_files}')
    #print(f'**keywords= {keywords}')
    log = SimpleLogger(log) if isinstance(log, str) else log
    bdf_filenames: list[Path] = get_bdf_filenames_to_run(
        bdf_filename_dirname, extensions, recursive=recursive,
        process_all=process_all)

    bdf_filenames_temp = []
    for bdf_filename in bdf_filenames:
        #print(f'bdf_filename={bdf_filename}')
        if bdf_filename in skip_files_path:
            log.warning(f'skipping bdf_filename={bdf_filename}')
            continue
        bdf_filenames_temp.append(bdf_filename)
    bdf_filenames = bdf_filenames_temp

    bdf_filenames_str = [str(bdf_filename) for bdf_filename in bdf_filenames]
    bdf_filenames_str_short = [os.path.basename(bdf_filename) for bdf_filename in bdf_filenames]
    # if show_folder:
    msg = '\n - '.join(bdf_filenames_str)
    # else:
    #     msg = '\n - '.join(bdf_filenames_str_short)
    #log.info(f'bdf_filenames = {bdf_filenames_str}')
    log.info(f'bdf_filenames:\n - {msg}')

    # widthcases = len(str(len(bdf_filenames))) + 1
    # msg2 = ''
    for ifile, bdf_filename, bdf_filenames_str_short in zip(count(), bdf_filenames, bdf_filenames_str_short):
        assert bdf_filename.exists(), print_bad_path(bdf_filename)
    #     if debug:
    #         ifile_str = f'{ifile+1}:'
    #         if show_folder:
    #             msg2 += f'{ifile_str:{widthcases}s} {bdf_filename}\n'
    #         else:
    #             msg2 += f'{ifile_str:{widthcases}s} {bdf_filenames_str_short}\n'
    # if msg2:
    #     log.debug(f'bdf_filenames:\n{msg2}')

    # if out_filename:
    #     with open(out_filename, 'w') as out_file:
    #         for bdf_filename in bdf_filenames:
    #             out_file.write(f'{str(bdf_filename)}\n')
    nfiles, all_call_args = run_jobs_by_filenames(
        bdf_filenames, nastran_exe, keywords, log,
        process_all=process_all, cleanup=cleanup,
        run=run, debug=debug)
    assert isinstance(nfiles, int), nfiles

    if out_filename:
        _write_outfile(out_filename, all_call_args)
    return nfiles


def _write_outfile(out_filename: PathLike,
                   all_call_args: list[list[str]]) -> None:
    # with open(out_filename, 'w') as out_file:
    #     for bdf_filename in bdf_filenames:
    #         out_file.write(f'{str(bdf_filename)}\n')
    with open(out_filename, 'w') as out_file:
        for call_args in all_call_args:
            assert len(call_args) > 0, call_args

            # write the output arg
            out = ''
            for arg in call_args:
                if os.path.exists(arg):
                    arg = os.path.relpath(arg, start='.')
                if isinstance(arg, str) and ' ' in arg:
                    out += f'{arg!r} '
                else:
                    out += f'{arg} '
            out = out.strip()
            #out = str(all_call_args)[1:-1]
            #out_file.write(f'{str(bdf_filename)}\n')
            out_file.write(f'{out}\n')


def run_jobs_by_filenames(bdf_filenames: list[PathLike],
                          nastran_exe: PathLike,
                          keywords: list[str],
                          log: SimpleLogger,
                          process_all: bool=False,
                          cleanup: bool=True,
                          run: bool=True,
                          debug: bool=False) -> tuple[int, list[list[str]]]:
    nfiles = len(bdf_filenames)
    eta = 'N/A'
    eta_next = 'N/A'
    t_run_min = 0.
    t_est_min = 0.
    t_est_hr = 0.
    t0 = time.time()

    all_call_args = []
    msg = f'{nfiles}/{nfiles}=100%:'
    nmsg = len(msg)
    is_keywords_list = False
    if isinstance(keywords, list) and len(keywords) > 0:
        keywords0 = keywords[0]
        is_keywords_list = isinstance(keywords0, list)
    if is_keywords_list:
        #print(f'keywords0 = {keywords0}')
        assert len(keywords) == len(bdf_filenames), f'keywords={keywords} \nbdf_filenames={bdf_filenames}'

    allow_op2_skip = not process_all
    for ifile, bdf_filename in enumerate(bdf_filenames):
        keywordsi = keywords[ifile] if is_keywords_list else keywords
        # print(f'keywords[{ifile}] = {keywordsi}')

        if not os.path.exists(bdf_filename):
            log.warning(f'skipping {str(bdf_filename)!r} because {bdf_filename!r} doesnt exist')
            continue

        base = os.path.splitext(str(bdf_filename))[0]
        op2_filename = base + '.op2'
        if os.path.exists(op2_filename) and allow_op2_skip:
            log.warning(f'skipping {str(bdf_filename)!r} because {op2_filename!r} already exists')
            continue

        #nfiles_remaining0 = nfiles - ifile
        nfiles_remaining1 = nfiles - (ifile + 1)
        percent0 = ifile / nfiles * 100
        percent1 = (ifile + 1) / nfiles * 100

        log.debug(f'ETA:{eta}; time remaining: {t_est_min:.0f} min = {t_est_hr:.1f} hr; '
                  f'time/run={t_run_min:.1f} min; ETA next:{eta_next}')
        msg0 = f'{ifile+1}/{nfiles}={percent0:.0f}%:'
        log.info(f'running  {msg0:<{nmsg}} {str(bdf_filename)}')
        return_code, call_args = run_nastran(bdf_filename, nastran_cmd=nastran_exe,
                                             keywords=keywordsi, cleanup=cleanup, run=run,
                                             debug=debug, log=log)
        msg1 = f'{ifile+1}/{nfiles}={percent1:.0f}%:'
        log.debug(f'finished {msg1:<{nmsg}} {str(bdf_filename)}; return_code={return_code}')
        #time.sleep(5)

        dt = time.time() - t0
        t_run_min = dt / (ifile + 1) / 60.
        t_est_sec = dt * nfiles_remaining1 / (ifile + 1)
        t_est_min = t_est_sec / 60.
        t_est_hr = t_est_min / 60.
        now = datetime.datetime.now()
        new = now + datetime.timedelta(minutes=t_est_min)
        nexti = now + datetime.timedelta(minutes=t_run_min)
        #print(f't_est_total(s) = {t_est_min*60:.6g}')
        #print(f't_est_next(s)  = {t_run_min*60:.6g}')
        eta = new.strftime("%Y-%m-%d %I:%M %p")  # '2025-01-29 05:30 PM'
        eta_next = nexti.strftime("%Y-%m-%d %I:%M %p")  # '2025-01-29 05:30 PM'
        all_call_args.append(call_args)
    log.info('done')
    return nfiles, all_call_args


if __name__ == '__main__':
    cmd_line_run_jobs()
