import os
import sys
import time
import shlex
import shutil
import datetime
import getpass
import subprocess
import argparse
from pathlib import Path
from typing import Optional

from cpylog import SimpleLogger, FileLogger
import pyNastran
from pyNastran.utils import print_bad_path, PathLike
#from pyNastran.bdf.test.run_jobs import run_jobs_by_filenames

def cmd_line_host_jobs(argv=None, quiet: bool=False) -> int:
    """
    host_jobs dirname
    host_jobs dirname1 dirname2  # TODO: add me
    bdf host_jobs . --test
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
    parser.add_argument('host_dirname', nargs='+', help='path to Nastran filename')
    #parser.add_argument('-o', '--overwrite', default=False, help='overwrite files')
    #parser.add_argument('-x', '--exe', default='nastran', help='path to Nastran execuable')
    #parser.add_argument('-c', '--cleanup', action='store_true', help='cleanup the junk output files (log, f04, plt)')
    #parser.add_argument('--args', help='additional arguments')
    parser.add_argument('--test', action='store_false', help='skip run the jobs')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                        version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    show_folder = False #not args.nofolder
    run = args.test
    debug = True  # args.debug
    #print(args)
    #nastran_args = args.args
    #if nastran_args is None or len(nastran_args) == 0:
        #keywords = []
    #else:
        #keywords = nastran_args.split()

    host_dirname = args.host_dirname

    #nastran_exe = args.exe
    # if '.bat' in nastran_exe or '.exe' in nastran_exe:
    #     nastran_exe = Path(nastran_exe)
    #     assert nastran_exe.exists(), print_bad_path(nastran_exe)
    #     assert nastran_exe.is_file(), nastran_exe

    #cleanup = args.cleanup
    #extensions = ['.dat', '.bdf']

    level = 'warning' if quiet else 'debug'
    #out_filename = '' if args.outfile is None else args.outfile
    username = getpass.getuser().lower()
    hosting_filenames = []
    for dirname in host_dirname:
        hosting_filename = os.path.join(dirname, f'hosting_{username}.txt')
        hosting_filenames.append(hosting_filename)
        with open(hosting_filename, 'w') as hosting_file:
            pass
    try:
        nfiles = host_jobs(
            host_dirname,
            #extensions=extensions,
            #cleanup=cleanup,
            #show_folder=show_folder,
            run=run,
            #out_filename=out_filename,
            debug=debug, log=level)
    except:
        for hosting_filename in hosting_filenames:
            remove_file(hosting_filename)
        nfiles = 0

    for hosting_filename in hosting_filenames:
        remove_file(hosting_filename)
    return nfiles


def host_jobs(host_dirnames: PathLike,  # | list[PathLike],
              #extensions: str | list[str],
              cleanup: bool=True,
              #recursive: bool=False,
              #show_folder: bool=True,
              run: bool=True,
              out_filename: str='',
              debug: bool=False,
              log: SimpleLogger | str='debug') -> int:
    """
    jobs_[mycompy].txt
      - a.bdf
      - b.bdf
    _jobs_[mycompy].txt
      - completed jobs
    exe_[mycompy].txt
      - # comment
      - nastran: C:\bin\nastran.exe

    TODO: Arguments (e.g., scr=yes old=no is not supported)
    TODO: Only supports nastran
    TODO: support logging per job
    """
    host_dirname = Path(host_dirnames[0])
    assert isinstance(host_dirname, PathLike), host_dirname
    user = getpass.getuser().lower()
    exe_paths_filename = host_dirname / f'exe_{user}.txt'

    tag = f'jobs_{user}'
    print(f'tag = {tag}')

    i = 0
    exe_paths_dict = {}
    exe_path_timestamp = ''
    while 1:
        print('-' * 20)
        print(f'i = {i}')
        #-------------------------------------------------------------
        # update the exe files
        try:
            exe_paths_dict, exe_path_timestamp = load_exe_paths(
                exe_paths_filename, exe_paths_dict, exe_path_timestamp)
        except:
            continue

        #-------------------------------------------------------------
        # get the list of files
        files_to_check = os.listdir(host_dirname)
        files_to_run = get_files_to_run(files_to_check, tag)
        print(f'files_to_check = {files_to_check}')

        for fname in files_to_run:
            print(f'fname = {fname!r}')

            log_filename = host_dirname / f'{fname}.log'
            log = FileLogger(
                 level='debug', #nlevels: int=1,
                 filename=log_filename,
                 mode='w',
                 include_stream=True)
            check_paths(exe_paths_dict, log)

            log.info(f'fname to run = {fname}')
            run_filename1 = os.path.join(host_dirname, fname)
            run_filename2 = os.path.join(host_dirname, '_' + fname)
            command_line_args_list = load_input_filenames(run_filename1)

            for iarg, args in enumerate(command_line_args_list):
                log.info(f'  args[{iarg}] = {args}')
            log.debug(f'  rename {fname} -> _{fname}')

            remove_file(run_filename2)
            if run:
                rename_file(run_filename1, run_filename2)
            nfiles = run_jobs_by_filenames(
                command_line_args_list, exe_paths_dict, log, #cleanup=cleanup,
                run=run)
            #asdf
        time.sleep(2)
        i += 1
        # if i == 10:
        #     break

def get_files_to_run(files_to_check: list[PathLike],
                     tag: str,) -> list[PathLike]:
    files_to_run = []
    for fname in files_to_check:
        if fname.startswith(tag) and fname.endswith('.txt'):
            files_to_run.append(fname)
    return files_to_run

def check_paths(exe_paths_dict: dict[str, PathLike], log: SimpleLogger) -> None:
    for key, path in exe_paths_dict.items():
        if not os.path.exists(path):
            log.warning(f'{str(path)!r} does not exist')
        else:
            log.info(f'{str(path)!r} exists')


def remove_file(filename: PathLike) -> None:
    if not os.path.exists(filename):
        return
    try:
        os.remove(filename)
    except:
        print(f'cant remove {filename!r}')
        raise
    return


def rename_file(filename1: PathLike, filename2: PathLike) -> None:
    try:
        os.rename(filename1, filename2)
    except:
        print(f'cant rename {filename1!r} to {filename2!r}')
        raise
    return


def load_input_filenames(run_filename1: PathLike) -> list[list[str]]:
    with open(run_filename1, 'r') as file_obj:
        lines = file_obj.readlines()
    lines = [line.strip() for line in lines]

    input_filenames = []
    input_args = []
    for line in lines:
        #print(line)
        # "nastran 'a.bdf'      old=no"
        # -> ['nastran', 'a.bdf', 'old=no']
        # "nastran 'file 2.bdf' old=no"
        # -> ['nastran', 'file 2.bdf', 'old=no']
        args = shlex.split(line)
        for arg in args:
            base, ext = os.path.splitext(arg)
            ext = ext.lower()
            if ext == '' or ext == '.exe':
                continue
            #print(f'  base={base!r} ext={ext!r}')
        #print(f'args = {args}')
        input_args.append(args)
    return input_args


def load_exe_paths(exe_paths_filename: PathLike,
                   exe_paths_dict_old: dict[str, str],
                   exe_path_timestamp: str) -> tuple[dict[str, str], str]:
    """
    Loads the paths to the different exe's (e.g., nastran, python).

    Parameters
    ----------
    exe_paths_filename: Path
        # comment
        key: path
    exe_paths_dict_old: dict[str, str]
        mapping of key (e.g., nastran) to the path
    exe_path_timestamp: str
        the latest loaded time of the file

    Returns
    -------

    """
    print(f'exe_path_timestamp={exe_path_timestamp!r} check')
    # if exe_path_timestamp == '':
    #     pass
    # else:
    #     print(f'exe_path_timestamp={exe_path_timestamp!r} exists...')

    exe_paths_dict = {}
    if not os.path.exists(exe_paths_filename):
        print(f'{str(exe_paths_filename)} does not exist')
        return exe_paths_dict, ''

    modified_time = time.ctime(os.path.getmtime(exe_paths_filename))
    print(f'modified_time = {modified_time}; type={type(modified_time)}')
    if modified_time == exe_path_timestamp:
        return exe_paths_dict_old, exe_path_timestamp

    #print("created: %s" % time.ctime(os.path.getctime(file)))
    with open(exe_paths_filename, 'r') as file_obj:
        lines = file_obj.readlines()
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            continue
        # nastran: path
        name, path = line.split(':', 1)
        exe_paths_dict[name] = path.strip()
    print(f'exe_paths_dict = {exe_paths_dict}')
    return exe_paths_dict, modified_time


def run_jobs_by_filenames(command_line_args: list[list[str]],
                          exe_paths_dict: dict[str, str],
                          log: SimpleLogger,
                          #cleanup: bool=True,
                          run: bool=True,
                          debug: bool=False) -> int:
    #log.debug(f'exe_paths_dict = {exe_paths_dict}')
    log.info('-'*80)
    nfiles = len(command_line_args)
    eta = 'N/A'
    eta_next = 'N/A'
    t_run_min = 0.
    t_est_min = 0.
    t_est_hr = 0.
    t0 = time.time()

    msg = f'{nfiles}/{nfiles}=100%:'
    nmsg = len(msg)
    #assert run is False, run
    for ifile, call_args in enumerate(command_line_args):
        # update nastran -> C:\bin\nastran.bat
        call0 = call_args[0]
        call0_new = exe_paths_dict.get(call0, call0)
        if call0 != call0_new:
            call_args[0] = call0_new

        #nfiles_remaining0 = nfiles - ifile
        nfiles_remaining1 = nfiles - (ifile + 1)
        percent0 = ifile / nfiles * 100
        percent1 = (ifile + 1) / nfiles * 100

        log.debug(f'ETA:{eta}; time remaining: {t_est_min:.0f} min = {t_est_hr:.1f} hr; time/run={t_run_min:.1f} min; ETA next:{eta_next}')
        msg0 = f'{ifile+1}/{nfiles}={percent0:.0f}%:'
        log.info(f'running  {msg0:<{nmsg}} {str(call_args)}')

        return_code = None
        if run:
            return_code = subprocess.call(call_args)
        msg1 = f'{ifile+1}/{nfiles}={percent1:.0f}%:'
        log.debug(f'finished {msg1:{nmsg}} {str(call_args)}; return_code={return_code}')
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
    log.info('done')
    return nfiles


def main():
    dirname = '.'
    args = ['bdf', 'host_jobs', dirname]
    cmd_line_host_jobs(argv=args)


if __name__ == '__main__':
    main()
