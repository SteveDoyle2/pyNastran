import os
import shutil
import warnings
import subprocess
from typing import Optional
from pyNastran.utils import PathLike, print_bad_path
from cpylog import SimpleLogger

KEYWORDS_CHECK = {
    'old': ['yes', 'no'],     # throw away old results
    'bat': ['yes', 'no'],     # linux only
    'news': ['yes', 'no', 'auto'],    # doesn't seem to work, but shrinks f06 size
    'notify': ['yes', 'no'],  # make nastran beep
    'fpe': ['yes', 'no'],  # ???
    'msgbell': ['yes', 'no'],  # ???
    'scr': ['yes', 'no', 'mini', 'post'],
    'endian': ['big', 'little'],
    #'mem=15gb'
    #'parallel=4'
}
# auth = 28000@LocalHost
# memory = 0.45*physical
# smem = 20.0X
# buffpool = 20.0X
# ishellpath = %NXN_BASE%\bin
# diag = 8
# program=FEMAP


def run_nastran(bdf_filename: PathLike,
                nastran_cmd: PathLike='nastran',
                keywords: Optional[str | list[str] | dict[str, str]]=None,
                run: bool=True,
                run_if_op2_exists: bool=True,
                run_in_bdf_dir: bool=True,
                cleanup: bool=False,
                debug: bool=False,
                log: Optional[SimpleLogger]=None) -> tuple[Optional[int], list[str]]:
    """
    Call a nastran subprocess with the given filename

    Parameters
    ----------
    bdf_filename: string
        Filename of the Nastran .bdf file
    nastran_cmd: str; default='nastran'
        path to Nastran executable or a command on the PATH
        (e.g., 'nastran' finds nastran.bat/nastran.exe on Windows)
    keywords: str/dict/list of strings, optional
        Default keywords are `'scr=yes'`, `'bat=no'`, `'old=no'`, and `'news=no'`
    run: bool; default=True
        lets you disable actually running Nastran to test out code/get the call arguments
    run_in_bdf_dir: bool; default=True
        True  : outputs (e.g., *.f06) will go to the input BDF directory (default)
        False : outputs (e.g., *.f06) will go to the current working directory
    cleanup: bool; default=False
        removes the *.asg, *.asm, *.log, *.f04, *.mon1, *.mon2 and *.plt files

    Returns
    -------
    return_code : int
        the nastran flag
    cmd_args : list[str]
        the nastran commands that go into subprocess

    Example
    -------
    keywords_str = 'scr=yes old=no mem=1024mb'
    keywords_list = ['scr=yes', 'old=no', 'mem=1024mb']
    keywords_dict = {
        'scr' : 'yes',
        'old' : 'no',
        'mem' : '1024mb',
    }
    run_nastran(bdf_filename, keywords_str)
    run_nastran(bdf_filename, keywords_list)
    run_nastran(bdf_filename, keywords_dict)

    """
    keywords_list = _get_keywords_list(
        keywords=keywords, log=log)
    assert os.path.exists(bdf_filename), print_bad_path(bdf_filename)

    # resolve everything before any chdir so relative paths stay valid
    bdf_filename_str = os.path.abspath(bdf_filename)
    bdf_directory = os.path.dirname(bdf_filename_str)
    base = os.path.splitext(os.path.basename(bdf_filename_str))[0]
    op2_filename = os.path.splitext(bdf_filename_str)[0] + '.op2'
    nastran_cmd_str = _resolve_nastran_cmd(nastran_cmd, required=run)

    # match the command line ("nastran old=no model.bdf" from the bdf dir)
    bdf_arg = os.path.basename(bdf_filename_str) if run_in_bdf_dir else bdf_filename_str
    call_args = [nastran_cmd_str, bdf_arg] + keywords_list
    if '=' in os.path.basename(bdf_filename_str):
        msg1 = f'Nastran cant run files with an = sign in them; {bdf_filename_str}'
        if log:
            log.debug(msg1)
        else:
            warnings.warn(msg1)
    if debug:
        msg2 = f'call_args = {call_args}'
        if log:
            log.debug(msg2)
        else:
            print(msg2)

    pwd = os.getcwd()
    switch_dir = run_in_bdf_dir
    if switch_dir:
        os.chdir(bdf_directory)

    return_code = None
    try:
        if run and (not os.path.exists(op2_filename) or run_if_op2_exists):
            return_code = subprocess.call(call_args)

        if run and cleanup:
            # outputs land in the current working directory
            exts = ['.asg', '.asm', '.f04', '.log', '.plt', '.mon1', '.mon2']
            for ext in exts:
                fname = base + ext
                if os.path.exists(fname):
                    print(f'removing {fname}')
                    os.remove(fname)
    finally:
        if switch_dir:
            os.chdir(pwd)
    return return_code, call_args


def _resolve_nastran_cmd(nastran_cmd: PathLike, required: bool=True) -> str:
    """
    Gets the full path to the Nastran executable.

    An explicit path is used as-is. Otherwise, the PATH is searched with
    shutil.which, which honors PATHEXT on Windows, so 'nastran' finds
    nastran.bat. subprocess with shell=False only appends '.exe', so
    passing the bare 'nastran' would raise FileNotFoundError.
    """
    nastran_cmd_str = str(nastran_cmd)
    if os.path.isfile(nastran_cmd_str):
        return os.path.abspath(nastran_cmd_str)

    exe = shutil.which(nastran_cmd_str)
    if exe is not None:
        return exe
    if required:
        raise FileNotFoundError(
            f'nastran_cmd={nastran_cmd_str!r} is not a file and was not found on the PATH')
    return nastran_cmd_str


def _get_keywords_list(keywords: Optional[str | list[str] |
                                          dict[str, str]]=None,
                       log: Optional[SimpleLogger]=None) -> list[str]:
    if keywords is None:
        keywords_list = ['scr=yes', 'bat=no', 'old=no', 'news=no', 'notify=no']  # 'mem=1024mb',
    elif isinstance(keywords, str):
        keywords_list = keywords.split()
    else:
        if isinstance(keywords, (list, tuple)):
            keywords_list = list(keywords)
        else:
            keywords_list = []
            for keyword, value in keywords.items():
                if value is None:
                    continue
                keywords_list.append(f'{keyword}={value}')

    #print(f'keywords_list = {keywords_list}')
    allowed_keywords = list(KEYWORDS_CHECK) + [
        'parallel', 'mem', 'sdirectory', 'sdir',
        'buffsize']  # CHEMIN={path}
    allowed_keywords.sort()
    for keyword_value in keywords_list:
        keyword, value = keyword_value.split('=', 1)
        keyword = keyword.lower()
        if keyword in {'parallel', 'buffsize'}:
            # parallel=number of processors (should check this...)
            # buffsize: default=32769   (size of buffer; larger=better by 2x)
            # authqueue: default=100000 (time in queue; seconds?)
            value_int = int(value)
        elif keyword == 'mem':
            assert value.upper().endswith(('MB', 'GB')), value
            num_str = value[:-2]
            num = float(num_str)
        elif keyword in {'sdirectory', 'sdir'}:  # Directory for scratch files
            path = value
            assert os.path.exists(path), print_bad_path(path)
        elif keyword in {'auth', 'authorize'}:
            pass
        elif keyword in KEYWORDS_CHECK:
            expected_values = KEYWORDS_CHECK[keyword]
            assert value.lower() in expected_values, f'keyword={keyword} value={value}; expected={expected_values}'
        else:
            msg = f'unsupported keyword; {keyword!r}={value!r}; allowed={allowed_keywords}'
            if log:
                log.warning(msg)
            else:
                warnings.warn(msg)
    return keywords_list
