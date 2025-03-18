import os
import warnings
import subprocess
from typing import Optional
from pyNastran.utils import PathLike


KEYWORDS_CHECK = {
    'old': ['yes', 'no'],
    'bat': ['yes', 'no'],
    'news': ['yes', 'no'],
    'notify': ['yes', 'no'],
    #'mem=15gb'
    #'parallel=4'
}
def run_nastran(bdf_filename: PathLike,
                nastran_cmd: PathLike='nastran',
                keywords: Optional[str | list[str] | dict[str, str]]=None,
                run: bool=True, run_in_bdf_dir: bool=True,
                cleanup: bool=False,
                debug: bool=False) -> tuple[Optional[int], list[str]]:
    """
    Call a nastran subprocess with the given filename

    Parameters
    ----------
    bdf_filename : string
        Filename of the Nastran .bdf file
    keywords : str/dict/list of strings, optional
        Default keywords are `'scr=yes'`, `'bat=no'`, `'old=no'`, and `'news=no'`
    run : bool; default=True
        let's you disable actually running Nastran to test out code/get the call arguments
    run_in_local_dir : bool; default=True
        True : output (e.g., *.f06) will go to the current working directory (default)
        False : outputs (e.g., *.f06) will go to the input BDF directory
    cleanup : bool; default=False
        remove the *.asg, *asm, *.log, *.f04, and *.plt files

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
    keywords_list = _get_keywords_list(keywords=keywords)

    pwd = os.getcwd()
    bdf_directory = os.path.dirname(bdf_filename)
    switch_dir = bdf_directory != '' and run_in_bdf_dir
    if switch_dir:
        os.chdir(bdf_directory)

    assert os.path.exists(bdf_filename), bdf_filename
    if nastran_cmd != 'nastran':
        assert os.path.exists(nastran_cmd), nastran_cmd
    nastran_cmd_str = str(nastran_cmd)
    call_args = [nastran_cmd_str, str(bdf_filename)] + keywords_list
    if '=' in str(bdf_filename):
        warnings.warn(f'Nastran cant run files with an = sign in them; {str(bdf_filename)}')
    if debug:
        print(f'call_args = {call_args}')
    return_code = None
    if run:
        return_code = subprocess.call(call_args)

    if run and cleanup:
        base = os.path.basename(bdf_filename)[0]
        fnames = [
            base + '.asg',
            base + '.asm',
            base + '.f04',
            base + '.log',
            base + '.plt',
            base + '.mon1',
            base + '.mon2',
        ]
        for fname in fnames:
            if os.path.exists(fname):
                print(f'removing {fname}')
                os.remove(fname)

    if switch_dir:
        os.chdir(pwd)
    return return_code, call_args

def _get_keywords_list(keywords: Optional[str | list[str] |
                                          dict[str, str]]=None) -> list[str]:
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

    for keyword_value in keywords_list:
        keyword, value = keyword_value.upper().split('=', 1)
        if keyword == 'parallel':
            value_int = int(value)
        elif keyword == 'mem':
            assert value.endswith(('MB', 'GB')), value
            num_str = value[:-2]
            num = float(num_str)
        elif keyword in KEYWORDS_CHECK:
            expected_values = KEYWORDS_CHECK[keyword]
            assert value.lower() in expected_values, f'keyword={keyword} value={value}; expected={expected_values}'
        else:
            warnings.warn(f'unsupported keyword; {keyword!r}={value!r}')

    return keywords_list
