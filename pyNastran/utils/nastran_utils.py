import os
import subprocess
from typing import List, Dict, Union, Optional, Tuple


def run_nastran(bdf_filename: str, nastran_cmd: str='nastran',
                keywords: Optional[Union[str, List[str], Dict[str, str]]]=None,
                run: bool=True, run_in_bdf_dir: bool=True) -> Tuple[Optional[int], List[str]]:
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

    Returns
    -------
    return_code : int
        the nastran flag
    cmd_args : List[str]
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
    call_args = [nastran_cmd, str(bdf_filename)] + keywords_list
    return_code = None
    if run:
        return_code = subprocess.call(call_args)

    if switch_dir:
        os.chdir(pwd)
    return return_code, call_args

def _get_keywords_list(keywords: Optional[Union[str,
                                                List[str],
                                                Dict[str, str]]]=None) -> List[str]:
    if keywords is None:
        keywords_list = ['scr=yes', 'bat=no', 'old=no', 'news=no'] # 'mem=1024mb',
    elif isinstance(keywords, str):
        keywords_list = keywords.split()
    else:
        if isinstance(keywords, (list, tuple)):
            keywords_list = keywords
        else:
            keywords_list = []
            for keyword, value in keywords.items():
                if value is None:
                    continue
                keywords_list.append(f'{keyword}={value}')
    return keywords_list
