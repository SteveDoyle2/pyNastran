import os
import subprocess
from typing import List, Dict, Union, Optional, Any, Tuple


def run_nastran(bdf_filename, nastran_cmd='nastran', keywords=None,
                run=True, run_in_bdf_dir=True):
    # type: (str, str, Optional[Union[List[str], Dict[str, str]]], bool, bool) -> Tuple(int, List[str])
    """
    Call a nastran subprocess with the given filename

    Parameters
    ----------
    bdf_filename : string
        Filename of the Nastran .bdf file
    keywords : dict/list of strings, optional
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
    """
    if keywords is None:
        keywords_list = ['scr=yes', 'bat=no', 'old=no','news=no'] # 'mem=1024mb',
    else:
        if isinstance(keywords, (list, tuple)):
            keywords_list = keywords
        else:
            keywords_list = []
            for keyword, value in keywords.items():
                if value is None:
                    continue
                keywords_list.append('%s=%s' % (keyword, value))

    pwd = os.getcwd()
    bdf_directory = os.path.dirname(bdf_filename)
    if run_in_bdf_dir:
        os.chdir(bdf_directory)
    call_args = [nastran_cmd, bdf_filename] + keywords_list
    return_code = None
    if run:
        return_code = subprocess.call(call_args)

    if run_in_bdf_dir:
        os.chdir(pwd)
    return return_code, call_args
