import subprocess
from six import iteritems
import typing


def run_nastran(fname, nastran_cmd='nastran', keywords=None, run=True):
    # type: (str, Dict[str, str], bool) -> float
    """
    Call a nastran subprocess with the given filename

    Parameters
    -----------
    fname : string
        Filename of the Nastran .bdf file
    keywords : dict/list of strings, optional
        Default keywords are `'scr=yes'`, `'bat=no'`, `'old=no'`, and `'news=no'`

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
            for keyword, value in iteritems(keywords):
                if value is None:
                    continue
                keywords_list.append('%s=%s' % (keyword, value))

    call_args = [nastran_cmd, fname] + keywords_list
    return_code = None
    if run:
        return_code = subprocess.call(call_args)
    return return_code, call_args

