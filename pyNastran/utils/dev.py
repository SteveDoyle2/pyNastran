import os
from typing import List, Any

import numpy as np

#from pyNastran.utils import object_attributes


def get_files_of_type(dirname: str, extension: str='.txt',
                      max_size: float=100., limit_file: str='no_dig.txt') -> List[str]:
    """
    Gets the list of all the files with a given extension in the specified directory

    Parameters
    ----------
    dirname : str
        the directory name
    extension : str; default='.txt'
        list of filetypes to get
    max_size : float; default=100.0
        size in MB for max file size
    limit_file : str; default=no_dig.txt
        the presence of this file indicates no folder digging
        should be done on this folder

    Returns
    -------
    files : List[str]
        list of all the files with a given extension in the specified directory

    """
    if not os.path.exists(dirname):
        return []

    filenames2 = []  # type: List[str]
    filenames = os.listdir(dirname)
    allow_digging = True
    if limit_file in filenames:
        allow_digging = False
    for filenamei in filenames:
        filename = os.path.join(dirname, filenamei)
        if os.path.isdir(filename):
            if allow_digging:
                filenames2 += get_files_of_type(filename, extension, max_size)
                #assert len(filenames2) > 0, dirnamei
            else:
                print('no digging in filename=%s; dirname=%s' % (filename, dirname))
        elif (os.path.isfile(filename) and
              os.path.splitext(filenamei)[1].endswith(extension) and
              os.path.getsize(filename) / 1048576. <= max_size):
            filenames2.append(filename)
    return filenames2


def list_print(lst: List[Any], float_fmt: str='%-4.2f') -> str:
    """
    Prints a list or numpy array in an abbreviated format.
    Supported element types: None, string, numbers. Useful for debugging.

    Parameters
    ----------
    lst : list / numpy array
        the value to print

    Returns
    -------
    msg : str
        the clean string representation of the object
    """
    def _print(val):
        if val is None or isinstance(val, str):
            return str(val)
        if isinstance(val, float):
            return float_fmt % val
        try:
            return '%g' % val
        except TypeError:
            print("parameter = %r" % val)
            raise

    try:
        # TODO: remove try block and fix bug in OP2 code or add a warning message
        if len(lst) == 0:
            return '[]'

        if isinstance(lst, (np.ndarray)) and lst.ndim == 2:
            row, col = lst.shape
            return (
                "["+",\n ".join(["["+",".join(
                    [float_fmt % lst[i, j]
                     for j in range(col)]) + "]" for i in range(row)])+"]")
        return "[" + ", ".join([_print(a) for a in lst]) + "]"
    except: # not a list
        return _print(lst)
