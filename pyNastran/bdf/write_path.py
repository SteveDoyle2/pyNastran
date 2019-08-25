"""
Defines following useful methods:
  - write_include(filename, is_windows=True)

"""
import sys
from typing import Tuple
from pathlib import PurePosixPath, PureWindowsPath


def write_include(filename: str, is_windows: bool=None) -> str:
    """
    Writes a bdf INCLUDE file line given an imported filename.

    Parameters
    ----------
    filename : str
        the filename to write
    is_windows : bool; default=None
        True/False : Windows has a special format for writing INCLUDE
            files, so the format for a BDF that will run on Linux and
            Windows is different.
        None : Check the platform

    For a model that will run on Linux:

    ..code-block:: python

      fname = r'/opt/NASA/test1/test2/test3/
      test4/formats/pynastran_v0.6/pyNastran/bdf/model.inc'
      write_include(fname, is_windows=False)

    We want:

    ..code-block:: python

      INCLUDE /opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/
              pyNastran/bdf/model.inc

    """
    is_windows = is_windows if is_windows is not None else sys.platform in ['win32']

    msg = "INCLUDE '"  # len=8
    #nmax = 72 - 8 # 64

    if is_windows:
        marker = '\\'
    else:
        marker = '/'

    sline = _split_path(filename, is_windows)
    #print('sline =', sline)
    nsline = len(sline)
    if len(filename) > 52: # 62
        pth = ''
        for isline, pathi in enumerate(sline):
            if pathi == '/':  # /home/etc -> [/, home, etc]
                pth += '%s' % marker
            else:
                pth += '%s%s' % (pathi, marker)
            if len(pth) > 52:
                if isline == nsline - 1: # if this is the last one...
                    break
                    #pth = pth.rstrip(marker)
                    #msg += pth
                    #pth = '\n'
                else:
                    pth += '\n        '
                    msg += pth
                    pth = ''
    else:
        pth = marker.join(sline)
    out = msg + pth.rstrip('\n ' + marker) + "'\n"
    return out


def _split_path(abspath: str, is_windows: bool) -> Tuple[str, ...]:
    """
    Takes a path and splits it into the various components.

    This is a helper method for ``write_include``

    """
    if is_windows:
        parts = PureWindowsPath(abspath).parts
    else:
        parts = PurePosixPath(abspath).parts
    return parts
