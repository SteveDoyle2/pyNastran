"""
Defines following useful methods:
  - write_include(filename, is_windows=True)

"""
import sys
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

    assert isinstance(filename, str), f'filename={filename} is not a string'
    if is_windows:
        marker = '\\'
    else:
        marker = '/'

    sline = _split_path(filename, is_windows)

    if len(filename) <= 62:
        # short path
        #
        # 62; 72-10=62
        # we need space for:
        #  - INCLUDE (9)
        #  - single quote at the end
        pth = marker.join(sline)
        pth2 = pth.rstrip('\n ' + marker)
        if not is_windows and pth2.startswith(marker):
            ## '//opt/work/fem.bdf' -> '/opt/work/fem.bdf'
            pth2 = marker + pth2.lstrip(marker)
        out = "INCLUDE '" + pth2 + "'\n"
    else:
        # multiline path
        pth = "INCLUDE '"
        all_paths = []
        for isline, pathi in enumerate(sline):
            if pathi == '/':  # /home/etc -> [/, home, etc]
                next_term = marker
            else:
                next_term = '%s%s' % (pathi, marker)

            if len(pth + next_term) < 71:  # we need space for the quote
                pth += next_term
            else:
                pth += '\n'
                all_paths.append(pth)
                pth = '        ' + next_term
                assert len(pth) < 71, pth

        if len(pth):
            all_paths.append(pth)
        pth = ''.join(all_paths).rstrip('\n \\')
        out = pth.rstrip('\n ' + marker) + "'\n"
    return out


def _split_path(abspath: str, is_windows: bool) -> tuple[str, ...]:
    """
    Takes a path and splits it into the various components.

    This is a helper method for ``write_include``

    """
    if is_windows:
        parts = PureWindowsPath(abspath).parts
    else:
        parts = PurePosixPath(abspath).parts
    return parts
