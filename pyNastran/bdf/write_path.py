"""
Defines following useful methods:
  - write_include
"""
from __future__ import print_function
import os


def write_include(filename, is_windows=True):
    """
    Writes a bdf INCLUDE file line given an imported filename.

    Parameters
    ----------
    filename : str
        the filename to write
    is_windows : bool; default=True
        Windows has a special format for writing INCLUDE files, so the
        format for a BDF that will run on Linux and Windows is
        different.  We could check the platform, but since you might
        need to change platforms, it's an option

    For a model that will run on Linux:

    ..code-block:: python

      fname = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/model.inc'
      write_include(fname, is_windows=False)

    We want:

    ..code-block:: python

      INCLUDE /opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/
              pyNastran/bdf/model.inc
    """
    msg = 'INCLUDE '  # len=8
    #nmax = 72 - 8 # 64

    if is_windows:
        marker = '\\'
    else:
        marker = '/'

    sline = _split_path(filename)
    #print("sline = %s" % sline)
    if len(filename) > 52: # 62
        pth = ''
        for p in sline:
            if p == '/':
                pth += '%s' % marker
            else:
                pth += '%s%s' % (p, marker)
            if len(pth) > 52:
                pth += '\n        '
                msg += pth
                pth = ''
        pth = pth.rstrip(' /\\\n')
    else:
        pth = marker.join(sline)
    return msg + pth


def _split_path(abspath):
    """
    Takes a path and splits it into the various components.

    This is a helper method for write_include
    """
    path = abspath

    basepaths = []
    npaths_old = -1
    npaths_new = 0
    dpaths = 1
    while dpaths:
        npaths_old = len(basepaths)
        basepath = os.path.basename(path)
        if basepath:
            basepaths.append(basepath)
            path = os.path.dirname(path)
        npaths_new = len(basepaths)
        dpaths = npaths_new - npaths_old
    if path:
        basepaths.append(path)

    basepaths.reverse()
    return basepaths


def main():
    include_name = r'C:\NASA\formats\pynastran_v0.6\pyNastran\bdf\writePath.py'
    print('%s' % write_include(include_name, is_windows=True))

    include_name = r'/opt/NASA/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print('%s' % write_include(include_name, is_windows=False))

    include_name = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print('%s' % write_include(include_name, is_windows=False))

    include_name = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print('%s' % write_include(include_name, is_windows=True))

if __name__ == '__main__':  ## pragma: no cover
    main()
