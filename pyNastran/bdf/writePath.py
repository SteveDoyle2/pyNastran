import os

def split_path(abspath):
    """
    Takes a path and splits it into the various components
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

def write_include(filename, is_windows=True):
    """
    Writes a bdf INCLUDE file line given an imported filename.

    :param filename: the filename to write
    :param is_windows: Windows has a special format for writing INCLUDE files
                       so the format for a BDF that will run on Linux and
                       Windows is different.

    Example (Linux)
    ---------------
    fname = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    write_include(fname, is_windows=False)

    returns...
    INCLUDE /opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/
            pyNastran/bdf/writePath.py
    """
    msg = 'INCLUDE '  # len=8
    #nmax = 72 - 8 # 64

    if is_windows:
        marker = '\\'
    else:
        marker = '/'

    sline = split_path(filename)
    print "sline =", sline
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


if __name__ == '__main__':  ## pragma: no cover
    fname = r'C:\NASA\formats\pynastran_v0.6\pyNastran\bdf\writePath.py'
    print '%s' % write_include(fname, is_windows=True)

    fname = r'/opt/NASA/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print '%s' % write_include(fname, is_windows=False)

    fname = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print '%s' % write_include(fname, is_windows=False)

    fname = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print '%s' % write_include(fname, is_windows=True)