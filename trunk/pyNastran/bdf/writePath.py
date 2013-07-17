import os

def split_path(abspath):
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

    #print "basepaths = ",basepaths
    #print "path = ",path
    #pths = []
    basepaths.reverse()
    #print "split_path = ",basepaths
    return basepaths

def write_include(filename, is_windows=True):
    msg = 'INCLUDE '  # len=8
    nmax = 72 - 8
    
    if is_windows:
        marker = '\\'
    else:
        marker = '/'
    
    print "abspath = ", filename, len(filename)
    sline = split_path(filename)

    #print "**", sline
    if len(filename) > 52:#  62
        pth = ''
        for p in sline:
            if p == '/':
                pth += '%s' % marker
            else:
                pth += '%s%s' % (p, marker)
        pth = pth.rstrip(marker)
        #pth = marker.join(sline)
    else:
        #asdf
        #print abspath
        pth = marker.join(sline)

    return pth

if __name__ == '__main__':  ## pragma: no cover
    fname = r'C:\NASA\m4\formats\pynastran_v0.6\pyNastran\bdf\writePath.py'
    print '%r' % write_include(fname, is_windows=True)

    fname = r'/opt/NASA/m4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
    print '%r' % write_include(fname, is_windows=False)