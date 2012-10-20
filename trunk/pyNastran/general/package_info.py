import os

def pkginfo(name):
    module = __import__(name)
    path = module.__path__[0]
    #print os.listdir(path)
    #print "path = ",path
    try:
        version = module.__version__
    except:
        if name == 'PyQt4':
            import PyQt4.QtCore
            version = PyQt4.QtCore.PYQT_VERSION_STR
        elif name == 'h5py':
            import h5py.version
            version = h5py.version.version
        else:
            version = None
    print "name=%s path=%s version=%s" %(name,path,version)
    return version  

if __name__ == '__main__':
    #pkginfo('numpy')
    #pkginfo('pyNastran')
    #pkginfo('scipy')
    #pkginfo('matplotlib')
    pkginfo('sphinx')
    pkginfo('vtk')  # doesn't work
    #pkginfo('h5py')
    #pkginfo('sympy')
    #pkginfo('setuptools')
    #pkginfo('PyQt4')
    #pkginfo('wx')

