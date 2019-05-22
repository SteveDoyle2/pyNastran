"""helper for setup.py"""
import os
import sys

# features in packages used by pyNastran
# numpy
#  - 1.12 min for 3.6
#  - 1.13: adds axis support to unique
#  - 1.14: adds encoding support to savetxt (unused)
#  - 1.14: adds proper writing of np.savetxt for open file objects
#          (used for unicode savetxt using with statement) in Python 3.6
#  - 1.15: min for Python 3.7? I guess 1.14 is fine for a requirement...
# scipy:
#  - 0.19:   min for Python 3.6
#  - 0.19:   min for Python 3.7?; last 0.x release
# matplotlib:
#  - 2.0: adds C0/C1 colors (use colors from default colormap);
#    min for Python 3.6
#  - 2.1: adds plt.subplots support (untested?)
#  - 2.2: min for Python 3.7

# the packages that change requirements based on python version
REQS = {
    '3.6' : {
        # package : (less than check, required)
        'numpy' : ('1.14', '>=1.14'),
        'scipy' : ('0.19', '>=0.19'),
        'matplotlib' : ('2.0', '>=2.0'),
    },
    '3.7' : {
        'numpy' : ('1.14', '>=1.14'),
        'scipy' : ('1.0', '>=1.0'),
        'matplotlib' : ('2.2', '>=2.2'),  # 2.2.4 adds Python 3.7 support
    },
}

def check_python_version():
    """verifies the python version"""
    imajor, minor1, minor2 = sys.version_info[:3]
    if sys.version_info < (3, 6, 0):  # 3.7.1 used
        sys.exit('Upgrade your Python to 3.6+; version=(%s.%s.%s)' % (
            imajor, minor1, minor2))


def int_version(name, version):
    """splits the version into a tuple of integers"""
    sversion = version.split('-')[0]
    #numpy
    #scipy
    #matplotlib
    #qtpy
    #vtk
    #cpylog
    #pyNastran
    if 'rc' not in name:
        # it's gotta be something...
        # matplotlib3.1rc1
        sversion = sversion.split('rc')[0]

    try:
        return [int(val) for val in sversion.split('.')]
    except ValueError:
        raise SyntaxError('cannot determine version for %s %s' % (name, sversion))


def str_version(version):
    """converts a tuple of intergers to a version number"""
    return '.'.join(str(versioni) for versioni in version)


def get_package_requirements(is_gui=True, python_version=None):
    """gets the requirements for setup.py"""
    if python_version is None:
        python_version = '%s.%s' % sys.version_info[:2]

    if python_version not in REQS:
        python_version = '3.7'
    vreqs = REQS[python_version]

    all_reqs = {}
    packages = []


    is_dev = (
        'TRAVIS' in os.environ or
        'APPVEYOR' in os.environ or
        'READTHEDOCS' in os.environ
    )
    is_travis = 'TRAVIS' in os.environ
    is_rtd = 'READTHEDOCS' in os.environ

    #if is_dev or is_gui:
        #try:
            #import vtk
            #vtk_version = '.'.join(vtk.VTK_VERSION.split('.'))
            #all_reqs['vtk'] = vtk_version
            #if vtk_version < '7.0.0':
                #print("vtk.VTK_VERSION = %r < '7.0.0'" % vtk.VTK_VERSION)
                #py_packages.append('vtk >= 7.0.0')
        #except ImportError:
            #py_packages.append('vtk >= 7.0.0')  # 8.x used

    py_packages = []

    if is_rtd:
        py_packages.append('numpy')
    else:
        version_check, required_version = vreqs['scipy']
        try:
            import numpy as np
            sver = np.lib.NumpyVersion(np.__version__)
            iver = int_version('numpy', sver.version)
            all_reqs['numpy'] = sver.version

            iversion_check = int_version('numpy', version_check)
            #srequired_version = int_version('numpy', required_version)
            #print('numpy %r %r' % (sver, iversion_check))
            if iver < iversion_check:
                print("numpy.__version__ = %r < %s" % (np.__version__, version_check))
                py_packages.append('numpy %s' % required_version)
                all_reqs['numpy'] = version_check
                py_packages.append('numpy %s' % required_version) # was 1.11; ,<1.13.0
        except ImportError:
            all_reqs['numpy'] = required_version
            py_packages.append('numpy %s' % required_version) # ,<1.13.0; 1.15.1 used

    if is_rtd:
        py_packages.append('scipy')
    else:
        version_check, required_version = vreqs['scipy']
        try:
            import scipy
            sver = scipy.version.short_version
            iver = int_version('scipy', sver)
            all_reqs['scipy'] = sver

            iversion_check = int_version('scipy', version_check)
            #srequired_version = int_version('scipy', required_version)
            print('scipy %r %r' % (sver, iversion_check))
            if iver < iversion_check:
                print("scipy.version.short_version = %r < %r" % (
                    scipy.version.short_version, version_check))
                all_reqs['scipy'] = required_version
                py_packages.append('scipy %s' % required_version)
        except ImportError:
            py_packages.append('scipy %s' % required_version)  # 1.1.0 used

    try:
        import six
        iver = int_version('six', six.__version__)
        all_reqs['six'] = str_version(iver)
        if iver < [1, 9, 0]:
            print("six.__version__ = %r < '1.9.0'" % six.__version__)
            all_reqs['six'] = '>= 1.9.0'
            py_packages.append('six >= 1.9.0')
    except ImportError:
        py_packages.append('six >= 1.9.0')  # 1.12.0 used

    if is_gui:
        version_check, required_version = vreqs['matplotlib']
        try:
            import matplotlib
            iver = int_version('matplotlib', matplotlib.__version__)
            all_reqs['matplotlib'] = str_version(iver)
            iversion_check = int_version('matplotlib', version_check)
            if iver < iversion_check:
                print("matplotlib.__version__ = %r < %r" % (matplotlib.__version__, version_check))
                    #matplotlib.__version__, str_version(iversion_check)))
                all_reqs['matplotlib'] = required_version # '>= 2.1.2'
                py_packages.append('matplotlib %s' % required_version)
        except ImportError:
            py_packages.append('matplotlib %s' % required_version)  # was 2.1.2; 2.2.3 used


    try:
        import cpylog
        iver = int_version('cpylog', cpylog.__version__)
        all_reqs['cpylog'] = str_version(iver)
        if iver <= [1, 0, 2]:
            print("cpylog.__version__ = %r != '1.0.2'" % cpylog.__version__)
            all_reqs['cpylog'] = '>= 1.0.2'
            py_packages.append('cpylog >= 1.0.2')
    except ImportError:
        py_packages.append('cpylog >= 1.0.2')  # 1.0.2 used


    try:
        import docopt
        iver = int_version('docopt', docopt.__version__)
        all_reqs['docopt'] = str_version(iver)
        if iver != [0, 6, 2]:
        #if docopt.__version__ != '0.6.2':
            print("docopt.__version__ = %r != '0.6.2'" % docopt.__version__)
            all_reqs['docopt'] = '== 0.6.2'
            py_packages.append('docopt == 0.6.2')
    except ImportError:
        py_packages.append('docopt == 0.6.2')  # 0.6.2 used

    if is_rtd:
        pass
    elif is_gui:
        try:
            import qtpy
            iver = int_version('qtpy', qtpy.__version__)
            all_reqs['qtpy'] = str_version(iver)
            if iver < [1, 4, 0]:
                print("qtpy.__version__ = %r < '1.4.0'" % qtpy.__version__)
                all_reqs['qtpy'] = '>= 1.4.0'
                py_packages.append('qtpy >= 1.4.0')
        except ImportError:
            py_packages.append('qtpy >= 1.4.0')  # 1.5.0 used

    if is_rtd:
        pass
    elif is_gui:
        try:
            import imageio
            if imageio.__version__ < '2.2.0':
                #print("imageio.version = %r < '2.2.0'" % imageio.__version__)
                all_reqs['imageio'] = '>= 2.2.0'
                py_packages.append('imageio >= 2.2.0')
            else:
                all_reqs['imageio'] = imageio.__version__
        except ImportError:
            py_packages.append('imageio >= 2.2.0')


    is_windows = 'nt' in os.name
    if is_travis and not is_windows:
        py_packages.append('python-coveralls')
        #py_packages.append('codecov')
        #py_packages.append('coverage')

    install_requires = py_packages
    #print(all_reqs)
    print('install_requires =', install_requires)
    return all_reqs, install_requires
