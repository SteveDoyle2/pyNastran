"""helper for setup.py"""
import os
import sys

PY2 = False
PY3 = True
if sys.version_info < (3, 0):
    PY2 = True
    PY3 = False

# features in packages used by pyNastran
# numpy 1.14
#  - adds encoding support to savetxt
# scipy 0.18.1:
#  - fixed kdtree used by nodal equivalencing
# matplotlib 2.1:
#  - adds plt.subplots support used by shear/moment/torque plotter

# the packages that change requirements based on python version
REQS = {
    '2.7' : {
        # package : (less than check, required)
        'numpy' : ('1.11', '>=1.12,<1.17'), #  1.12???
        'scipy' : ('0.18.1', '>=0.18.1,<1.13'),
        'matplotlib' : ('2.1', '>=2.1,<3'),
    },
    '3.6' : {
        'numpy' : ('1.12', '>=1.12'),  # 1.12 adds Python 3.6 support
        'scipy' : ('0.19', '>=0.19'),  # 0.19 adds Python 3.6 support
        'matplotlib' : ('2.1', '>=2.1'),  # guess
    },
    '3.7' : {
        'numpy' : ('1.12', '>=1.15'),
        'scipy' : ('0.19', '>=0.19'), #  1.12???
        'matplotlib' : ('2.2', '>=2.2'),  # 2.2.4 adds Python 3.7 support
    },
}

def check_python_version():
    """verifies the python version"""
    imajor, minor1, minor2 = sys.version_info[:3]
    if sys.version_info < (2, 7, 7):  # 2.7.15 used
        # makes sure we don't get the following bug:
        #   Issue #19099: The struct module now supports Unicode format strings.
        sys.exit('Upgrade your Python to >= 2.7.7 or 3.6+; version=(%s.%s.%s)' % (
            imajor, minor1, minor2))

    if PY3:
        if sys.version_info < (3, 6, 0):  # 3.7.1 used
            sys.exit('Upgrade your Python to >= 2.7.7 or 3.6+; version=(%s.%s.%s)' % (
                imajor, minor1, minor2))


def int_version(version):
    """splits the version into a tuple of integers"""
    sversion = version.split('-')[0]
    return [int(val) for val in sversion.split('.')]


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
    py2_packages = []
    py3_packages = []


    is_dev = (
        'TRAVIS' in os.environ or
        'APPVEYOR' in os.environ or
        'READTHEDOCS' in os.environ
    )
    is_travis = 'TRAVIS' in os.environ
    is_rtd = 'READTHEDOCS' in os.environ

    if is_dev or is_gui:
        try:
            import vtk
            vtk_version = '.'.join(vtk.VTK_VERSION.split('.'))
            all_reqs['vtk'] = vtk_version
            if vtk_version < '7.0.0':
                print("vtk.VTK_VERSION = %r < '7.0.0'" % vtk.VTK_VERSION)
                py2_packages.append('vtk >= 7.0.0')
        except ImportError:
            py2_packages.append('vtk >= 7.0.0')  # 8.x used

        py2_packages += [
            ##'dill'
        ]

        py3_packages += [
            #'pillow >= 2.7.0',
            ##'dill'
        ]

    py_packages = []

    if is_rtd:
        py_packages.append('numpy')
    else:
        try:
            import numpy as np
            sver = np.lib.NumpyVersion(np.__version__)
            iver = int_version(sver.version)
            all_reqs['numpy'] = sver.version

            version_check, required_version = vreqs['scipy']
            iversion_check = int_version(version_check)
            #srequired_version = int_version(required_version)
            #print('numpy %r %r' % (sver, iversion_check))
            if iver < iversion_check:
                print("numpy.__version__ = %r < %s" % (np.__version__, version_check))
                py_packages.append('numpy %s' % version_check)
                all_reqs['numpy'] = version_check
                py_packages.append('numpy %s' % version_check) # was 1.11; ,<1.13.0
        except ImportError:
            all_reqs['numpy'] = version_check
            py_packages.append('numpy %s' % version_check) # ,<1.13.0; 1.15.1 used

    if is_rtd:
        py_packages.append('scipy')
    else:
        try:
            import scipy
            sver = scipy.version.short_version
            iver = int_version(sver)
            # numpy:
            #  - 1.12 min for 3.6
            #  - 1.15 min for 3.7
            # scipy:
            #  - 0.19 min for 3.6
            #  - ???  min for 3.7
            all_reqs['scipy'] = sver

            version_check, required_version = vreqs['scipy']
            iversion_check = int_version(version_check)
            #srequired_version = int_version(required_version)
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
        iver = int_version(six.__version__)
        all_reqs['six'] = str_version(iver)
        if iver < [1, 9, 0]:
            print("six.__version__ = %r < '1.9.0'" % six.__version__)
            all_reqs['six'] = '>= 1.9.0'
            py_packages.append('six >= 1.9.0')
    except ImportError:
        py_packages.append('six >= 1.9.0')  # 1.12.0 used

    if is_gui:
        try:
            import matplotlib
            iver = int_version(matplotlib.__version__)
            all_reqs['matplotlib'] = str_version(sver)
            version_check, required_version = vreqs['matplotlib']
            iversion_check = int_version(version_check)
            if iver < iversion_check:
                print("matplotlib.__version__ = %r < %r" % (matplotlib.__version__, version_check))
                    #matplotlib.__version__, str_version(iversion_check)))
                all_reqs['matplotlib'] = required_version # '>= 2.1.2'
                py_packages.append('matplotlib %s' % required_version)
        except ImportError:
            py_packages.append('matplotlib %s' % required_version)  # was 2.1.2; 2.2.3 used


    try:
        import cpylog
        sver = int_version(cpylog.__version__)
        all_reqs['cpylog'] = str_version(sver)
        if sver <= [1, 0, 2]:
            print("cpylog.__version__ = %r != '1.0.2'" % cpylog.__version__)
            all_reqs['cpylog'] = '>= 1.0.2'
            py_packages.append('cpylog >= 1.0.2')
    except ImportError:
        py_packages.append('cpylog >= 1.0.2')  # 1.0.2 used


    try:
        import docopt
        iver = int_version(docopt.__version__)
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
            iver = int_version(qtpy.__version__)
            all_reqs['qtpy'] = str_version(iver)
            if iver < [1, 4, 0]:
                print("qtpy.__version__ = %r < '1.4.0'" % qtpy.__version__)
                all_reqs['qtpy'] = '>= 1.4.0'
                py_packages.append('qtpy >= 1.4.0')
        except ImportError:
            py_packages.append('qtpy >= 1.4.0')  # 1.5.0 used

    if PY2:
        try:
            import typing
        except ImportError:
            # PY2
            all_reqs['typing'] = '>= 3.6.4'
            py_packages.append('typing >= 3.6.4')  # 3.6.6 used

        try:
            import pathlib2
        except ImportError:
            all_reqs['pathlib2'] = '>= 2.3.0'
            py_packages.append('pathlib2 >= 2.3.0')  # 2.3.2 used

        try:
            import scandir
            iver = int_version(scandir.__version__)
            all_reqs['scandir'] = str_version(iver)
            if sver < [1, 7, 0]:
                print("scandir.__version__ = %r < '1.7.0'" % scandir.__version__)
                all_reqs['scandir'] = '>= 1.7.0'
                py_packages.append('scandir >= 1.7.0')
        except ImportError:
            py_packages.append('scandir >= 1.7.0')  # 1.9.0 used

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

    #py_packages = [
    #    'numpy >= 1.9.2',
    #    'scipy >= 0.16.0, scipy < 0.18.0',
    #]

    is_windows = 'nt' in os.name
    if is_travis and not is_windows:
        py_packages.append('python-coveralls')
        #py_packages.append('codecov')
        #py_packages.append('coverage')

    install_requires = py_packages + [
        # -*- Extra requirements: -*-
        ##'cython',
    ] + py2_packages + py3_packages,
    #print(all_reqs)
    #print(install_requires)
    return all_reqs, install_requires
