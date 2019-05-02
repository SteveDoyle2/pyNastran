"""helper for setup.py"""
import os
import sys

PY2 = False
PY3 = True
if sys.version_info < (3, 0):
    PY2 = True
    PY3 = False

def check_python_version():
    """verifies the python version"""
    imajor, minor1, minor2 = sys.version_info[:3]
    if sys.version_info < (2, 7, 7):  # 2.7.15 used
        # makes sure we don't get the following bug:
        #   Issue #19099: The struct module now supports Unicode format strings.
        sys.exit('Upgrade your Python to >= 2.7.7 or 3.5+; version=(%s.%s.%s)' % (
            imajor, minor1, minor2))

    if PY3:
        if sys.version_info < (3, 6, 0):  # 3.7.1 used
            sys.exit('Upgrade your Python to >= 2.7.7 or 3.6+; version=(%s.%s.%s)' % (
                imajor, minor1, minor2))


def get_package_requirements(is_gui=True):
    """gets the requirements for setup.py"""
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
            ver = np.lib.NumpyVersion(np.__version__)
            all_reqs['numpy'] = ver.version
            if ver < '1.11.0':
                print("np.__version__ = %r < '1.11.0'" % np.__version__)
                py_packages.append('numpy >= 1.11.0')
                all_reqs['numpy'] = '>= 1.11.0'
                py_packages.append('numpy >= 1.11.0') # ,<1.13.0
        except ImportError:
            all_reqs['numpy'] = '>= 1.11.0'
            py_packages.append('numpy >= 1.11.0') # ,<1.13.0; 1.15.1 used

    if is_rtd:
        py_packages.append('scipy')
    else:
        try:
            import scipy
            ver = scipy.version.short_version
            all_reqs['scipy'] = ver
            if ver < '1.0.0':
                print("scipy.version.short_version = %r < '1.0.0'" % scipy.version.short_version)
                all_reqs['scipy'] = '>= 1.0.0'
                py_packages.append('scipy >= 1.0.0')
        except ImportError:
            py_packages.append('scipy >= 1.0.0')  # 1.1.0 used

    try:
        import six
        sver = [int(val) for val in six.__version__.split('-')[0].split('.')]
        all_reqs['six'] = '.'.join((str(val) for val in sver))
        if sver < [1, 11, 0]:
            print("six.__version__ = %r < '1.11.0'" % six.__version__)
            all_reqs['six'] = '>= 1.11.0'
            py_packages.append('six >= 1.11.0')
    except ImportError:
        py_packages.append('six >= 1.11.0')  # 1.12.0 used

    if is_gui:
        try:
            import matplotlib
            sver = [int(val) for val in matplotlib.__version__.split('-')[0].split('.')]
            all_reqs['matplotlib'] = '.'.join((str(val) for val in sver))
            if sver < [2, 1, 2]:
                print("matplotlib.__version__ = %r < '2.1.2'" % matplotlib.__version__)
                all_reqs['matplotlib'] = '>= 2.1.2'
                py_packages.append('matplotlib >= 2.1.2')
        except ImportError:
            py_packages.append('matplotlib >= 2.1.2')  # 2.2.3 used


    try:
        import cpylog
        sver = [int(val) for val in cpylog.__version__.split('-')[0].split('.')]
        all_reqs['cpylog'] = '.'.join((str(val) for val in sver))
        if sver <= [1, 0, 2]:
            print("cpylog.__version__ = %r != '1.0.2'" % cpylog.__version__)
            all_reqs['cpylog'] = '>= 1.0.2'
            py_packages.append('cpylog >= 1.0.2')
    except ImportError:
        py_packages.append('cpylog >= 1.0.2')  # 1.0.2 used


    try:
        import docopt
        sver = [int(val) for val in docopt.__version__.split('-')[0].split('.')]
        all_reqs['docopt'] = '.'.join((str(val) for val in sver))
        if sver != [0, 6, 2]:
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
            sver = [int(val) for val in qtpy.__version__.split('-')[0].split('.')]
            all_reqs['qtpy'] = '.'.join((str(val) for val in sver))
            if sver < [1, 4, 0]:
                print("qtpy.__version__ = %r < '1.4.0'" % qtpy.__version__)
                all_reqs['qtpy'] = '>= 1.4.0'
                py_packages.append('qtpy >= 1.4.0')
        except ImportError:
            py_packages.append('qtpy >= 1.4.0')  # 1.5.0 used

    try:
        import typing
    except ImportError:
        # PY2
        all_reqs['typing'] = '>= 3.6.4'
        py_packages.append('typing >= 3.6.4')  # 3.6.6 used


    if PY2:
        try:
            import pathlib2
        except ImportError:
            all_reqs['pathlib2'] = '>= 2.3.0'
            py_packages.append('pathlib2 >= 2.3.0')  # 2.3.2 used

        try:
            import scandir
            sver = [int(val) for val in scandir.__version__.split('-')[0].split('.')]
            all_reqs['scandir'] = '.'.join((str(val) for val in sver))
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
    return all_reqs, install_requires
