#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

PY2 = False
if sys.version_info < (3, 0):
    PY2 = True
if sys.version_info < (2, 7, 7):  # 2.7.13 used
    imajor, minor1, minor2 = sys.version_info[:3]
    # makes sure we don't get the following bug:
    #   Issue #19099: The struct module now supports Unicode format strings.
    sys.exit('Upgrade your Python to >= 2.7.7; version=(%s.%s.%s)' % (imajor, minor1, minor2))

import pyNastran
packages = find_packages()+['gui/icons/*.*']
#print "packages = ",packages
#sys.exit()

py2_packages = []
py_packages = []

try:
    import numpy as np
    ver = np.lib.NumpyVersion(np.__version__)
    if ver < '1.11.0':
        print("np.__version__ = %r < '1.11.0'" % np.__version__)
        py_packages.append('numpy >= 1.11.0,<1.13.0')
except ImportError:
    py_packages.append('numpy >= 1.11.0,<1.13.0')

try:
    import scipy
    ver = scipy.version.short_version
    if ver < '0.18.1':
        print("scipy.version.short_version = %r < '0.18.1'" % scipy.version.short_version)
        py_packages.append('scipy >= 0.17.0')
except ImportError:
    py_packages.append('scipy >= 0.18.1')  # 0.18.1 used

try:
    import six
    sver = [int(val) for val in six.__version__.split('-')[0].split('.')]
    if sver < [1, 10, 0]:
        print("six.__version__ = %r < '1.10.0'" % six.__version__)
        py_packages.append('six >= 1.10.0')
except ImportError:
    py_packages.append('six >= 1.10.0')  # 1.10.0 used


#try:
#    import matplotlib
#    sver = [int(val) for val in matplotlib.__version__.split('-')[0].split('.')]
#    if sver < [1, 5, 1]:
#        print("matplotlib.__version__ = %r < '1.5.1'" % six.__version__)
#        py_packages.append('matplotlib >= 1.5.1')
#except ImportError:
#    py_packages.append('matplotlib >= 1.5.1')


try:
    import docopt
    sver = [int(val) for val in docopt.__version__.split('-')[0].split('.')]
    if sver != [0, 6, 2]:
    #if docopt.__version__ != '0.6.2':
        print("docopt.__version__ = %r != '0.6.2'" % docopt.__version__)
        py_packages.append('docopt == 0.6.2')
except ImportError:
    py_packages.append('docopt == 0.6.2')  # 0.6.2 used


try:
    import typing
except ImportError:
    py_packages.append('typing >= 3.6.1')  # 3.6.1 used


if PY2:
    try:
        import pathlib2
    except ImportError:
        py_packages.append('pathlib2 >= 2.2.0')  # 2.2.0 used

    try:
        import scandir
        sver = [int(val) for val in scandir.__version__.split('-')[0].split('.')]
        if sver < [1, 4, 0]:
            print("scandir.__version__ = %r < '1.4.0'" % scandir.__version__)
            py_packages.append('scandir >= 1.4.0')
    except ImportError:
        py_packages.append('scandir >= 1.4.0')  # 1.4.0 used

#py_packages = [
#    'numpy >= 1.9.2',
#    'scipy >= 0.16.0, scipy < 0.18.0',
#]

install_requires = py_packages + [
    # -*- Extra requirements: -*-
    #'docopt == 0.6.2',
    ##'matplotlib >= 1.3.0',
    #'six >= 1.9.0',
    ##'cython',
] + py2_packages,


# set up all icons
icon_path = os.path.join('pyNastran', 'gui', 'icons')
icon_files = os.listdir(icon_path)
icon_files2 = []
for icon_file in icon_files:
    if icon_file.endswith('.png'):
        icon_files2.append(os.path.join(icon_path, icon_file))

exclude_words = [
    'pyNastran.dev.bdf_vectorized', 'pyNastran.dev.bdf_vectorized.cards',
    'pyNastran.f06.dev',
    'pyNastran.op2.dev', 'pyNastran.op2.dev.original',
    'pyNastran.converters.dev', 'pyNastran.xdb',]
packages = find_packages(exclude=['ez_setup', 'examples', 'tests'] + exclude_words)
for exclude_word in exclude_words:
    packages = [package for package in packages if exclude_word not in package]
#print(packages, len(packages)) # 83

setup(
    name='pyNastran',
    version=pyNastran.__version__,
    description=pyNastran.__desc__,
    long_description="""\
""",
    classifiers=[
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        ], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords='',
    python_requires='>2.7.6, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
    author=pyNastran.__author__,
    author_email=pyNastran.__email__,
    url=pyNastran.__website__,
    license=pyNastran.__license__,
    packages=packages,
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requires,
    #{'': ['license.txt']}
    #package_data={'': ['*.png']},
    #data_files=[(icon_path, icon_files2)],
    package_data={
        # https://pythonhosted.org/setuptools/setuptools.html#including-data-files
        # If any package contains *.png files, include them:
        '': ['*.png'],
        #'mypkg': ['data/*.dat'],
    },
    entry_points={
        'console_scripts': [
            'run_nastran_double_precision = pyNastran.bdf.test.run_nastran_double_precision:cmd_line',
            'test_bdf  = pyNastran.bdf.test.test_bdf:main',
            'test_op2  = pyNastran.op2.test.test_op2:main',
            'test_op4  = pyNastran.op4.test.test_op4:main',
            'test_abaqus = pyNastran.converters.abaqus.test_abaqus:main',
            'test_pynastrangui = pyNastran.gui.test.test_gui:main',
            #'test_f06  = pyNastran.f06.test.test_f06:main',

            'format_converter = pyNastran.converters.type_converter:main',
            'pyNastranGUI = pyNastran.gui.gui:cmd_line',
            'bdf = pyNastran.bdf.mesh_utils.utils:cmd_line',
            'f06 = pyNastran.f06.utils:cmd_line',

            #'pyNastranv = pyNastran.dev.bdf_vectorized.solver.solver:main',
            #'test_bdfv = pyNastran.dev.bdf_vectorized.test.test_bdf_vectorized2:main',
            #'nastran_to_code_aster = pyNastran.converters.dev.code_aster.nastran_to_code_aster:main',
        ]# + py2_gui_scripts
    },
    test_suite='pyNastran.all_tests',
)

