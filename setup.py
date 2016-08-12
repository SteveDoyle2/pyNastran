#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

PY2 = False
if sys.version_info < (3, 0):
    PY2 = True
if sys.version_info < (2, 7, 7):
    imajor, minor1, minor2 = sys.version_info[:3]
    # makes sure we don't get the following bug:
    #   Issue #19099: The struct module now supports Unicode format strings.
    sys.exit('Upgrade your Python to >= 2.7.7; version=(%s.%s.%s)' % (imajor, minor1, minor2))

import pyNastran
packages = find_packages()+['gui/icons/*.*']
#print "packages = ",packages
#sys.exit()

#extra = {}
#if sys.version_info >= (3,):
    #extra['use_2to3'] = True
    #extra['convert_2to3_doctests'] = ['src/your/module/README.txt']  # what does this do?
    #extra['use_2to3_fixers'] = ['your.fixers']

py2_gui_scripts = []
py2_packages = []
py3_gui_scripts = []
py3_packages = []

if sys.version_info <= (3,):
    py2_gui_scripts = ['pyNastranGUI = pyNastran.gui.gui:main',]

    try:
        import vtk
        vtk_version = '.'.join(vtk.VTK_VERSION.split('.'))
        if vtk_version < '5.10.1':
            print("vtk.VTK_VERSION = %r < '5.10.1'" % vtk.VTK_VERSION)
            py2_packages.append('vtk >= 5.10.1')
    except ImportError:
        py2_packages.append('vtk >= 5.10.1')

    #try:  # is this still required?
        #import PIL
        #if PIL.VERSION < '2.7.0':
            #print("PIL.version = %r < '2.7.0'" % PIL.VERSION)
            #py2_packages.append('pillow >= 2.7.0')
    #except ImportError:
        #py2_packages.append('pillow >= 2.7.0')

    py2_packages += [
        #'vtk >= 5.10.1',
        #'pillow >= 2.7.0',
        ##'dill'
        ##'wx >= 2.8.12.0',
    ]
else:
    py3_gui_scripts = ['pyNastranGUI = pyNastran.gui.gui:main',]

#    try:
#        import vtk
#        vtk_version = '.'.join(vtk.VTK_VERSION.split('.'))
#        if vtk_version < '7.0.0':
#            print("vtk.VTK_VERSION = %r < '7.0.0'" % vtk.VTK_VERSION)
#            py3_packages.append('vtk >= 7.0.0')
#    except ImportError:
#        py3_packages.append('vtk >= 7.0.0')

    #try:  # is this still required?
        #import PIL
        #if PIL.VERSION < '2.7.0':
            #print("PIL.version = %r < '2.7.0'" % PIL.VERSION)
            #py3_packages.append('pillow >= 2.7.0')
    #except ImportError:
        #py3_packages.append('pillow >= 2.7.0')

    py3_packages += [
        #'vtk >= 5.10.0',
        #'pillow >= 2.7.0',
        ##'dill'
        ##'wx >= 2.8.12.0',
    ]

py_packages = []

try:
    import numpy as np
    ver = np.lib.NumpyVersion(np.__version__)
    if ver < '1.9.2':
        print("np.__version__ = %r < '1.9.2'" % np.__version__)
        py_packages.append('numpy >= 1.9.2')
except ImportError:
    py_packages.append('numpy >= 1.9.2')

try:
    import scipy
    ver = scipy.version.short_version
    if ver < '0.16.0':
        print("scipy.version.short_version = %r < '0.16.0'" % scipy.version.short_version)
        py_packages.append('scipy >= 0.16.0')
except ImportError:
    py_packages.append('scipy >= 0.16.0')

try:
    import six
    sver = [int(val) for val in six.__version__.split('-')[0].split('.')]
    if sver < [1, 9, 0]:
        print("six.__version__ = %r < '1.9.0'" % six.__version__)
        py_packages.append('six >= 1.9.0')
except ImportError:
    py_packages.append('six >= 1.9.0')

try:
    import docopt
    sver = [int(val) for val in docopt.__version__.split('-')[0].split('.')]
    if sver != [0, 6, 2]:
    #if docopt.__version__ != '0.6.2':
        print("docopt.__version__ = %r != '0.6.2'" % docopt.__version__)
        py_packages.append('docopt == 0.6.2')
except ImportError:
    py_packages.append('docopt == 0.6.2')


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
] + py2_packages + py3_packages,


# set up all icons
icon_path = os.path.join('pyNastran', 'gui', 'icons')
icon_files = os.listdir(icon_path)
icon_files2 = []
for icon_file in icon_files:
    if icon_file.endswith('.png'):
        icon_files2.append(os.path.join(icon_path, icon_file))

setup(
    name='pyNastran',
    version=pyNastran.__version__,
    description=pyNastran.__desc__,
    long_description="""\
""",
    classifiers=[
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        ], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords='',
    author=pyNastran.__author__,
    author_email=pyNastran.__email__,
    url=pyNastran.__website__,
    license=pyNastran.__license__,
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
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
            'test_f06  = pyNastran.f06.test.test_f06:main',

            #'test_bdf2 = pyNastran.bdf.dev_vectorized.test.test_bdf:main',
            #'pyNastran2 = pyNastran.bdf.dev_vectorized.solver.solver:main',
            #'nastranToCodeAster = pyNastran.converters.toCodeAster:main',
            'format_converter = pyNastran.converters.type_converter:main',
            'bdf = pyNastran.bdf.mesh_utils.utils:cmd_line',
        ] + py2_gui_scripts + py3_gui_scripts
    },
    test_suite='pyNastran.all_tests',
)

