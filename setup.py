#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

PY2 = False
if sys.version_info < (3, 0):
    PY2 = True
if sys.version_info < (2, 7, 7):
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
if sys.version_info <= (3,):
    py2_gui_scripts = ['pyNastranGUI = pyNastran.gui.gui:main',]
    py2_packages += [
        'vtk >= 5.10.0',
        'pillow >= 2.7.0',
        #'wx >= 2.8.12.0',
    ]

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
        ], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords='',
    author=pyNastran.__author__,
    author_email=pyNastran.__email__,
    url=pyNastran.__website__,
    license=pyNastran.__license__,
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    include_package_data=True,
    zip_safe=False,
    install_requires=[

        # -*- Extra requirements: -*-
        'numpy >= 1.9.0',
        'scipy >= 0.15.0',
        'docopt == 0.6.2',
        #'matplotlib >= 1.3.0',
        #'vtk >= 5.10.0',   # should be installed, but not for the docs
        'pillow >= 2.7.0',
        'six >= 1.8.0',
        #'wx >= 2.8.12.0',
        #'cython',
    ],# + py2_packages,
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
        ]
    },
    test_suite='pyNastran.all_tests',
)

