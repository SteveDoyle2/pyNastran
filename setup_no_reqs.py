#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

PY2 = False
PY3 = True
if sys.version_info < (3, 0):
    PY2 = True
    PY3 = False

imajor, minor1, minor2 = sys.version_info[:3]
if sys.version_info < (2, 7, 7):  # 2.7.15 used
    # makes sure we don't get the following bug:
    #   Issue #19099: The struct module now supports Unicode format strings.
    sys.exit('Upgrade your Python to >= 2.7.7 or 3.5+; version=(%s.%s.%s)' % (imajor, minor1, minor2))

if PY3:
    if sys.version_info < (3, 5, 0):  # 3.7.1 used
        sys.exit('Upgrade your Python to >= 2.7.7 or 3.5+; version=(%s.%s.%s)' % (imajor, minor1, minor2))
import pyNastran
packages = find_packages()+['gui/icons/*.*']

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

            'format_converter = pyNastran.converters.type_converter:main',
            'pyNastranGUI = pyNastran.gui.gui:cmd_line',
            'bdf = pyNastran.bdf.mesh_utils.utils:cmd_line',
            'f06 = pyNastran.f06.utils:cmd_line',

            #'pyNastranv = pyNastran.dev.bdf_vectorized.solver.solver:main',
            #'test_bdfv = pyNastran.dev.bdf_vectorized.test.test_bdf_vectorized2:main',
            'test_bdfv = pyNastran.dev.bdf_vectorized2.test.test_bdf:main',
            #'nastran_to_code_aster = pyNastran.converters.dev.code_aster.nastran_to_code_aster:main',
        ]
    },
    test_suite='pyNastran.all_tests',
)

