#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

import pyNastran
from packages import (check_python_version, get_package_requirements,
                      update_version_file, PYTHON_REQUIRES,
                      LONG_DESCRIPTION, CLASSIFIERS, EXCLUDE_WORDS)

add_vtk_qt = True
bdist = False
if 'bdist_wheel' in sys.argv:
    add_vtk_qt = False
    bdist = True
    assert '\r' not in LONG_DESCRIPTION, LONG_DESCRIPTION

check_python_version()
install_requires = get_package_requirements(is_gui=True, add_vtk_qt=add_vtk_qt, bdist=bdist)

packages = find_packages() + ['gui/icons/*.*']
#print("packages = %s" % packages)

# set up all icons
icon_path = os.path.join('pyNastran', 'gui', 'icons')
icon_files = os.listdir(icon_path)
icon_files2 = []
for icon_file in icon_files:
    if icon_file.endswith('.png'):
        icon_files2.append(os.path.join(icon_path, icon_file))

packages = find_packages(exclude=['ez_setup', 'examples', 'tests'] + EXCLUDE_WORDS)
for exclude_word in EXCLUDE_WORDS:
    packages = [package for package in packages if exclude_word not in package]
#print(packages, len(packages)) # 83
update_version_file()

setup(
    name='pyNastran',
    version=pyNastran.__version__,
    description=pyNastran.__desc__,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/x-rst',
    classifiers=CLASSIFIERS,
    keywords='',
    python_requires=PYTHON_REQUIRES,
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
            #'run_nastran_double_precision = pyNastran.bdf.test.run_nastran_double_precision:cmd_line',
            'test_bdf  = pyNastran.bdf.test.test_bdf:main',
            'test_op2  = pyNastran.op2.test.test_op2:main',
            'test_op4  = pyNastran.op4.test.test_op4:main',
            #'test_abaqus = pyNastran.converters.abaqus.test_abaqus:main',
            'test_pynastrangui = pyNastran.gui.test.test_gui:main',

            'format_converter = pyNastran.converters.format_converter:cmd_line_format_converter',
            'pyNastranGUI = pyNastran.gui.gui:cmd_line',
            'bdf = pyNastran.bdf.mesh_utils.utils:cmd_line',
            'f06 = pyNastran.f06.utils:cmd_line',
        ]
    },
    test_suite='pyNastran.all_tests',
)
