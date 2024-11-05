#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

import pyNastran
from packages import check_python_version, update_version_file

add_vtk_qt = True
bdist = False
if 'bdist_wheel' in sys.argv:
    add_vtk_qt = False
    bdist = True

check_python_version()
packages = find_packages() + ['gui/icons/*.*']
#print("packages = %s" % packages)

# set up all icons
icon_path = os.path.join('pyNastran', 'gui', 'icons')
icon_files = os.listdir(icon_path)
icon_files2 = []
for icon_file in icon_files:
    if icon_file.endswith('.png'):
        icon_files2.append(os.path.join(icon_path, icon_file))

# packages = find_packages(exclude=['ez_setup', 'examples', 'tests'] + EXCLUDE_WORDS)
# for exclude_word in EXCLUDE_WORDS:
#    packages = [package for package in packages if exclude_word not in package]
#print(packages, len(packages)) # 83
#update_version_file()

setup(
    name='pyNastran',
    version=pyNastran.__version__,
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
)
