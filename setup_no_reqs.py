#!/usr/bin/env python
import os
from setuptools import setup, find_packages

import pyNastran
from packages import (check_python_version, get_package_requirements,
                      update_version_file, PYTHON_REQUIRES,
                      LONG_DESCRIPTION, CLASSIFIERS, EXCLUDE_WORDS,
                      CONSOLE_SCRIPTS)


check_python_version()
all_reqs, install_requires = get_package_requirements(is_gui=True)

packages = find_packages() + ['gui/icons/*.*']

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

#revision = get_git_revision_short_hash()
#__version__ = '1.3.0+%s' % revision
#__releaseDate__ = '2019/6/xx'
#__releaseDate2__ = 'JUNE xx, 2019'

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
        'console_scripts': CONSOLE_SCRIPTS,
    },
    test_suite='pyNastran.all_tests',
)
print()
for package in install_requires:
    print('did not install %s' % package)
