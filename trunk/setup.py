#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

import pyNastran
packages = find_packages()+['gui/icons/*.*']
#print "packages = ",packages
#sys.exit()

extra = {}
if sys.version_info >= (3,):
    extra['use_2to3'] = True
    #extra['convert_2to3_doctests'] = ['src/your/module/README.txt']  # what does this do?
    #extra['use_2to3_fixers'] = ['your.fixers']

setup(name='pyNastran',
      version=pyNastran.__version__,
      description=pyNastran.__desc__,
      long_description="""\
""",
      classifiers=[
                   'Natural Language :: English',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 3'], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
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
          'numpy >= 1.3.0',
          'scipy >= 0.6.0',
#          'matplotlib >= 1.1.0',
#          'vtk >= 5.8.0',
#          'wx >= 2.8.12.0',
#          'cython',
      ],
      entry_points = {
        'console_scripts': [
                            'test_bdf = pyNastran.bdf.test.test_bdf:main',
                            'test_op2 = pyNastran.op2.test.test_op2:main',
                            'test_f06 = pyNastran.f06.test.test_f06:main',
                            'pyNastranGUI = pyNastran.gui.gui:main',
                           #'nastranToCodeAster = pyNastran.converters.toCodeAster:main',
                           ]
      },
      test_suite = 'pyNastran.all_tests',
      )
