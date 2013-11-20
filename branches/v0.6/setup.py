#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

import pyNastran
packages = find_packages()+['gui/icons/*.*']
#print "packages = ",packages
#sys.exit()

#extra = {}
#if sys.version_info >= (3,):
    #extra['use_2to3'] = True
    #extra['convert_2to3_doctests'] = ['src/your/module/README.txt']  # what does this do?
    #extra['use_2to3_fixers'] = ['your.fixers']

py2_scripts = []
if sys.version_info <= (3,):
    py2_scripts = ['pyNastranGUI = pyNastran.gui.gui:main',
                   'pyNastranGUI_wx = pyNastran.gui.gui_wx:main',
                   'pyNastranGUI_qt = pyNastran.gui.gui_qt:main',]

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
          'docopt == 0.6.1',
          'numpy >= 1.8.0',
          'scipy >= 0.6.0',
#          'matplotlib >= 1.1.0',
          'vtk >= 5.8.0',
          'pillow >= 1.1.7',
#          'wx >= 2.8.12.0',
#          'cython',
      ],#{'': ['license.txt']}
      #package_data={'': ['*.png']},
      data_files=[('pyNastran/gui/icons',
        [
            'pyNastran/gui/icons/topen.png',
            'pyNastran/gui/icons/tbdf.png',
            'pyNastran/gui/icons/top2.png',
            #'pyNastran/gui/icons/topen.png',
            'pyNastran/gui/icons/twireframe.png',
            'pyNastran/gui/icons/tsolid.png',
            #'pyNastran/gui/icons/tsave.png',
            #'pyNastran/gui/icons/tsaveas.png',
            'pyNastran/gui/icons/texit.png',
            'pyNastran/gui/icons/tcamera.png',
            #'pitz/pitztypes/tracpitz.py.sample',
        ])],

      entry_points = {
        'console_scripts': [
                            'test_bdf = pyNastran.bdf.test.test_bdf:main',
                            'test_op2 = pyNastran.op2.test.test_op2:main',
                            'test_f06 = pyNastran.f06.test.test_f06:main',
                           #'nastranToCodeAster = pyNastran.converters.toCodeAster:main',
                           ] + py2_scripts
      },
      test_suite = 'pyNastran.all_tests',
      )
