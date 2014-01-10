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
    py2_gui_scripts = ['pyNastranGUI = pyNastran.gui.gui:main',
                       'pyNastranGUI_wx = pyNastran.gui.gui_wx:main',
                       'pyNastranGUI_qt = pyNastran.gui.gui_qt:main',]

# set up all icons
icon_path = os.path.join('pyNastran', 'gui', 'icons')
icon_files = os.listdir(icon_path)
icon_files2 = []
for icon_file in icon_files:
    if icon_file.endswith('.png'):
        icon_files2.append(os.path.join(icon_path, icon_file))

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
          'numpy >= 1.6.2',
          'scipy >= 0.6.0',
#          'matplotlib >= 1.1.0',
          'vtk >= 5.8.0',
          'pillow >= 1.1.7',
#          'wx >= 2.8.12.0',
#          'cython',
      ],#{'': ['license.txt']}
      #package_data={'': ['*.png']},
      data_files=[tuple([icon_path] + icon_files2)],

      entry_points = {
        'console_scripts': [
                            'run_nastran_double_precision = pyNastran.bdf.test.run_nastran_double_precision:cmd_line',
                            'test_bdf  = pyNastran.bdf.test.test_bdf:main',
                            'test_op2  = pyNastran.op2.test.test_op2:main',
                            'test_f06  = pyNastran.f06.test.test_f06:main',
                            'pyNastran  = pyNastran.applications.pyNastranSolver:main',

                            'test_bdf2 = pyNastran.bdf2.test.test_bdf:main',
                            'pyNastran2 = pyNastran.bdf2.solver.solver:main',
                           #'nastranToCodeAster = pyNastran.converters.toCodeAster:main',
                           ] + py2_gui_scripts
      },
      test_suite = 'pyNastran.all_tests',
      )
