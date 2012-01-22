import os
import sys
from setuptools import setup, find_packages

import pyNastran

setup(name='pyNastran',
      version=pyNastran.__version__,
      description=pyNastran.__desc__,
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Steven P. Doyle',
      author_email='mesheb82@gmail.com',
      url=pyNastran.__website__,
      license=pyNastran.__license__,
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
#          'numpy >= 1.3.0',
#          'scipy >= 0.6.0',
#          'vtk'
#          'wx'
      ],
      entry_points = {
        'console_scripts': [
                            'test_bdf = pyNastran.bdf.test.test_bdf:main',
                            'test_op2 = pyNastran.op2.test.test_op2:main',
                           ]
      },
      )
