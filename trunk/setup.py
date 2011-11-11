from setuptools import setup, find_packages
import sys, os

version = '0.2dev'

setup(name='pyNastran',
      version=version,
      description="General Purpose nastran file reader/writer",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='mesheb82',
      author_email='',
      url='',
      license='GPL V3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
#          'numpy >= 1.3.0',
#          'scipy >= 0.6.0',  #Last built with 0.7.1
#          'matplotlib >= 0.99.1',
      ],
      entry_points = {
      },
      )
