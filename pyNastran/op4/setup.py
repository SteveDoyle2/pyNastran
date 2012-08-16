""" Example of building a module with a Cython file. See the distutils
and numpy distutils documentations for more info:
http://docs.scipy.org/doc/numpy/reference/distutils.html
"""
# Author: Gael Varoquaux
# License: BSD

import numpy
from Cython.Distutils import build_ext


def configuration(parent_package='', top_path=None):
    """ Function used to build our configuration.
    """
    from numpy.distutils.misc_util import Configuration

    # The configuration object that hold information on all the files
    # to be built.
    config = Configuration('', parent_package, top_path)
    config.add_extension('cop4',
                         sources=['op4.pyx'],
                         # libraries=['m'],
                         depends=['op4.c'],
                         extra_compile_args=["-g"],
                         extra_link_args=["-g"],
                         include_dirs=[numpy.get_include()])
    return config


if __name__ == '__main__':
    # Retrieve the parameters of our local configuration
    params = configuration(top_path='').todict()

    # Override the C-extension building so that it knows about '.pyx'
    # Cython files
    params['cmdclass'] = dict(build_ext=build_ext)

    # Call the actual building/packaging function (see distutils docs)
    from numpy.distutils.core import setup
    setup(**params)
