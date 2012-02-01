from distutils.core      import setup
from distutils.extension import Extension
from Cython.Distutils    import build_ext

from numpy import get_include

ext = Extension(
        "loadop4",
        ['loadop4.pyx', '_loadop4.c'],
        language='c',
        include_dirs=[get_include()],
        )

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext],
)
