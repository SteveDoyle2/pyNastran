"""helper for setup.py"""
import os
import sys
from typing import List

EXCLUDE_WORDS = [
    'pyNastran.f06.dev',
    'pyNastran.op2.dev', 'pyNastran.op2.dev.original',
    'pyNastran.converters.dev',

    'pyNastran.dev',
    'pyNastran.dev.xdb',
    'pyNastran.dev.bdf_vectorized', 'pyNastran.dev.bdf_vectorized.cards',
    'pyNastran.dev.bdf_vectorized3', 'pyNastran.dev.op2_vectorized3',
]

CONSOLE_SCRIPTS = [
    #'run_nastran_double_precision = pyNastran.bdf.test.run_nastran_double_precision:cmd_line',
    'test_bdf  = pyNastran.bdf.test.test_bdf:main',
    'test_op2  = pyNastran.op2.test.test_op2:main',
    'test_op4  = pyNastran.op4.test.test_op4:main',

    'pyNastranGUI = pyNastran.gui.gui:cmd_line',
    'bdf = pyNastran.bdf.mesh_utils.utils:cmd_line',
    'f06 = pyNastran.f06.utils:cmd_line',
    'format_converter = pyNastran.converters.format_converter:cmd_line_format_converter',
    'abaqus_to_nastran = pyNastran.converters.abaqus.abaqus_to_nastran:cmd_abaqus_to_nastran',
    'test_pynastrangui = pyNastran.gui.test.test_gui:main',

    'test_bdfv = pyNastran.dev.bdf_vectorized3.test.test_bdf:main',
    #'pyNastranv = pyNastran.dev.bdf_vectorized.solver.solver:main',
    #'test_bdfv2 = pyNastran.dev.bdf_vectorized2.test.test_bdf:main',
    #'test_abaqus = pyNastran.converters.abaqus.test_abaqus:main',
    #'nastran_to_code_aster = pyNastran.converters.dev.code_aster.nastran_to_code_aster:main',
]


def check_python_version() -> None:
    """verifies the python version"""
    imajor, minor1, minor2 = sys.version_info[:3]
    if sys.version_info < (3, 10, 0):  # 3.10.4 used
        sys.exit('Upgrade your Python to 3.10+; version=(%s.%s.%s)' % (
            imajor, minor1, minor2))

def int_version(name: str, version: str) -> List[int]:
    """splits the version into a tuple of integers"""
    sversion = version.split('-')[0]
    #numpy
    #scipy
    #matplotlib
    #qtpy
    #vtk
    #cpylog
    #pyNastran
    if 'rc' not in name:
        # it's gotta be something...
        # matplotlib3.1rc1
        sversion = sversion.split('rc')[0]

    try:
        return [int(val) for val in sversion.split('.')]
    except ValueError:
        raise SyntaxError('cannot determine version for %s %s' % (name, sversion))


def str_version(version: str) -> str:
    """converts a tuple of integers to a version number"""
    return '.'.join(str(versioni) for versioni in version)


def update_version_file():
    """
    Creates the version.py file with the github string
    to lock down the version when the user the following
    on the dev version:
    >>> python setup.py install

    instead of:
    >>> python setup.py develop

    This is handy
    """
    import pyNastran
    args = sys.argv

    # 'setup.py install'        will not quit out
    # 'setup.py install --user' will not quit out
    # 'pip install .'           ???
    if 'develop' in args:
        return

    #if 'install' in args and len(args) >= 3:   # don't care about dev
        # args =  ['setup.py', 'install', '--user']
        # args =  ['setup.py', 'install', '-home=<path>']
        # args =  ['setup.py', 'install', '--prefix=<path>']
        #
        #pass
    pkg_path = pyNastran.__path__[0]
    init_filename = os.path.join(pkg_path, '__init__.py')
    version_filename = os.path.join(pkg_path, 'version.py')
    with open(version_filename, 'w') as version_file:
        version_file.write(f'__version__ = {pyNastran.__version__!r}\n')
        version_file.write(f'__releaseDate__ = {pyNastran.__releaseDate__!r}\n')
        version_file.write(f'__releaseDate2__ = {pyNastran.__releaseDate2__!r}\n')

    with open(init_filename, 'r') as init_file:
        data = init_file.read()

    data2 = data.replace('is_installed = False', 'is_installed = True')
    with open(init_filename, 'w') as init_file_out:
        data = init_file_out.write(data2)
    #__version__ = '1.3.0+%s' % revision
