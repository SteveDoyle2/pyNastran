"""helper for setup.py"""
import os
import sys
from typing import List

# Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
CLASSIFIERS = [
    'Natural Language :: English',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
]

PYTHON_REQUIRES = '>=3.10'

EXCLUDE_WORDS = [
    'pyNastran.f06.dev',
    'pyNastran.op2.dev', 'pyNastran.op2.dev.original',
    'pyNastran.converters.dev',

    'pyNastran.dev',
    'pyNastran.dev.xdb',
    'pyNastran.dev.bdf_vectorized', 'pyNastran.dev.bdf_vectorized.cards',
]


# features in packages used by pyNastran
# numpy
#  - 1.12 min for 3.6
#  - 1.13: adds axis support to unique
#  - 1.14: adds encoding support to savetxt (unused)
#  - 1.14: adds proper writing of np.savetxt for open file objects
#          (used for unicode savetxt using with statement) in Python 3.6
#  - 1.15: min for Python 3.7? I guess 1.14 is fine for a requirement...
#  - 1.19.4 : buggy
# scipy:
#  - 0.18.1: fixed kdtree used by nodal equivalencing; min for Python 2.7
#  - 0.19:   min for Python 3.6
#  - 0.19:   min for Python 3.7?; last 0.x release
# matplotlib:
#  - 1.5: min for Python 2.7; last 1.x release
#  - 2.0: adds C0/C1 colors (use colors from default colormap);
#    min for Python 3.6
#  - 2.1: adds plt.subplots support (untested?)
#  - 2.2: min for Python 3.7

# the packages that change requirements based on python version
REQS = {
    '3.8' : {
        'numpy' : ('1.17', '>=1.17.3,!=1.19.4'),
        'scipy' : ('1.3.2', '>=1.3.2'),
        'matplotlib' : ('3.1.2', '>=3.1.2'),
    },
    '3.9' : {
        'numpy' : ('1.19.3', '>=1.19.3,!=1.19.4'),
        'scipy' : ('1.5.4', '>=1.5.4'),
        'matplotlib' : ('3.4.0', '>=3.4.0'),
    },
    '3.10' : {
        'numpy' : ('1.21.2', '>=1.21.2'),
        'scipy' : ('1.7.1', '>=1.7.1'),
        'matplotlib' : ('3.5.0', '>=3.5.0'),
    },
}
MAX_VERSION = '3.10'

def check_python_version() -> None:
    """verifies the python version"""
    imajor, minor1, minor2 = sys.version_info[:3]
    if sys.version_info < (3, 9, 0):  # 3.9.4 used
        sys.exit('Upgrade your Python to 3.9+; version=(%s.%s.%s)' % (
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


def get_package_requirements(is_gui: bool=True, add_vtk_qt: bool=True,
                             python_version: str=None, bdist: bool=False):
    """
    Gets the requirements for setup.py

    Parameters
    ----------
    is_gui: bool; default=True
        add matplotlib, qtpy, pillow, imageio
        not vtk or pyqt/pyside because they're harder to install
    add_vtk_qt : bool; default=True
        does nothing?
    python_version: str; default=None -> sys.version_info
        allows us to get dynamic requirements
    bdist: bool; default=False
        loosen the requirements on numpy, scipy, etc.

    Returns
    -------
    all_reqs : Dict[name, version]
        name : str
            the name of the requirement (e.g., 'numpy')
        version: List[int]
            the required version number -> {'numpy': [1, 2, 3],}
    install_requires : List[str]
        the requirements -> ['numpy>=1.2.3',]

    """
    if python_version is None:
        python_version = '%s.%s' % sys.version_info[:2]

    if python_version not in REQS:
        python_version = MAX_VERSION
    vreqs = REQS[python_version]

    all_reqs = {}

    is_continuous_integration = (
        #'APPVEYOR' in os.environ or
        #'READTHEDOCS' in os.environ
        'TRAVIS' in os.environ or
        'TRAVIS_PYTHON_VERSION' in os.environ or
        'GITHUB_ACTOR' in os.environ
    )
    #user_name = getpass.getuser()
    #user_name not in ['travis']

    is_rtd = 'READTHEDOCS' in os.environ

    #if is_dev or is_gui:
        #try:
            #import vtk
            #vtk_version = '.'.join(vtk.VTK_VERSION.split('.'))
            #all_reqs['vtk'] = vtk_version
            #if vtk_version < '7.0.0':
                #print("vtk.VTK_VERSION = %r < '7.0.0'" % vtk.VTK_VERSION)
                #install_requires.append('vtk >= 7.0.0')
        #except ImportError:
            #install_requires.append('vtk >= 7.0.0')  # 8.x used

    install_requires = []

    found_numpy = False
    if is_rtd:
        install_requires.append('numpy')
    else:
        version_check, required_version = vreqs['numpy']
        if bdist:
            all_reqs['numpy'] = required_version
            install_requires.append('numpy %s' % required_version) # 1.22.3 used
        else:
            found_numpy = _add_numpy(version_check, required_version,
                                     all_reqs, install_requires)

    if is_rtd:
        install_requires.append('scipy')
    else:
        version_check, required_version = vreqs['scipy']
        if bdist:
            all_reqs['scipy'] = required_version
            install_requires.append('scipy %s' % required_version)  # 1.8.0 used
        else:
            _add_scipy(version_check, required_version,
                       all_reqs, install_requires)

    if is_gui:
        version_check, required_version = vreqs['matplotlib']
        if bdist:
            all_reqs['matplotlib'] = required_version
            install_requires.append('matplotlib %s' % required_version)  # 3.2.0 used
        else:
            _add_matplotlib(version_check, required_version,
                            all_reqs, install_requires)

    required_version_str = '1.4.0'
    if bdist:
        all_reqs['cpylog'] = f'>= {required_version_str}'
        install_requires.append(f'cpylog >= {required_version_str}')  # 1.4.0 used
    else:
        _add_cpylog(required_version_str, all_reqs, install_requires)

    if not is_rtd:
        _add_typish(all_reqs, install_requires)
        #_add_nptyping(all_reqs, install_requires)

    if bdist:
        all_reqs['docopt-ng'] = '>= 0.8.1'
        install_requires.append('docopt-ng >= 0.8.1')  # 0.8.1 used
    else:
        _add_docopt(all_reqs, install_requires)

    if is_rtd:
        pass
    elif is_gui:
        _add_qtpy(all_reqs, install_requires)
        _add_pillow(all_reqs, install_requires)
        _add_imageio(found_numpy, all_reqs, install_requires)


    #is_windows = 'nt' in os.name
    if is_continuous_integration:
        #install_requires.append('coverage>=4.4.2')
        #install_requires.append('python-coveralls>=2.9')
        #install_requires.append('coveralls>=1.7')
        install_requires.append('codecov')
        #install_requires.append('coverage')

    #print(all_reqs)
    print('install_requires =', install_requires)
    return all_reqs, install_requires

def _add_docopt(all_reqs, install_requires):
    required_version_str = '0.8.1'
    try:
        import docopt
        iver = int_version('docopt', docopt.__version__)
        all_reqs['docopt-ng'] = str_version(iver)
        if iver < [0, 8, 1]:
            print(f'docopt.__version__ = {docopt.__version__!r} < {required_version_str!r}')
            all_reqs['docopt-ng'] = f'>= {required_version_str}'
            install_requires.append(f'docopt-ng >= {required_version_str}')
    except ImportError:
        all_reqs['docopt-ng'] = f'>= {required_version_str}'
        install_requires.append(f'docopt-ng >= {required_version_str}')  # 0.8.1 used

def _add_numpy(version_check, required_version,
               all_reqs, install_requires):
    found_numpy = False
    try:
        import numpy as np
        found_numpy = True
    except RuntimeError:
        raise RuntimeError(f'numpy=1.19.4 is buggy; install a different version')
    except ImportError:
        all_reqs['numpy'] = required_version
        install_requires.append(f'numpy {required_version}') # 1.22.3 used

    if found_numpy:
        sver = np.lib.NumpyVersion(np.__version__)
        iver = int_version('numpy', sver.version)
        all_reqs['numpy'] = sver.version

        iversion_check = int_version('numpy', version_check)
        #srequired_version = int_version('numpy', required_version)
        #print('numpy %r %r' % (sver, iversion_check))
        if iver < iversion_check:
            print("numpy.__version__ = %r < %s" % (np.__version__, version_check))
            install_requires.append(f'numpy {required_version}')
            all_reqs['numpy'] = version_check
            install_requires.append(f'numpy {required_version}')
    return found_numpy

def _add_scipy(version_check, required_version,
               all_reqs, install_requires):
    try:
        import scipy
        sver = scipy.version.short_version
        iver = int_version('scipy', sver)
        all_reqs['scipy'] = sver

        iversion_check = int_version('scipy', version_check)
        #srequired_version = int_version('scipy', required_version)
        #print('scipy %r %r' % (sver, iversion_check))
        #print(iver, iversion_check)
        if iver < iversion_check:
            print("scipy.version.short_version = %r < %r" % (
                scipy.version.short_version, version_check))
            all_reqs['scipy'] = required_version
            install_requires.append(f'scipy {required_version}')
    except ImportError:
        all_reqs['scipy'] = required_version
        install_requires.append(f'scipy {required_version}')  # 1.8.0 used

def _add_matplotlib(version_check, required_version,
                    all_reqs, install_requires):
    try:
        import matplotlib
        iver = int_version('matplotlib', matplotlib.__version__)
        all_reqs['matplotlib'] = str_version(iver)
        iversion_check = int_version('matplotlib', version_check)
        if iver < iversion_check:
            print("matplotlib.__version__ = %r < %r" % (matplotlib.__version__, version_check))
                #matplotlib.__version__, str_version(iversion_check)))
            all_reqs['matplotlib'] = required_version
            install_requires.append(f'matplotlib {required_version}')
    except ImportError:
        all_reqs['matplotlib'] = required_version
        install_requires.append(f'matplotlib {required_version}')  # 3.5.1 used

def _add_cpylog(required_version_str, all_reqs, install_requires):
    try:
        import cpylog
        iver = int_version('cpylog', cpylog.__version__)
        all_reqs['cpylog'] = str_version(iver)
        if iver < [1, 4, 0]:
            print(f"cpylog.__version__ = {cpylog.__version__!r} != {required_version_str!r}")
            all_reqs['cpylog'] = f'>= {required_version_str}'
            install_requires.append(f'cpylog >= {required_version_str}')
    except ImportError:
        all_reqs['cpylog'] = f'>= {required_version_str}'
        install_requires.append(f'cpylog >= {required_version_str}')  # 1.4.0 used

def _add_typish(all_reqs, install_requires):
    """
    actual requirement somewhere between 1.6.0 and 1.9.1
    1.5.3 fails
    1.6.0 installs (does it work?)
    1.7.0 according to nptyping
    1.9.1 installs
    """
    required_version_str = '1.7.0'
    try:
        import typish
        iver = int_version('typish', typish.__version__)
        all_reqs['typish'] = str_version(iver)
        if iver < [1, 7, 0]:
            print(f"typish.__version__ = {typish.__version__!r} != {required_version_str!r}")
            all_reqs['typish'] = f'>= {required_version_str}'
            install_requires.append(f'typish >= {required_version_str}')
    except ImportError:
        all_reqs['typish'] = f'>= {required_version_str}'
        install_requires.append(f'typish >= {required_version_str}')  # 1.3.1 used


def _add_nptyping(all_reqs, install_requires):
    required_version_str = '>= 1.0.1, !=1.1.0'
    try:
        import nptyping
        iver = int_version('nptyping', nptyping.__version__)
        all_reqs['nptyping'] = str_version(iver)
        if iver < [1, 0, 1] or iver == [1, 1, 0]:
            print(f"nptyping.__version__ = {nptyping.__version__!r} not {required_version_str!r}")
            all_reqs['nptyping'] = required_version_str
            install_requires.append(f'nptyping {required_version_str}')
    except ImportError:
        all_reqs['nptyping'] = '>= 1.0.1'
        install_requires.append('nptyping >= 1.0.1, !=1.1.0')  # 1.0.1 used

def _add_qtpy(all_reqs, install_requires):
    try:
        import qtpy
        iver = int_version('qtpy', qtpy.__version__)
        all_reqs['qtpy'] = str_version(iver)
        if iver < [1, 4, 0]:
            print("qtpy.__version__ = %r < '1.4.0'" % qtpy.__version__)
            all_reqs['qtpy'] = '>= 1.4.0'
            install_requires.append('qtpy >= 1.4.0')
    #except ImportError:  # also PythonQtError
    except: # ImportError, PythonQtError
        # how do you import PythonQtError?
        all_reqs['qtpy'] = '>= 1.4.0'
        install_requires.append('qtpy >= 1.4.0')  # 1.9.0 used


def _add_pillow(all_reqs, install_requires):
    try:
        import PIL
        iver = int_version('pillow', PIL.__version__)
        all_reqs['pillow'] = str_version(iver)

        # at least 5.2.0, but not 7.1.0
        if iver > [7, 1, 0]:
            pass
            #print("pillow.__version__ = %r > '5.2.0', !='7.1.0" % PIL.__version__)
            #all_reqs['pillow'] = 'pillow >=5.2.0, !7.1.0'
            #install_requires.append('pillow >= 5.2.0, !=7.1.0')
        elif iver == [7, 1, 0]:
            print("pillow.__version__ = %r > '5.2.0', !='7.1.0" % PIL.__version__)
            all_reqs['pillow'] = 'pillow >=5.2.0, !=7.1.0'
            install_requires.append('pillow >= 5.2.0, !=7.1.0')
        elif iver < [5, 2, 0]:
            print("pillow.__version__ = %r > '5.2.0', !=7.1.0" % PIL.__version__)
            all_reqs['pillow'] = 'pillow >=5.2.0, !=7.1.0'
            install_requires.append('pillow >= 5.2.0, !=7.1.0')
    except ImportError:
        #print('missing pillow...')
        all_reqs['pillow'] = 'pillow >=5.2.0, !=7.1.0'
        install_requires.append('pillow >= 5.2.0, !=7.1.0')  # 7.1.1 used


def _add_imageio(found_numpy, all_reqs, install_requires):
    # pip messes up dependency resolution if you use newer versions of numpy
    imageio_str_ver = '>= 2.2.0, <3'

    try:
        import imageio
        imageio_int_ver = int_version('imageio', imageio.__version__)
        if imageio_int_ver > [3, ]:
            print("imageio.version = %r > 3.0'" % imageio.__version__)
            all_reqs['imageio'] = imageio_str_ver
            install_requires.append(f'imageio {imageio_str_ver}')
        elif imageio_int_ver < [2, 2, 0]:
            print("imageio.version = %r < '2.2.0'" % imageio.__version__)
            all_reqs['imageio'] = imageio_str_ver
            install_requires.append(f'imageio {imageio_str_ver}')
        else:
            all_reqs['imageio'] = imageio.__version__
    except ImportError:
        all_reqs['imageio'] = imageio_str_ver
        install_requires.append(f'imageio {imageio_str_ver}')  # 2.6.1 used

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


def cat_files(*filenames, encoding='utf8', sep='\n'):
    """Get the long description from the relevant file"""
    # thanks to harold!
    buf = []
    for filename in filenames:
        with open(filename, encoding=encoding) as file_obj:
            buf.append(file_obj.read())
    return sep.join(buf)

LONG_DESCRIPTION = cat_files('README.rst')
#assert '\r' not in LONG_DESCRIPTION, LONG_DESCRIPTION
#for i, line in enumerate(LONG_DESCRIPTION.split('\n')):
    #print(f'%4i: %s' % (i, line))
    # print(LONG_DESCRIPTION)
