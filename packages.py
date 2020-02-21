"""helper for setup.py"""
import os
import sys

# features in packages used by pyNastran
# numpy
#  - 1.12 min for 3.6
#  - 1.13: adds axis support to unique
#  - 1.14: adds encoding support to savetxt (unused)
#  - 1.14: adds proper writing of np.savetxt for open file objects
#          (used for unicode savetxt using with statement) in Python 3.6
#  - 1.15: min for Python 3.7? I guess 1.14 is fine for a requirement...
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
    '3.7' : {
        'numpy' : ('1.14', '>=1.14'),
        'scipy' : ('1.0', '>=1.0'),
        'matplotlib' : ('2.2', '>=2.2'),  # 2.2.4 adds Python 3.7 support
    },
    '3.8' : {  # TODO: not updated
        'numpy' : ('1.14', '>=1.14'),
        'scipy' : ('1.0', '>=1.0'),
        'matplotlib' : ('2.2', '>=2.2'),  # 2.2.4 adds Python 3.7 support
    },
}

def check_python_version():
    """verifies the python version"""
    imajor, minor1, minor2 = sys.version_info[:3]
    if sys.version_info < (3, 7, 0):  # 3.7.4 used
        sys.exit('Upgrade your Python to 3.7+; version=(%s.%s.%s)' % (
            imajor, minor1, minor2))


def int_version(name, version):
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


def str_version(version):
    """converts a tuple of intergers to a version number"""
    return '.'.join(str(versioni) for versioni in version)


def get_package_requirements(is_gui=True, add_vtk_qt=True, python_version=None):
    """gets the requirements for setup.py"""
    if python_version is None:
        python_version = '%s.%s' % sys.version_info[:2]

    if python_version not in REQS:
        python_version = '3.7'
    vreqs = REQS[python_version]

    all_reqs = {}

    #is_dev = (
        #'TRAVIS' in os.environ or
        #'APPVEYOR' in os.environ or
        #'READTHEDOCS' in os.environ
    #)
    is_travis = 'TRAVIS' in os.environ or 'TRAVIS_PYTHON_VERSION' in os.environ
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

    if is_rtd:
        install_requires.append('numpy')
    else:
        version_check, required_version = vreqs['numpy']
        try:
            import numpy as np
            sver = np.lib.NumpyVersion(np.__version__)
            iver = int_version('numpy', sver.version)
            all_reqs['numpy'] = sver.version

            iversion_check = int_version('numpy', version_check)
            #srequired_version = int_version('numpy', required_version)
            #print('numpy %r %r' % (sver, iversion_check))
            if iver < iversion_check:
                print("numpy.__version__ = %r < %s" % (np.__version__, version_check))
                install_requires.append('numpy %s' % required_version)
                all_reqs['numpy'] = version_check
                install_requires.append('numpy %s' % required_version)
        except ImportError:
            all_reqs['numpy'] = required_version
            install_requires.append('numpy %s' % required_version) # 1.18.1 used

    if is_rtd:
        install_requires.append('scipy')
    else:
        version_check, required_version = vreqs['scipy']
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
                install_requires.append('scipy %s' % required_version)
        except ImportError:
            all_reqs['scipy'] = required_version
            install_requires.append('scipy %s' % required_version)  # 1.4.1 used

    if is_gui:
        version_check, required_version = vreqs['matplotlib']
        try:
            import matplotlib
            iver = int_version('matplotlib', matplotlib.__version__)
            all_reqs['matplotlib'] = str_version(iver)
            iversion_check = int_version('matplotlib', version_check)
            if iver < iversion_check:
                print("matplotlib.__version__ = %r < %r" % (matplotlib.__version__, version_check))
                    #matplotlib.__version__, str_version(iversion_check)))
                all_reqs['matplotlib'] = required_version
                install_requires.append('matplotlib %s' % required_version)
        except ImportError:
            all_reqs['matplotlib'] = required_version
            install_requires.append('matplotlib %s' % required_version)  # 3.2.0 used


    required_version_str = '1.3.1'
    try:
        import cpylog
        iver = int_version('cpylog', cpylog.__version__)
        all_reqs['cpylog'] = str_version(iver)
        if iver < [1, 3, 1]:
            print(f"cpylog.__version__ = {cpylog.__version__!r} != {required_version_str!r}")
            all_reqs['cpylog'] = f'>= {required_version_str}'
            install_requires.append(f'cpylog >= {required_version_str}')
    except ImportError:
        all_reqs['cpylog'] = f'>= {required_version_str}'
        install_requires.append(f'cpylog >= {required_version_str}')  # 1.3.1 used

    try:
        import nptyping
        #iver = int_version('nptyping', nptyping.__version__)
        #all_reqs['nptyping'] = str_version(iver)
        #if iver < [0, 3, 1]:
            #print("nptyping.__version__ = %r < '0.3.1'" % nptyping.__version__)
        all_reqs['nptyping'] = '>= 0.3.1'
        install_requires.append('nptyping >= 0.3.1')
    except ImportError:
        all_reqs['nptyping'] = '>= 0.3.1'
        install_requires.append('nptyping >= 0.3.1')  # 0.3.1 used

    try:
        import docopt
        iver = int_version('docopt', docopt.__version__)
        all_reqs['docopt-ng'] = str_version(iver)
        if iver < [0, 7, 0]:
            print("docopt.__version__ = %r < '0.7.0'" % docopt.__version__)
            all_reqs['docopt-ng'] = '>= 0.7.0'
            install_requires.append('docopt-ng >= 0.7.0')
    except ImportError:
        all_reqs['docopt-ng'] = '>= 0.7.0'
        install_requires.append('docopt-ng >= 0.7.0')  # 0.7.2 used

    if is_rtd:
        pass
    elif is_gui:
        try:
            import qtpy
            iver = int_version('qtpy', qtpy.__version__)
            all_reqs['qtpy'] = str_version(iver)
            if iver < [1, 4, 0]:
                print("qtpy.__version__ = %r < '1.4.0'" % qtpy.__version__)
                all_reqs['qtpy'] = '>= 1.4.0'
                install_requires.append('qtpy >= 1.4.0')
        except ImportError:
            all_reqs['qtpy'] = '>= 1.4.0'
            install_requires.append('qtpy >= 1.4.0')  # 1.9.0 used

        try:
            import PIL
            iver = int_version('pillow', PIL.__version__)
            all_reqs['pillow'] = str_version(iver)
            if iver < [1, 4, 0]:
                print("pillow.__version__ = %r > '1.4.0'" % PIL.__version__)
                all_reqs['pillow'] = 'pillow >=1.4.0'
                install_requires.append('pillow >= 1.4.0')
        except ImportError:
            #print('missing pillow...')
            all_reqs['pillow'] = 'pillow >=1.4.0'
            install_requires.append('pillow >= 1.4.0')  # 6.2.0 used

        try:
            import imageio
            if imageio.__version__ < '2.2.0':
                print("imageio.version = %r < '2.2.0'" % imageio.__version__)
                all_reqs['imageio'] = '>= 2.2.0'
                install_requires.append('imageio >= 2.2.0')
            else:
                all_reqs['imageio'] = imageio.__version__
        except ImportError:
            all_reqs['imageio'] = '>= 2.2.0'
            install_requires.append('imageio >= 2.2.0')  # 2.6.1 used


    #is_windows = 'nt' in os.name
    if is_travis:
        #install_requires.append('coverage>=4.4.2')
        #install_requires.append('python-coveralls>=2.9')
        #install_requires.append('coveralls>=1.7')
        install_requires.append('codecov')
        #install_requires.append('coverage')

    #print(all_reqs)
    print('install_requires =', install_requires)
    return all_reqs, install_requires

def update_version_file():
    import pyNastran
    if 'install' not in sys.argv:
        return

    pkg_path = pyNastran.__path__[0]
    init_filename = os.path.join(pkg_path, '__init__.py')
    version_filename = os.path.join(pkg_path, 'version.py')
    with open(version_filename, 'w') as version_file:
        version_file.write(f'__version__ = {pyNastran.__version__!r}\n')
        version_file.write(f'__releaseDate__ = {pyNastran.__releaseDate__!r}\n')

    with open(init_filename, 'r') as init_file:
        data = init_file.read()

    data2 = data.replace('is_installed = False', 'is_installed = True')
    with open(init_filename, 'w') as init_file_out:
        data = init_file_out.write(data2)

    #__version__ = '1.3.0+%s' % revision
