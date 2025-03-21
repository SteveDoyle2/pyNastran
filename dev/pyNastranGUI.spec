# -*- PyInstaller input file -*-
# -*- mode: python           -*-
import os
import sys
try:
    import pyInstaller
    pyInstaller_path = [os.path.dirname(pyinstaller.__file__)]
    print("pyInstaller_path = %r" % pyInstaller_path)
except ImportError:
    pyInstaller_path = []
    #asdf


#pyInstaller_path = [r'C:\NASA\dev\pyinstaller']
#pyInstaller_path = [r'C:\NASA\m4\formats\git\pyNastran\env\Lib\site-packages\PyInstaller']
print(f'sys.prefix = {sys.prefix}')
#print(f'pyInstaller_path = {pyInstaller_path}')

IS_H5PY = False
INCLUDE_BDFV = False
INCLUDE_PANDAS = False
INCLUDE_MATPLOTLIB = True

DEBUG = False
IS_RELEASE = True
USE_TODAY = True
IS_ONEDIR = False # True=onedir; False=onefile

MAKE_SPLASH = True # requires tk
BUILD_EXE = True
BUILD_COLLECT = False


# get pyNastran location
print('sys.version_info.major =', sys.version_info.major)
pkg_path = os.path.abspath(os.path.join('.', '..', 'pyNastran'))
init_filename = os.path.join(pkg_path, '__init__.py')
assert os.path.exists(init_filename), init_filename

#--------------------------------------------------------------------------
def get_git_revision_short_hash_date() -> tuple[str, tuple[str, str, str]]:
    """determines the git revision; only works if the packages was checked
    out using git"""
    import subprocess
    #ghash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])

    # independent of pyNastran location as long as there is a git folder
    #   what about if you use setup_user.py install?
    #   what about if you don't have git?
    # can raise a subprocess.CalledProcessError, which means the return code != 0
    ghash_bytes = subprocess.check_output(['git', 'describe', '--always'],
                                          cwd=pkg_path)
    ghash = ghash_bytes.decode('utf-8').rstrip()

    #git show -s --format=%cd --date=short
    #Mon Mar 11 23:00:39 2024 -0700
    #
    #git show --no-patch --format=%ai
    #2024-03-11 16:37:20 -0700
    date_bytes = subprocess.check_output(
        ['git', 'show', '-s', '--format=%cd', '--date=short'],
        cwd=pkg_path)
    date = date_bytes.decode('utf-8').rstrip()
    year, month, day = date.split('-')
    month = month.lstrip('0')
    day = day.lstrip('0')
    date_out = (year, month, day)
    # 1.5.0+dev.0eccfa918
    return 'dev.%s' % ghash, date_out

revision, (year, month, day) = get_git_revision_short_hash_date()

#--------------------------------------------------------------------------

#env_path = r'C:\NASA\m4\formats\git\pyNastran\env'
env_path = r'C:\env312'
site_packages_path = os.path.join(env_path, 'Lib', 'site-packages')
#vtk_path = os.path.join(site_packages_path, 'vtkmodules')
vtk_lib_path = os.path.join(site_packages_path, 'vtk.libs')
if not os.path.exists(site_packages_path):
    from pyNastran.utils import print_bad_path
    msg = print_bad_path(site_packages_path)
    raise RuntimeError(msg)

#-------------------------------------------------------------------------
## this block gets/sets the version so it doesnt use git

#import pyNastran
#pyNastran.is_pynastrangui_exe = True
#pyNastran.__version__ = '%r' % pyNastran.__version__
#pyNastran.__releaseDate__ =  '%r' % pyNastran.__releaseDate__
#pyNastran.__releaseDate2__ = '%r' % datei.strftime('%d %B %Y')

import os
import sys
import shutil
import datetime

# getting pyNastran version without using the __init__.py file
# because we need to hack on it
import subprocess
try:
    ghash = subprocess.check_output(['git', 'describe', '--always'], cwd=os.path.dirname(init_filename))
    ghash = ghash.decode('utf-8').rstrip()
except:
    # git isn't installed
    ghash = 'no.checksum.error'

# hacking on the __init___.py file to set:
#   is_pynastrangui_exe = True
init_filename_bkp = init_filename + '.bkp'
shutil.copyfile(init_filename, init_filename_bkp)

print('opening %s to hack the is_pynastrangui_exe flag' % init_filename)
with open(init_filename, 'r') as init_file:
    lines = init_file.readlines()

# update the init.py with:
# - is_pynastrangui_exe = True
#
# and get the version
#
lines2 = []
for line in lines:
    if 'is_pynastrangui_exe = False' in line:
        line = 'is_pynastrangui_exe = True\n'
    elif '__version__ = ' in line:
        # __version__ = '1.1.0+%s' % revision
        if "'" not in line:
            assert IS_RELEASE
            continue
        version_fmt = line.split("'")[1]  # '1.1.0+%s'
    elif 'import' not in line and '__releaseDate2__' in line:
        __releaseDate2__ = line.split("'")[1]
    elif 'import' not in line and '__releaseDate__' in line:
        #raise RuntimeError(line)
        __releaseDate__ = line.split("'")[1]
    #elif 'import' not in line and 'revision' in line and 'def' not in line and '"""' not in line:
        #raise RuntimeError(line)
        #revision = line.split("'")[1]
    elif 'import' not in line and '__version_release__' in line:
        __version_release__ = line.split("'")[1]
    lines2.append(line)

assert len(lines2) > 8, lines2
with open(init_filename, 'w') as init_file:
    for line in lines2:
        init_file.write(line)

#version = '1.1.0+dev.%s' % ghash
version = version_fmt
if not IS_RELEASE:
    #print(version_fmt, ghash)
    # hack because version was '{version}.{revision}+dev.%s' or close enough
    version_fmt = '1.5.0+dev.%s'
    version = version_fmt % ghash

# write the version
version_filename = os.path.join(pkg_path, 'version.py')
with open(version_filename, 'w') as version_file:
    version_file.write('# -*- coding: utf-8 -*-\n')
    version_file.write('# this file is autogenerated by pyInstaller\n')
    version_file.write('__version_release__ = f%r\n' % __version_release__)
    version_file.write('revision = %r\n' % revision)
    version_file.write('__version__ = f%r\n' % version)
    if USE_TODAY:
        datei = datetime.date.today()
        version_file.write('__releaseDate__ = %r\n' % str(datei))  # 2016-2-5
        version_file.write('__releaseDate2__ = %r\n' % datei.strftime('%d %B %Y')) # 5 Feb 2016
    else:
        version_file.write('__releaseDate__ = %r\n' % __releaseDate__)
        version_file.write('__releaseDate2__ = %r\n' % __releaseDate2__)

#-------------------------------------------------------------------------

#a1 = os.path.join(pkg_path, 'bdf','bdf.py')
#a2 = os.path.join(pkg_path, 'op2','op2.py')
#a3 = os.path.join(pkg_path, 'f06','f06.py')
a4 = os.path.join(pkg_path, 'gui','gui.py')
a5 = os.path.join(pkg_path, 'version.py')

#analyze_files = [a1, a2, a3, a4]
analyze_files = [a4, a5]

icon_path = os.path.join(pkg_path, 'gui', 'icons')

icon_names = os.listdir(icon_path)
icons = []
for icon_name in icon_names:
    i1 = os.path.join('icons', icon_name)   # where to put the icon in the exe
    i2 = os.path.join(icon_path, icon_name) # where the icon is now
    icon = (i1, i2, 'DATA')
    icons.append(icon)

#icon_main = os.path.join(iconPath,'guiDemo_128.ico')
icon_main = 'logo.ico' # pyInstaller doesn't like paths in the ico
                       # it also won't handle png files
assert os.path.exists(icon_main), '%s doesnt exist' % icon_main

#-------------------------------------------------------------------------
# main code

#from PyInstaller import compat
#from os import listdir

##mkldir = compat.base_prefix + "/Lib/site-packages/numpy/core"
#logger = logging.getLogger(__name__)
#logger.info("MKL installed as part of numpy, importing that!")
#binaries = [(mkldir + "/" + mkl, '') for mkl in listdir(mkldir) if mkl.startswith('mkl_')]

#python_path = 'F:\Anaconda'
python_path = os.path.dirname(sys.executable)

mkl_dlls = [
    os.path.join(python_path, 'evns', 'py35', 'Library', 'bin', 'mkl_def3.dll')
    # others?
]

has_mkl_dlls = False
if mkl_dlls:
    mkl_dll_base = os.path.basename(mkl_dlls[0])
    #assert os.path.exists(mkl_dll), '%s doesnt exist' % mkl_dll
    has_mkl_dlls = os.path.exists(mkl_dll_base)

binaries = []
if sys.platform == 'win32':
    binaries = [
        ('msvcp100.dll', 'C:\\Windows\\System32\\msvcp100.dll', 'BINARY'),
        ('msvcr100.dll', 'C:\\Windows\\System32\\msvcr100.dll', 'BINARY'),
    ]
    if has_mkl_dlls:
        for mkl_dll in mkl_dlls:
            binaries.append(
                (mkl_dll_base, mkl_dll, 'BINARY')
            )
ver = sys.version_info
python_version = '%i%i' % (ver.major, ver.minor)
print('python_version', python_version)
print('python_path', python_path)

if 'Anaconda' not in python_path:
    python_version_dllname = 'python%s.dll' % (python_version)  # python37.dll
    binaries.append(
        ('python%s.dll' % python_version, os.path.join(python_path, python_version_dllname), 'BINARY'),
    )


INCLUDE_VTK_DLLS = False
vtk_imports = [
    'vtk',
    'vtkmodules',
    'vtkmodules.all',
    'vtkmodules.vtkCommonDataModel',
    'vtkmodules.vtkCommonExecutionModel',
    #'vtkmodules.vtkCommonCore',
    'vtkmodules.vtkCommonMath',
    'vtkmodules.vtkCommonTransforms',
    'vtkmodules.vtkFiltersCore',
    'vtkmodules.vtkFiltersExtraction',
    'vtkmodules.vtkFiltersGeneral',
    'vtkmodules.vtkFiltersGeometry',
    'vtkmodules.vtkFiltersPoints',
    'vtkmodules.vtkFiltersSources',
    'vtkmodules.vtkIOXMLParser',
    'vtkmodules.vtkRenderingCore',
    'vtkmodules.vtkRenderingUI',
    'vtkmodules.vtkRenderingLOD',
    'vtkmodules.vtkRenderingCore',
    'vtkmodules.vtkRenderingAnnotation',
    'vtkmodules.vtkInteractionStyle',
    'vtkmodules.vtkInteractionWidgets',
    'vtkmodules.vtkIOImage',
    'vtkmodules.vtkRenderingLabel',
    'vtkmodules.vtkFiltersHybrid',
    'vtkmodules.util.vtkConstants',
]

if INCLUDE_VTK_DLLS:
    # I never actually figured how to make vtk 9.3 build...it's hard
    # vtk.libs is a new directory
    # I think there's a "glew" DLL missing
    vtk_version = '9.3'
    assert os.path.exists(vtk_lib_path), vtk_lib_path
    pyInstaller_path.append(vtk_lib_path)

    vtk_lib_names = os.listdir(vtk_lib_path)
    vtk_name_to_lib_name = {
        name.split('-')[0] : name for name in vtk_lib_names
        if name.startswith('vtk') and name.endswith('.dll')
    }
    print(vtk_name_to_lib_name)
    for mod_name_all in vtk_imports:
        mod_name_split = mod_name_all.split('.')
        if len(mod_name_split) == 1:
            continue
        mod_name = mod_name_split[-1]

        print('mod', mod_name_all, mod_name)
        #vtkFiltersProgrammable.cp39-win_amd64.pyd
        #vtkCommonPython-9.2.dll
        #pyd_name = f'{mod_name}.cp{pythonversion}}-win_amd64.pyd'
        #dll_name = f'{mod_name}-{vtk_version}.dll'
        #dll_source_path = os.path.join(vtk_path, dll_name)
        if mod_name not in vtk_name_to_lib_name:
            print(f"vtk mod_name {mod_name} doesn't exist?")
            continue
        dll_name = vtk_name_to_lib_name[mod_name]
        dll_source_path = os.path.join(vtk_lib_path, dll_name)
        assert os.path.exists(dll_source_path), dll_source_path
        binary = (f'vtkmodules/{dll_name}', dll_source_path, 'BINARY')
        binaries.append(binary)
        print(binary)

    vtk_dll_modules = [
        # modules
        'vtkFiltersHybrid',
        'vtkFiltersPoints',
        'vtkFiltersSources',
        'vtkInteractionStyle',
        'vtkInteractionWidgets',
        'vtkIOImage',
        'vtkIOXMLParser',
        'vtkRenderingAnnotation',
        'vtkRenderingCore',
        'vtkRenderingLabel',
        'vtkRenderingLOD',
        'vtkRenderingUI',
        # lib
        'vtkCommonColor',
        'vtkFiltersModeling',
        'vtkFiltersTexture',
        'vtkfreetype',
        'vtkglew',
        'vtkImagingColor',
        'vtkImagingCore',
        'vtkImagingGeneral',
        'vtkImagingHybrid',
        'vtkImagingSources',
        'vtkjpeg',
        'vtkmetaio',
        'vtkpng',
        'vtkRenderingContext2D',
        'vtkRenderingOpenGL2',
        'vtktiff',
    ]
    if 0:
        vtk_dll_modules = list(set(vtk_dll_modules))
        for dll_base in vtk_dll_modules:
            dll_name = vtk_name_to_lib_name[dll_base]
            dll_source_path = os.path.join(vtk_lib_path, dll_name)
            binary1 = (f'vtkmodules/{dll_name}', dll_source_path, 'BINARY')
            binary2 = (f'vtk.libs/{dll_name}', dll_source_path, 'BINARY')
            binaries.append(binary1)
            binaries.append(binary2)
            print(binary)
    for dll_name in vtk_name_to_lib_name.values():
        dll_source_path = os.path.join(vtk_lib_path, dll_name)
        assert os.path.exists(dll_source_path), dll_source_path
        binary1 = (f'vtkmodules/{dll_name}', dll_source_path, 'BINARY')
        binary2 = (f'vtk.libs/{dll_name}', dll_source_path, 'BINARY')
        binary3 = (dll_name, dll_source_path, 'BINARY')
        binaries.append(binary1)
        binaries.append(binary2)
        binaries.append(binary3)
        #print(binary)
    #asdf

pathex = pyInstaller_path + [
    python_path,
    os.path.join(python_path, 'Lib'),
    os.path.join(python_path, 'Lib', 'site-packages'),

    pkg_path,
    os.path.join(pkg_path, 'bdf'),
    os.path.join(pkg_path, 'op2'),
    os.path.join(pkg_path, 'f06'),
    os.path.join(pkg_path, 'gui'),
]

pandas_imports = [
    'pandas', 'pandas._libs', 'pandas._libs.algos', 'pandas._libs.calendar',
    'pandas._libs.collections', 'pandas._libs.csv', 'pandas._libs.datetime',
    'pandas._libs.dateutil', 'pandas._libs.decimal', 'pandas._libs.distutils',
    'pandas._libs.groupby', 'pandas._libs.hashing', 'pandas._libs.hashtable',
    'pandas._libs.index', 'pandas._libs.interval', 'pandas._libs.join',
    'pandas._libs.json', 'pandas._libs.lib', 'pandas._libs.locale',
    'pandas._libs.numbers', 'pandas._libs.numpy', 'pandas._libs.operator',
    'pandas._libs.os', 'pandas._libs.pandas', 'pandas._libs.parsers',
    'pandas._libs.period', 'pandas._libs.pytz', 'pandas._libs.random',
    'pandas._libs.re', 'pandas._libs.reshape', 'pandas._libs.sparse',
    'pandas._libs.sys', 'pandas._libs.testing', 'pandas._libs.thread',
    'pandas._libs.time', 'pandas._libs.tslib', 'pandas._libs.warnings',
    'pandas._libs.window', 'pandas._version', 'pandas.api', 'pandas.api.types',
    'pandas.api.types.pandas', 'pandas.compat', 'pandas.compat.StringIO',
    'pandas.compat.__builtin__', 'pandas.compat.cPickle',
    'pandas.compat.cStringIO', 'pandas.compat.chainmap',
    'pandas.compat.chainmap_impl', 'pandas.compat.collections',
    'pandas.compat.copy', 'pandas.compat.dateutil', 'pandas.compat.distutils',
    'pandas.compat.functools', 'pandas.compat.httplib', 'pandas.compat.inspect',
    'pandas.compat.itertools', 'pandas.compat.numpy',
    'pandas.compat.numpy.distutils', 'pandas.compat.numpy.function',
    'pandas.compat.numpy.numpy', 'pandas.compat.numpy.pandas',
    'pandas.compat.numpy.re', 'pandas.compat.openpyxl_compat',
    'pandas.compat.pandas', 'pandas.compat.pickle',
    'pandas.compat.pickle_compat', 'pandas.compat.re', 'pandas.compat.struct',
    'pandas.compat.sys', 'pandas.compat.thread', 'pandas.compat.types',
    'pandas.compat.unicodedata', 'pandas.core', 'pandas.core.__future__',
    'pandas.core.algorithms', 'pandas.core.api', 'pandas.core.base',
    'pandas.core.bottleneck', 'pandas.core.categorical', 'pandas.core.codecs',
    'pandas.core.collections', 'pandas.core.common', 'pandas.core.computation',
    'pandas.core.computation.abc', 'pandas.core.computation.align',
    'pandas.core.computation.api', 'pandas.core.computation.ast',
    'pandas.core.computation.common', 'pandas.core.computation.datetime',
    'pandas.core.computation.distutils', 'pandas.core.computation.engines',
    'pandas.core.computation.eval', 'pandas.core.computation.expr',
    'pandas.core.computation.expressions', 'pandas.core.computation.functools',
    'pandas.core.computation.inspect', 'pandas.core.computation.itertools',
    'pandas.core.computation.numexpr', 'pandas.core.computation.numpy',
    'pandas.core.computation.operator', 'pandas.core.computation.ops',
    'pandas.core.computation.pandas', 'pandas.core.computation.pprint',
    'pandas.core.computation.pytables', 'pandas.core.computation.scope',
    'pandas.core.computation.struct', 'pandas.core.computation.sys',
    'pandas.core.computation.tokenize', 'pandas.core.computation.warnings',
    'pandas.core.config', 'pandas.core.config_init', 'pandas.core.contextlib',
    'pandas.core.copy', 'pandas.core.datetime', 'pandas.core.distutils',
    'pandas.core.dtypes', 'pandas.core.dtypes.api', 'pandas.core.dtypes.cast',
    'pandas.core.dtypes.collections', 'pandas.core.dtypes.common',
    'pandas.core.dtypes.concat', 'pandas.core.dtypes.datetime',
    'pandas.core.dtypes.dtypes', 'pandas.core.dtypes.generic',
    'pandas.core.dtypes.inference', 'pandas.core.dtypes.missing',
    'pandas.core.dtypes.numbers', 'pandas.core.dtypes.numpy',
    'pandas.core.dtypes.pandas', 'pandas.core.dtypes.re', 'pandas.core.dtypes.sys',
    'pandas.core.dtypes.warnings', 'pandas.core.frame', 'pandas.core.functools',
    'pandas.core.gc', 'pandas.core.generic', 'pandas.core.groupby',
    'pandas.core.index', 'pandas.core.indexes', 'pandas.core.indexes.__future__',
    'pandas.core.indexes.accessors', 'pandas.core.indexes.api',
    'pandas.core.indexes.base', 'pandas.core.indexes.category',
    'pandas.core.indexes.datetime', 'pandas.core.indexes.datetimelike',
    'pandas.core.indexes.datetimes', 'pandas.core.indexes.frozen',
    'pandas.core.indexes.functools', 'pandas.core.indexes.interval',
    'pandas.core.indexes.multi', 'pandas.core.indexes.numeric',
    'pandas.core.indexes.numpy', 'pandas.core.indexes.operator',
    'pandas.core.indexes.pandas', 'pandas.core.indexes.period',
    'pandas.core.indexes.range', 'pandas.core.indexes.sys',
    'pandas.core.indexes.timedeltas', 'pandas.core.indexes.warnings',
    'pandas.core.indexing', 'pandas.core.internals', 'pandas.core.itertools',
    'pandas.core.json', 'pandas.core.keyword', 'pandas.core.missing',
    'pandas.core.nanops', 'pandas.core.numpy', 'pandas.core.operator',
    'pandas.core.ops', 'pandas.core.pandas', 'pandas.core.panel',
    'pandas.core.panel4d', 'pandas.core.panelnd', 'pandas.core.re',
    'pandas.core.resample', 'pandas.core.reshape', 'pandas.core.reshape.api',
    'pandas.core.reshape.concat', 'pandas.core.reshape.copy',
    'pandas.core.reshape.itertools', 'pandas.core.reshape.merge',
    'pandas.core.reshape.numpy', 'pandas.core.reshape.pandas',
    'pandas.core.reshape.pivot', 'pandas.core.reshape.re', 'pandas.core.reshape.reshape',
    'pandas.core.reshape.string', 'pandas.core.reshape.tile', 'pandas.core.reshape.util',
    'pandas.core.reshape.warnings', 'pandas.core.series', 'pandas.core.sorting',
    'pandas.core.sparse', 'pandas.core.sparse.__future__', 'pandas.core.sparse.api',
    'pandas.core.sparse.array', 'pandas.core.sparse.frame', 'pandas.core.sparse.list',
    'pandas.core.sparse.numpy', 'pandas.core.sparse.pandas',
    'pandas.core.sparse.scipy_sparse', 'pandas.core.sparse.series',
    'pandas.core.sparse.warnings', 'pandas.core.strings', 'pandas.core.sys',
    'pandas.core.textwrap', 'pandas.core.tokenize', 'pandas.core.tools',
    'pandas.core.tools.collections', 'pandas.core.tools.datetime',
    'pandas.core.tools.datetimes', 'pandas.core.tools.dateutil',
    'pandas.core.tools.numeric', 'pandas.core.tools.numpy', 'pandas.core.tools.pandas',
    'pandas.core.tools.timedeltas', 'pandas.core.types', 'pandas.core.util',
    'pandas.core.util.hashing', 'pandas.core.util.itertools', 'pandas.core.util.numpy',
    'pandas.core.util.pandas', 'pandas.core.warnings', 'pandas.core.weakref',
    'pandas.core.window', 'pandas.core.xlsxwriter', 'pandas.datetime', 'pandas.errors',
    'pandas.errors.pandas', 'pandas.io', 'pandas.io.__future__', 'pandas.io.abc',
    'pandas.io.api', 'pandas.io.clipboards', 'pandas.io.codecs',
    'pandas.io.collections', 'pandas.io.common', 'pandas.io.contextlib',
    'pandas.io.copy', 'pandas.io.csv', 'pandas.io.date_converters',
    'pandas.io.datetime', 'pandas.io.dateutil', 'pandas.io.distutils',
    'pandas.io.excel', 'pandas.io.feather_format', 'pandas.io.formats',
    'pandas.io.formats.__future__', 'pandas.io.formats.common',
    'pandas.io.formats.console', 'pandas.io.formats.csv', 'pandas.io.formats.distutils',
    'pandas.io.formats.format', 'pandas.io.formats.itertools',
    'pandas.io.formats.locale', 'pandas.io.formats.numpy', 'pandas.io.formats.os',
    'pandas.io.formats.pandas', 'pandas.io.formats.printing', 'pandas.io.formats.shutil',
    'pandas.io.formats.sys', 'pandas.io.formats.terminal', 'pandas.io.formats.textwrap',
    'pandas.io.functools', 'pandas.io.gbq', 'pandas.io.html', 'pandas.io.httplib',
    'pandas.io.itertools', 'pandas.io.json', 'pandas.io.json.collections',
    'pandas.io.json.copy', 'pandas.io.json.json', 'pandas.io.json.normalize',
    'pandas.io.json.numpy', 'pandas.io.json.os', 'pandas.io.json.pandas',
    'pandas.io.json.table_schema', 'pandas.io.mmap', 'pandas.io.msgpack',
    'pandas.io.msgpack._packer', 'pandas.io.msgpack._unpacker',
    'pandas.io.msgpack._version', 'pandas.io.msgpack.collections',
    'pandas.io.msgpack.exceptions', 'pandas.io.msgpack.os',
    'pandas.io.msgpack.pandas', 'pandas.io.numbers', 'pandas.io.numpy', 'pandas.io.os',
    'pandas.io.packers', 'pandas.io.pandas', 'pandas.io.parsers', 'pandas.io.pickle',
    'pandas.io.py', 'pandas.io.pytables', 'pandas.io.re', 'pandas.io.sas',
    'pandas.io.sas.pandas', 'pandas.io.sas.sasreader', 'pandas.io.sql', 'pandas.io.stata',
    'pandas.io.struct', 'pandas.io.sys', 'pandas.io.textwrap', 'pandas.io.time',
    'pandas.io.urllib', 'pandas.io.urllib2', 'pandas.io.urlparse', 'pandas.io.warnings',
    'pandas.io.zlib', 'pandas.json', 'pandas.pandas', 'pandas.plotting',
    'pandas.plotting.__future__', 'pandas.plotting._compat', 'pandas.plotting._converter',
    'pandas.plotting._core', 'pandas.plotting._misc', 'pandas.plotting._style',
    'pandas.plotting._tools', 'pandas.plotting.collections', 'pandas.plotting.contextlib',
    'pandas.plotting.cycler', 'pandas.plotting.datetime', 'pandas.plotting.dateutil',
    'pandas.plotting.distutils', 'pandas.plotting.math', 'pandas.plotting.matplotlib',
    'pandas.plotting.numpy', 'pandas.plotting.pandas', 'pandas.plotting.re',
    'pandas.plotting.warnings', 'pandas.stats', 'pandas.stats.__future__',
    'pandas.stats.api', 'pandas.stats.moments', 'pandas.stats.numpy', 'pandas.stats.pandas',
    'pandas.stats.warnings', 'pandas.sys', 'pandas.testing', 'pandas.tools',
    'pandas.tools.pandas', 'pandas.tools.plotting', 'pandas.tools.sys',
    'pandas.tools.warnings', 'pandas.tseries', 'pandas.tseries.api',
    'pandas.tseries.datetime', 'pandas.tseries.dateutil', 'pandas.tseries.frequencies',
    'pandas.tseries.functools', 'pandas.tseries.numpy', 'pandas.tseries.offsets',
    'pandas.tseries.operator', 'pandas.tseries.pandas', 'pandas.tseries.pytz',
    'pandas.tseries.re', 'pandas.tseries.warnings', 'pandas.util',
    'pandas.util.__future__', 'pandas.util._decorators', 'pandas.util._depr_module',
    'pandas.util._move', 'pandas.util._print_versions', 'pandas.util._tester',
    'pandas.util._validators', 'pandas.util.codecs', 'pandas.util.contextlib',
    'pandas.util.datetime', 'pandas.util.distutils', 'pandas.util.functools',
    'pandas.util.importlib', 'pandas.util.inspect', 'pandas.util.locale', 'pandas.util.numpy',
    'pandas.util.os', 'pandas.util.pandas', 'pandas.util.platform', 'pandas.util.pytest',
    'pandas.util.re', 'pandas.util.string', 'pandas.util.struct', 'pandas.util.subprocess',
    'pandas.util.sys', 'pandas.util.tempfile', 'pandas.util.testing', 'pandas.util.textwrap',
    'pandas.util.traceback', 'pandas.util.types', 'pandas.util.warnings', 'pandas.warnings',
]
h5py_imports = [
    'h5py', 'h5py._conv', 'h5py._errors', 'h5py._hl', 'h5py._hl.attrs',
    'h5py._hl.base', 'h5py._hl.compat', 'h5py._hl.dataset', 'h5py._hl.datatype',
    'h5py._hl.files', 'h5py._hl.filters', 'h5py._hl.group', 'h5py._hl.os',
    'h5py._hl.selections', 'h5py._hl.selections2', 'h5py._hl.six',
    'h5py._hl.sys', 'h5py._objects', 'h5py._proxy', 'h5py.defs',
    'h5py.functools', 'h5py.gc', 'h5py.h5', 'h5py.h5a', 'h5py.h5ac', 'h5py.h5d',
    'h5py.h5ds', 'h5py.h5f', 'h5py.h5fd', 'h5py.h5g', 'h5py.h5i', 'h5py.h5l',
    'h5py.h5o', 'h5py.h5p', 'h5py.h5py', 'h5py.h5r', 'h5py.h5s', 'h5py.h5t',
    'h5py.h5z', 'h5py.highlevel', 'h5py.numpy', 'h5py.operator', 'h5py.platform',
    'h5py.sys', 'h5py.tests', 'h5py.tests.__future__', 'h5py.tests.common',
    'h5py.tests.hl', 'h5py.tests.hl.test_attribute_create',
    'h5py.tests.hl.test_dataset_getitem', 'h5py.tests.hl.test_dataset_swmr',
    'h5py.tests.hl.test_datatype', 'h5py.tests.hl.test_dims_dimensionproxy',
    'h5py.tests.hl.test_file', 'h5py.tests.hl.test_threads', 'h5py.tests.old',
    'h5py.tests.old.test_attrs', 'h5py.tests.old.test_attrs_data',
    'h5py.tests.old.test_base', 'h5py.tests.old.test_dataset',
    'h5py.tests.old.test_datatype', 'h5py.tests.old.test_dimension_scales',
    'h5py.tests.old.test_file', 'h5py.tests.old.test_file_image',
    'h5py.tests.old.test_group', 'h5py.tests.old.test_h5',
    'h5py.tests.old.test_h5f', 'h5py.tests.old.test_h5p',
    'h5py.tests.old.test_h5t', 'h5py.tests.old.test_objects',
    'h5py.tests.old.test_selections', 'h5py.tests.old.test_slicing',
    'h5py.tests.sys', 'h5py.utils', 'h5py.version', 'h5py.warnings',
    'h5py.weakref',
]

matplotlib_imports = [
    'matplotlib', 'cycler', 'contourpy', 'dateutil',
]

excludes = [
    #'unittest',   # why cant I remove this?
    'wx', 'nose', 'Tkinter',

    'distutils', 'distutils.distutils', 'distutils.errors', 'distutils.os',
    'distutils.re', 'distutils.string', 'distutils.sys',
    'distutils.sysconfig', 'distutils.types', 'distutils.version',

    'beautifulsoup4', 'bitarray', 'bottleneck', 'bzip2', 'cdecimal',
    'cffi', 'comtypes', 'conda-build', 'configobj', 'console_shortcut',
    'cryptography', 'cython', 'dask', 'docutils', 'fastcache', 'flask',
    'freetype', 'funcsigs', 'greenlet', 'grin',

    'idna',
    'ipaddress', 'ipython-notebook', 'ipython-qtconsole',
    'ipython_genutils', 'itsdangerous', 'jedi', 'jinja2', 'jpeg',
    'jsonschema', 'jupyter_client', 'jupyter_core', 'launcher',
    'libsodium', 'markupsafe', 'mistune', 'multipledispatch',
    'nbformat', 'nltk', 'node-webkit', 'nose', 'patsy', 'pickleshare',
    'ply', 'pyasn1', 'pycosat', 'pycparser', 'pycrypto', 'pycurl',
    'pyflakes', 'pyopenssl', 'pyparsing', 'pyreadline',
    'python-dateutil', 'rope', 'scikit-image', 'simplegeneric',
    'singledispatch', 'sockjs-tornado', 'ssl_match_hostname',
    'statsmodels', 'sympy', 'tk', 'toolz', 'ujson', 'unicodecsv',
    'xlrd', 'xlwt', 'zeromq', 'zlib', 'alabaster',
    'anaconda-client', 'appdirs', 'astroid', 'astroid', 'astropy'
    'babel', 'backports_abc', 'blackwidow', 'blaze-core', 'bokeh',
    'boto', 'clyent', 'coverage',
    'curl', 'cytoolz', 'datashape', 'decorator', 'freeimage',
    'gl2ps', 'oce', 'pythonocc-core', 'tbb', 'enum34', 'et_xmlfile',
    'futures', 'gevent', 'gevent-websocket', 'hdf5', 'ipykernel',
    'ipython', 'ipywidgets', 'jdcal', 'jupyter', 'jupyter_console',
    'lazy-object-proxy', 'libtiff', 'llvmlite', 'logilab-common',
    'lxml', 'matplotlib', 'menuinst', 'MeshPy',
    'msvc_runtime', 'nbconvert', 'networkx', 'notebook', 'numba',

    'numexpr', 'numexpr.__config__', 'numexpr.__future__',
    'numexpr._winreg', 'numexpr.cpuinfo', 'numexpr.expressions',
    'numexpr.inspect', 'numexpr.interpreter', 'numexpr.necompiler',
    'numexpr.numexpr', 'numexpr.numpy', 'numexpr.operator', 'numexpr.os',
    'numexpr.pkg_resources', 'numexpr.platform', 'numexpr.re',
    'numexpr.subprocess', 'numexpr.sys', 'numexpr.tests',
    'numexpr.tests.numexpr', 'numexpr.tests.test_numexpr',
    'numexpr.threading', 'numexpr.types', 'numexpr.utils',
    'numexpr.version', 'numexpr.warnings',

    'numpydoc', 'odo', 'openmdao', 'openpyxl', 'openssl',

    '_pytest', '_pytest.__future__', '_pytest._code', '_pytest._code.__future__',
    '_pytest._code._ast', '_pytest._code._py2traceback', '_pytest._code.bisect',
    '_pytest._code.code', '_pytest._code.inspect', '_pytest._code.py', '_pytest._code.re',
    '_pytest._code.source', '_pytest._code.sys', '_pytest._code.tokenize',
    '_pytest._code.types', '_pytest._code.weakref', '_pytest._pluggy', '_pytest._pytest',
    '_pytest.argparse', '_pytest.assertion', '_pytest.assertion._ast',
    '_pytest.assertion._pytest', '_pytest.assertion.ast', '_pytest.assertion.collections',
    '_pytest.assertion.errno', '_pytest.assertion.fnmatch', '_pytest.assertion.imp',
    '_pytest.assertion.itertools', '_pytest.assertion.marshal', '_pytest.assertion.os',
    '_pytest.assertion.pprint', '_pytest.assertion.py', '_pytest.assertion.re',
    '_pytest.assertion.rewrite', '_pytest.assertion.struct', '_pytest.assertion.sys',
    '_pytest.assertion.types', '_pytest.assertion.util', '_pytest.bdb',
    '_pytest.cacheprovider', '_pytest.capture', '_pytest.codecs', '_pytest.collections',
    '_pytest.compat', '_pytest.config', '_pytest.contextlib', '_pytest.debugging',
    '_pytest.doctest', '_pytest.enum', '_pytest.fixtures', '_pytest.fnmatch',
    '_pytest.freeze_support', '_pytest.functools', '_pytest.helpconfig', '_pytest.hookspec',
    '_pytest.inspect', '_pytest.itertools', '_pytest.json', '_pytest.junitxml',
    '_pytest.main', '_pytest.mark', '_pytest.math', '_pytest.monkeypatch', '_pytest.nose',
    '_pytest.os', '_pytest.pastebin', '_pytest.platform', '_pytest.py', '_pytest.pytest',
    '_pytest.python', '_pytest.re', '_pytest.recwarn', '_pytest.resultlog', '_pytest.runner',
    '_pytest.setuponly', '_pytest.setupplan', '_pytest.shlex', '_pytest.skipping',
    '_pytest.sys', '_pytest.tempfile', '_pytest.terminal', '_pytest.time', '_pytest.tmpdir',
    '_pytest.traceback', '_pytest.types', '_pytest.unittest', '_pytest.vendored_packages',
    '_pytest.vendored_packages.inspect', '_pytest.vendored_packages.pluggy',
    '_pytest.vendored_packages.sys', '_pytest.warnings',

    'path.py', 'pep8', 'pi', 'pip',
    'psutil', 'py', 'pylint', 'pyside','pyside2', 'pytest',
    'pytools', 'pytz', 'pyyaml', 'pyzmq',
    'qtconsole', 'requests', 'ruamel_yaml', 'RunSnakeRun',
    'scikit-learn', 'setuptools',
    'snowballstemmer', 'sphinx', 'sphinx_rtd_theme', 'spyder',
    'spyder-app', 'sqlalchemy', 'sqlitedict', 'SquareMap', 'tornado',
    'traitlets', 'werkzeug', 'wheel',
    'wrapt', 'wxpython',

    'xlsxwriter', 'xlsxwriter.StringIO', 'xlsxwriter.app', 'xlsxwriter.chart',
    'xlsxwriter.chart_area', 'xlsxwriter.chart_bar', 'xlsxwriter.chart_column',
    'xlsxwriter.chart_doughnut', 'xlsxwriter.chart_line', 'xlsxwriter.chart_pie',
    'xlsxwriter.chart_radar', 'xlsxwriter.chart_scatter', 'xlsxwriter.chart_stock',
    'xlsxwriter.chartsheet', 'xlsxwriter.codecs', 'xlsxwriter.collections',
    'xlsxwriter.comments', 'xlsxwriter.compatibility', 'xlsxwriter.contenttypes',
    'xlsxwriter.copy', 'xlsxwriter.core', 'xlsxwriter.custom', 'xlsxwriter.datetime',
    'xlsxwriter.decimal', 'xlsxwriter.drawing', 'xlsxwriter.format',
    'xlsxwriter.fractions', 'xlsxwriter.io', 'xlsxwriter.operator', 'xlsxwriter.os',
    'xlsxwriter.packager', 'xlsxwriter.re', 'xlsxwriter.relationships',
    'xlsxwriter.shape', 'xlsxwriter.sharedstrings', 'xlsxwriter.shutil',
    'xlsxwriter.stat', 'xlsxwriter.struct', 'xlsxwriter.styles', 'xlsxwriter.sys',
    'xlsxwriter.table', 'xlsxwriter.tempfile', 'xlsxwriter.theme',
    'xlsxwriter.utility', 'xlsxwriter.vml', 'xlsxwriter.warnings',
    'xlsxwriter.workbook', 'xlsxwriter.worksheet', 'xlsxwriter.xmlwriter',
    'xlsxwriter.zipfile',

    'xlwings',

    # not required...strange...
    'conda', 'conda-env', 'pywin32', 'python', 'vs2008_runtime',
    #'pyqt', 'anaconda',

    # things we are using
    'libpng',
    #'sip', 'colorama', 'numpy', 'pillow', 'qt',
    #'vtk', 'six', 'mkl', 'mkl-service',
]
if INCLUDE_PANDAS:
    excludes.extend(pandas_imports)
if not INCLUDE_BDFV:
    excludes.append('pytables')
if not INCLUDE_MATPLOTLIB:
    excludes.extend(matplotlib_imports)
#excludes = ['vtkmodules.all', ]
excludes = []

qt_imports = [
    'PyQt5.QtOpenGL', 'PyQt5.Qsci',
    'PyQt5.QtPrintSupport',  # needed for Qsci
]
# 9.1
#'vtkmodules', 'vtkmodules.all', 'vtkmodules.qt', 'vtkmodules.util', 'vtkmodules.util.misc', 'vtkmodules.util.numpy_support',
#'vtkmodules.util.vtkConstants', 'vtkmodules.util.vtkVariant',
#'vtkmodules.vtkAcceleratorsVTKmCore', 'vtkmodules.vtkAcceleratorsVTKmDataModel', 'vtkmodules.vtkAcceleratorsVTKmFilters', 'vtkmodules.vtkChartsCore', 'vtkmodules.vtkCommonColor', 'vtkmodules.vtkCommonComputationalGeometry', 'vtkmodules.vtkCommonCore', 'vtkmodules.vtkCommonDataModel', 'vtkmodules.vtkCommonExecutionModel', 'vtkmodules.vtkCommonMath', 'vtkmodules.vtkCommonMisc', 'vtkmodules.vtkCommonPython', 'vtkmodules.vtkCommonSystem', 'vtkmodules.vtkCommonTransforms', 'vtkmodules.vtkDomainsChemistry', 'vtkmodules.vtkDomainsChemistryOpenGL2', 'vtkmodules.vtkFiltersAMR', 'vtkmodules.vtkFiltersCore', 'vtkmodules.vtkFiltersExtraction', 'vtkmodules.vtkFiltersFlowPaths', 'vtkmodules.vtkFiltersGeneral', 'vtkmodules.vtkFiltersGeneric', 'vtkmodules.vtkFiltersGeometry', 'vtkmodules.vtkFiltersHybrid', 'vtkmodules.vtkFiltersHyperTree', 'vtkmodules.vtkFiltersImaging', 'vtkmodules.vtkFiltersModeling', 'vtkmodules.vtkFiltersParallel', 'vtkmodules.vtkFiltersParallelDIY2', 'vtkmodules.vtkFiltersParallelImaging', 'vtkmodules.vtkFiltersParallelStatistics', 'vtkmodules.vtkFiltersPoints', 'vtkmodules.vtkFiltersProgrammable', 'vtkmodules.vtkFiltersPython', 'vtkmodules.vtkFiltersSMP', 'vtkmodules.vtkFiltersSelection', 'vtkmodules.vtkFiltersSources', 'vtkmodules.vtkFiltersStatistics', 'vtkmodules.vtkFiltersTexture', 'vtkmodules.vtkFiltersTopology', 'vtkmodules.vtkFiltersVerdict', 'vtkmodules.vtkGeovisCore', 'vtkmodules.vtkIOAMR', 'vtkmodules.vtkIOAsynchronous', 'vtkmodules.vtkIOCGNSReader', 'vtkmodules.vtkIOCONVERGECFD', 'vtkmodules.vtkIOChemistry', 'vtkmodules.vtkIOCityGML', 'vtkmodules.vtkIOCore', 'vtkmodules.vtkIOEnSight', 'vtkmodules.vtkIOExodus', 'vtkmodules.vtkIOExport', 'vtkmodules.vtkIOExportGL2PS', 'vtkmodules.vtkIOExportPDF', 'vtkmodules.vtkIOGeoJSON', 'vtkmodules.vtkIOGeometry', 'vtkmodules.vtkIOH5Rage', 'vtkmodules.vtkIOH5part', 'vtkmodules.vtkIOHDF', 'vtkmodules.vtkIOIOSS', 'vtkmodules.vtkIOImage', 'vtkmodules.vtkIOImport', 'vtkmodules.vtkIOInfovis', 'vtkmodules.vtkIOLSDyna', 'vtkmodules.vtkIOLegacy', 'vtkmodules.vtkIOMINC', 'vtkmodules.vtkIOMotionFX', 'vtkmodules.vtkIOMovie', 'vtkmodules.vtkIONetCDF', 'vtkmodules.vtkIOOMF', 'vtkmodules.vtkIOOggTheora', 'vtkmodules.vtkIOPIO', 'vtkmodules.vtkIOPLY', 'vtkmodules.vtkIOParallel', 'vtkmodules.vtkIOParallelExodus', 'vtkmodules.vtkIOParallelLSDyna', 'vtkmodules.vtkIOParallelXML', 'vtkmodules.vtkIOSQL', 'vtkmodules.vtkIOSegY', 'vtkmodules.vtkIOTRUCHAS', 'vtkmodules.vtkIOTecplotTable', 'vtkmodules.vtkIOVPIC', 'vtkmodules.vtkIOVeraOut', 'vtkmodules.vtkIOVideo', 'vtkmodules.vtkIOXML', 'vtkmodules.vtkIOXMLParser', 'vtkmodules.vtkIOXdmf2', 'vtkmodules.vtkImagingColor', 'vtkmodules.vtkImagingCore', 'vtkmodules.vtkImagingFourier', 'vtkmodules.vtkImagingGeneral', 'vtkmodules.vtkImagingHybrid', 'vtkmodules.vtkImagingMath', 'vtkmodules.vtkImagingMorphological', 'vtkmodules.vtkImagingOpenGL2', 'vtkmodules.vtkImagingSources', 'vtkmodules.vtkImagingStatistics', 'vtkmodules.vtkImagingStencil', 'vtkmodules.vtkInfovisCore', 'vtkmodules.vtkInfovisLayout', 'vtkmodules.vtkInteractionImage', 'vtkmodules.vtkInteractionStyle', 'vtkmodules.vtkInteractionWidgets', 'vtkmodules.vtkParallelCore', 'vtkmodules.vtkPythonContext2D', 'vtkmodules.vtkRenderingAnnotation', 'vtkmodules.vtkRenderingContext2D', 'vtkmodules.vtkRenderingContextOpenGL2', 'vtkmodules.vtkRenderingCore', 'vtkmodules.vtkRenderingExternal', 'vtkmodules.vtkRenderingFreeType', 'vtkmodules.vtkRenderingGL2PSOpenGL2', 'vtkmodules.vtkRenderingImage', 'vtkmodules.vtkRenderingLICOpenGL2', 'vtkmodules.vtkRenderingLOD', 'vtkmodules.vtkRenderingLabel', 'vtkmodules.vtkRenderingMatplotlib', 'vtkmodules.vtkRenderingOpenGL2', 'vtkmodules.vtkRenderingParallel', 'vtkmodules.vtkRenderingSceneGraph', 'vtkmodules.vtkRenderingUI', 'vtkmodules.vtkRenderingVR', 'vtkmodules.vtkRenderingVolume', 'vtkmodules.vtkRenderingVolumeAMR', 'vtkmodules.vtkRenderingVolumeOpenGL2', 'vtkmodules.vtkRenderingVtkJS', 'vtkmodules.vtkTestingRendering', 'vtkmodules.vtkViewsContext2D', 'vtkmodules.vtkViewsCore', 'vtkmodules.vtkViewsInfovis', 'vtkmodules.vtkWebCore', 'vtkmodules.vtkWebGLExporter'

vtk_imports_all = [  # 9.3
    'vtk',
    'vtkmodules', 'vtkmodules.all',
    'vtkmodules.qt', 'vtkmodules.util', 'vtkmodules.util.misc', 'vtkmodules.util.numpy_support',
    'vtkmodules.util.vtkConstants', 'vtkmodules.util.vtkVariant', 'vtkmodules.vtkAcceleratorsVTKmCore',
    'vtkmodules.vtkAcceleratorsVTKmDataModel', 'vtkmodules.vtkAcceleratorsVTKmFilters',
    'vtkmodules.vtkChartsCore', 'vtkmodules.vtkCommonColor', 'vtkmodules.vtkCommonComputationalGeometry',
    'vtkmodules.vtkCommonCore', 'vtkmodules.vtkCommonDataModel', 'vtkmodules.vtkCommonExecutionModel',
    'vtkmodules.vtkCommonMath', 'vtkmodules.vtkCommonMisc', 'vtkmodules.vtkCommonPython', 'vtkmodules.vtkCommonSystem',
    'vtkmodules.vtkCommonTransforms', 'vtkmodules.vtkDomainsChemistry', 'vtkmodules.vtkDomainsChemistryOpenGL2',
    'vtkmodules.vtkFiltersAMR', 'vtkmodules.vtkFiltersCellGrid', 'vtkmodules.vtkFiltersCore',
    'vtkmodules.vtkFiltersExtraction', 'vtkmodules.vtkFiltersFlowPaths', 'vtkmodules.vtkFiltersGeneral',
    'vtkmodules.vtkFiltersGeneric', 'vtkmodules.vtkFiltersGeometry', 'vtkmodules.vtkFiltersGeometryPreview',
    'vtkmodules.vtkFiltersHybrid', 'vtkmodules.vtkFiltersHyperTree', 'vtkmodules.vtkFiltersImaging',
    'vtkmodules.vtkFiltersModeling', 'vtkmodules.vtkFiltersParallel', 'vtkmodules.vtkFiltersParallelDIY2',
    'vtkmodules.vtkFiltersParallelImaging', 'vtkmodules.vtkFiltersParallelStatistics',
    'vtkmodules.vtkFiltersPoints', 'vtkmodules.vtkFiltersProgrammable', 'vtkmodules.vtkFiltersPython',
    'vtkmodules.vtkFiltersReduction', 'vtkmodules.vtkFiltersSMP', 'vtkmodules.vtkFiltersSelection',
    'vtkmodules.vtkFiltersSources', 'vtkmodules.vtkFiltersStatistics', 'vtkmodules.vtkFiltersTensor',
    'vtkmodules.vtkFiltersTexture', 'vtkmodules.vtkFiltersTopology', 'vtkmodules.vtkFiltersVerdict',
    'vtkmodules.vtkGeovisCore', 'vtkmodules.vtkIOAMR', 'vtkmodules.vtkIOAsynchronous',
    'vtkmodules.vtkIOCGNSReader', 'vtkmodules.vtkIOCONVERGECFD', 'vtkmodules.vtkIOCellGrid',
    'vtkmodules.vtkIOCesium3DTiles', 'vtkmodules.vtkIOChemistry', 'vtkmodules.vtkIOCityGML',
    'vtkmodules.vtkIOCore', 'vtkmodules.vtkIOEnSight', 'vtkmodules.vtkIOExodus', 'vtkmodules.vtkIOExport',
    'vtkmodules.vtkIOExportGL2PS', 'vtkmodules.vtkIOExportPDF', 'vtkmodules.vtkIOFLUENTCFF',
    'vtkmodules.vtkIOGeoJSON', 'vtkmodules.vtkIOGeometry', 'vtkmodules.vtkIOH5Rage',
    'vtkmodules.vtkIOH5part', 'vtkmodules.vtkIOHDF', 'vtkmodules.vtkIOIOSS', 'vtkmodules.vtkIOImage',
    'vtkmodules.vtkIOImport', 'vtkmodules.vtkIOInfovis', 'vtkmodules.vtkIOLSDyna', 'vtkmodules.vtkIOLegacy',
    'vtkmodules.vtkIOMINC', 'vtkmodules.vtkIOMotionFX', 'vtkmodules.vtkIOMovie', 'vtkmodules.vtkIONetCDF',
    'vtkmodules.vtkIOOMF', 'vtkmodules.vtkIOOggTheora', 'vtkmodules.vtkIOPIO', 'vtkmodules.vtkIOPLY',
    'vtkmodules.vtkIOParallel', 'vtkmodules.vtkIOParallelExodus', 'vtkmodules.vtkIOParallelLSDyna',
    'vtkmodules.vtkIOParallelXML', 'vtkmodules.vtkIOSQL', 'vtkmodules.vtkIOSegY', 'vtkmodules.vtkIOTRUCHAS',
    'vtkmodules.vtkIOTecplotTable', 'vtkmodules.vtkIOVPIC', 'vtkmodules.vtkIOVeraOut', 'vtkmodules.vtkIOVideo',
    'vtkmodules.vtkIOXML', 'vtkmodules.vtkIOXMLParser', 'vtkmodules.vtkIOXdmf2', 'vtkmodules.vtkImagingColor',
    'vtkmodules.vtkImagingCore', 'vtkmodules.vtkImagingFourier', 'vtkmodules.vtkImagingGeneral',
    'vtkmodules.vtkImagingHybrid', 'vtkmodules.vtkImagingMath', 'vtkmodules.vtkImagingMorphological',
    'vtkmodules.vtkImagingOpenGL2', 'vtkmodules.vtkImagingSources', 'vtkmodules.vtkImagingStatistics',
    'vtkmodules.vtkImagingStencil', 'vtkmodules.vtkInfovisCore', 'vtkmodules.vtkInfovisLayout',
    'vtkmodules.vtkInteractionImage', 'vtkmodules.vtkInteractionStyle', 'vtkmodules.vtkInteractionWidgets',
    'vtkmodules.vtkParallelCore', 'vtkmodules.vtkPythonContext2D', 'vtkmodules.vtkRenderingAnnotation',
    'vtkmodules.vtkRenderingCellGrid', 'vtkmodules.vtkRenderingContext2D', 'vtkmodules.vtkRenderingContextOpenGL2',
    'vtkmodules.vtkRenderingCore', 'vtkmodules.vtkRenderingExternal', 'vtkmodules.vtkRenderingFreeType',
    'vtkmodules.vtkRenderingGL2PSOpenGL2', 'vtkmodules.vtkRenderingHyperTreeGrid', 'vtkmodules.vtkRenderingImage',
    'vtkmodules.vtkRenderingLICOpenGL2', 'vtkmodules.vtkRenderingLOD', 'vtkmodules.vtkRenderingLabel',
    'vtkmodules.vtkRenderingMatplotlib', 'vtkmodules.vtkRenderingOpenGL2', 'vtkmodules.vtkRenderingParallel',
    'vtkmodules.vtkRenderingSceneGraph', 'vtkmodules.vtkRenderingUI', 'vtkmodules.vtkRenderingVR',
    'vtkmodules.vtkRenderingVolume', 'vtkmodules.vtkRenderingVolumeAMR', 'vtkmodules.vtkRenderingVolumeOpenGL2',
    'vtkmodules.vtkRenderingVtkJS', 'vtkmodules.vtkTestingRendering', 'vtkmodules.vtkViewsContext2D',
    'vtkmodules.vtkViewsCore', 'vtkmodules.vtkViewsInfovis', 'vtkmodules.vtkWebCore', 'vtkmodules.vtkWebGLExporter',
    # tested - 9.4.1
    # Unable to find a valid OpenGL 3.2 or later implementation
    'vtkmodules.util.data_model', 'vtkmodules.util.execution_model',
]

import vtk
vtk_version_str = vtk.VTK_VERSION
print(f'vtk_version_str = {vtk.VTK_VERSION!r}') # '9.3.1'

if vtk_version_str >= '9.4.0':
    raise RuntimeError('use vtk<9.4')

# can we get rid of scipy.optimize?
hiddenimports = [
    #'vtk', 'vtk.vtkCommonPythonSIP', 'vtk.vtkFilteringPythonSIP',
    #'vtk.vtkRenderingPythonSIP',
    'scipy._lib.messagestream', # 'scipy',
    'pygments',
] + qt_imports + vtk_imports_all
if INCLUDE_BDFV:
    hiddenimports += ['pytables'] + pandas_imports


if not IS_H5PY:
    excludes += h5py_imports

excludes = None
a = Analysis(analyze_files,
             pathex=pathex,
             excludes=excludes,
             hiddenimports=hiddenimports,
             hookspath=None,
             noarchive=IS_ONEDIR)
pyz = PYZ(a.pure)

png_filename = os.path.join('..', 'pyNastran', 'gui', 'images', 'bwb_stress_black.png')
assert os.path.exists(png_filename), png_filename

#print("help(EXE) = \n", help(EXE))
if BUILD_EXE:
    exe_binaries = a.binaries + binaries + icons
    if DEBUG:
        # v: add a debug message each time a module is initialized
        # u: unbuffered stdio
        exe_binaries += [('v', None, 'OPTION'), ('u', None, 'OPTION')]
    if MAKE_SPLASH:
        splash = Splash(png_filename,
                        binaries=a.binaries,
                        datas=a.datas,
                        text_pos=(10, 50),
                        text_size=14,  # was 12
                        text_color='white')  # was black
        exe = EXE(pyz,
                  a.scripts,
                  splash,
                  splash.binaries,
                  exe_binaries, # was never active..
                  a.zipfiles,
                  a.datas,
                  exclude_binaries=IS_ONEDIR,
                  name=os.path.join('build\\pyi.win32\\pyNastranGUI', 'pyNastranGUI.exe'),
                  debug=DEBUG,
                  strip=None, #
                  disable_windowed_traceback=False,
                  #upx=True,  # some compression thing; uses upx-dir defined somewhere...
                  icon=icon_main,
                  console=True )
    else:
        exe = EXE(pyz,
                  a.scripts,
                  exe_binaries,
                  a.zipfiles,
                  a.datas,
                  exclude_binaries=IS_ONEDIR,
                  name=os.path.join('build\\pyi.win32\\pyNastranGUI', 'pyNastranGUI.exe'),
                  debug=DEBUG,
                  disable_windowed_traceback=False,
                  strip=None, #
                  #upx=True,  # some compression thing; uses upx-dir defined somewhere...
                  icon=icon_main,
                  console=True )

print('*'*80)
#print("help(COLLECT) = \n",help(COLLECT))
if BUILD_COLLECT:
    coll = COLLECT(exe,
                   a.binaries + binaries + icons,
                   a.zipfiles,
                   a.datas,
                   exclude_binaries=1,
                   #icon=icon_main,
                   strip=None,
                   upx=True,
                   name=os.path.join('dist', 'pyNastranGUI'))

#-------------------------------------------------------------------------
# fix the __init__.py file

shutil.copyfile(init_filename_bkp, init_filename)
os.remove(init_filename_bkp)
