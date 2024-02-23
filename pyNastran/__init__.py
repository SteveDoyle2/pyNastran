# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
is_installed = True

import sys

# this is still a requirement, but disabling it so readthedocs works
if sys.version_info < (3, 9):  # pragma: no cover
    IMAJOR, MINOR1, MINOR2 = sys.version_info[:3]
    raise ImportError('Upgrade your Python to >= 3.9.0; version=(%s.%s.%s)' % (
        IMAJOR, MINOR1, MINOR2))

__version_release__ = '1.4.0'

# only for release; 1.4.0
__version__ = __version_release__
__releaseDate__ = '2024/2/22'
__releaseDate2__ = 'FEBRUARY 22, 2024'

__author__ = 'Steven Doyle'
__email__ = 'mesheb82@gmail.com'
__desc__ = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__license__ = 'BSD-3'
__copyright__ = f'Copyright {__license__}; 2011-2024'
__pyside_copyright__ = 'Copyright LGPLv3 - pySide'
__pyqt_copyright__ = 'Copyright GPLv3 - PyQt'
__website__ = 'https://github.com/SteveDoyle2/pyNastran'

DEV = 'dev' in __version__
if DEV:
    __docs__ = 'https://pynastran-git.readthedocs.io/en/latest/quick_start/index.html'
else:
    # 1.3
    # we don't do separate doc releases for 1.3 vs 1.3.1
    __docs__ = f'https://pynastran-git.readthedocs.io/en/{__version__[:3]}/quick_start/index.html'

__issue__ = 'https://github.com/SteveDoyle2/pyNastran/issues'
__discussion_forum__ = 'https://groups.google.com/forum/#!forum/pynastran-discuss'

is_release = True  ## True=turns on skipping of tables that aren't supported
