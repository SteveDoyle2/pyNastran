import os
import sys

# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
if is_pynastrangui_exe:
    # pyInstaller
    from pyNastran.version import __version__, __releaseDate__
else:
    __version__ = '1.0.0'
    __releaseDate__ = '2017/5/25'
    __releaseDate2__ = 'MAY 25, 2017'

__author__  = 'Steven Doyle, Saullo Castro, hurlei, Paul Blelloch, Nikita Kalutsky'
__email__   = 'mesheb82@gmail.com'
__desc__    = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__license__     = 'BSD-3'
__copyright__   = 'Copyright %s; 2011-2017' % __license__
__pyside_copyright__ = 'Copyright LGPLv3 - pySide'
__pyqt_copyright__ = 'Copyright GPLv3 - PyQt'
__website__     = 'https://github.com/SteveDoyle2/pyNastran'

is_release = True  ## turns on skipping of tables that aren't supported
