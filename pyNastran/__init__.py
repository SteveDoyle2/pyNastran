import os
import sys

# this variable is automatically set by the .spec file; should be False
is_pynastrangui_exe = False
__version__ = '0.8.0'
__releaseDate__ = '2016/8/21'
__releaseDate2__ = 'AUGUST 21, 2016'

__author__  = 'Steven Doyle, Al Danials, Marcin Gasiorek, hurlei, saullocastro, Paul Blelloch, Nikita Kalutsky'
__email__   = 'mesheb82@gmail.com'
__desc__    = 'Nastran BDF/F06/OP2/OP4 File reader/editor/writer/viewer'
__license__     = 'LGPLv3'
__copyright__   = 'Copyright %s; 2011-2016' % __license__
__pyqt_copyright__ = 'Copyright GPLv3; 2011-2016'
__website__     = 'https://github.com/SteveDoyle2/pyNastran'

is_release = True  ## turns on skipping of tables that aren't supported
